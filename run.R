#------------------------------------------------------------------------------#
#                   Code to simulate data and fit models
#------------------------------------------------------------------------------#

library(multivar) # subgrouping feature available on GitHub (https://github.com/zackfisher/multivar)
library(gimme)
library(igraph)

# ------------------------------------------------- #
#  Generate conditions
# ------------------------------------------------- #
k     <- c(50,100)        # number of individuals
d     <- c(10)            # number of variables
n     <- c(50,100)        # time series length
s     <- c(2,3)           # number of subgroups
iter  <- c(1:50)          # number of replications
bal   <- c("bal","unbal") # balanced or unbalanced subgroup membership

cond  <- expand.grid(iter,k,d,n,s,bal) 

ncond <- c(rep(3,50),rep(4,50),rep(5,50),rep(6,50),
           rep(9,50),rep(10,50),rep(11,50),rep(12,50),
           rep(15,50),rep(16,50),rep(17,50),rep(18,50),
           rep(21,50),rep(22,50),rep(23,50),rep(24,50))

cond  <- cbind(cond, ncond) 

# ------------------------------------------------- #
#  Simulate data and fit models
# ------------------------------------------------- #
reps <- nrow(cond)

for (i in 1:reps){
  
  # select condition
  iter <- as.numeric(cond[i,1]) 
  k    <- as.numeric(cond[i,2]) 
  d    <- as.numeric(cond[i,3]) 
  n    <- as.numeric(cond[i,4])
  s    <- as.numeric(cond[i,5]) 
  bal  <- as.character(cond[i,6])
  c    <- as.numeric(cond[i,7])
  
  # set base filename
  base_filename  <- paste0("_cond",c,"_k",k,"_d",d,"_n",n,"_s",s,"_",bal,"_",iter,".RDS")
  
  # set up true subgroup assignment vector
  if (bal == "bal" & s == 2){
    subgroup <- c(rep(1,k/2), rep(2,k/2))
  } else if (bal == "bal" & s == 3 & k == 50){
    subgroup <- c(rep(1,16), rep(2,17), rep(3,17))
  } else if (bal == "bal" & s == 3 & k == 100){
    subgroup <- c(rep(1,33), rep(2,33), rep(3,34))
  } else if (bal == "unbal" & s == 2){
    subgroup <- c(rep(1,k*0.3), rep(2,k*0.7))
  } else if (bal == "unbal" & s == 3){
    subgroup <- c(rep(1,k*0.2), rep(2,k*0.2), rep(3,k*0.6))
  }
  
  # check to make sure subgroup assignments are correct
  if(length(subgroup) != k) stop("Length of subgroup does not equal k")
  if(length(unique(subgroup)) != s) stop("Number of subgroups does not equal s")
  
  # set proportion of nonzero paths
  p_com <- d/d^2 # proportion of common paths (AR coefs)
  p_sub <- 0.05  # proportion of subgroup paths
  p_ind <- 0.05  # proportion of individual paths 
  
  # ---- simulate data ----
  set.seed(iter)
  
  nonstationary <- TRUE
  
  while(nonstationary){
    
    # set up true matrices
    true_com <- matrix(0, d, d, byrow = TRUE)
    true_sub <- array(0, dim = c(d, d, k))
    true_ind <- array(0, dim = c(d, d, k))
    
    # set up elements for subgrouping
    diag_pos <- 1 + 0:(d - 1) * (d + 1)
    com_pos  <- sample(x = diag_pos, size = p_com*d^2)
    
    sub_pos  <- sample(x = setdiff(1:d^2, c(diag_pos, com_pos)), size = p_sub*s*d^2)
    sub_pos  <- split(sub_pos, factor(sort(rank(sub_pos) %% s)))
    sub_pos  <- lapply(1:k, function(i){sub_pos[[subgroup[i]]]})
    
    ind_pos  <- lapply(seq_along(1:k), function(i){
      sample(x = setdiff(1:d^2, c(diag_pos, com_pos, unlist(sub_pos))), size = p_ind*d^2, replace = F)
    })
    
    # fill matrices
    true_com[com_pos] <- runif(length(com_pos), 0, 1)
    
    true_sub <- lapply(seq_along(1:k), function(j){
      true_sub[,,j][sub_pos[[j]]] <- runif(length(sub_pos[[j]]), 0, 1)
      true_sub[,,j]
    })
    
    true_ind <- lapply(seq_along(1:k), function(j){
      true_ind[,,j][ind_pos[[j]]] <- runif(length(ind_pos[[j]]), 0, 1)
      true_ind[,,j]
    })
    
    mats <- lapply(seq_along(true_ind), function(j){
      true_com + true_sub[[j]] + true_ind[[j]]
    })
    
    # check eigenvalues
    nonstationary <- any(unlist(lapply(mats, function(x){
      max_eig <- max(abs(eigen(x)$values))
      if(max_eig > .99){
        TRUE
      } else {
        FALSE
      }
    })))
    
  }
  
  # simulate time series
  sim  <- multivar::multivar_sim(
    k = k,   # number of individuals
    d = d,   # number of variables
    n = n,   # time series length
    mat_total = mats,
    sigma = diag(d) # innovation covariance matrix
  )
  
  sim$mat_com         <- true_com
  sim$mat_ind_unique  <- true_ind
  sim$mat_sub_unique  <- true_sub
  sim$mat_ind_final   <- mats
  
  saveRDS(sim, file.path("data", paste0("data", base_filename)))
  
  # ---- fit standard multivar ----
  fit_nosub <- constructModel(data = sim$data, nlambda1 = 10, nlambda2 = 10, pendiag = FALSE)
  fit_nosub <- cv.multivar(fit_nosub)
  
  saveRDS(fit_nosub, file.path("results", paste0("res_nosub", base_filename)))
  
  # ---- fit confirmatory subgrouping multivar ----
  fit_conf <- constructModel(data = sim$data, nlambda1 = 10, nlambda2 = 10, ntau = 10,
                             subgroup_membership = subgroup, pendiag = FALSE)
  fit_conf <- cv.multivar(fit_conf)
  
  saveRDS(fit_conf, file.path("results", paste0("res_conf", base_filename)))
  
  # ---- fit subgrouping multivar ----
  fit_walkunique <- constructModel(data = sim$data, nlambda1 = 10, nlambda2 = 10, ntau = 10,
                                   subgroup = TRUE, pendiag = FALSE)
  fit_walkunique <- cv.multivar(fit_walkunique)
  
  saveRDS(fit_walkunique, file.path("results", paste0("res_walkunique", base_filename)))
  
  # ---- fit sgimme ----
  sgimme_fit <- gimme(data= sim$data, subgroup = TRUE, VAR = TRUE, sub_feature = "lagged")
  
  saveRDS(sgimme_fit, file.path("results", paste0("res_sgimme", base_filename)))
  
}

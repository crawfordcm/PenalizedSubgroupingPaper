#------------------------------------------------------------------------------#
#                     Code to analyze simulated data
#------------------------------------------------------------------------------#

library("SimDesign")
library("caret")
library("mclust")

source("calculate_measures.R")

# ------------------------------------------------- #
#  Analyze
# ------------------------------------------------- #
kj       <- c(50,100) 
dm       <- c(10)     
nl       <- c(50,100)
su       <- c(2,3)
datasets <- 1:50
mods     <- c("nosub","conf","sgimme","sgcvar","walkunique")
cond     <- c(1:24)
bal      <- c("bal","unbal")

cnt <- 1
res <- list()

for(g in c(1:length(mods))){
  for(c in c(1:length(cond))){
    for(j in c(1:length(kj))){
      for(m in c(1:length(dm))){
        for(l in c(1:length(nl))){
          for(u in c(1:length(su))){
            for(b in c(1:length(bal))){
              for(i in c(1:length(datasets))){
                
                # select condition
                k    <- kj[j]
                d    <- dm[m]
                n    <- nl[l]
                s    <- su[u]
                co   <- cond[c]
                ba   <- bal[b]
                it   <- datasets[i]
                mod  <- mods[g]
                
                # set base filename
                base_filename  <- paste0("_cond",co,"_k",k,"_d",d,"_n",n,"_s",s,"_",ba,"_",it,".RDS")
                
                # setup true subgroup membership vector
                if (ba == "bal" & s == 2){
                  subgroup <- c(rep(1, k/2), rep(2, k/2))
                } else if (ba == "bal" & s == 3 & k == 50){
                  subgroup <- c(rep(1, 16), rep(2, 17), rep(3, 17))
                } else if (ba == "bal" & s == 3 & k == 100){
                  subgroup <- c(rep(1, 33), rep(2, 33), rep(3, 34))
                } else if (ba == "unbal" & s == 2){
                  subgroup <- c(rep(1, k*0.3), rep(2, k*0.7))
                } else if (ba == "unbal" & s == 3){
                  subgroup <- c(rep(1, k*0.2), rep(2, k*0.2), rep(3, k*0.6))
                }
                
                # read in simulated data and results
                if(file.exists(file.path("results", paste0("res_", mod, base_filename)))){
                  
                  sim <- readRDS(file.path("data", paste0("data", base_filename)))
                  fit <- readRDS(file.path("results", paste0("res_", mod, base_filename)))
                  
                  # extract true transition matrices
                  true_mats <- sim$mat_ind_final
                  
                  # extract estimated matrices
                  if(mod == "sgimme"){
                    
                    est_mats <- fit$path_est_mats
                    est_sub  <- fit$fit$sub_membership
                    
                  } else if(mod == "sgcvar"){
                    
                    est_mats <- fit$est_mats
                    est_sub  <- fit$recovered_subgroups
                    
                  } else{
                    
                    est_mats <- fit$mats$total
                    est_sub  <- fit$obj@subgroup
                    
                  }
                  
                  # true subgroup membership
                  true_sub <- subgroup
                  
                  # compute summary statistics
                  res[[cnt]] <- as.data.frame(do.call("rbind", calculate_measures(
                    true_mats = true_mats, 
                    est_mats = est_mats, 
                    dataset = i, 
                    n = n,
                    d = d,
                    k = k,
                    s = s,
                    alabel = ba,
                    mods = mod,
                    cond = co,
                    true_sub = true_sub,
                    est_sub = est_sub,
                    pendiag = FALSE
                  )))
                  
                } else {
                  
                  # do nothing
                  
                }
                
                cnt <- cnt + 1
                
              }
            }
          }
        }
      }
    }
  }
}

# ------------------------------------------------- #
#  Summarize results
# ------------------------------------------------- #
# combine results
df <- as.data.frame(do.call("rbind",res))

# make columns numeric
num_columns <- c("dataset", "n", "d", "k", "s","mcc", "Accuracy", "Kappa", "AccuracyLower", 
                 "AccuracyUpper", "AccuracyNull", "AccuracyPValue", "McnemarPValue", 
                 "Sensitivity", "Specificity", "Pos Pred Value", "Neg Pred Value", 
                 "Precision", "Recall", "F1", "Prevalence", "Detection Rate", 
                 "Detection Prevalence", "Balanced Accuracy", "biasabs","bias","biasabs_nonzero",
                 "biasabs_zero","bias_nonzero","bias_zero", "biasrel","ari","fail","rmse","sub_num")

df[, num_columns] <- lapply(num_columns, function(x) as.numeric(df[[x]]))

# save full results
saveRDS(df, file.path("files", "results_all.RDS"))

# summarize results across conditions
dt_table <- data.table::data.table(df)[, list(
  Nreps = .N,
  sens = mean(Sensitivity),
  sd_sens = sd(Sensitivity),
  spec = mean(Specificity),
  sd_spec = sd(Specificity),
  mcc = mean(mcc),
  sd_mcc = sd(mcc),
  biasabs = mean(biasabs),
  sd_biasabs = sd(biasabs),
  rmse = mean(rmse),
  sd_rmse = sd(rmse),
  ari = mean(ari),
  sd_ari = sd(ari),
  sub_num = mean(sub_num),
  fail = sum(fail),
  failperc = (sum(fail)/.N)*100
), by = c("label","mods","n","d","k","s")]

# save summary 
saveRDS(dt_table, file.path("files", "results_summary.RDS"))




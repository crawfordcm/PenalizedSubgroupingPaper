# ------------------------------------------------- #
# Summary statistics function
# ------------------------------------------------- #
calculate_measures <- function(true_mats, # list of true matrices
                               est_mats,  # list of estimated matrices
                               dataset,   # dataset number
                               n,         # time series length
                               d,         # number of variables 
                               k,         # number of individuals
                               s,         # number of subgroups
                               alabel,    # balanced/unbalanced membership
                               mods,      # model
                               cond,      # condition number
                               true_sub,  # true subgroup membership
                               est_sub,   # estimated subgroup membership
                               pendiag = TRUE){ # penalize AR coefs
  
  summary_stats <- lapply(seq_along(true_mats), function(j){
    
    # only look at lagged dynamics
    if(mods == "sgimme"){
      est_mats[[j]] <- est_mats[[j]][,1:d]
    }
    
    # don't count AR coefs in model recovery metrics if not penalized 
    if(!pendiag){
      
      data1 <- c(est_mats[[j]][lower.tri(est_mats[[j]])]) != 0
      data2 <- c(est_mats[[j]][upper.tri(est_mats[[j]])]) != 0
      data <- as.factor(c(data1, data2))
      
      reference1 <- c(true_mats[[j]][lower.tri(true_mats[[j]])]) != 0
      reference2 <- c(true_mats[[j]][upper.tri(true_mats[[j]])]) != 0
      reference  <- as.factor(c(reference1, reference2))
      
    } else {
      
      data <- as.factor(c(est_mats[[j]]) != 0)
      reference <- as.factor(c(true_mats[[j]]) != 0)
    }
    
    
    if(length(levels(data)) == 1){
      fail <- 1
    } else{
      fail <- 0
    }
    
    confuseMat <- caret::confusionMatrix(
      data = data, 
      reference = reference,
      positive = "TRUE"
    )
    
    sens <- confuseMat$byClass["Sensitivity"]
    spec <- confuseMat$byClass["Specificity"]
    
    mcc <- mltools::mcc(actuals = reference, preds = data)
    
    biasrel <- SimDesign::bias(
      estimate = c(est_mats[[j]])[c(true_mats[[j]])!=0],
      parameter = c(true_mats[[j]])[c(true_mats[[j]])!=0],
      type = "relative",
      abs = TRUE,
      percent = TRUE,
      unname = FALSE
    )
    
    biasabs <- SimDesign::bias(
      estimate = c(est_mats[[j]]),
      parameter = c(true_mats[[j]]),
      type = "bias",
      abs = TRUE,
      percent = FALSE,
      unname = FALSE
    )
    
    bias <- SimDesign::bias(
      estimate = c(est_mats[[j]]),
      parameter = c(true_mats[[j]]),
      type = "bias",
      abs = FALSE,
      percent = FALSE,
      unname = FALSE
    )
    
    biasabs_zero <- SimDesign::bias(
      estimate = c(est_mats[[j]])[c(true_mats[[j]])==0],
      parameter = c(true_mats[[j]])[c(true_mats[[j]])==0],
      type = "bias",
      abs = TRUE,
      percent = FALSE,
      unname = FALSE
    )
    
    biasabs_nonzero <- SimDesign::bias(
      estimate = c(est_mats[[j]])[c(true_mats[[j]])!=0],
      parameter = c(true_mats[[j]])[c(true_mats[[j]])!=0],
      type = "bias",
      abs = TRUE,
      percent = FALSE,
      unname = FALSE
    )
    
    bias_zero <- SimDesign::bias(
      estimate = c(est_mats[[j]])[c(true_mats[[j]])==0],
      parameter = c(true_mats[[j]])[c(true_mats[[j]])==0],
      type = "bias",
      abs = FALSE,
      percent = FALSE,
      unname = FALSE
    )
    
    bias_nonzero <- SimDesign::bias(
      estimate = c(est_mats[[j]])[c(true_mats[[j]])!=0],
      parameter = c(true_mats[[j]])[c(true_mats[[j]])!=0],
      type = "bias",
      abs = FALSE,
      percent = FALSE,
      unname = FALSE
    )
    
    rmse <- SimDesign::RMSE(
      estimate = c(est_mats[[j]]),
      parameter = c(true_mats[[j]]),
      type = "RMSE",
      MSE = FALSE,
      percent = FALSE,
      unname = FALSE
    )
    
    ari <- adjustedRandIndex(est_sub, true_sub)
    
    sub_num <- length(unique(est_sub))
    
    c(  "dataset" = dataset,
        "n" = n,
        "d" = d,
        "k" = k,
        "s" = s,
        "mods" = mods,
        "mcc" = mcc,
         confuseMat$overall,
         confuseMat$byClass, 
        "bias" = bias, 
        "biasabs" = biasabs, 
        "biasabs_zero" = biasabs_zero, 
        "biasabs_nonzero" = biasabs_nonzero, 
        "bias_zero" = bias_zero, 
        "bias_nonzero" = bias_nonzero, 
        "biasrel"= biasrel,
        "rmse" = rmse,
        "label" = alabel,
        "ari" = ari,
        "sub_num" = sub_num,
        "fail" = fail
        )
    
  })
  
  return(summary_stats)
  
}



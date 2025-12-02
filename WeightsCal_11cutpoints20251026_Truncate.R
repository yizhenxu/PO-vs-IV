datapath = "/uufs/chpc.utah.edu/common/home/u6050801/Weighted_Survival/"

for(cc in c(0.01, 0.05, 0.1)){
  
  if(cc==0.01) cname = "p01"
  if(cc==0.05) cname = "p05"
  if(cc==0.1) cname = "p1"
  
  setwd("/uufs/chpc.utah.edu/common/home/u6050801/Weighted_Survival/")
  
  ### mITT
  #main      <- readRDS("~/Library/CloudStorage/Box-Box/Work/metrc/poiv/Test2/Progs/Table1/Weighted_New_Assigned/RDS/main.rds")
  
  main <- readRDS("./AData.rds") #all(mainITT==AData,na.rm=T) TRUE
  
  intvars   <- c("insStatus", "fixation_retained", "n_debridements", 
                 "age", "vr12_q1_exvg", "bmi", "fixationToDebridementb", "index_hosp_durationb", 
                 "gram_positive","gram_neg","renal_disease_indicator", "subabuse_indicator", "diabetes","Male","fracture")
  
  n=dim(main)[1]
  main$fractureC[is.na(main$fractureC)]=0
  main$fracture = main$fractureC
  #main$fracture[is.na(main$fractureC)]=0
  #main$fracture[main$fractureB==0 & main$fractureC==0]=0
  #main$college[is.na(main$college)]=0
  main$fixationToDebridementb = as.integer((main$fixationToDebridement<=median(main$fixationToDebridement)))
  main$index_hosp_durationb = as.integer((main$index_hosp_duration<=median(main$index_hosp_duration)))
  
  main$tx_assigned = as.integer(main$treatAssignEng == "PO")
  table(main$tx_assigned)
  
  lgform    <- paste( "tx_assigned ~", paste( intvars, sep=" ", collapse=" + ") )
  
  library(rstanarm)
  lg <- stan_glm(
    as.formula(lgform),
    data = main,
    family = binomial(link = "logit"),
    prior_intercept = normal(0, 10),
    QR = TRUE,
    refresh = 0,
    # for speed of example only
    chains = 2, iter = 10000, seed = 123
  )
  
  f <- t(posterior_epred(lg)) # 233 x 2000
  f[ main[,"tx_assigned"] == 0 ,] <- 1 - f[ main[,"tx_assigned"] == 0 ,]
  cat( paste( "----- min f =",min(f),"   max f =",max(f),"   min wt =", min(1/f), "   max wt =", max(1/f),"-----\n"))
  #plot(lg, plotfun = "trace")
  f   <-  pmax(f,cc)
  f   <-  pmin(f,1-cc)
  wt_mat<- 1/f
  round(range(wt_mat),2 )# 1.002915 147.445379
  
  wt_mat = cbind(main$study_id, wt_mat)
  
  saveRDS(wt_mat, file = paste0(datapath, "weights/mITT_Weights20251026_11_",cname,".Rds"))
  #######################################################################
  #######################################################################
  #######################################################################
  #rm(list = ls())
  gc()
  
  ### PP (As treated)
  main <- readRDS("./AData.rds") 
  
  intvars   <- c("insStatus", "fixation_retained", "n_debridements", 
                 "age", "vr12_q1_exvg", "bmi", "fixationToDebridementb", "index_hosp_durationb", 
                 "gram_positive","gram_neg","renal_disease_indicator", "subabuse_indicator", "diabetes","Male","fracture")
  
  main$tx = as.integer(main$ppgroup5 %in% c("PO/PO Only", "IV/PO Only"))
  table(main$tx)
  
  n=dim(main)[1]
  main$fractureC[is.na(main$fractureC)]=0
  main$fracture = main$fractureC
  #main$fracture = rep(1,n)
  #main$fracture[is.na(main$fractureB)]=0
  #main$fracture[main$fractureB==0 & main$fractureC==0]=0
  main$fixationToDebridementb = as.integer((main$fixationToDebridement<=median(main$fixationToDebridement)))
  main$index_hosp_durationb = as.integer((main$index_hosp_duration<=median(main$index_hosp_duration)))
  
  lgform    <- paste( "tx ~", paste( intvars, sep=" ", collapse=" + ") )
  
  library(rstanarm)
  lg <- stan_glm(
    as.formula(lgform),
    data = main,
    family = binomial(link = "logit"),
    prior_intercept = normal(0, 10),
    QR = TRUE,
    refresh = 0,
    # for speed of example only
    chains = 2, iter = 10000, seed = 321
  )
  
  f <- t(posterior_epred(lg)) # 205 x 2000
  f[ main[,"tx"] == 0 ,] <- 1 - f[ main[,"tx"] == 0 ,]
  cat( paste( "----- min f =",min(f),"   max f =",max(f),"   min wt =", min(1/f), "   max wt =", max(1/f),"-----\n"))
  #plot(lg, plotfun = "trace")
  f   <-  pmax(f,cc)
  f   <-  pmin(f,1-cc)
  wt_mat<- 1/f
  round(range(wt_mat),2 )# 1.00 392.68 for PP main, 1.0 234.4 for ITT main
  
  wt_mat = cbind(main$study_id, wt_mat)
  
  saveRDS(wt_mat, file = paste0(datapath, "weights/PP_Weights20251026_11_",cname,".Rds"))
  
  
  #######################################################################
  #######################################################################
  #######################################################################
  #rm(list = ls())
  gc()
  
  ### PP no crossover
  main <- readRDS("./AData.rds") 
  
  intvars   <- c("insStatus", "fixation_retained", "n_debridements", 
                 "age", "vr12_q1_exvg", "bmi", "fixationToDebridementb", "index_hosp_durationb", 
                 "gram_positive","gram_neg","renal_disease_indicator", "subabuse_indicator", "diabetes","Male","fracture")
  
  main$tx = as.integer(main$ppgroup5 %in% c("PO/PO Only", "IV/PO Only"))
  table(main$tx)
  
  main$tx_assigned = as.integer(main$treatAssignEng == "PO")
  table(main$tx_assigned)
  
  table(main$tx_assigned, main$tx)
  #table(df$A, main$treatAssignEng)
  #table(df$A, main$treatAssign)
  
  main = main[main$tx_assigned == main$tx,] # 233 -> 214
  table(main$tx)
  
  n=dim(main)[1]
  main$fractureC[is.na(main$fractureC)]=0
  main$fracture = main$fractureC
  #main$fracture = rep(1,n)
  #main$fracture[is.na(main$fractureB)]=0
  #main$fracture[main$fractureB==0 & main$fractureC==0]=0
  main$fixationToDebridementb = as.integer((main$fixationToDebridement<=median(main$fixationToDebridement)))
  main$index_hosp_durationb = as.integer((main$index_hosp_duration<=median(main$index_hosp_duration)))
  
  lgform    <- paste( "tx ~", paste( intvars, sep=" ", collapse=" + ") )
  
  library(rstanarm)
  lg <- stan_glm(
    as.formula(lgform),
    data = main,
    family = binomial(link = "logit"),
    prior_intercept = normal(0, 10),
    QR = TRUE,
    refresh = 0,
    # for speed of example only
    chains = 2, iter = 10000, seed = 321
  )
  
  f <- t(posterior_epred(lg)) # 205 x 2000
  f[ main[,"tx"] == 0 ,] <- 1 - f[ main[,"tx"] == 0 ,]
  cat( paste( "----- min f =",min(f),"   max f =",max(f),"   min wt =", min(1/f), "   max wt =", max(1/f),"-----\n"))
  #plot(lg, plotfun = "trace")
  f   <-  pmax(f,cc)
  f   <-  pmin(f,1-cc)
  wt_mat<- 1/f
  round(range(wt_mat),2 )# 1.00 392.68 for PP main, 1.0 234.4 for ITT main
  
  wt_mat = cbind(main$study_id, wt_mat)
  
  saveRDS(wt_mat, file = paste0(datapath, "weights/PP_nocs_Weights20251026_11_",cname,".Rds"))
  
}






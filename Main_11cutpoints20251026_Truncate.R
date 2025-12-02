


oncluster = 1
#remove "cutpoints"
#update "SamplePosteriorDepGamma.R" such that
#things are aligh with the single influential feature and he want sot bve outside the screen

rerun = 0
posterior_idx = seq(2000, 10000, length=1001)

time.max = 365; K = 11
cutpoints = seq(0, time.max, length.out = (K+1))[-1]

library(survival)
#install.packages("BayesSurvival")
library(BayesSurvival)
library(ggplot2)
library(dplyr)
library(coda)

if(oncluster){
  #datapath = "/users/yxu2/delete/"
  datapath = "/uufs/chpc.utah.edu/common/home/u6050801/Weighted_Survival/"
} else {
  datapath = "~/Library/CloudStorage/GoogleDrive-xuyizhen00@gmail.com/My Drive/My Drive/Desktop/2023 fall/dan/nonparametric survival/For Paper/K11/"
}



ReshapeData1 <- function(df, time = "time", event = "event",
                         cutpoints){
  #Reshape the data
  #Input: dataframe with columns 'time' and 'event'
  #Output: number of failures and total exposure time per interval
  
  K = length(cutpoints)
  
  df.small <- data.frame(time = df[[time]], event = df[[event]])
  
  df.long <- survSplit( Surv(time, event) ~ 1, data = df.small,
                        cut = c(0, cutpoints),
                        id = "id", episode = "interval")
  
  df.long$interval <- (df.long$interval - 1) #interval count starts at 2 otherwise
  
  df.long$delta <- (df.long$time - df.long$tstart) #time spent per interval
  
  #We will continue with the following data:
  failures <- aggregate(event ~ interval, df.long, sum)$event
  exposures <- aggregate(delta ~ interval, df.long, sum)$delta
  
  #If nobody is under followup anymore during the final interval(s),
  #they will be discarded in the aggregate() - step. We augment the
  #data with zeroes for those intervals.
  indices.no.followup <- which( 1:K > max(df.long$interval) )
  
  failures[indices.no.followup] <- 0
  exposures[indices.no.followup] <- 0
  
  return(list(failures = failures, exposures = exposures))
}

ReshapeData_weighted <- function(df, weight = paste0("wt",1:npost), time = "time", event = "event",
                                 cutpoints){
  #Reshape the data
  #Input: dataframe with columns 'time' and 'event'
  #Output: number of failures and total exposure time per interval
  
  K = length(cutpoints)
  
  npost = length(weight)
  wt_mat =  df[, weight, drop = FALSE]
  colnames(wt_mat) = paste0("weight",1:npost)
  
  df.small <- data.frame(time = df[[time]], event = df[[event]])
  df.small = cbind(df.small, wt_mat)
  
  df.long <- survSplit( Surv(time, event) ~ ., data = df.small,
                        cut = c(0, cutpoints),
                        id = "id", episode = "interval")
  
  df.long$interval <- (df.long$interval - 1) #interval count starts at 2 otherwise
  
  df.long$delta <- (df.long$time - df.long$tstart) #time spent per interval
  
  #We will continue with the following data:
  library(data.table); setDT(df.long)
  tmp = lapply(1:npost, function(k) df.long[, sum(event*get(paste0("weight",k))), by=interval]$V1)
  failures = do.call(rbind, tmp) # npost x K
  rm(tmp)
  tmp = lapply(1:npost, function(k) df.long[, sum(delta*get(paste0("weight",k))), by=interval]$V1)
  exposures = do.call(rbind, tmp) # npost x K
  
  #If nobody is under followup anymore during the final interval(s),
  #they will be discarded in the aggregate() - step. We augment the
  #data with zeroes for those intervals.
  indices.no.followup <- which( 1:K > max(df.long$interval) )
  
  failures[,indices.no.followup] <- 0
  exposures[,indices.no.followup] <- 0
  
  return(list(failures = failures, exposures = exposures))
}


SamplePosteriorDepGamma_weighted <- function(failures, exposures, alpha.dep = 1, alpha0.dep = 1.5, beta0.dep = 1){
  
  K <- ncol(failures)
  N <- nrow(failures)
  
  epsilon <- 0.001 #to prevent an initial 'current' estimate of zero
  samples <- matrix(0, nrow = N, ncol = K)
  
  samples[1, ] <- failures[1,]/(exposures[1,] + 1) + epsilon #initialization
  
  if(K == 1){
    for(i in 1:N){
      samples[i, 1] <- rgamma(N, shape = (failures[i,1] + alpha0.dep), rate = (exposures[i,1] + beta0.dep))
    }
    #samples[ , 1] <- rgamma(N, shape = (failures[1] + alpha0.dep), rate = (exposures[1] + beta0.dep))
  } #end if(K == 1)
  
  
  if(K == 2){
    for(i in 2:N){
      samples[i, 1] <- MCMCDepGammaFirst(samples[(i-1), 1], samples[(i-1), 2], failures[i,1], exposures[i,1], alpha.dep, alpha0.dep, beta0.dep)
      
      samples[i, K] <- rgamma(1, failures[i,K] + alpha.dep, alpha.dep/samples[i, (K-1)] + exposures[i,K]) #final interval, direct sampling possible
    }
  } #end if(K == 2)
  
  if(K > 2){
    for(i in 2:N){
      samples[i, 1] <- MCMCDepGammaFirst(samples[(i-1), 1], samples[(i-1), 2], failures[i,1], exposures[i,1], alpha.dep, alpha0.dep, beta0.dep)
      
      for(k in 2:(K-1)){
        samples[i, k] <- MCMCDepGammaIntermediate(samples[(i-1), k], samples[i, (k-1)], samples[(i-1), (k+1)], failures[i,k], exposures[i,k], alpha.dep)
      }
      
      samples[i, K] <- rgamma(1, failures[i,K] + alpha.dep, alpha.dep/samples[i, (K-1)] + exposures[i,K]) #final interval, direct sampling possible
    }
  } #end if(K > 2)
  
  return(samples)
  
}


SamplePosteriorIndepGamma_weighted <- function(failures, exposures, alpha.indep = 1.5, beta.indep = 1){
  
  K <- ncol(failures)
  N <- nrow(failures)
  
  samples <- matrix(0, nrow = N, ncol = K)
  
  for(k in 1:K){
    for(i in 1:N){
      samples[i, k] <- rgamma(N, shape = (failures[i,k] + alpha.indep), rate = (exposures[i,k] + beta.indep))
    }
  }
  
  return(samples)
  
}


summ = function(x){
  rn = 1
  sx = sort(x)
  n = length(x)
  rr = paste0(round(mean(sx)*100, rn), " (", round(sx[round(n*0.025)]*100, rn), ', ',  round(sx[round(n*0.975)]*100, rn), ")")
  return(rr)
}


summ_med = function(x){
  rn = 2
  sx = sort(x)
  n = length(x)
  rr = paste0(round(median(sx), rn), " (", round(sx[round(n*0.025)], rn), ', ',  round(sx[round(n*0.975)], rn), ")")
  return(rr)
}

library(data.table)
library(ggplot2)

# --- helpers (same as before) ---
to_step <- function(x, y) {
  xs <- rep(x, each = 2)
  ys <- rep(y, each = 2)
  data.table(time = xs[-1], value = ys[-length(ys)])
}

summarize_surv_steps <- function(cut_vec, H_mat, label) {
  times <- c(0, cut_vec, max(cut_vec)+1)         # add t=0
  H_ext <- cbind(0, H_mat)       # H(0)=0
  S_draws <- exp(-H_ext)         # survival draws
  
  S_mean  <- colMeans(S_draws)
  #S_lower <- apply(S_draws, 2, quantile, probs = 0.025)
  #S_upper <- apply(S_draws, 2, quantile, probs = 0.975)
  
  eval.vec = seq(0, (max(cut_vec) - 0.0001), length = 10*K)
  surv.pm <- approxfun(times, S_mean)
  surv.pm.at.grid <- surv.pm(eval.vec)
  #surv.pm.at.grid = S_mean
  surv.radius <- RadiusCredibleSet(S_draws, S_mean, alpha = 0.05)
  
  S_lower = pmax(0, surv.pm.at.grid - surv.radius)
  S_upper = pmin(1, surv.pm.at.grid + surv.radius)
  
  step_mean  <- to_step(eval.vec, surv.pm.at.grid)[, group := label]
  step_lower <- to_step(eval.vec, S_lower)
  step_upper <- to_step(eval.vec, S_upper)
  
  step_band <- data.table(time = step_mean$time,
                          lower = step_lower$value,
                          upper = step_upper$value,
                          group = label)
  
  list(mean = step_mean, band = step_band)
}

# --- main: survival plot ---
plot_survival <- function(pcuts, pH, labels = c("IV", "Oral"),
                          colors = c("red", "blue"),
                          max_time = 365, title_here) {
  # Summaries per group
  surv_summaries <- lapply(seq_along(pcuts), function(i) {
    summarize_surv_steps(pcuts[[i]], pH[[i]], labels[i])
  })
  
  step_mean <- rbindlist(lapply(surv_summaries, `[[`, "mean"))
  step_band <- rbindlist(lapply(surv_summaries, `[[`, "band"))
  
  step_mean[, group := factor(group, levels = labels)]
  step_band[, group := factor(group, levels = labels)]
  
  ggplot() +
    geom_ribbon(data = step_band,
                aes(x = time, ymin = lower, ymax = upper, fill = group),
                alpha = 0.25) +
    geom_line(data = step_mean,
              aes(x = time, y = value, color = group),
              linewidth = 1) +
    scale_color_manual(values = setNames(colors, labels)) +
    scale_fill_manual(values  = setNames(colors, labels)) +
    coord_cartesian(xlim = c(0, max_time), ylim = c(0, 1)) +
    labs(x = "Time", y = "Survival Probability",
         color = NULL, fill = NULL,
         title = title_here) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
}

# --- new: CDF plot (F(t) = 1 - S(t)) ---
summarize_cdf_steps <- function(cut_vec, H_mat, label) {
  times <- c(0, cut_vec, max(cut_vec)+1)
  H_ext <- cbind(0, H_mat)
  S_draws <- exp(-H_ext)
  F_draws <- 1 - S_draws
  
  F_mean  <- colMeans(F_draws)
  
  eval.vec = seq(0, (max(cut_vec) - 0.0001), length = 10*K)
  surv.pm <- approxfun(times, F_mean)
  surv.pm.at.grid <- surv.pm(eval.vec)
  #surv.pm.at.grid = S_mean
  surv.radius <- RadiusCredibleSet(F_draws, F_mean, alpha = 0.05)
  
  F_lower = pmax(0, surv.pm.at.grid - surv.radius)
  F_upper = pmin(1, surv.pm.at.grid + surv.radius)
  
  step_mean  <- to_step(eval.vec, surv.pm.at.grid)[, group := label]
  step_lower <- to_step(eval.vec, F_lower)
  step_upper <- to_step(eval.vec, F_upper)
  
  step_band <- data.table(time = step_mean$time,
                          lower = step_lower$value,
                          upper = step_upper$value,
                          group = label)
  
  list(mean = step_mean, band = step_band)
}

plot_cdf <- function(pcuts, pH, labels = c("IV", "Oral"),
                     colors = c("red", "blue"),
                     max_time = 365, title_here) {
  cdf_summaries <- lapply(seq_along(pcuts), function(i) {
    summarize_cdf_steps(pcuts[[i]], pH[[i]], labels[i])
  })
  
  step_mean <- rbindlist(lapply(cdf_summaries, `[[`, "mean"))
  step_band <- rbindlist(lapply(cdf_summaries, `[[`, "band"))
  
  step_mean[, group := factor(group, levels = labels)]
  step_band[, group := factor(group, levels = labels)]
  
  ggplot() +
    geom_ribbon(data = step_band,
                aes(x = time, ymin = lower, ymax = upper, fill = group),
                alpha = 0.25) +
    geom_line(data = step_mean,
              aes(x = time, y = value, color = group),
              linewidth = 1) +
    scale_color_manual(values = setNames(colors, labels)) +
    scale_fill_manual(values  = setNames(colors, labels)) +
    coord_cartesian(xlim = c(0, max_time), ylim = c(0, 1)) +
    labs(x = "Time", y = "Cumulative Incidence",
         color = NULL, fill = NULL,
         title = title_here) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
}

truncate_upper_col <- function(M, prob = 0.99) {
  stopifnot(is.matrix(M))
  p <- apply(M, 2, quantile, probs = prob, na.rm = TRUE, type = 7)
  pmin(M, rep(p, each = nrow(M)))
}
# wt_mat_trunc <- truncate_upper_col(wt_mat, 0.99)

####################################################################################################

for(cc in c(0.01, 0.05, 0.1)){
  gc()
  if(cc==0.01) cname = "p01"
  if(cc==0.05) cname = "p05"
  if(cc==0.1) cname = "p1"
  
  SurvPlots = cdfPlots = vector("list", 6)
  
  ### ITT
  
  AData <- readRDS(paste0(datapath, "AData.rds"))
  
  library(dplyr)
  df_raw <- readRDS(paste0(datapath, "Yizhen_Data.Rds"))
  df_raw = df_raw %>%
    left_join(AData %>% select(study_id, maxFUTimeAllVisits), by = "study_id")
  #all(df_raw$daysToReinfection<= df_raw$maxFUTimeAllVisits, na.rm=T)# TRUE, all events <= censoring
  #df_raw$max_futime_all_visits[df_raw$max_futime_all_visits == 0] = 1
  #all(df_raw$max_futime_all_visits == pmin(df_raw$maxFUTimeAllVisits, 365, na.rm = TRUE)) # TRUE, created censoring = min(true censoring, 365)
  #range(df_raw$maxFUTimeAllVisits) # 1, 1631
  df <- data.frame(id = as.numeric(as.factor(df_raw$study_id)), study_id = df_raw$study_id)
  df$eventtime <- pmin(df_raw$maxFUTimeAllVisits, df_raw$daysToReinfection, na.rm = TRUE)
  df$eventtime[is.infinite(df$eventtime)] <- NA  
  df$status <- 1*(df_raw$adjReinfectionResult == "Y")
  df$A <- 1*(df_raw$treatment_group=="Oral")
  
  if(oncluster){
    wt_mat = readRDS(paste0(datapath, "weights/mITT_Weights20251026_11_",cname,".Rds"))
    ##wt_mat <- truncate_upper_col(wt_mat, 0.99)
    #wt_mat <- readRDS("/users/yxu2/delete/mITT_Weights20251026_11_",cname,".Rds")
  } else {
    wt_mat <- readRDS(paste0(datapath,"weights/mITT_Weights20251026_11_",cname,".Rds"))
  }
  
  
  
  npost <- ncol(wt_mat) - 1
  colnames(wt_mat) = c("study_id", paste0("wt", 1:npost) )
  wt_mat = as.data.frame(wt_mat)
  
  #df = cbind(df, wt_mat)
  library(dplyr)
  df = df %>% left_join(wt_mat, by = "study_id")
  
  N = npost;
  alpha = 0.05; alpha.dep = 1; alpha0.dep = 1.5; beta0.dep = 1; alpha.indep = 1.5; beta.indep = 1;
  surv.factor = 10; surv.epsilon = 0.0000000001
  
  time = "eventtime"
  event = "status"
  weight = paste0("wt",1:npost)
  
  tabshow = matrix(NA, ncol = 3, nrow = 2)
  colnames(tabshow) = c("Oral", "IV", "Oral - IV")
  rownames(tabshow) = c("Unadjusted", "Adjusted")
  
  tabgeweke = matrix(NA, ncol = 3, nrow = 2)
  colnames(tabgeweke) = c("Oral", "IV", "Oral - IV")
  rownames(tabgeweke) = c("Unadjusted", "Adjusted")
  
  df_stored = df
  
  ### unadjusted
  set.seed(2025)
  post.samples.unadj = vector("list",2)
  pcuts = pH = vector("list",2)
  
  a = 0 # IV
  dff = df[df$A == a,]
  #cutpoints =  sort(unique(dff$eventtime))#sort(unique(c(df_raw$daysToReinfection, df_raw$maxFUTimeAllVisits)))
  #interval.length = diff(c(0,cutpoints))
  data.summary <- ReshapeData1(dff,  time, event, cutpoints)
  fail <- data.summary$failures
  exp <- data.summary$exposures
  post.samples.unadj[[a+1]] <- SamplePosteriorDepGamma(fail, exp, 10000, alpha.dep, alpha0.dep, beta0.dep)
  #cumhaz <- t(apply(post.samples.unadj[[a+1]]%*%diag(interval.length), 1, cumsum)) #divide by K to account for interval length
  cumhaz <- t(apply(post.samples.unadj[[a+1]] * time.max/K, 1, cumsum)) #divide by K to account for interval length
  tmp0 = 1 - exp(-cumhaz[,which(cutpoints == 365)])
  tabshow[1,2] = summ(tmp0[posterior_idx])
  tabgeweke[1,2] = coda::geweke.diag(mcmc(tmp0[posterior_idx]))$z
  pcuts[[a+1]] = cutpoints; pH[[a+1]] = cumhaz[posterior_idx,]
  
  rm(cumhaz)
  
  a = 1 # Oral
  dff = df[df$A == a,]
  #cutpoints =  sort(unique(dff$eventtime))#sort(unique(c(df_raw$daysToReinfection, df_raw$maxFUTimeAllVisits)))
  #interval.length = diff(c(0,cutpoints))
  data.summary <- ReshapeData1(dff,  time, event, cutpoints)
  fail <- data.summary$failures
  exp <- data.summary$exposures
  post.samples.unadj[[a+1]] <- SamplePosteriorDepGamma(fail, exp, 10000, alpha.dep, alpha0.dep, beta0.dep)
  cumhaz <- t(apply(post.samples.unadj[[a+1]] * time.max/K, 1, cumsum)) #divide by K to account for interval length
  tmp1 = 1 - exp(-cumhaz[,which(cutpoints == 365)])
  tabshow[1,1] = summ(tmp1[posterior_idx])
  tabgeweke[1,1] = geweke.diag(mcmc(tmp1[posterior_idx]))$z
  pcuts[[a+1]] = cutpoints; pH[[a+1]] = cumhaz[posterior_idx,]
  
  rm(cumhaz)
  
  tabshow[1,3] = summ(tmp1[posterior_idx] - tmp0[posterior_idx]) # Oral - IV
  tabgeweke[1,3] = geweke.diag(mcmc(tmp1[posterior_idx] - tmp0[posterior_idx]))$z
  
  mean(tmp1[posterior_idx] < tmp0[posterior_idx]) # 52.6%
  
  ### unadjusted plot
  # survival/cdf plots
  SurvPlots[[1]] = plot_survival(pcuts, pH, labels = c("IV","Oral"),
                                 colors = c("red","blue"), max_time = 365, title_here = "mITT Unadjusted")
  cdfPlots[[1]] = plot_cdf(pcuts, pH, labels = c("IV","Oral"),
                           colors = c("red","blue"), max_time = 365, title_here = "mITT Unadjusted")
  
  
  # Combine into a single long data frame with custom labels
  plot_df <- data.frame(
    value = c(tmp0, tmp1),
    group = factor(rep(c("IV", "Oral"), each = length(tmp0)))
  )
  
  # Plot overlaying densities
  p1 = ggplot(plot_df, aes(x = value, color = group, fill = group)) +
    geom_density(aes(y = ..scaled..), alpha = 0.3) +
    labs(x = "Value", y = "Density", title = "mITT Unadjusted") +
    theme_minimal()+
    coord_cartesian(xlim = c(0.2, 0.6), ylim = c(0, 1))
  
  plot_df <- data.frame(difftmp = tmp1 - tmp0)
  
  p7 = ggplot(plot_df, aes(x = difftmp)) +
    geom_density( aes(y = ..scaled..), color = "black") +
    labs(x = "Difference", y = "Density", title = "mITT Unadjusted") +
    theme_minimal()+
    coord_cartesian(xlim = c(-0.2, 0.3))
  
  
  ### adjusted
  set.seed(2025)
  if(rerun){
    load(paste0(datapath, "mod_post_mITT_Truncate_", cname,".RData"))
    #load(paste0(datapath, "mod_post.RData"))
  } else {
    post.samples = vector("list",2)
  }
  
  pcuts = pH = vector("list",2)
  
  a = 0
  dff = df[df$A == a,]
  #cutpoints =  sort(unique(dff$eventtime))#sort(unique(c(df_raw$daysToReinfection, df_raw$maxFUTimeAllVisits)))
  #interval.length = diff(c(0,cutpoints))
  if(!rerun){
    data.summary <- ReshapeData_weighted(dff, weight, time, event, cutpoints)
    fail <- data.summary$failures
    exp <- data.summary$exposures
    post.samples[[a+1]] <- SamplePosteriorDepGamma_weighted(fail, exp, alpha.dep, alpha0.dep, beta0.dep)
  }
  cumhaz <- t(apply(post.samples[[a+1]] * time.max/K, 1, cumsum)) #divide by K to account for interval length
  tmp0 = 1 - exp(-cumhaz[,which(cutpoints == 365)])
  tabshow[2,2] = summ(tmp0[posterior_idx])
  tabgeweke[2,2] = geweke.diag(mcmc(tmp0[posterior_idx]))$z
  pcuts[[a+1]] = cutpoints; pH[[a+1]] = cumhaz[posterior_idx,]
  
  rm(cumhaz)
  
  a = 1
  dff = df[df$A == a,]
  #cutpoints =  sort(unique(dff$eventtime))#sort(unique(c(df_raw$daysToReinfection, df_raw$maxFUTimeAllVisits)))
  #interval.length = diff(c(0,cutpoints))
  if(!rerun){
    data.summary <- ReshapeData_weighted(dff, weight, time, event, cutpoints)
    fail <- data.summary$failures
    exp <- data.summary$exposures
    post.samples[[a+1]] <- SamplePosteriorDepGamma_weighted(fail, exp, alpha.dep, alpha0.dep, beta0.dep)
  }
  cumhaz <- t(apply(post.samples[[a+1]] * time.max/K, 1, cumsum)) #divide by K to account for interval length
  tmp1 = 1 - exp(-cumhaz[,which(cutpoints == 365)])
  tabshow[2,1] = summ(tmp1[posterior_idx])
  tabgeweke[2,1] = geweke.diag(mcmc(tmp1[posterior_idx]))$z
  pcuts[[a+1]] = cutpoints; pH[[a+1]] = cumhaz[posterior_idx,]
  
  rm(cumhaz)
  
  tabshow[2,3] = summ(tmp1[posterior_idx] - tmp0[posterior_idx])
  tabgeweke[2,3] = geweke.diag(mcmc(tmp1[posterior_idx] - tmp0[posterior_idx]))$z
  
  if(!rerun){
    save(post.samples, file = paste0(datapath, "mod_post_mITT_Truncate_", cname,".RData"))
  }
  
  ### adjusted plot
  # survival/cdf plots
  SurvPlots[[2]] = plot_survival(pcuts, pH, labels = c("IV","Oral"),
                                 colors = c("red","blue"), max_time = 365, title_here = "mITT Adjusted")
  cdfPlots[[2]] = plot_cdf(pcuts, pH, labels = c("IV","Oral"),
                           colors = c("red","blue"), max_time = 365, title_here = "mITT Adjusted")
  
  # Combine into a single long data frame with custom labels
  plot_df <- data.frame(
    value = c(tmp0, tmp1),
    group = factor(rep(c("IV", "Oral"), each = length(tmp0)))
  )
  
  # Plot overlaying densities
  p2 = ggplot(plot_df, aes(x = value, color = group, fill = group)) +
    geom_density(aes(y = ..scaled..), alpha = 0.3) +
    labs(x = "Value", y = "Density", title = "mITT Adjusted") +
    theme_minimal()+
    coord_cartesian(xlim = c(0.2, 0.6), ylim = c(0, 1))
  
  plot_df <- data.frame(difftmp = tmp1 - tmp0)
  
  p8 = ggplot(plot_df, aes(x = difftmp)) +
    geom_density( aes(y = ..scaled..), color = "black") +
    labs(x = "Difference", y = "Density", title = "mITT Adjusted") +
    theme_minimal()+
    coord_cartesian(xlim = c(-0.2, 0.3))
  
  noquote(tabshow)
  tabmITT = copy(tabshow) # 2169
  tabgeweke_mITT = copy(tabgeweke)
  ####################################################################################################
  ####################################################################################################
  ####################################################################################################
  ### As-treated
  
  library(BayesSurvival)
  
  #all(df$study_id == df_raw$study_id) # TRUE
  #df$A <- 1*(df_raw$finalPPGroup=="PO only") # from Yizhen_Data20250424.Rds
  df = df[,1:4]
  df$A <- 1 - AData$TrtActual
  
  if(oncluster){
    wt_mat =  readRDS(paste0(datapath, "weights/PP_Weights20251026_11_",cname,".Rds"))
    #wt_mat <- readRDS("/users/yxu2/delete/PP_Weights20250426_11_",cname,".Rds")
  } else {
    wt_mat <- readRDS(paste0(datapath,"weights/PP_Weights20251026_11_",cname,".Rds"))
  }
  
  npost <- ncol(wt_mat) - 1
  colnames(wt_mat) = c("study_id", paste0("wt", 1:npost) )
  wt_mat = as.data.frame(wt_mat)
  
  library(dplyr)
  df = df %>% left_join(wt_mat, by = "study_id")
  
  N = npost;
  alpha = 0.05; alpha.dep = 1; alpha0.dep = 1.5; beta0.dep = 1; alpha.indep = 1.5; beta.indep = 1;
  surv.factor = 10; surv.epsilon = 0.0000000001
  
  time = "eventtime"
  event = "status"
  weight = paste0("wt",1:npost)
  
  tabshow = matrix(NA, ncol = 3, nrow = 2)
  colnames(tabshow) = c("Oral", "IV", "Oral - IV")
  rownames(tabshow) = c("Unadjusted", "Adjusted")
  
  tabgeweke = matrix(NA, ncol = 3, nrow = 2)
  colnames(tabgeweke) = c("Oral", "IV", "Oral - IV")
  rownames(tabgeweke) = c("Unadjusted", "Adjusted")
  
  df_stored = df
  
  ### unadjusted
  set.seed(2025)
  post.samples.unadj = vector("list",2)
  
  pcuts = pH = vector("list",2)
  
  a = 0 # IV
  dff = df[df$A == a,]
  data.summary <- ReshapeData1(dff,  time, event, cutpoints)
  fail <- data.summary$failures
  exp <- data.summary$exposures
  post.samples.unadj[[a+1]] <- SamplePosteriorDepGamma(fail, exp, 10000, alpha.dep, alpha0.dep, beta0.dep)
  cumhaz <- t(apply(post.samples.unadj[[a+1]] * time.max/K, 1, cumsum)) #divide by K to account for interval length
  tmp0 = 1 - exp(-cumhaz[,which(cutpoints == 365)])
  tabshow[1,2] = summ(tmp0[posterior_idx])
  tabgeweke[1,2] = geweke.diag(mcmc(tmp0[posterior_idx]))$z
  pcuts[[a+1]] = cutpoints; pH[[a+1]] = cumhaz[posterior_idx,]
  
  rm(cumhaz)
  
  a = 1 # Oral
  dff = df[df$A == a,]
  data.summary <- ReshapeData1(dff,  time, event, cutpoints)
  fail <- data.summary$failures
  exp <- data.summary$exposures
  post.samples.unadj[[a+1]] <- SamplePosteriorDepGamma(fail, exp, 10000, alpha.dep, alpha0.dep, beta0.dep)
  cumhaz <- t(apply(post.samples.unadj[[a+1]] * time.max/K, 1, cumsum)) #divide by K to account for interval length
  tmp1 = 1 - exp(-cumhaz[,which(cutpoints == 365)])
  tabshow[1,1] = summ(tmp1[posterior_idx])
  tabgeweke[1,1] = geweke.diag(mcmc(tmp1[posterior_idx]))$z
  pcuts[[a+1]] = cutpoints; pH[[a+1]] = cumhaz[posterior_idx,]
  
  rm(cumhaz)
  
  tabshow[1,3] = summ(tmp1[posterior_idx] - tmp0[posterior_idx]) # Oral - IV
  tabgeweke[1,3] = geweke.diag(mcmc(tmp1[posterior_idx] - tmp0[posterior_idx]))$z
  
  mean(tmp1[posterior_idx] < tmp0[posterior_idx]) # 33.7%
  
  ### unadjusted plot
  # survival/cdf plots
  SurvPlots[[3]] = plot_survival(pcuts, pH, labels = c("IV","Oral"),
                                 colors = c("red","blue"), max_time = 365, title_here = "As-treated Unadjusted")
  cdfPlots[[3]] = plot_cdf(pcuts, pH, labels = c("IV","Oral"),
                           colors = c("red","blue"), max_time = 365, title_here = "As-treated Unadjusted")
  
  # Combine into a single long data frame with custom labels
  plot_df <- data.frame(
    value = c(tmp0, tmp1),
    group = factor(rep(c("IV", "Oral"), each = length(tmp0)))
  )
  
  # Plot overlaying densities
  p3 = ggplot(plot_df, aes(x = value, color = group, fill = group)) +
    geom_density(aes(y = ..scaled..), alpha = 0.3) +
    labs(x = "Value", y = "Density", title = "As-treated Unadjusted") +
    theme_minimal()+
    coord_cartesian(xlim = c(0.2, 0.6), ylim = c(0, 1))
  
  plot_df <- data.frame(difftmp = tmp1 - tmp0)
  
  p9 = ggplot(plot_df, aes(x = difftmp)) +
    geom_density( aes(y = ..scaled..), color = "black") +
    labs(x = "Difference", y = "Density", title = "As-treated Unadjusted") +
    theme_minimal()+
    coord_cartesian(xlim = c(-0.2, 0.3))
  
  
  ### adjusted
  set.seed(2025)
  if(rerun){
    load(paste0(datapath, "mod_post_PP_Truncate_", cname,".RData"))
  } else {
    post.samples = vector("list",2)
  }
  
  pcuts = pH = vector("list",2)
  
  a = 0
  dff = df[df$A == a,]
  if(!rerun){
    data.summary <- ReshapeData_weighted(dff, weight, time, event, cutpoints)
    fail <- data.summary$failures
    exp <- data.summary$exposures
    post.samples[[a+1]] <- SamplePosteriorDepGamma_weighted(fail, exp, alpha.dep, alpha0.dep, beta0.dep)
  }
  cumhaz <- t(apply(post.samples[[a+1]] * time.max/K, 1, cumsum)) #divide by K to account for interval length
  tmp0 = 1 - exp(-cumhaz[,which(cutpoints == 365)])
  tabshow[2,2] = summ(tmp0[posterior_idx])
  tabgeweke[2,2] = geweke.diag(mcmc(tmp0[posterior_idx]))$z
  pcuts[[a+1]] = cutpoints; pH[[a+1]] = cumhaz[posterior_idx,]
  
  rm(cumhaz)
  
  a = 1 #PO
  dff = df[df$A == a,]
  if(!rerun){
    data.summary <- ReshapeData_weighted(dff, weight, time, event, cutpoints)
    fail <- data.summary$failures
    exp <- data.summary$exposures
    post.samples[[a+1]] <- SamplePosteriorDepGamma_weighted(fail, exp, alpha.dep, alpha0.dep, beta0.dep)
  }
  cumhaz <- t(apply(post.samples[[a+1]] * time.max/K, 1, cumsum)) #divide by K to account for interval length
  tmp1 = 1 - exp(-cumhaz[,which(cutpoints == 365)])
  tabshow[2,1] = summ(tmp1[posterior_idx])
  tabgeweke[2,1] = geweke.diag(mcmc(tmp1[posterior_idx]))$z
  pcuts[[a+1]] = cutpoints; pH[[a+1]] = cumhaz[posterior_idx,]
  
  rm(cumhaz)
  
  tabshow[2,3] = summ(tmp1[posterior_idx] - tmp0[posterior_idx])
  tabgeweke[2,3] = geweke.diag(mcmc(tmp1[posterior_idx] - tmp0[posterior_idx]))$z
  
  mean(tmp1[posterior_idx] < tmp0[posterior_idx]) # 5.3%
  
  if(!rerun){
    save(post.samples, file = paste0(datapath, "mod_post_PP_Truncate_", cname,".RData")) # 2161
  }
  
  ### adjusted plot
  # survival/cdf plots
  SurvPlots[[4]] = plot_survival(pcuts, pH, labels = c("IV","Oral"),
                                 colors = c("red","blue"), max_time = 365, title_here = "As-treated Adjusted")
  cdfPlots[[4]] = plot_cdf(pcuts, pH, labels = c("IV","Oral"),
                           colors = c("red","blue"), max_time = 365, title_here = "As-treated Adjusted")
  
  # Combine into a single long data frame with custom labels
  plot_df <- data.frame(
    value = c(tmp0, tmp1),
    group = factor(rep(c("IV", "Oral"), each = length(tmp0)))
  )
  
  # Plot overlaying densities
  p4 = ggplot(plot_df, aes(x = value, color = group, fill = group)) +
    geom_density(aes(y = ..scaled..), alpha = 0.3) +
    labs(x = "Value", y = "Density", title = "As-treated Adjusted") +
    theme_minimal()+
    coord_cartesian(xlim = c(0.2, 0.6), ylim = c(0, 1))
  
  plot_df <- data.frame(difftmp = tmp1 - tmp0)
  
  p10 = ggplot(plot_df, aes(x = difftmp)) +
    geom_density( aes(y = ..scaled..), color = "black") +
    labs(x = "Difference", y = "Density", title = "As-treated Adjusted") +
    theme_minimal()+
    coord_cartesian(xlim = c(-0.2, 0.3))
  
  noquote(tabshow)
  tabPP = copy(tabshow)
  tabgeweke_PP = copy(tabgeweke)
  
  ####################################################################################################
  ####################################################################################################
  ####################################################################################################
  ### PP no crossover
  
  library(BayesSurvival)
  
  #all(df$study_id == df_raw$study_id) # TRUE
  #df$A <- 1*(df_raw$finalPPGroup=="PO only") # from Yizhen_Data20250424.Rds
  df = df[,1:4]
  df$A <- 1 - AData$TrtActual
  assignedA <- 1*(df_raw$treatment_group=="Oral")
  df = df[assignedA == df$A,]
  table(df$A)
  
  if(oncluster){
    wt_mat =  readRDS(paste0(datapath, "weights/PP_nocs_Weights20251026_11_",cname,".Rds"))
    #wt_mat <- readRDS("/users/yxu2/delete/PP_nocs_Weights20251026_11_",cname,".Rds")
  } else {
    wt_mat <- readRDS(paste0(datapath,"weights/PP_nocs_Weights20251026_11_",cname,".Rds"))
  }
  
  npost <- ncol(wt_mat) - 1
  colnames(wt_mat) = c("study_id", paste0("wt", 1:npost) )
  wt_mat = as.data.frame(wt_mat)
  
  library(dplyr)
  df = df %>% left_join(wt_mat, by = "study_id")
  
  N = npost;
  alpha = 0.05; alpha.dep = 1; alpha0.dep = 1.5; beta0.dep = 1; alpha.indep = 1.5; beta.indep = 1;
  surv.factor = 10; surv.epsilon = 0.0000000001
  
  time = "eventtime"
  event = "status"
  weight = paste0("wt",1:npost)
  
  tabshow = matrix(NA, ncol = 3, nrow = 2)
  colnames(tabshow) = c("Oral", "IV", "Oral - IV")
  rownames(tabshow) = c("Unadjusted", "Adjusted")
  
  tabgeweke = matrix(NA, ncol = 3, nrow = 2)
  colnames(tabgeweke) = c("Oral", "IV", "Oral - IV")
  rownames(tabgeweke) = c("Unadjusted", "Adjusted")
  
  ### unadjusted
  set.seed(2025)
  post.samples.unadj = vector("list",2)
  
  pcuts = pH = vector("list",2)
  
  a = 0 # IV
  dff = df[df$A == a,]
  data.summary <- ReshapeData1(dff,  time, event, cutpoints)
  fail <- data.summary$failures
  exp <- data.summary$exposures
  post.samples.unadj[[a+1]] <- SamplePosteriorDepGamma(fail, exp, 10000, alpha.dep, alpha0.dep, beta0.dep)
  cumhaz <- t(apply(post.samples.unadj[[a+1]] * time.max/K, 1, cumsum)) #divide by K to account for interval length
  tmp0 = 1 - exp(-cumhaz[,which(cutpoints == 365)])
  tabshow[1,2] = summ(tmp0[posterior_idx])
  tabgeweke[1,2] = geweke.diag(mcmc(tmp0[posterior_idx]))$z
  pcuts[[a+1]] = cutpoints; pH[[a+1]] = cumhaz[posterior_idx,]
  
  rm(cumhaz)
  
  a = 1 # Oral
  dff = df[df$A == a,]
  data.summary <- ReshapeData1(dff,  time, event, cutpoints)
  fail <- data.summary$failures
  exp <- data.summary$exposures
  post.samples.unadj[[a+1]] <- SamplePosteriorDepGamma(fail, exp, 10000, alpha.dep, alpha0.dep, beta0.dep)
  cumhaz <- t(apply(post.samples.unadj[[a+1]] * time.max/K, 1, cumsum)) #divide by K to account for interval length
  tmp1 = 1 - exp(-cumhaz[,which(cutpoints == 365)])
  tabshow[1,1] = summ(tmp1[posterior_idx])
  tabgeweke[1,1] = geweke.diag(mcmc(tmp1[posterior_idx]))$z
  pcuts[[a+1]] = cutpoints; pH[[a+1]] = cumhaz[posterior_idx,]
  
  rm(cumhaz)
  
  tabshow[1,3] = summ(tmp1[posterior_idx] - tmp0[posterior_idx]) # Oral - IV
  tabgeweke[1,3] = geweke.diag(mcmc(tmp1[posterior_idx] - tmp0[posterior_idx]))$z
  
  mean(tmp1[posterior_idx] < tmp0[posterior_idx]) # 41.96%
  
  ### unadjusted plot
  # survival/cdf plots
  SurvPlots[[5]] = plot_survival(pcuts, pH, labels = c("IV","Oral"),
                                 colors = c("red","blue"), max_time = 365, title_here = "Per Protocol Unadjusted")
  cdfPlots[[5]] = plot_cdf(pcuts, pH, labels = c("IV","Oral"),
                           colors = c("red","blue"), max_time = 365, title_here = "Per Protocol Unadjusted")
  
  # Combine into a single long data frame with custom labels
  plot_df <- data.frame(
    value = c(tmp0, tmp1),
    group = factor(rep(c("IV", "Oral"), each = length(tmp0)))
  )
  
  # Plot overlaying densities
  p5 = ggplot(plot_df, aes(x = value, color = group, fill = group)) +
    geom_density(aes(y = ..scaled..), alpha = 0.3) +
    labs(x = "Value", y = "Density", title = "Per Protocol Unadjusted") +
    theme_minimal()+
    coord_cartesian(xlim = c(0.2, 0.6), ylim = c(0, 1))
  
  plot_df <- data.frame(difftmp = tmp1 - tmp0)
  
  p11 = ggplot(plot_df, aes(x = difftmp)) +
    geom_density( aes(y = ..scaled..), color = "black") +
    labs(x = "Difference", y = "Density", title = "Per Protocol Unadjusted") +
    theme_minimal()+
    coord_cartesian(xlim = c(-0.2, 0.3))
  
  
  ### adjusted
  set.seed(2025)
  if(rerun){
    load(paste0(datapath, "mod_post_PP_nocs_Truncate_", cname,".RData"))
  } else {
    post.samples = vector("list",2)
  }
  
  pcuts = pH = vector("list",2)
  
  a = 0
  dff = df[df$A == a,]
  if(!rerun){
    data.summary <- ReshapeData_weighted(dff, weight, time, event, cutpoints)
    fail <- data.summary$failures
    exp <- data.summary$exposures
    post.samples[[a+1]] <- SamplePosteriorDepGamma_weighted(fail, exp, alpha.dep, alpha0.dep, beta0.dep)
  }
  cumhaz <- t(apply(post.samples[[a+1]] * time.max/K, 1, cumsum)) #divide by K to account for interval length
  tmp0 = 1 - exp(-cumhaz[,which(cutpoints == 365)])
  tabshow[2,2] = summ(tmp0[posterior_idx])
  tabgeweke[2,2] = geweke.diag(mcmc(tmp0[posterior_idx]))$z
  pcuts[[a+1]] = cutpoints; pH[[a+1]] = cumhaz[posterior_idx,]
  
  rm(cumhaz)
  
  a = 1
  dff = df[df$A == a,]
  if(!rerun){
    data.summary <- ReshapeData_weighted(dff, weight, time, event, cutpoints)
    fail <- data.summary$failures
    exp <- data.summary$exposures
    post.samples[[a+1]] <- SamplePosteriorDepGamma_weighted(fail, exp, alpha.dep, alpha0.dep, beta0.dep)
  }
  cumhaz <- t(apply(post.samples[[a+1]] * time.max/K, 1, cumsum)) #divide by K to account for interval length
  tmp1 = 1 - exp(-cumhaz[,which(cutpoints == 365)])
  tabshow[2,1] = summ(tmp1[posterior_idx])
  tabgeweke[2,1] = geweke.diag(mcmc(tmp1[posterior_idx]))$z
  pcuts[[a+1]] = cutpoints; pH[[a+1]] = cumhaz[posterior_idx,]
  
  rm(cumhaz)
  
  tabshow[2,3] = summ(tmp1[posterior_idx] - tmp0[posterior_idx])
  tabgeweke[2,3] = geweke.diag(mcmc(tmp1[posterior_idx] - tmp0[posterior_idx]))$z
  
  if(!rerun){
    save(post.samples, file = paste0(datapath, "mod_post_PP_nocs_Truncate_", cname,".RData"))
  }
  
  ### adjusted plot
  # survival/cdf plots
  SurvPlots[[6]] = plot_survival(pcuts, pH, labels = c("IV","Oral"),
                                 colors = c("red","blue"), max_time = 365, title_here = "Per Protocol Adjusted")
  cdfPlots[[6]] = plot_cdf(pcuts, pH, labels = c("IV","Oral"),
                           colors = c("red","blue"), max_time = 365, title_here = "Per Protocol Adjusted")
  
  # Combine into a single long data frame with custom labels
  plot_df <- data.frame(
    value = c(tmp0, tmp1),
    group = factor(rep(c("IV", "Oral"), each = length(tmp0)))
  )
  
  # Plot overlaying densities
  p6 = ggplot(plot_df, aes(x = value, color = group, fill = group)) +
    geom_density(aes(y = ..scaled..), alpha = 0.3) +
    labs(x = "Value", y = "Density", title = "Per Protocol Adjusted") +
    theme_minimal()+
    coord_cartesian(xlim = c(0.2, 0.6), ylim = c(0, 1))
  
  plot_df <- data.frame(difftmp = tmp1 - tmp0)
  
  p12 = ggplot(plot_df, aes(x = difftmp)) +
    geom_density( aes(y = ..scaled..), color = "black") +
    labs(x = "Difference", y = "Density", title = "Per Protocol Adjusted") +
    theme_minimal()+
    coord_cartesian(xlim = c(-0.2, 0.3))
  
  noquote(tabshow)
  tabPPnocs = copy(tabshow)
  tabgeweke_PPnocs = copy(tabgeweke)
  
  ####################################################################################################
  ####################################################################################################
  ####################################################################################################
  ### Plot pdf
  library(ggplot2)
  library(cowplot)
  
  ## ----- Row 1: with shared legend -----
  
  # remove x labels and legends
  p1c <- p1 + theme(axis.title.x = element_blank(),
                    legend.position = "none")
  p2c <- p2 + theme(axis.title.x = element_blank(),
                    legend.position = "none")
  p3c <- p3 + theme(axis.title.x = element_blank(),
                    legend.position = "none")
  p4c <- p4 + theme(axis.title.x = element_blank(),
                    legend.position = "none")
  p5c <- p5 + theme(axis.title.x = element_blank(),
                    legend.position = "none")
  p6c <- p6 + theme(axis.title.x = element_blank(),
                    legend.position = "none")
  
  # arrange p1 and p2
  row1 <- plot_grid(p1c, p2c, p5c, p6c, p3c, p4c, ncol = 6, align = "hv")
  
  # extract shared legend
  shared_legend <- get_legend(p1 + theme(legend.position = "right"))
  
  # combine plots + legend
  row1_with_legend <- plot_grid(row1, shared_legend, ncol = 2, rel_widths = c(1, 0.2))
  
  # add shared x axis label
  row1_final <- add_sub(row1_with_legend, "Risk of Reinfection in 1 Year", vjust = -0.5)
  
  
  ## ----- Row 2: no legend, but keep blank space -----
  
  
  p7c <- p7 + theme(axis.title.x = element_blank(),
                    legend.position = "none")
  p8c <- p8 + theme(axis.title.x = element_blank(),
                    legend.position = "none")
  p9c <- p9 + theme(axis.title.x = element_blank(),
                    legend.position = "none")
  p10c <- p10 + theme(axis.title.x = element_blank(),
                      legend.position = "none")
  p11c <- p11 + theme(axis.title.x = element_blank(),
                      legend.position = "none")
  p12c <- p12 + theme(axis.title.x = element_blank(),
                      legend.position = "none")
  row2 <- plot_grid(p7c, p8c, p11c, p12c, p9c, p10c,  ncol = 6, align = "hv")
  
  # add an empty placeholder where legend would be
  blank <- ggplot() + theme_void()
  
  # combine plots + blank legend column
  row2_with_blank <- plot_grid(row2, blank, ncol = 2, rel_widths = c(1, 0.2))
  
  # add shared x axis label
  row2_final <- add_sub(row2_with_blank, "Difference between Oral and IV", vjust = -0.5)
  
  
  ## ----- Stack rows -----
  final <- plot_grid(row1_final, row2_final, ncol = 1, rel_heights = c(1, 1))
  
  #ggdraw(final)
  ggsave(paste0(datapath,"final_plot_", cname,".pdf"), final, width = 15, height = 6)
  
  
  ####################################################################################################
  ####################################################################################################
  ####################################################################################################
  ###  Table into csv
  library(dplyr)
  library(readr)
  
  df_mITT <- as.data.frame(tabmITT) %>%
    tibble::rownames_to_column("Model") %>%
    mutate(source = "mITT")
  
  df_PP <- as.data.frame(tabPP) %>%
    tibble::rownames_to_column("Model") %>%
    mutate(source = "As-treated")
  
  df_PPnocs <- as.data.frame(tabPPnocs) %>%
    tibble::rownames_to_column("Model") %>%
    mutate(source = "PP")
  
  gw_mITT <- as.data.frame(round(tabgeweke_mITT,2)) %>%
    tibble::rownames_to_column("Model") %>%
    mutate(source = "mITT")
  
  gw_PP <- as.data.frame(round(tabgeweke_PP,2)) %>%
    tibble::rownames_to_column("Model") %>%
    mutate(source = "As-treated")
  
  gw_PPnocs <- as.data.frame(round(tabgeweke_PPnocs,2)) %>%
    tibble::rownames_to_column("Model") %>%
    mutate(source = "PP")
  
  #A rule of thumb is that if |Z| > 1.96, the chain hasn't converged well (at 95% level).
  cv_mITT <- as.data.frame(ifelse(abs(tabgeweke_mITT) > 1.96, "No", "Yes")) %>%
    tibble::rownames_to_column("Model") %>%
    mutate(source = "mITT")
  cv_PP <- as.data.frame(ifelse(abs(tabgeweke_PP) > 1.96, "No", "Yes")) %>%
    tibble::rownames_to_column("Model") %>%
    mutate(source = "As-treated")
  cv_PPnocs <- as.data.frame(ifelse(abs(tabgeweke_PPnocs) > 1.96, "No", "Yes")) %>%
    tibble::rownames_to_column("Model") %>%
    mutate(source = "PP")
  
  # Combine into one table
  tab_combined <- rbind(c("Reinfection at Day 365", rep("", 4)),
                        df_mITT, df_PPnocs, df_PP,
                        c("Geweke Z Scores", rep("", 4)),
                        gw_mITT, gw_PPnocs, gw_PP,
                        c("Convergence", rep("", 4)),
                        cv_mITT, cv_PPnocs, cv_PP)
  
  
  # save to CSV
  write_csv(tab_combined, paste0(datapath, "tab_combined_", cname,".csv"))
  
  ####################################################################################################
  ####################################################################################################
  ####################################################################################################
  ### plot survival and cdf plots
  library(patchwork)
  
  # Arrange SurvPlots into 2x2 grid
  surv_grid <- (SurvPlots[[1]] | SurvPlots[[2]])/
    (SurvPlots[[5]] | SurvPlots[[6]]) /
    (SurvPlots[[3]] | SurvPlots[[4]])
  
  # Arrange cdfPlots into 2x2 grid
  cdf_grid <- (cdfPlots[[1]] | cdfPlots[[2]]) /
    (cdfPlots[[5]] | cdfPlots[[6]])/
    (cdfPlots[[3]] | cdfPlots[[4]])
  
  # Save both into a PDF
  pdf(paste0(datapath, "surv_and_cdf_band_", cname,".pdf"), width = 10, height = 9)  # adjust size as needed
  print(surv_grid)
  print(cdf_grid)
  dev.off()
}





if(0){
  p01 = read.csv(paste0(datapath, "tab_combined_", "p01",".csv"))
  p05 = read.csv(paste0(datapath, "tab_combined_", "p05",".csv"))
  p1 = read.csv(paste0(datapath, "tab_combined_", "p1",".csv"))
}



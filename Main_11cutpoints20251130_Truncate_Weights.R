library(survival)
library(BayesSurvival)
library(ggplot2)
library(dplyr)
library(data.table)
library(tableone)
library(Exact)
library(MetBrewer)
library(zoo)
library(coda)
library(patchwork)


oncluster = 1
rerun = 1
writeOutput <- F

posterior_idx = seq(2000, 10000, length=1001)

time.max = 365; K = 11 # Why aren't we just calculating K = (n/log(n))^.5 for each sample (~5 cutpoints) as built into the main function of the package?
cutpoints = seq(0, time.max, length.out = (K+1))[-1]

if(oncluster){
  #datapath = "/users/yxu2/delete/"
  datapath = "/uufs/chpc.utah.edu/common/home/u6050801/Weighted_Survival/"
} else {
  datapath = "~/Library/CloudStorage/GoogleDrive-xuyizhen00@gmail.com/My Drive/My Drive/Desktop/2023 fall/dan/nonparametric survival/For Paper/K11/"
}
setwd(datapath)
source("theme_min.r")
###############################################
### functions

# mat: matrix with
#   - column 1 = original subject ID (numeric, e.g. 1025–5001)
#   - columns 2:K = posterior draws
summarise_weights <- function(mat, analysis_label) {
  
  subj_raw  <- as.numeric(mat[, 1])          # original IDs (1025, 1092, 5001, ...)
  draws_mat_raw <- mat[, -1, drop = FALSE]       # posterior draws only
  draws_mat <- apply(draws_mat_raw, 2, function(x) x/sum(x) * 100)
  # Map original IDs to consecutive integers 1, 2, 3, ... in numeric order
  # Example: sort(unique(subj_raw)) = c(1025, 1032, 1092, 5001)
  #          -> new IDs = c(1, 2, 3, 4)
  subj_sorted <- sort(unique(subj_raw))
  subj_int    <- match(subj_raw, subj_sorted)
  
  means <- rowMeans(draws_mat)
  
  qs <- t(apply(draws_mat, 1, quantile, probs = c(0.025, 0.975)))
  L <- qs[, 1]
  U <- qs[, 2]
  
  data.frame(
    SubjectIndex   = subj_int,               # 1, 2, 3, ...
    SubjectOrigID  = subj_raw,              # 1025, 1092, 5001, ...
    M        = means,
    L        = L,
    U        = U,
    analysis = analysis_label
  )
}


## function to make cumulative–weight summaries
summarise_cumweights <- function(mat, analysis_label) {
  
  # column 1 = original subject ID (e.g. 1025–5001)
  subj_raw       <- as.numeric(mat[, 1])
  draws_mat_raw  <- mat[, -1, drop = FALSE]  # posterior draws only
  
  # order subjects by their original numeric ID
  ord            <- order(subj_raw)
  subj_sorted    <- subj_raw[ord]
  draws_ord      <- draws_mat_raw[ord, , drop = FALSE]
  
  # normalize each posterior draw so weights sum to 1 (then convert to %)
  draws_prop     <- apply(draws_ord, 2, function(x) x / sum(x))
  # cumulative sums across subjects for each posterior draw
  cum_mat        <- apply(draws_prop, 2, cumsum) * 100  # now 0–100%
  
  # posterior mean and 95% CrI of cumulative % at each subject rank
  M  <- rowMeans(cum_mat)
  qs <- t(apply(cum_mat, 1, quantile, probs = c(0.025, 0.975)))
  L  <- qs[, 1]
  U  <- qs[, 2]
  
  data.frame(
    SubjectRank  = seq_along(subj_sorted),  # 1, 2, 3, ...
    SubjectOrigID = subj_sorted,            # original IDs if you want them
    M = M,
    L = L,
    U = U,
    analysis = analysis_label,
    stringsAsFactors = FALSE
  )
}

###############################################
for(cc in c(0, 0.01, 0.05, 0.1)){
  gc()
  if(cc==0) cname = ""
  if(cc==0.01) cname = "p01"
  if(cc==0.05) cname = "p05"
  if(cc==0.1) cname = "p1"
  
  if(cc == 0){
    wt_mat1 <- readRDS(paste0(datapath, "weights/mITT_Weights20250912_11.Rds"))
    wt_mat2 <- readRDS(paste0(datapath, "weights/PP_Weights20250426_11.Rds"))
    wt_mat3 <- readRDS(paste0(datapath, "weights/PP_nocs_Weights20250912_11.Rds"))
  } else {
    wt_mat1 = readRDS(paste0(datapath, "weights/mITT_Weights20251026_11_",cname,".Rds"))# ITT
    wt_mat2 =  readRDS(paste0(datapath, "weights/PP_Weights20251026_11_",cname,".Rds"))# RCT
    wt_mat3 =  readRDS(paste0(datapath, "weights/PP_nocs_Weights20251026_11_",cname,".Rds")) # RCT compliers
    
  }
  
 
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  ###########
  ### weight bar plots
  pd_ITT          <- summarise_weights(wt_mat1, "ITT")
  pd_RCT          <- summarise_weights(wt_mat2, "RCT")
  pd_RCT_complier <- summarise_weights(wt_mat3, "RCT-compliers")
  
  pd <- bind_rows(pd_ITT, pd_RCT, pd_RCT_complier)
  
  # If you want SubjectIndex as factor for plotting:
  pd$Subject <- factor(pd$SubjectIndex, levels = sort(unique(pd$SubjectIndex)))
  
  library(ggplot2)
  
  # use numeric index for x
  pd$Subject_idx <- pd$SubjectIndex
  
  pos <- position_dodge(width = 0.7)
  
  p <- ggplot(pd, aes(x = Subject_idx, y = M, colour = analysis)) +
    geom_linerange(aes(ymin = L, ymax = U), position = pos, linewidth = 0.3) +
    geom_point(position = pos, size = 0.9) +
    labs(x = "Subject", y = "Weight %", colour = "analysis") +
    scale_colour_manual(
      values = c(
        "ITT"            = "#009E73",
        "RCT-compliers"  = "#D55E00",
        "RCT"            = "#000000"
      )
    ) +
    theme_minimal(base_size = 9) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "right",
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background  = element_rect(fill = "white", colour = NA),
      text            = element_text(colour = "black"),
      axis.text       = element_text(colour = "black"),
      axis.title      = element_text(colour = "black"),
      legend.text     = element_text(colour = "black"),
      legend.title    = element_text(colour = "black"))+ ylim(c(0, 7.5))

  
  # save as a long, short figure like your example
  if(cc == 0){
    ggsave("subject_weights_by_analysis.png",
           p, width = 13, height = 2, dpi = 300, bg = "white" )
  } else {
    ggsave(paste0("subject_weights_by_analysis_",cname,".png"),
           p, width = 13, height = 2, dpi = 300, bg = "white" )
  }
  
  ###########
  ### cumulative weights
  ## build pd_new from the three matrices
  pd_ITT          <- summarise_cumweights(wt_mat1, "ITT")
  pd_RCT          <- summarise_cumweights(wt_mat2, "RCT")
  pd_RCT_complier <- summarise_cumweights(wt_mat3, "RCT-compliers")
  
  pd_new <- bind_rows(pd_ITT, pd_RCT_complier, pd_RCT)
  
  ## plot: cumulative percent of total weight
  p_cum <- ggplot(pd_new,
                  aes(x = SubjectRank, y = M,
                      colour = analysis, fill = analysis)) +
    geom_ribbon(aes(ymin = L, ymax = U), alpha = 0.3, colour = NA) +
    geom_line(linewidth = 0.8) +
    scale_y_continuous(
      name   = "Percent of total weight",
      limits = c(0, 100),
      breaks = seq(0, 100, by = 25),
      labels = function(x) paste0(x, "%")
    ) +
    scale_x_continuous(name = "Subjects") +
    scale_colour_manual(
      values = c(
        "ITT"            = "#009E73",
        "RCT-compliers"  = "#D55E00",
        "RCT"            = "#000000"
      )
    ) +
    scale_fill_manual(
      values = c(
        "ITT"            = "#009E73",
        "RCT-compliers"  = "#D55E00",
        "RCT"            = "#000000"
      )
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      legend.title     = element_text(),
      legend.position  = "right"
    )
  
  if(cc == 0){
    ggsave("cum_weights_plot.png", p_cum,
           width = 8, height = 5, dpi = 300, bg = "white")
  } else {
    ggsave(paste0("cum_weights_plot_",cname,".png"), p_cum,
           width = 8, height = 5, dpi = 300, bg = "white")
  }
  
}
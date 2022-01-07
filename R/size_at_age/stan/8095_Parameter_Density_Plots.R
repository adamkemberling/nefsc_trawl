
# This code creates a to illustrate the 90% and 60% credible intervals of parameter estimates from a hierarchical von Bertalanffy growth function (VBGF) following Midway et al. 2015

##### READ AND PREPARE DATA #####
haddock <- readRDS(file = "Data/haddock_bio.rds")

# Load model
load("Stan models/VBGF_Stan_models_2022_01_03.RData") # May need to change

# Get parameters
model <- mod_list[[2]] # Hierarchical model

######################################################
# Get Linf, K, and t0
draws <- as.data.frame(model)
cohort.id <- data.frame(cohort = unique(haddock$yearclass), cohort_no = as.numeric(unique(haddock$yearclass)))
design_mat <- matrix(c(rep(1, nrow(haddock))), nrow = nrow(haddock))


################
# Midway plot ##
################
cols <- c("#007FFF","#FF7F00") # Colors for the VBGF lines
cols.points <- c("#80bfff", "#ffc180")
par(mar=c(3.7 , 0.5 , .25 , .125) +0.1, tcl=-.25 , mgp=c(2.5,  .7 ,  0) ,  oma=c(0 , 0, 0 , 0), cex = 0.75)


###### Plot Loo density #######
plot(NA,NA, main = NA, xlab =expression(italic("L"[paste(infinity,sep="")])~(cm)), yaxt = "n", ylab= NA, xlim = c(50, 100), ylim = c(-44.5,6.5), cex.lab = 1.25)
for ( j in 1:nrow(cohort.id)){
  
  # Get parameters
  log_linf_re <- draws[,paste("log_linf_re[",j,"]", sep = "")] # Subset the mcmc chains
  betas_linf_re <- draws[,paste("B_log_linf[",1,"]", sep = "")] # Subset the mcmc chains
  
  dat.sub <- exp(design_mat %*% betas_linf_re + log_linf_re) # Linf on natural scale
  quantiles_025_975 <- quantile(dat.sub,c(0.025, 0.975))
  quantiles_10_90 <- quantile(dat.sub,c(0.1, 0.9))
  quantiles_99 <- quantile(dat.sub,c(0.99))
  
  spacing.ind <- 7-(j + (.35)) # Spacing index for plots
  
  # Plot the credible intervals
  segments(quantiles_025_975[1],spacing.ind,quantiles_025_975[2],spacing.ind, lwd = 3 , col = cols.points[1])
  segments(quantiles_10_90[1],spacing.ind,quantiles_10_90[2],spacing.ind, lwd = 5 , col = cols[1])
  
  # Plot the median estimate
  med <- median(dat.sub)
  segments(med,spacing.ind-.05,med,spacing.ind+.05, lwd = 3 , col = 1)
  
  # Label lines with cohort year
  text(quantiles_99, spacing.ind, paste(cohort.id[j,1],sep=" ") , cex = 0.75) 
}


###### plot K density ######
plot(NA,NA, main = NA, xlab =expression(italic("k")~(yr^-1)), yaxt = "n", ylab= NA, xlim = c(.12, 0.45), ylim = c(-44.5,6.5), cex.lab = 1.25)
for ( j in 1:nrow(cohort.id)){
  
  # Get parameters
  log_k_re <- draws[,paste("log_k_re[",j,"]", sep = "")] # Subset the mcmc chains
  betas_k_re <- draws[,paste("B_log_k[",1,"]", sep = "")] # Subset the mcmc chains
  
  dat.sub <- exp(design_mat %*% betas_k_re + log_k_re) # K on natural scale
  quantiles_025_975 <- quantile(dat.sub,c(0.025, 0.975))
  quantiles_10_90 <- quantile(dat.sub,c(0.1, 0.9))
  quantiles_99 <- quantile(dat.sub,c(0.99))
  
  spacing.ind <- 7-(j + (.35)) # Spacing index for plots
  
  # Plot the credible intervals
  segments(quantiles_025_975[1],spacing.ind,quantiles_025_975[2],spacing.ind, lwd = 3 , col = cols.points[1])
  segments(quantiles_10_90[1],spacing.ind,quantiles_10_90[2],spacing.ind, lwd = 5 , col = cols[1])
  
  # Plot the median estimate
  med <- median(dat.sub)
  segments(med,spacing.ind-.05,med,spacing.ind+.05, lwd = 3 , col = 1)
  
  # Label lines with cohort year
  text(quantiles_99, spacing.ind, paste(cohort.id[j,1],sep=" ") , cex = 0.75) 
}


###### plot t0 density ######
plot(NA,NA, main = NA, xlab =expression(italic(t[paste(0,"s,r",sep="")])~(yr)), yaxt = "n", ylab= NA, xlim = c(-2.25,-.05), ylim = c(-44.5,6.5), cex.lab = 1.25)
for ( j in 1:nrow(cohort.id)){
  
  # Get parameters
  t0_re <- draws[,paste("t0_re[",j,"]", sep = "")] # Subset the mcmc chains
  B_t0 <- draws[,paste("B_t0[",1,"]", sep = "")] # Subset the mcmc chains
  
  dat.sub <- design_mat %*% B_t0 + t0_re # on natural scale
  quantiles_025_975 <- quantile(dat.sub,c(0.025, 0.975))
  quantiles_10_90 <- quantile(dat.sub,c(0.1, 0.9))
  quantiles_99 <- quantile(dat.sub,c(0.99))
  
  spacing.ind <- 7-(j + (.35)) # Spacing index for plots
  
  # Plot the credible intervals
  segments(quantiles_025_975[1],spacing.ind,quantiles_025_975[2],spacing.ind, lwd = 3 , col = cols.points[1])
  segments(quantiles_10_90[1],spacing.ind,quantiles_10_90[2],spacing.ind, lwd = 5 , col = cols[1])
  
  # Plot the median estimate
  med <- median(dat.sub)
  segments(med,spacing.ind-.05,med,spacing.ind+.05, lwd = 3 , col = 1)
  
  # Label lines with cohort year
  text(quantiles_99, spacing.ind, paste(cohort.id[j,1],sep=" ") , cex = 0.75) 
}


###### plot Gallucci OMEGA density ######
plot(NA,NA, main = NA, xlab =expression(omega~(mm~yr^-1)), yaxt = "n", ylab= NA, xlim = c(5,30), ylim = c(-44.5,6.5), cex.lab = 1.25)
for ( j in 1:nrow(cohort.id)){
  
  # Get parameters
  log_linf_re <- draws[,paste("log_linf_re[",j,"]", sep = "")] # Subset the mcmc chains
  betas_linf_re <- draws[,paste("B_log_linf[",1,"]", sep = "")] # Subset the mcmc chains
  
  # Get parameters
  log_k_re <- draws[,paste("log_k_re[",j,"]", sep = "")] # Subset the mcmc chains
  betas_k_re <- draws[,paste("B_log_k[",1,"]", sep = "")] # Subset the mcmc chains
  
  dat.sub <- exp(design_mat %*% betas_k_re + log_k_re) * exp(design_mat %*% betas_linf_re + log_linf_re) # omega on natural scale
  quantiles_025_975 <- quantile(dat.sub,c(0.025, 0.975))
  quantiles_10_90 <- quantile(dat.sub,c(0.1, 0.9))
  quantiles_99 <- quantile(dat.sub,c(0.99))
  
  spacing.ind <- 7-(j + (.35)) # Spacing index for plots
  
  # Plot the credible intervals
  segments(quantiles_025_975[1],spacing.ind,quantiles_025_975[2],spacing.ind, lwd = 3 , col = cols.points[1])
  segments(quantiles_10_90[1],spacing.ind,quantiles_10_90[2],spacing.ind, lwd = 5 , col = cols[1])
  
  # Plot the median estimate
  med <- median(dat.sub)
  segments(med,spacing.ind-.05,med,spacing.ind+.05, lwd = 3 , col = 1)
  
  # Label lines with cohort year
  text(quantiles_99, spacing.ind, paste(cohort.id[j,1],sep=" ") , cex = 0.75) 
}

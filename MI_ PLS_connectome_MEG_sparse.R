#Global connectome measures
#Setup
library(mixOmics)
library(psych)
library(Matrix)
library(gtools)
library(ggplot2)
library(MCMCpack)


#Load data
mydata_fa <- readRDS('U:/Documents/ACE_Data/Thesis_Analysis/MI_FA/FA_results/mydata_fa_scaled.RData')

##### Start Questions, Outcomes, Cognition setup (full sample) ######
wd <- "U:/Documents/ACE_Data/Thesis_Analysis/MI_PLS"
setwd(wd)
questions <- NULL
cognitive <- NULL
academic <- NULL
behaviour <- NULL
for(imp in 1:5){
  questions[[imp]] <- mydata_fa[[imp]][,1:22] 
  cognitive[[imp]] <- mydata_fa[[imp]][,23:26] 
  academic[[imp]] <- mydata_fa[[imp]][,27:28] 
  behaviour[[imp]] <- mydata_fa[[imp]][,29:30] 
}
##### End Questions, Outcomes, Cognition setup (full sample) ######

##### Start MRI Connectome setup (MRI sample) ######
wd <- setwd("U:/Documents/ACE_Data/Thesis_Analysis/MI_PLS_MRI")

#Load data
graph_measures <- read.csv('graph_measures_FA.csv', header = TRUE)
node_degree <- read.csv('node_degree_FA.csv', header = FALSE)
node_strength <- read.csv('node_strength_FA.csv', header = FALSE)
node_efficiency <- read.csv('node_efficiency_FA.csv', header = FALSE)
tracts <- read.csv('FA_values_by_tract_reshaped.csv', header = TRUE)
half_nets <- read.csv('half_nets_FA.csv', header = FALSE) #Unravelled half adjacency matrix here (SubjectXconnection)
labels_85 <- read.csv('labels_85.txt',header = FALSE)
unravelled_labels_85 <- read.csv('unravelled_labels_85.csv')
colnames(half_nets) <- t(unravelled_labels_85[,2])
colnames(node_degree) <- t(labels_85)
colnames(node_strength) <- t(labels_85)

#Scale the brain data 
graph_measures[,2:ncol(graph_measures)] <- scale(graph_measures[,2:ncol(graph_measures)], center=TRUE, scale=TRUE)
node_degree <- scale(node_degree, center=TRUE, scale=TRUE)
node_strength <- scale(node_strength, center=TRUE, scale=TRUE)
tracts[,3:ncol(tracts)] <- scale(tracts[,3:ncol(tracts)], center=TRUE, scale=TRUE)
half_nets <- scale(half_nets, center=TRUE, scale=TRUE)

#Get MRI ID's
MRI.ID_FA <- read.csv('MRI.ID_for_FA_data.csv', header = TRUE)
MRI.ID_brain <- data.frame(graph_measures$MRI.ID)
colnames(MRI.ID_brain) <- 'MRI.ID'
node_degree <- cbind(MRI.ID_brain, data.frame(node_degree))
colnames(node_degree)[1] <- 'MRI.ID'
node_strength <- cbind(MRI.ID_brain, data.frame(node_strength))
colnames(node_strength)[1] <- 'MRI.ID'
half_nets <- cbind(MRI.ID_brain, data.frame(half_nets))
colnames(half_nets)[1] <- 'MRI.ID'
#Remove zero columns
# na_node_strength <- which(colSums(!is.na(node_strength)) == 0) #If needed to remove zero columns
# na_half_nets <- which(colSums(!is.na(half_nets)) == 0) #If needed to remove zero columns
# node_strength <- node_strength[,colSums(!is.na(node_strength)) > 0] #If needed to remove zero columns
# half_nets <- half_nets[,colSums(!is.na(half_nets)) > 0] #If needed to remove zero columns
#node_degree <- node_degree[,colSums(!is.na(node_degree)) > 0] #If needed to remove zero columns

#Remove participants not in MRI ID list
mydata_fa_MRI <- NULL
for(imp in 1:5){
  mydata_fa_MRI[[imp]] <- cbind(MRI.ID_FA,mydata_fa[[imp]])
  mydata_fa_MRI[[imp]] <- merge(mydata_fa_MRI[[imp]],MRI.ID_brain,by='MRI.ID')
  mydata_fa_MRI[[imp]]$MRI.ID <- NULL
}
graph_measures <- merge(graph_measures,MRI.ID_FA,by='MRI.ID')
graph_measures$MRI.ID <- NULL
node_degree <- merge(node_degree,MRI.ID_FA,by='MRI.ID')
node_degree$MRI.ID <- NULL
node_strength <- merge(node_strength,MRI.ID_FA,by='MRI.ID')
node_strength$MRI.ID <- NULL
tracts <- merge(tracts,MRI.ID_brain,by='MRI.ID') #First merge by brain ID to remove high motion scans
tracts <- merge(tracts,MRI.ID_FA,by='MRI.ID')
tracts$MRI.ID <- NULL
tracts$ID <- NULL
half_nets <- merge(half_nets,MRI.ID_FA,by='MRI.ID')
half_nets$MRI.ID <- NULL

#Set up data
questions <- NULL
cognitive <- NULL
academic <- NULL
behaviour <- NULL
for(imp in 1:5){
  questions[[imp]] <- mydata_fa_MRI[[imp]][,1:22] 
  cognitive[[imp]] <- mydata_fa_MRI[[imp]][,23:26] 
  academic[[imp]] <- mydata_fa_MRI[[imp]][,27:28] 
  behaviour[[imp]] <- mydata_fa_MRI[[imp]][,29:30] 
}

wd <- "U:/Documents/ACE_Data/Thesis_Analysis/MI_PLS_MRI/results_non_sparse/"
setwd(wd)
##### End MRI Connectome setup (MRI sample) ######

##### Start MEG Connectome setup (MEG sample) ######
#Getting MEG connectome
wd <- "U:/Documents/ACE_Data/Thesis_Analysis/MI_PLS_MEG/"
setwd(wd)

#Load connectome data
#Select frequency range and correlation matrix type
freq_range <- 'beta_nets_LASSO' #'alpha_nets', 'beta_nets', 'beta_nets_LASSO'
corr_type <- 'envPartialCorrelationRegularized_z'

#If you want to threshold it 
use_mask <- FALSE
mask <- data.frame(read.csv(paste(wd, freq_range, '/', 'groupEnvPartialCorrelation_z.csv', sep='' ), header=FALSE))
z_thresh <- 3.75848587174052
mask[mask<z_thresh] <- 0
mask[mask>=z_thresh] <- 1

data_dir <- paste(wd, freq_range, '/', corr_type, '/', sep='' )
setwd(wd)
ID_brain <- data.frame(read.csv('IDs_for_MEG_over_5min.csv', header = TRUE)[,1])
colnames(ID_brain) <- 'ID'
k <- 1
for(i in ID_brain$ID){
  fname <- paste(data_dir, 'nets_bAdspm12_', i, '_mc_rest.csv', sep='')
  nets <- data.frame(read.csv(fname, header=FALSE))
  #nets[nets<z_thresh] <- NA
  if(use_mask == TRUE){
    nets <- mask*nets
  }
  if(k ==1){
    all_nets <- nets[lower.tri(nets)]
  }
  if(k>1){
    all_nets <- rbind(all_nets,nets[lower.tri(nets)])
  }
  k <- k+1
}
half_nets <- scale(all_nets,  center=TRUE, scale=TRUE)

#Other data
mydata_fa <- readRDS('U:/Documents/ACE_Data/Thesis_Analysis/MI_FA/FA_results/mydata_fa_scaled.RData')
labels_68 <- read.csv('labels_MEG.txt',header = FALSE)
unravelled_labels_68 <- read.csv('unravelled_labels_68.csv')
colnames(half_nets) <- t(unravelled_labels_68[,2])

#Match on IDs, removing participants not in ID list
ID_FA <- read.csv('IDs_for_FA_scores_MEG.csv', header = TRUE)
mydata_fa_MEG <- NULL
for(imp in 1:5){
  mydata_fa_MEG[[imp]] <- cbind(ID_FA,mydata_fa[[imp]])
  mydata_fa_MEG[[imp]] <- merge(mydata_fa_MEG[[imp]],ID_brain,by='ID')
  mydata_fa_MEG[[imp]]$ID <- NULL
}
half_nets <- cbind(ID_brain,half_nets)
half_nets <- merge(half_nets,ID_FA,by='ID')
half_nets$ID <- NULL

#Set up data
questions <- NULL
cognitive <- NULL
academic <- NULL
behaviour <- NULL
for(imp in 1:5){
  questions[[imp]] <- mydata_fa_MEG[[imp]][,1:22] 
  cognitive[[imp]] <- mydata_fa_MEG[[imp]][,23:26] 
  academic[[imp]] <- mydata_fa_MEG[[imp]][,27:28] 
  behaviour[[imp]] <- mydata_fa_MEG[[imp]][,29:30] 
}

wd <- ("U:/Documents/ACE_Data/Thesis_Analysis/MI_PLS_MEG/results_non_sparse/")
setwd(wd)

##### End MEG Connectome setup (MEG sample) ######

#### Run PLS Permutation ####
x_orig <- NULL
y_orig <- NULL
orig_fit <- NULL
orig_cov_LV <- NULL
orig_cor_LV <- NULL
nperm = 1000
permuted_cov_LV <- list()
permuted_cor_LV <- list()
permuted_cov_LV_proc <- list()
permuted_cor_LV_proc <- list()
permuted_LV <- list()
perm_fit <- NULL
sample_for_perm <- NULL
pval <- list()
n_CI <- round((0.05*(nperm+1)) - 1,0) #maximum nunber of permuted SV's that can be greater than origional and it still be significant
CI_upper <- list()
F <- NULL
G <- NULL
H <- NULL

signs_X <- NULL
signs_Y <- NULL

for(contrast in c(1)){
  for(imp in 1:5){
    if(contrast==1){#PLS connection strength - academic
      x_orig[[imp]] <- half_nets
      A <- paste(freq_range, '_', corr_type, sep='')
      y_orig[[imp]] <- academic[[imp]]
      B <- 'Academic'
      ncomp=2
    }
    if(contrast==2){#PLS connection strength - behaviour
      x_orig[[imp]] <- half_nets
      A <- paste(freq_range, '_', corr_type, sep='')
      y_orig[[imp]] <- behaviour[[imp]]
      B <- 'Behaviour'
      ncomp=2
    }
    if(contrast==3){#PLS connection strength- questions
      x_orig[[imp]] <- half_nets
      A <- paste(freq_range, '_', corr_type, sep='')
      y_orig[[imp]] <- questions[[imp]]
      B <- 'Questions'
      ncomp=ncol(questions[[imp]])
    }
    
    #sparse PLS
    orig_fit[[imp]] <- spls(x_orig[[imp]], y_orig[[imp]],  mode = "canonical", ncomp =ncomp, keepX=c(100,100), keepY = c(2,2), scale=TRUE)
    
    if(imp >1){
      orig_fit[[imp]]$loadings$X <- procrustes(as.matrix(orig_fit[[imp]]$loadings$X), as.matrix(orig_fit[[1]]$loadings$X))$X.new
      orig_fit[[imp]]$loadings$Y <- procrustes(as.matrix(orig_fit[[imp]]$loadings$Y), as.matrix(orig_fit[[1]]$loadings$Y))$X.new
      orig_fit[[imp]]$variates$X <- procrustes(as.matrix(orig_fit[[imp]]$variates$X), as.matrix(orig_fit[[1]]$variates$X))$X.new
      orig_fit[[imp]]$variates$Y <- procrustes(as.matrix(orig_fit[[imp]]$variates$Y), as.matrix(orig_fit[[1]]$variates$Y))$X.new
    }
    orig_cov_LV[[imp]] <- diag(cov(orig_fit[[imp]]$variates$X, orig_fit[[imp]]$variates$Y)) #Covariance between LV = Singular Value
    orig_cor_LV[[imp]] <- diag(cor(orig_fit[[imp]]$variates$X, orig_fit[[imp]]$variates$Y)) #Correlation between LV 
    
    #Fit PLS model to permuted data to get p values
    for (i in 1:nperm){
      if(i == 1){
        print(A)
        print(B)
        print(imp)
      }
      if(i == seq(0,nperm,50)){
        print(i)
      }
      if(imp ==1){
        sample_for_perm[[i]] <- sample(nrow(x_orig[[imp]])) #to enseure same sample is pulled out for each imputed set
      }
      x_perm <- x_orig[[imp]][sample_for_perm[[i]],] #Permuted data
      perm_fit[[imp]] <- pls(x_perm, y_orig[[imp]],  mode = "canonical", ncomp = ncomp, scale=TRUE)
      if(i ==1){
        permuted_cov_LV[[imp]] <- diag(cov(perm_fit[[imp]]$variates$X, perm_fit[[imp]]$variates$Y))
        permuted_cor_LV[[imp]] <- diag(cor(perm_fit[[imp]]$variates$X, perm_fit[[imp]]$variates$Y))
      }
      else{
        x_perm_orig <- x_orig[[1]][sample_for_perm[[i]],] #Permuted data
        perm_fit_orig <- pls(x_perm_orig, y_orig[[1]],  mode = "canonical", ncomp = ncomp, scale=TRUE)
        perm_fit[[imp]]$variates$X <- procrustes(as.matrix(perm_fit[[imp]]$variates$X), as.matrix(perm_fit_orig$variates$X))$X.new
        perm_fit[[imp]]$variates$Y <- procrustes(as.matrix(perm_fit[[imp]]$variates$Y), as.matrix(perm_fit_orig$variates$Y))$X.new
        permuted_cov_LV[[imp]] <- rbind(permuted_cov_LV[[imp]],diag(cov(perm_fit[[imp]]$variates$X, perm_fit[[imp]]$variates$Y)))
        permuted_cor_LV[[imp]] <- rbind(permuted_cor_LV[[imp]],diag(cor(perm_fit[[imp]]$variates$X, perm_fit[[imp]]$variates$Y)))
      }
    }
    for(comp in 1:ncomp){
      if(comp ==1){
        pval[[imp]] <- (sum(permuted_cov_LV[[imp]][,comp]>orig_cov_LV[[imp]][comp])+1)/(nperm + 1)
        CI_upper[[imp]] <- data.frame(sort(permuted_cov_LV[[imp]][,comp], decreasing = TRUE))[n_CI,1]
      }
      else {
        pval[[imp]] <- rbind(pval[[imp]],(sum(permuted_cov_LV[[imp]][,comp]>orig_cov_LV[[imp]][comp])+1)/(nperm + 1))
        CI_upper[[imp]] <- rbind(CI_upper[[imp]],data.frame(sort(permuted_cov_LV[[imp]][,comp], decreasing = TRUE))[n_CI,1])
      }
    }
    pval[[imp]] <- round(pval[[imp]], 5)
    pval[[imp]] <- data.frame(cbind(1:ncomp, pval[[imp]]))
    
    if(imp==1){
      F <- pval[[imp]]
      G <- CI_upper[[imp]]
      H <- orig_cov_LV[[imp]]
      
    }
    if(imp>1){
      F <- F + pval[[imp]]
      G <- G + CI_upper[[imp]]
      H <- H + orig_cov_LV[[imp]]
    }
  }
  
  pvalues <- NULL
  G <- NULL
  H <- NULL
  J <- NULL
  for(imp in 1:5){
    if(imp==1){
      pvalues <- pval[[imp]][,2]
      G <- CI_upper[[imp]]
      H <- orig_cov_LV[[imp]]
      J <- orig_cor_LV[[imp]]
      
      
    }
    if(imp>1){
      pvalues <-cbind(pvalues,pval[[imp]][,2])
      G <- G + CI_upper[[imp]]
      H <- H+ orig_cov_LV[[imp]]
      J <- J + orig_cor_LV[[imp]]
    }
  }
  #Use Rubin's Z tranform method to pool p values
  p_pooled <- NULL 
  for(comp in 1:ncomp){
    z <- qnorm(1-pvalues[comp,])  # transform to z-scale
    zmean <- mean(z)
    imp_var <- sum((zmean-z)^2)/(5-1)
    total_var <- 1 + (1 + (1/5))*imp_var
    p_pooled[comp] <- 1-pnorm( zmean / sqrt(total_var)) # average and transform back
  }
  
  #Save results 
  flip <- 1
  #flip <- sign(orig_fit[[1]]$loadings$Y[1,1])
  exp_var <- cbind(data.frame(round(orig_fit[[1]]$explained_variance$X,4)),data.frame(round(orig_fit[[1]]$explained_variance$Y,4)))
  First_LV <- cbind(data.frame(flip*orig_fit[[1]]$variates$X[,1]),data.frame(flip*orig_fit[[1]]$variates$Y[,1]))
  X_1st_loadings <- data.frame(flip*orig_fit[[1]]$loadings$X[,2])
  Y_1st_loadings <- data.frame(flip*orig_fit[[1]]$loadings$Y[,2])
  for(imp in 2:5){
    exp_var <- cbind(exp_var, data.frame(round(orig_fit[[imp]]$explained_variance$X,4)),data.frame(round(orig_fit[[imp]]$explained_variance$Y,4)))
    First_LV <- cbind(First_LV,data.frame(flip*orig_fit[[imp]]$variates$X[,1]),data.frame(flip*orig_fit[[imp]]$variates$Y[,1]))
    X_1st_loadings <- cbind(X_1st_loadings,flip*orig_fit[[imp]]$loadings$X[,2])
    Y_1st_loadings <- cbind(Y_1st_loadings,flip*orig_fit[[imp]]$loadings$Y[,2])
  }
  colnames_all <- rep(c(A,B),5)
  colnames(exp_var) <- colnames_all
  colnames(First_LV) <- colnames_all
  colnames(X_1st_loadings) <- rep(A,5)
  colnames(Y_1st_loadings) <- rep(B,5)
  #Save tables
  fname <- paste(wd, '/' , A, '_', B,'_exp_var.csv', sep='' )
  write.csv(exp_var, file = fname)
  fname <- paste(wd, '/' ,A, '_', B,'_1st_LV.csv', sep='' )
  write.csv(First_LV, file = fname)
  fname <- paste(wd, '/' ,A, '_', B,  '_X_2nd_loadings.csv', sep='' )
  write.csv(X_1st_loadings, file = fname)
  fname <- paste(wd, '/' , A, '_', B, '_Y_2nd_loadings.csv', sep='' )
  write.csv(Y_1st_loadings, file = fname)
  av_pval <- cbind(c(1:ncomp), data.frame(p_pooled))
  av_CI_upper <- G/5
  av_orig_cov_LV <- H/5
  
  fname <- paste(wd, '/' , A, '_', B,'_pvals.csv', sep='' )
  df <- cbind(av_pval,av_CI_upper)
  colnames(df) <- c('p','SV_CI')
  write.csv(df, file = fname)
  
  av_pval_group <- data.frame(rep(1,nrow(av_pval))) #To get rule of thumb reliable loadings
  av_pval_group[which(av_pval[,2] > 0.05),] <- 0
  av_pval <- cbind(av_pval, av_pval_group)
  colnames(av_pval) <- c('X1', 'X2', 'significant')
  av_pval$significant <- as.factor(av_pval$significant)
  
  ggplot(av_pval, aes(X1,X2, fill=significant)) + 
    geom_bar(stat="identity") + 
    xlab("Component") + 
    ylab("P Value") +
    scale_x_continuous(breaks = seq(1, ncomp, by = 1))+
    scale_fill_manual(values=c( "#999999", "#ec6c20")) +
    theme_minimal(base_size=15)+
    theme(axis.text.x = element_text(size=10, hjust=1), axis.text.y = element_text(size=10, hjust=1))+
    guides(fill=FALSE)
  fname <- paste(wd, '/' , A, '_', B,'_pvals.png', sep='' )
  ggsave(fname, width = 10, height = 10, units = "cm")
  #Plot components with CI
  CI_lower = rep(0, ncomp)
  ggplot(data.frame(av_orig_cov_LV), aes(av_orig_cov_LV)) + 
    geom_line(aes(x=1:ncomp,y=av_orig_cov_LV), colour="#ec6c20") + 
    geom_ribbon(aes(x=1:ncomp,ymin=CI_lower, ymax=av_CI_upper), alpha=0.2)  + 
    xlab("Component") + 
    ylab("Covariance between LVs")+
    scale_x_continuous(breaks = seq(1, ncomp, by = 1)) +
    theme_minimal(base_size=15)+
    theme(axis.text.x = element_text(size=10, hjust=1), axis.text.y = element_text(size=10, hjust=1))
  fname <- paste(wd, '/' , A, '_', B,'_SVs.png', sep='' )
  ggsave(fname, width = 20, height = 10, units = "cm")
  
}

###########Bootstrap #######################
library(boot)

cols <- c('0' = '#999999','1'='#ec6c20') #set colours so that sig are orange
ncomp = 1 #Number components to extract
comp=1 #Which of these components to bootstrap
nboot = 500
orig_fit <- NULL
perm_fit <- NULL
X_boot_result <- NULL
Y_boot_result <- NULL
X_SE <- NULL
Y_SE <- NULL
x_orig <- NULL
y_orig <- NULL
for(contrast in c(1)){
  for(imp in 1:5){
    if(contrast==1){#PLS connection strength - academic
      x_orig[[imp]] <- half_nets
      A <- paste(freq_range, '_', corr_type, sep='')
      y_orig[[imp]] <- academic[[imp]]
      B <- 'Academic'
    }
    if(contrast==2){#PLS connection strength - behaviour
      x_orig[[imp]] <- half_nets
      A <- paste(freq_range, '_', corr_type, sep='')
      y_orig[[imp]] <- behaviour[[imp]]
      B <- 'Behaviour'
    }
    if(contrast==3){#PLS connection strength- questions
      x_orig[[imp]] <- half_nets
      A <- paste(freq_range, '_', corr_type, sep='')
      y_orig[[imp]] <- questions[[imp]]
      B <- 'Questions'
    }
    
    print(A)
    print(B)
    print(imp)
    data <- cbind(x_orig[[imp]], y_orig[[imp]])
    orig_fit[[imp]] <- pls(x_orig[[imp]], y_orig[[imp]],  mode = "canonical", ncomp =ncomp, scale=TRUE)
    X_boot <- function(d,j){
      d2 <- d[j,]
      rownames(d2)<- NULL
      perm_fit<- pls(d2[,1:ncol(x_orig[[imp]])], d2[,(ncol(x_orig[[imp]])+1):ncol(d2)],  mode = "canonical", ncomp =ncomp, scale=TRUE)
      proc <- NULL
      proc$X <- procrustes(perm_fit$loadings$X,orig_fit[[1]]$loadings$X)
      return(proc$X$X.new[,comp])
    }
    Y_boot <- function(d,j){
      d2 <- d[j,]
      rownames(d2)<- NULL
      perm_fit <- pls(d2[,1:ncol(x_orig[[imp]])], d2[,(ncol(x_orig[[imp]])+1):ncol(d2)],  mode = "canonical", ncomp =ncomp, scale=TRUE)
      proc <- NULL
      proc$Y <- procrustes(perm_fit$loadings$Y,orig_fit[[1]]$loadings$Y)
      return(proc$Y$X.new[,comp])
    }
    X_boot_result[[imp]] <- boot(data, X_boot, R=nboot)
    Y_boot_result[[imp]] <- boot(data, Y_boot, R=nboot)
    X_SE[[imp]] <- diag(sqrt(var(X_boot_result[[imp]]$t)))
    Y_SE[[imp]] <- diag(sqrt(var(Y_boot_result[[imp]]$t)))
  }
  G <- NULL
  H <- NULL
  I <- NULL
  J <- NULL
  for(imp in 1:5){
    if(imp==1){
      G <- X_boot_result[[imp]]$t0
      H <- Y_boot_result[[imp]]$t0
      I <- (X_SE[[imp]])^2
      J <- (Y_SE[[imp]])^2
      
    }
    if(imp>1){
      G <- G  + X_boot_result[[imp]]$t0
      H <- H + Y_boot_result[[imp]]$t0
      I <- I + (X_SE[[imp]])^2
      J <- J + (Y_SE[[imp]])^2 
    }
  }
  av_X_boot_result <- G/5
  av_Y_boot_result <- H/5
  av_X_SE <- sqrt(I/5) #Note this is only the within imputation standard error
  av_Y_SE <- sqrt(J/5)
  
  #update the standard error
  F <- NULL
  G <- NULL
  for(imp in 1:5){
    if(imp==1){
      F <- (X_boot_result[[imp]]$t0-av_X_boot_result)^2
      G <- (Y_boot_result[[imp]]$t0-av_Y_boot_result)^2
      
    }
    if(imp>1){
      F <- F + (X_boot_result[[imp]]$t0-av_X_boot_result)^2 
      G <- G + (Y_boot_result[[imp]]$t0-av_Y_boot_result)^2
      
    }
  }
  pooled_X_SE <- sqrt((av_X_SE)^2 + ((1+ (1/m))*(1/(m-1))*F)) #Remember var = SE^2
  pooled_Y_SE <- sqrt((av_Y_SE)^2 + ((1+ (1/m))*(1/(m-1))*G)) #Remember var = SE^2
  
  X_ratio <- data.frame(abs(av_X_boot_result)/pooled_X_SE)
  X_ratio_group <- data.frame(rep(1,nrow(X_ratio))) #To get rule of thumb reliable loadings
  X_ratio_group[which(X_ratio < 2),] <- 0
  Y_ratio <- data.frame(abs(av_Y_boot_result)/pooled_Y_SE)
  Y_ratio_group <- data.frame(rep(1,nrow(Y_ratio))) #To get rule of thumb reliable loadings
  Y_ratio_group[which(Y_ratio < 2),] <- 0
  for_brain_plot <- av_X_boot_result*X_ratio_group
  
  flip <- sign(orig_fit[[1]]$loadings$Y[1,1])
  av_X_boot_result <- flip*av_X_boot_result
  av_Y_boot_result <- flip*av_Y_boot_result
  X_load <- NULL
  X_load <- data.frame(round(av_X_boot_result, 2),av_X_boot_result-pooled_X_SE, av_X_boot_result+pooled_X_SE, X_ratio_group,colnames(x_orig[[imp]]))
  colnames(X_load) <- c('loading', 'SE_low', 'SE_high','reliable', 'name')
  fname <- paste(wd,'/', A, '_', B, '_', comp, '_boot_X_loadings.csv', sep='' )
  write.csv(X_load, file = fname)
  X_load <- X_load[order(-abs(X_load$loading)),]
  X_load$name <- factor(X_load$name, levels = X_load$name[order(-abs(X_load$loading))])
  X_load$reliable <- factor(X_load$reliable)
  
  ggplot(data = X_load,
         aes(x = name, y = loading,  fill=reliable)) +
    geom_bar(stat = 'identity', position = 'identity', width = 0.9) +
    #geom_errorbar(aes(ymin=-1*CI_low,ymax=-1*CI_high)) +
    geom_errorbar(aes(ymin=SE_low,ymax=SE_high)) +
    scale_fill_manual(values=cols) +
    theme_minimal(base_size=15)+
    theme(axis.text.x = element_text(angle = 90, size=10, hjust=1,vjust=0.4), axis.title.x = element_blank()) +
    guides(fill=FALSE)
  
  fname <- paste(wd,'/', A, '_', B,'_X_loadings.png', sep='' )
  ggsave(fname, width = 30, height = 20, units = "cm")
  
  Y_load <- data.frame(round(av_Y_boot_result, 2),av_Y_boot_result-pooled_Y_SE, av_Y_boot_result+pooled_Y_SE, Y_ratio_group,colnames(y_orig[[imp]]))
  colnames(Y_load) <- c('loading', 'SE_low', 'SE_high','reliable', 'name')
  fname <- paste(wd,'/', A, '_', B, '_', comp,'_boot_Y_loadings.csv', sep='' )
  write.csv(Y_load, file = fname)
  Y_load <- Y_load[order(-abs(Y_load$loading)),]
  Y_load$name <- factor(Y_load$name, levels = Y_load$name[order(-abs(Y_load$loading))])
  Y_load$reliable <- factor(Y_load$reliable)
  
  ggplot(data = Y_load,
         aes(x = name, y = loading,  fill=reliable)) +
    geom_bar(stat = 'identity', position = 'identity', width = 0.9) +
    #geom_errorbar(aes(ymin=-1*CI_low,ymax=-1*CI_high)) +
    geom_errorbar(aes(ymin=SE_low,ymax=SE_high)) +
    scale_fill_manual(values=cols) +
    theme_minimal(base_size=15)+
    theme(axis.text.x = element_text(angle = 90, size=10, hjust=1,vjust=0.4), axis.title.x = element_blank())+
    guides(fill=FALSE)
  
  fname <- paste(wd,'/', A, '_', B, '_', comp,'_Y_loadings.png', sep='' )
  ggsave(fname, width = 8, height = 10, units = "cm")
  
  X_load$table <- rep(A, nrow(X_load))
  Y_load$table <- rep(B, nrow(Y_load))
  
  
  X_load <- X_load[1:120,]
  
  all_load <- rbind(X_load, Y_load)
  ymin <- min(all_load$SE_low)
  if(ymin>0){
    ymin <- 0
  }
  ggplot(data = all_load,
         aes(x = name, y = loading,  fill=reliable)) +
    geom_bar(stat = 'identity', position = 'dodge', width = 0.9) +
    facet_grid(~table, space='free', scale="free")+
    #geom_errorbar(aes(ymin=-1*CI_low,ymax=-1*CI_high)) +
    geom_errorbar(aes(ymin=SE_low,ymax=SE_high)) +
    scale_fill_manual(values=cols) +
    ylim(ymin, max(all_load$SE_high)) +
    #ggtitle("Barchart of Loadings For First Component") + 
    theme_bw(base_size=15)+
    theme(axis.text.x = element_text(angle = 90, size=10, hjust=1,vjust=0.4), axis.title.x = element_blank()) +
    guides(fill=FALSE)+
    theme(strip.background = element_blank(),
          strip.text.x = element_blank())
  
  fname <- paste(wd,'/', A, '_', B, '_', comp,'_X_and_Y_loadings.png', sep='' )
  ggsave(fname, width = 40, height = 20, units = "cm")
  
}





#Plot loadings later
X_load <- read_csv("U:/Documents/ACE_Data/PLS/PLS_results/Node_degree_Academic_boot_X_loadings.csv")
X_load$type <- c(rep('Cognitive',4), rep('Environment', 20), rep('SES', 4)) #Questionnaire types

X_load <- X_load[order(-abs(X_load$loading)),]
X_load$name <- factor(X_load$name, levels = X_load$name[order(-abs(X_load$loading))])
C = -1
ggplot(data = X_load,
       aes(x = name, y = C*loading, fill=type)) +
  geom_bar(stat = 'identity',position = 'identity') +
  geom_errorbar(aes(ymin=C*CI_low,ymax=C*CI_high)) +
  #geom_errorbar(aes(ymin=SE_low,ymax=SE_high)) +
  ggtitle("Barchart of Loadings For First Component") +
  theme(axis.text.x = element_text(angle = 90, size=10), axis.title.x = element_blank())
fname <- paste('U:/Documents/ACE_Data/PLS/PLS_results/Questions_Behaviour_X_loadings.png', sep='' )
ggsave(fname, width = 30, height = 20, units = "cm")


########################### Other bits ############################################

#Plot SV vs permuted SV's
comp = 1 #Which component to plot
ggplot(data.frame(permuted_cov_LV[,comp]), aes(x=permuted_cov_LV[,comp])) + geom_density() + geom_vline(xintercept=orig_cov_LV[comp], color="red")

#Plot loadings
#Note, get labels for the connectome variables below
plotLoadings(orig_fit, comp=comp)

#These bits are useful for getting the unravelled network matrices back into adjacency matrix format
#Recreate adjacency matrix
X_load <- read.csv('Connection.strength_Questions_boot_X_loadings.csv', header = TRUE)
for_brain_plot <- X_load$loading*X_load$reliable
adj <- matrix(data=NA,nrow=68,ncol=68)
var <- for_brain_plot
adj[lower.tri(adj)] <- as.numeric(var) #Note, has to be lower so that it gets filled in in the right order
adj[is.na(adj)] <- 0

adj = adj + t(adj)
adj <- data.frame(adj)
rownames(adj) <- t(labels_68) 
colnames(adj) <- t(labels_68)
write.csv(adj,file = "adj_pls_questions_connection_strength.csv")

#labels for unravelled nets
adj$labs <- rownames(adj)
adj_melt <- na.omit(melt(adj, 'labs', variable_name='Var1'))
adj_melt$variable <- factor(adj_melt$variable, levels=rev(levels(adj_melt$variable)))
labels <- transform(adj_melt, BLOCKID = paste(labs, variable, sep = " "))
colnames(x_orig) <- t(labels[,4])

#Note, if you use PLS mode regression you can use this to tune the PLS using cross-validation:
tune.pls <- perf(orig_fit, validation = 'Mfold', folds = 5, criterion = 'all', progressBar = FALSE)

#Matrix algebra version of deriving singular values
#Maths can be found here: http://people.revoledu.com/kardi/tutorial/LinearAlgebra/SVD.html
cov_mat <-cov(as.matrix(x_orig), as.matrix(y_orig), use = "pairwise.complete.obs")
Orig_SVs <-t(as.matrix(orig_fit$loadings$X))%*%as.matrix(cov_mat)%*%as.matrix(orig_fit$loadings$Y)

x_orig <- x_orig[[1]]
y_orig <- y_orig[[1]]



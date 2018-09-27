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

wd <- "U:/Documents/ACE_Data/MI_FA/MI_PLS_MRI/results_non_sparse"
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


flip <- 1 #set to 1 or -1 to rotate factor loadings
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
x_orig <- y_orig <- NULL
signs_X <- NULL
signs_Y <- NULL

for(contrast in c(1)){
  for(imp in 1:5){
    if(contrast==1){  #PLS questionnaires- academic
      x_orig[[imp]] <- cbind(questions[[imp]])
      A <- 'Questions'
      y_orig[[imp]] <- academic[[imp]]
      B <- 'Academic' 
    }
    if(contrast==2){#PLS questionnaires- behaviour
      x_orig[[imp]] <- cbind(questions[[imp]])
      A <- 'Questions'
      y_orig[[imp]] <- behaviour[[imp]]
      B <- 'Behaviour'  
    }
    if(contrast==3){#PLS cog and outcome
      x_orig[[imp]] <- questions[[imp]]
      A <- 'Questions'
      y_orig[[imp]] <- cognitive[[imp]]
      B <- 'Cognitive'  
    }
    if(contrast==4){#PLS Cog and outcome
      x_orig[[imp]] <- cognitive[[imp]]
      A <- 'Cognitive'
      y_orig[[imp]] <- academic[[imp]]
      B <- 'Academic'  
    }
    if(contrast==5){#PLS cog and outcome
      x_orig[[imp]] <- cognitive[[imp]]
      A <- 'Cognitive'
      y_orig[[imp]] <- behaviour[[imp]]
      B <- 'Behaviour'  
    }
    ##### MRI #####
    if(contrast==10){#PLS global- questions
      x_orig[[imp]] <- questions[[imp]]
      A <- 'Questions'
      y_orig[[imp]] <- graph_measures
      B <- 'Global.connectome'
    } 
    if(contrast==11){#PLS node degree- questions
      x_orig[[imp]] <- questions[[imp]]
      A <- 'Questions'
      y_orig[[imp]] <- node_degree
      B <- 'Node.degree'
      is_na <- which(colSums(!is.na(y_orig[[imp]])) == 0) #If needed to remove zero columns
      y_orig[[imp]] <- y_orig[[imp]][,colSums(!is.na(y_orig[[imp]])) > 0] #If needed to remove zero columns
    }
    if(contrast==12){#PLS node_strength- questions
      x_orig[[imp]] <- node_strength
      A <- 'Node.strength'
      is_na <- which(colSums(!is.na(x_orig[[imp]])) == 0) #If needed to remove zero columns
      x_orig[[imp]] <- x_orig[[imp]][,colSums(!is.na(x_orig[[imp]])) > 0] #If needed to remove zero columns
      y_orig[[imp]] <- questions[[imp]]
      B <- 'Questions'
    }
    if(contrast==15){#PLS connection strength- questions
      x_orig[[imp]] <- half_nets
      A <- 'Connection.strength'
      is_na <- which(colSums(!is.na(x_orig[[imp]])) == 0) #If needed to remove zero columns
      x_orig[[imp]] <- x_orig[[imp]][,colSums(!is.na(x_orig[[imp]])) > 0] #If needed to remove zero columns
      
      y_orig[[imp]] <- questions[[imp]]
      B <- 'Questions'
    }

    #PLS global brain - academic
    if(contrast==7){
      x_orig[[imp]] <- graph_measures
      A <- 'Global.connectome'
      y_orig[[imp]] <- academic[[imp]]
      B <- 'Academic'
    } 
    if(contrast==8){  #PLS global brain- academic
      x_orig[[imp]] <- graph_measures
      A <- 'Global.connectome'
      y_orig[[imp]] <- behaviour[[imp]]
      B <- 'Behaviour'
    }
    if(contrast==9){#PLS tracts- questionnaires
      x_orig[[imp]] <- tracts
      A <- 'Tract.FA'
      is_na <- which(colSums(!is.na(x_orig[[imp]])) == 0) #If needed to remove zero columns
      x_orig[[imp]] <- x_orig[[imp]][,colSums(!is.na(x_orig[[imp]])) > 0] #If needed to remove zero columns
      y_orig[[imp]] <- academic[[imp]]
      B <- 'Academic'
    }
    if(contrast==5){#PLS node_strength- academic
      x_orig[[imp]] <- node_strength
      A <- 'Node.strength'
      is_na <- which(colSums(!is.na(x_orig[[imp]])) == 0) #If needed to remove zero columns
      x_orig[[imp]] <- x_orig[[imp]][,colSums(!is.na(x_orig[[imp]])) > 0] #If needed to remove zero columns
      
      y_orig[[imp]] <- academic[[imp]]
      B <- 'Academic'
    }
    if(contrast==6){#PLS node_strength- behaviour
      x_orig[[imp]] <- node_strength
      A <- 'Node.strength'
      is_na <- which(colSums(!is.na(x_orig[[imp]])) == 0) #If needed to remove zero columns
      x_orig[[imp]] <- x_orig[[imp]][,colSums(!is.na(x_orig[[imp]])) > 0] #If needed to remove zero columns
      
      y_orig[[imp]] <- behaviour[[imp]]
      B <- 'Behaviour'
    }
    if(contrast==7){#PLS connection strength - academic
      x_orig[[imp]] <- half_nets
      A <- 'Connection.strength'
      is_na <- which(colSums(!is.na(x_orig[[imp]])) == 0) #If needed to remove zero columns
      x_orig[[imp]] <- x_orig[[imp]][,colSums(!is.na(x_orig[[imp]])) > 0] #If needed to remove zero columns
      
      y_orig[[imp]] <- academic[[imp]]
      B <- 'Academic_rotated'
    }
    if(contrast==8){#PLS connection strength - behaviour
      x_orig[[imp]] <- half_nets
      A <- 'Connection.strength'
      is_na <- which(colSums(!is.na(x_orig[[imp]])) == 0) #If needed to remove zero columns
      x_orig[[imp]] <- x_orig[[imp]][,colSums(!is.na(x_orig[[imp]])) > 0] #If needed to remove zero columns
      
      y_orig[[imp]] <- behaviour[[imp]]
      B <- 'Behaviour_rotated'
    }
   
    if(contrast==16){#PLS global- cognitive
      x_orig[[imp]] <- graph_measures
      A <- 'Global.connectome'
      y_orig[[imp]] <- cognitive[[imp]]
      B <- 'cognitive'
    } 
    if(contrast==17){#PLS node_strength- cognitive
      x_orig[[imp]] <- node_strength
      A <- 'Node.strength'
      is_na <- which(colSums(!is.na(x_orig[[imp]])) == 0) #If needed to remove zero columns
      x_orig[[imp]] <- x_orig[[imp]][,colSums(!is.na(x_orig[[imp]])) > 0] #If needed to remove zero columns
      
      y_orig[[imp]] <- cognitive[[imp]]
      B <- 'cognitive'
    }
    if(contrast==18){#PLS connection strength- cognitive
      x_orig[[imp]] <- half_nets
      A <- 'Connection.strength'
      is_na <- which(colSums(!is.na(x_orig[[imp]])) == 0) #If needed to remove zero columns
      x_orig[[imp]] <- x_orig[[imp]][,colSums(!is.na(x_orig[[imp]])) > 0] #If needed to remove zero columns
      
      y_orig[[imp]] <- cognitive[[imp]]
      B <- 'cognitive'
    }

    #fit PLS
    ncomp=min(cbind(ncol(x_orig[[imp]]),ncol(y_orig[[imp]])))
    orig_fit[[imp]] <- pls(x_orig[[imp]], y_orig[[imp]],  mode = "canonical", ncomp =ncomp, scale=TRUE, all.outputs=TRUE)
    #Note, don't use procrustes at this point as we're comparing singular values an procrustes reduces these
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
        #Note, don't use procrustes at this point as we're comparing singular values an procrustes reduces these
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
  }
  #Get pvals and CI for permutation
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
    z <- qnorm(pvalues[comp,])  # transform to z-scale
    num <- mean(z)
    den <- sqrt(1 + var(z))
    p_pooled[comp] <- pnorm( num / den) # average and transform back
  }
  av_pval <- cbind(c(1,2), data.frame(p_pooled))
  av_CI_upper <- G/5
  av_orig_cov_LV <- H/5
  av_orig_cor_LV <- J/5
  
  
  fname <- paste(wd, '/' , A, '_', B,'_pvals_CI_and_av_cov_LVs.csv', sep='' )
  df <- cbind(av_pval,av_CI_upper, av_orig_cov_LV, av_orig_cor_LV)
  colnames(df) <- c('comp', 'p','Av.CI', 'Av.cov.LVs', 'Av.cor.LVs')
  write.csv(df, file = fname)
  
  fname <- paste(wd, '/' , A, '_', B,'_pvals_each_imputation.csv', sep='' )
  df <- cbind(c(1,2), pvalues)
  colnames(df) <- c('comp', 'imp 1', 'imp 2', 'imp 3', 'imp 4', 'imp 5')
  write.csv(df, file = fname)
  
  #Save results 
  exp_var <- cbind(data.frame(round(orig_fit[[1]]$explained_variance$X,4)),data.frame(round(orig_fit[[1]]$explained_variance$Y,4)))
  for(imp in 2:5){
    exp_var <- cbind(exp_var, data.frame(round(orig_fit[[imp]]$explained_variance$X,4)),data.frame(round(orig_fit[[imp]]$explained_variance$Y,4)))
  }
  colnames_all <- rep(c(A,B),5)
  colnames(exp_var) <- colnames_all
  fname <- paste(wd, '/' , A, '_', B,'_exp_var.csv', sep='' )
  write.csv(exp_var, file = fname)
  
  flip <- sign(orig_fit[[1]]$loadings$Y[1,1])
  for(comp in 1:ncomp){
    variate <- cbind(data.frame(flip*orig_fit[[1]]$variates$X[,comp]),data.frame(flip*orig_fit[[1]]$variates$Y[,comp]))
    for(imp in 2:5){
      variate <- cbind(variate,data.frame(flip*orig_fit[[imp]]$variates$X[,comp]),data.frame(flip*orig_fit[[imp]]$variates$Y[,comp]))
    }
    colnames(variate) <- colnames_all
    fname <- paste(wd, '/' ,A, '_', B,'_LV_', comp, '.csv', sep='' )
    write.csv(variate, file = fname)
  }
  
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
brain_vis <- FALSE
for(contrast in c(1)){
  for(comp in 1:ncomp){
    
    for(imp in 1:5){
      #PLS global brain - academic
      if(contrast==1){
        x_orig[[imp]] <- questions[[imp]]
        A <- 'Questions'
        y_orig[[imp]] <- academic[[imp]]
        B <- 'Academic'
      } 
      if(contrast==2){  #PLS global brain- academic
        x_orig[[imp]] <- graph_measures
        A <- 'Global.connectome'
        y_orig[[imp]] <- behaviour[[imp]]
        B <- 'Behaviour'
      }
      if(contrast==3){#PLS tracts- questionnaires
        x_orig[[imp]] <- tracts
        A <- 'Tract.FA'
        y_orig[[imp]] <- academic[[imp]]
        B <- 'Academic'
      }
      if(contrast==4){  #PLS node degree- behaviour
        x_orig[[imp]] <- node_degree
        A <- 'Node.degree'
        y_orig[[imp]] <- behaviour[[imp]]
        B <- 'Behaviour'
        
      }
      if(contrast==5){#PLS node_strength- academic
        x_orig[[imp]] <- node_strength
        A <- 'Node.strength'
        y_orig[[imp]] <- academic[[imp]]
        B <- 'Academic'
      }
      if(contrast==6){#PLS node_strength- behaviour
        x_orig[[imp]] <- node_strength
        A <- 'Node.strength'
        y_orig[[imp]] <- behaviour[[imp]]
        B <- 'Behaviour'
      }
      if(contrast==7){#PLS connection strength - academic
        x_orig[[imp]] <- half_nets
        A <- 'Connection.strength'
        y_orig[[imp]] <- academic[[imp]]
        B <- 'Academic'
        nnodes <- 85
        labels_brain <- labels_85
        
        brain_vis <- TRUE
        
      }
      if(contrast==8){#PLS connection strength - behaviour
        x_orig[[imp]] <- half_nets[sample_for_perm[[1]],]
        A <- 'Connection.strength'
        y_orig[[imp]] <- behaviour[[imp]]
        B <- 'Behaviour'
      }
      if(contrast==9){  #PLS global brain- SDQ
        x_orig[[imp]] <- graph_measures
        A <- 'Global.connectome'
        y_orig[[imp]] <- behaviour[[imp]][,2]
        B <- 'SDQ.only'
      }
      if(contrast==10){  #PLS global brain- SDQ
        x_orig[[imp]] <- node_degree
        A <- 'Node.degree'
        y_orig[[imp]] <- behaviour[[imp]][,2]
        B <- 'SDQ.only'
      }
      if(contrast==11){  #PLS global brain- SDQ
        x_orig[[imp]] <- node_strength
        A <- 'Node.strength'
        y_orig[[imp]] <- behaviour[[imp]][,2]
        B <- 'SDQ.only'
      }
      if(contrast==12){#PLS global- questions
        x_orig[[imp]] <- graph_measures
        A <- 'Global.connectome'
        y_orig[[imp]] <- questions[[imp]]
        B <- 'Questions'
      } 
      if(contrast==13){#PLS node degree- questions
        x_orig[[imp]] <- node_degree
        A <- 'Node.degree'
        y_orig[[imp]] <- questions[[imp]]
        B <- 'Questions'
      }
      if(contrast==14){#PLS node_strength- questions
        x_orig[[imp]] <- node_strength
        A <- 'Node.strength'
        y_orig[[imp]] <- questions[[imp]]
        B <- 'Questions'
      }
      if(contrast==15){#PLS connection strength- questions
        x_orig[[imp]] <- half_nets
        A <- 'Connection.strength'
        y_orig[[imp]] <- questions[[imp]]
        B <- 'Questions'
        nnodes <- 85
        labels_brain <- labels_85
        
        brain_vis <- TRUE
      }
      if(contrast==16){#PLS global- cognitive
        x_orig[[imp]] <- graph_measures
        A <- 'Global.connectome'
        y_orig[[imp]] <- cognitive[[imp]]
        B <- 'cognitive'
      } 
      if(contrast==17){#PLS node_strength- cognitive
        x_orig[[imp]] <- node_strength
        A <- 'Node.strength'
        y_orig[[imp]] <- cognitive[[imp]]
        B <- 'cognitive'
      }
      if(contrast==18){#PLS connection strength- cognitive
        x_orig[[imp]] <- half_nets
        A <- 'Connection.strength'
        y_orig[[imp]] <- cognitive[[imp]]
        B <- 'cognitive'
        nnodes <- 85
        brain_vis <- TRUE
      }
      if(contrast==19){#PLS enviornment only (no cog) and outcome
        x_orig[[imp]] <- questions[[imp]]
        A <- 'Questions'
        y_orig[[imp]] <- cognitive[[imp]]
        B <- 'Cognitive'  
        
      }
      if(contrast==20){#PLS enviornment only (no cog) and outcome
        x_orig[[imp]] <- questions[[imp]]
        A <- 'Questions_no_Cog'
        y_orig[[imp]] <- behaviour[[imp]]
        B <- 'Behaviour'  
        
      }
      if(contrast==21){#PLS Cog and outcome
        x_orig[[imp]] <- cognitive[[imp]]
        A <- 'Cognitive_no_ques'
        y_orig[[imp]] <- academic[[imp]]
        B <- 'Academic'  
        
      }
      if(contrast==22){#PLS cog and outcome
        x_orig[[imp]] <- cognitive[[imp]]
        A <- 'Cognitive_no_ques'
        y_orig[[imp]] <- behaviour[[imp]]
        B <- 'Behaviour'
        
      }
      if(contrast==24){
        x_orig[[imp]] <- half_nets
        A <- 'Connection.strength'
        y_orig[[imp]] <- questions[[imp]][,1:3]
        B <- 'SES'
      }
      
      print(A)
      print(B)
      print(imp)
      data <- cbind(x_orig[[imp]], y_orig[[imp]])
      orig_fit[[imp]] <- pls(x_orig[[imp]], y_orig[[imp]],  mode = "canonical", ncomp =ncomp, scale=TRUE)
      X_boot <- function(d,j){
        d2 <- d[j,]
        rownames(d2)<- NULL
        perm_fit<- pls(d2[,1:ncol(x_orig[[imp]])], d2[,(ncol(x_orig[[imp]])+1):ncol(d2)],  mode = "canonical", ncomp =ncomp, scale=TRUE, all.outputs=TRUE)
        proc <- NULL
        proc$X <- procrustes(as.matrix(perm_fit$loadings$X[,comp]),as.matrix(orig_fit[[1]]$loadings$X[,comp]))
        return(proc$X$X.new)
      }
      X_boot_star <- function(d,j){
        d2 <- d[j,]
        rownames(d2)<- NULL
        perm_fit<- pls(d2[,1:ncol(x_orig[[imp]])], d2[,(ncol(x_orig[[imp]])+1):ncol(d2)],  mode = "canonical", ncomp =ncomp, scale=TRUE, all.outputs=TRUE)
        proc <- NULL
        proc$X <- procrustes(as.matrix(perm_fit$loadings.star[[1]][,comp]),as.matrix(orig_fit[[1]]$loadings.star[[1]][,comp]))
        return(proc$X$X.new)
      }
      X_boot_cor <- function(d,j){
        d2 <- d[j,]
        rownames(d2)<- NULL
        perm_fit<- pls(d2[,1:ncol(x_orig[[imp]])], d2[,(ncol(x_orig[[imp]])+1):ncol(d2)],  mode = "canonical", ncomp =ncomp, scale=TRUE, all.outputs=TRUE)
        proc <- NULL
        proc$X <- procrustes(as.matrix(perm_fit$loadings.star[[2]][,comp]),as.matrix(orig_fit[[1]]$loadings.star[[2]][,comp]))
        perm_fit$variates$Y[,comp] <- perm_fit$Y %*% perm_fit$loadings.star[[2]][,comp])
        corr_X_LVY <- cor(perm_fit$X, perm_fit$variates$Y[,comp])
        return(proc$X$X.new)
      }
      Y_boot <- function(d,j){
        d2 <- d[j,]
        rownames(d2)<- NULL
        perm_fit <- pls(d2[,1:ncol(x_orig[[imp]])], d2[,(ncol(x_orig[[imp]])+1):ncol(d2)],  mode = "canonical", ncomp =ncomp, scale=TRUE, all.outputs=TRUE)
        proc <- NULL
        proc$Y <- procrustes(as.matrix(perm_fit$loadings$Y[,comp]),as.matrix(orig_fit[[1]]$loadings$Y[,comp]))
        return(proc$Y$X.new)
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
    pooled_X_SE <- sqrt((av_X_SE)^2 + ((1+ (1/5))*(1/(5-1))*F)) #Remember var = SE^2
    pooled_Y_SE <- sqrt((av_Y_SE)^2 + ((1+ (1/5))*(1/(5-1))*G)) #Remember var = SE^2
    
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
    X_load <- data.frame(round(av_X_boot_result, 5),av_X_boot_result-pooled_X_SE, av_X_boot_result+pooled_X_SE, X_ratio_group,colnames(x_orig[[imp]]))
    colnames(X_load) <- c('loading', 'SE_low', 'SE_high','reliable', 'name')
    fname <- paste(wd,'/', A, '_', B, '_', comp, '_boot_X_loadings.csv', sep='' )
    write.csv(X_load, file = fname)
    
    #Save brain adjacency matrix
    if(brain_vis == TRUE){
      adj <- matrix(data=NA,nrow=nnodes,ncol=nnodes)
      brain_load <- X_load
      brain_load$reliable[brain_load$loading==0] <- 0
      brain_load$reliable[is.na(brain_load$loading)] <- 0
      
      # # Add in the variables that were removed #Not needed if full data (all variables, including NA) are used for bootstrap
      # na_loads <- data.frame(matrix(data=0, nrow=length(is_na), ncol=ncol(brain_load)))
      # colnames(na_loads) <- colnames(brain_load)
      # full_loads <- rbind(brain_load, na_loads)
      # full_loads_sort <- full_loads[colnames(x_orig[[1]]),]
      
      var <- data.frame(brain_load$loading*brain_load$reliable)
      adj[lower.tri(adj)] <- as.numeric(var[,1]) #Note, has to be lower so that it gets filled in in the right order
      adj <- data.frame(adj)
      adj[is.na(adj)] <- 0
      adj = adj + t(adj)
      rownames(adj) <- t(labels_brain)
      colnames(adj) <- t(labels_brain)
      fname <- paste(wd,'/', A, '_', B, '_', comp, '_adj_reliable_loading.csv', sep='' )
      write.csv(adj, file = fname, row.names=FALSE)
      adj_bin <- sign(adj)
      adj_bin[adj_bin==1] <- 2
      adj_bin[adj_bin==-1] <- 1
      fname <- paste(wd,'/', A, '_', B, '_', comp, '_adj_reliable_loading_bin.csv', sep='' )
      write.csv(adj_bin, file = fname, row.names=FALSE)
      adj_pos <- adj
      adj_pos[adj_pos<0] <- 0
      fname <- paste(wd,'/', A, '_', B, '_', comp, '_adj_reliable_loading_positive.csv', sep='' )
      write.csv(adj_pos, file = fname, row.names=FALSE)
      adj_neg <- adj
      adj_neg[adj_neg>0] <- 0
      adj_neg <- abs(adj_neg)
      
      fname <- paste(wd,'/', A, '_', B, '_', comp, '_adj_reliable_loading_negative.csv', sep='' )
      write.csv(adj_neg, file = fname, row.names=FALSE)
    }
    X_load <- X_load[order(-abs(X_load$loading)),]
    X_load$name <- factor(X_load$name, levels = X_load$name[order(-abs(X_load$loading))])
    X_load$reliable <- factor(X_load$reliable)
    
    ggplot(data = X_load,
           aes(x = name, y = loading,  fill=reliable)) +
      geom_bar(stat = 'identity', position = 'identity', width = 0.9) +
      #geom_errorbar(aes(ymin=-1*CI_low,ymax=-1*CI_high)) +
      geom_errorbar(aes(ymin=SE_low,ymax=SE_high), size=0.4) +
      scale_fill_manual(values=cols) +
      theme_minimal(base_size=10)+
      theme(axis.text.x = element_text(angle = 90, size=10, hjust=1,vjust=0.4), axis.title.x = element_blank()) +
      guides(fill=FALSE)
    
    fname <- paste(wd,'/', A, '_', B, '_', comp,'_X_loadings.png', sep='' )
    ggsave(fname, width = 15, height = 10, units = "cm")
    
    Y_load <- data.frame(round(av_Y_boot_result, 5),av_Y_boot_result-pooled_Y_SE, av_Y_boot_result+pooled_Y_SE, Y_ratio_group,colnames(y_orig[[imp]]))
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
      geom_errorbar(aes(ymin=SE_low,ymax=SE_high), size=0.4) +
      scale_fill_manual(values=cols) +
      theme_minimal(base_size=10)+
      theme(axis.text.x = element_text(angle = 90, size=10, hjust=1,vjust=0.4), axis.title.x = element_blank())+
      guides(fill=FALSE)
    
    fname <- paste(wd,'/', A, '_', B, '_', comp,'_Y_loadings.png', sep='' )
    ggsave(fname, width = 15, height = 10, units = "cm")
    
    X_load$table <- rep(A, nrow(X_load))
    Y_load$table <- rep(B, nrow(Y_load))
    
    if(A =='Connection.strength'){
      X_load <- X_load[1:120,]
    }
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
      geom_errorbar(aes(ymin=SE_low,ymax=SE_high), size=0.4) +
      scale_fill_manual(values=cols) +
      ylim(ymin, max(all_load$SE_high)) +
      theme_bw(base_size=10)+
      theme(axis.text.x = element_text(angle = 90, size=10, hjust=1,vjust=0.4), axis.title.x = element_blank()) +
      guides(fill=FALSE)+
      theme(strip.background = element_blank(),
            strip.text.x = element_blank())
    
    fname <- paste(wd,'/', A, '_', B, '_', comp,'_X_and_Y_loadings.png', sep='' )
    ggsave(fname, width = 15, height = 10, units = "cm")
    
  }
}




















#Scatter plots for global measures
#Get average scores across imputations...
academic_av <- NULL
behaviour_av <- NULL
for(imp in 1:5){
  if(imp==1){
    academic_av <- academic[[imp]]
    behaviour_av <- behaviour[[imp]]
  }
  if(imp>1){
    academic_av <- academic_av  + academic[[imp]]
    behaviour_av <- behaviour_av  + behaviour[[imp]]
  }
}
academic_av <- academic_av/5
behaviour_av <- behaviour_av/5

#Plot academic
df <- cbind(academic_av, graph_measures)
df.melted <- melt(df, id=c('WJ.Maths', 'WJ.Reading'),variable.name = "Graph.measure", 
                  value.name = "Graph.measure.value")
df.melted <- melt(df.melted, id=c('Graph.measure', 'Graph.measure.value'), variable.name = "Outcome", 
                  value.name = "Outcome.value")
ggplot(data = df.melted,
       aes(x = Graph.measure.value, y = Outcome.value, colour=Outcome)) +
  geom_point(show.legend =FALSE) +
  stat_smooth(method=lm, show.legend =FALSE) +
  facet_grid(Outcome~ Graph.measure)+
  theme_bw(base_size=15)+
  theme(axis.text.x = element_text(size=10, hjust=1), axis.text.y = element_text(size=10, hjust=1), axis.title.x=element_blank(), axis.title.y=element_blank())+
  theme(strip.text.x = element_text(size=10),
        strip.text.y = element_text(size=10),
        panel.border = element_rect(color = "gray", fill = NA, size = 1),
        strip.background = element_rect(colour="gray", size = 1,fill='white'), legend.key = element_blank()) 

fname <- paste(wd, '/academic_graph_measures_lr.png', sep='' )
ggsave(fname, width = 20, height = 10, units = "cm")

#Plot behaviour
df <- cbind(behaviour_av, graph_measures)
df.melted <- melt(df, id=c('BRIEF', 'SDQ'),variable.name = "Graph.measure", 
                  value.name = "Graph.measure.value")
df.melted <- melt(df.melted, id=c('Graph.measure', 'Graph.measure.value'), variable.name = "Outcome", 
                  value.name = "Outcome.value")
ggplot(data = df.melted,
       aes(x = Graph.measure.value, y = Outcome.value, colour=Outcome)) +
  geom_point(show.legend =FALSE) +
  stat_smooth(method=lm, show.legend =FALSE) +
  facet_grid(Outcome~ Graph.measure)+
  theme_bw(base_size=15)+
  theme(axis.text.x = element_text(size=10, hjust=1), axis.text.y = element_text(size=10, hjust=1), axis.title.x=element_blank(), axis.title.y=element_blank())+
  theme(strip.text.x = element_text(size=10),
        strip.text.y = element_text(size=10),
        panel.border = element_rect(color = "gray", fill = NA, size = 1),
        strip.background = element_rect(colour="gray", size = 1,fill='white'), legend.key = element_blank()) 

fname <- paste(wd, '/behaviour_graph_measures_lr.png', sep='' )
ggsave(fname, width = 20, height = 10, units = "cm") 
 

fit <- lm(BRIEF ~ ., data=df)
summary(fit)



#Plot loadings later
Y_load <- read.csv("U:/Documents/ACE_Data/MI_FA/MI_connectome/sparse/Connection_strength/s_Connection.strength_Questions_2_boot_y_loadings.csv")
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
adj <- matrix(data=NA,nrow=85,ncol=85)
var <- for_brain_plot
adj[lower.tri(adj)] <- as.numeric(var) #Note, has to be lower so that it gets filled in in the right order
adj[is.na(adj)] <- 0

adj = adj + t(adj)
adj <- data.frame(adj)
rownames(adj) <- t(labels_85) 
colnames(adj) <- t(labels_85)
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



#Setup
library(mixOmics)
library(psych)
library(Matrix)
library(gtools)
library(ggplot2)
library(MCMCpack)
library(boot)


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

wd <- "U:/Documents/ACE_Data/Thesis_Analysis/MI_PLS_MRI/results_non_sparse/Questions-MRI-Outcomes"
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

x_orig <- NULL
y_orig <- NULL
m_orig <- NULL
data <- NULL
ncomp <- 2
comp <- 1
for(contrast in c(2)){
  for(imp in 1:5){
    #PLS brain- academic
    if(contrast==1){
      x_orig[[imp]] <- questions[[imp]][,1:3]
      A <- 'SES'
      m_orig[[imp]] <- questions[[imp]][,c(4:6,8:22)]
      B <- 'Questions_excl_hours'
      y_orig[[imp]] <- academic[[imp]]
      C <- 'Academic'

    } 
    if(contrast==2){  #PLS questionnaires- behaviour
      x_orig[[imp]] <- questions[[imp]][,1:3]
      A <- 'SES'
      m_orig[[imp]] <- questions[[imp]][,c(4:6,8:22)]
      B <- 'Questions_excl_hours'
      y_orig[[imp]] <- behaviour[[imp]]
      C <- 'Behaviour'

    }
    if(contrast==3){  #PLS questionnaires- cog- academic
      x_orig[[imp]] <- questions[[imp]]
      A <- 'Questions'
      m_orig[[imp]] <- cognitive[[imp]]
      B <- 'Cognitive'
      y_orig[[imp]] <- academic[[imp]]
      C <- 'Academic'

    } 
   if(contrast==4){  #PLS questionnaires- cog- behaviour
     x_orig[[imp]] <- questions[[imp]]
     A <- 'Questions'
     m_orig[[imp]] <- cognitive[[imp]]
     B <- 'Cognitive'
     y_orig[[imp]] <- behaviour[[imp]]
     C <- 'Behaviour'

   } 

    if(contrast==5){  #PLS cog- questionnaires - behaviour
      x_orig[[imp]] <- cognitive[[imp]]
      A <- 'Cognitive'
      m_orig[[imp]] <- questions[[imp]]
      B <- 'Questions'
      y_orig[[imp]] <- behaviour[[imp]]
      C <- 'Behaviour'
    } 
    
    ##### MRI #####
    if(contrast==10){#PLS questions- global- academic
      x_orig[[imp]] <- questions[[imp]]
      A <- 'Questions'
      m_orig[[imp]] <- graph_measures
      B <- 'Global.connectome'
      y_orig[[imp]] <- academic[[imp]]
      C <- 'Academic'

    } 
    if(contrast==11){#PLS questions- node degree- academic
      x_orig[[imp]] <- questions[[imp]]
      A <- 'Questions'
      m_orig[[imp]] <- node_degree
      B <- 'Node.degree'
      y_orig[[imp]] <- academic[[imp]]
      C <- 'Academic'
      
    }
    if(contrast==12){#PLS questions- node_strength- academic
      x_orig[[imp]] <- questions[[imp]]
      A <- 'Questions'
      m_orig[[imp]] <- node_strength
      B <- 'Node.strength'
      y_orig[[imp]] <- academic[[imp]]
      C <- 'Academic'
    }
    if(contrast==13){#PLS questions- connection strength- academic
      x_orig[[imp]] <- questions[[imp]]
      A <- 'Questions'
      m_orig[[imp]] <- half_nets
      B <- 'Connection.strength'
      nnodes <- 85
      labels_brain <- labels_85
      brain_vis <- TRUE
      y_orig[[imp]] <- academic[[imp]]
      C <- 'Academic'
    }
    if(contrast==14){#PLS questions- global- behaviour
      x_orig[[imp]] <- questions[[imp]]
      A <- 'Questions'
      m_orig[[imp]] <- graph_measures
      B <- 'Global.connectome'
      y_orig[[imp]] <- behaviour[[imp]]
      C <- 'Behaviour'
    } 
    if(contrast==15){#PLS questions- node degree- behaviour
      x_orig[[imp]] <- questions[[imp]]
      A <- 'Questions'
      m_orig[[imp]] <- node_degree
      B <- 'Node.degree'
      y_orig[[imp]] <- behaviour[[imp]]
      C <- 'Behaviour'
    }
    if(contrast==16){#PLS questions- node_strength- behaviour
      x_orig[[imp]] <- questions[[imp]]
      A <- 'Questions'
      m_orig[[imp]] <- node_strength
      B <- 'Node.strength'
      y_orig[[imp]] <- behaviour[[imp]]
      C <- 'Behaviour'
    }
    if(contrast==17){#PLS questions- connection strength- Behaviour
      x_orig[[imp]] <- questions[[imp]]
      A <- 'Questions'
      m_orig[[imp]] <- half_nets
      B <- 'Connection.strength'
      nnodes <- 85
      labels_brain <- labels_85
      brain_vis <- TRUE
      y_orig[[imp]] <- behaviour[[imp]]
      C <- 'Behaviour'
    }
    data[[imp]] = list(X=x_orig[[imp]] , M=m_orig[[imp]], Y = as.matrix(y_orig[[imp]]))
    
  }
  print(paste(A,B,C))
  
  #Initialise variables
  orig_fit <- NULL
  orig_cov_LV <- NULL
  orig_cor_LV <- NULL
  orig_ab <- NULL
  a <- b<- c<- c_med <- NULL
  #Fit the actual PLS model
  for(imp in 1:5){
    orig_fit[[imp]] <- block.pls(data[[imp]], indY = 3,mode = "canonical",scheme='horst', ncomp =ncomp, scale=TRUE, max.iter = 1000)
    
    orig_cor_LV[[imp]] <- diag(cor(orig_fit[[imp]]$variates$X, orig_fit[[imp]]$variates$Y)) #Correlation between LV 
    #Standardise the LVs and run linear models to calculate the indirect path coefficient
    orig_fit[[imp]]$variates$X <- scale(orig_fit[[imp]]$variates$X)
    orig_fit[[imp]]$variates$M <- scale(orig_fit[[imp]]$variates$M)
    orig_fit[[imp]]$variates$Y <- scale(orig_fit[[imp]]$variates$Y)
    
    a[[imp]]<-  cbind(lm(orig_fit[[imp]]$variates$M[,1] ~ orig_fit[[imp]]$variates$X[,1])$coefficients[2], lm(orig_fit[[imp]]$variates$M[,2] ~ orig_fit[[imp]]$variates$X[,2])$coefficients[2])
    b[[imp]]<-  cbind(lm(orig_fit[[imp]]$variates$Y[,1] ~ orig_fit[[imp]]$variates$X[,1]+orig_fit[[imp]]$variates$M[,1])$coefficients[3], lm(orig_fit[[imp]]$variates$Y[,2] ~ orig_fit[[imp]]$variates$X[,2]+orig_fit[[imp]]$variates$M[,2])$coefficients[3])
    c[[imp]]<-  cbind(lm(orig_fit[[imp]]$variates$Y[,1] ~ orig_fit[[imp]]$variates$X[,1])$coefficients[2], lm(orig_fit[[imp]]$variates$Y[,2] ~ orig_fit[[imp]]$variates$X[,2])$coefficients[2])
    c_med[[imp]]<-  cbind(lm(orig_fit[[imp]]$variates$Y[,1] ~ orig_fit[[imp]]$variates$X[,1]+orig_fit[[imp]]$variates$M[,1])$coefficients[2], lm(orig_fit[[imp]]$variates$Y[,2] ~ orig_fit[[imp]]$variates$X[,2]+orig_fit[[imp]]$variates$M[,2])$coefficients[2])
    orig_ab[[imp]]<- a[[imp]]*b[[imp]]
    
  }
  
  #Save results 
  exp_var <- cbind(data.frame(round(orig_fit[[1]]$explained_variance$X,4)),data.frame(round(orig_fit[[1]]$explained_variance$M,4)),data.frame(round(orig_fit[[1]]$explained_variance$X,4)))
  First_LV <- cbind(data.frame(orig_fit[[1]]$variates$X[,1]),data.frame(orig_fit[[1]]$variates$M[,1]), data.frame(orig_fit[[1]]$variates$Y[,1]))
  X_1st_loadings <- data.frame(orig_fit[[1]]$loadings$X[,1])
  M_1st_loadings <- data.frame(orig_fit[[1]]$loadings$M[,1])
  Y_1st_loadings <- data.frame(orig_fit[[1]]$loadings$Y[,1])
  path_model <- rbind(c[[1]], c_med[[1]], a[[1]], b[[1]], orig_ab[[1]])
  rownames(path_model) <- c(paste(A, '->', C ),paste(A, '->', C, 'with_mediator'),paste(A, '->', B ),paste(B, '->', C ),paste(A, '->', B,'->', C  ))
  
  for(imp in 2:5){
    exp_var <- cbind(exp_var, data.frame(round(orig_fit[[imp]]$explained_variance$X,4)),data.frame(round(orig_fit[[imp]]$explained_variance$M,4)),data.frame(round(orig_fit[[imp]]$explained_variance$Y,4)))
    First_LV <- cbind(First_LV,data.frame(orig_fit[[imp]]$variates$X[,1]),data.frame(orig_fit[[1]]$variates$M[,1]),data.frame(orig_fit[[imp]]$variates$Y[,1]))
    X_1st_loadings <- cbind(X_1st_loadings,orig_fit[[imp]]$loadings$X[,1])
    M_1st_loadings <- cbind(M_1st_loadings,orig_fit[[imp]]$loadings$M[,1])
    Y_1st_loadings <- cbind(Y_1st_loadings,orig_fit[[imp]]$loadings$Y[,1])
    path_model <- cbind(path_model,rbind(c[[imp]], c_med[[imp]], a[[imp]], b[[imp]], orig_ab[[imp]]))
    
  }
  colnames(exp_var) <- rep(c(A,C),5)
  colnames(First_LV) <- rep(c(A,B,C),5)
  colnames(X_1st_loadings) <- rep(A,5)
  colnames(M_1st_loadings) <- rep(B,5)
  colnames(Y_1st_loadings) <- rep(C,5)
  colnames(path_model) <- rep(c('comp1','comp2'),5)
  
  #Save tables
  fname <- paste(wd, '/' , A, '_', B,'_', C,'_exp_var.csv', sep='' )
  write.csv(exp_var, file = fname)
  fname <- paste(wd, '/' ,A, '_', B,'_', C,'_1st_LV.csv', sep='' )
  write.csv(First_LV, file = fname)
  fname <- paste(wd, '/' ,A, '_', B,'_', C,  '_X_1st_loadings.csv', sep='' )
  write.csv(X_1st_loadings, file = fname)
  fname <- paste(wd, '/' , A, '_', B,'_', C, '_Y_1st_loadings.csv', sep='' )
  write.csv(Y_1st_loadings, file = fname)
  fname <- paste(wd, '/' , A, '_', B,'_', C, '_M_1st_loadings.csv', sep='' )
  write.csv(M_1st_loadings, file = fname)
  fname <- paste(wd, '/' , A, '_', B,'_', C, '_path_model.csv', sep='' )
  write.csv(path_model, file = fname)
  
  #Permuted PLS to get significance of each component
  nperm = 1000
  permuted_cov_LV <- list()
  permuted_cor_LV <- list()
  permuted_cov_LV_proc <- list()
  permuted_cor_LV_proc <- list()
  permuted_LV <- list()
  perm_fit <- NULL
  sample_for_perm <- NULL
  perm_ab <- NULL
  pval <- list()
  n_CI <- round((0.05*(nperm+1)) - 1,0) #maximum nunber of permuted SV's that can be greater than origional and it still be significant
  CI_upper <- list()
  for(imp in 1:5){
    for (i in 1:nperm){
      if(imp ==1){
        sample_for_perm[[i]] <- sample(nrow(m_orig[[imp]])) #to enseure same sample is pulled out for each imputed set
      }
      m_perm <- m_orig[[imp]][sample_for_perm[[i]],] #Permuted data
      data_perm = list(X=x_orig[[imp]] , M=m_perm, Y = as.matrix(y_orig[[imp]]))
      perm_fit[[imp]] <- block.pls(data_perm, indY = 3,mode = "canonical",scheme='horst', ncomp =ncomp, scale=TRUE, max.iter = 1000)
      
      
      if(i ==1){
        permuted_cov_LV[[imp]] <- diag(cov(perm_fit[[imp]]$variates$X, perm_fit[[imp]]$variates$Y))
        permuted_cor_LV[[imp]] <- diag(cor(perm_fit[[imp]]$variates$X, perm_fit[[imp]]$variates$Y))
        perm_fit[[imp]]$variates$X <- scale(perm_fit[[imp]]$variates$X)
        perm_fit[[imp]]$variates$M <- scale(perm_fit[[imp]]$variates$M)
        perm_fit[[imp]]$variates$Y <- scale(perm_fit[[imp]]$variates$Y)
        a[[imp]]<-  cbind(lm(perm_fit[[imp]]$variates$M[,1] ~ perm_fit[[imp]]$variates$X[,1])$coefficients[2], lm(perm_fit[[imp]]$variates$M[,2] ~ perm_fit[[imp]]$variates$X[,2])$coefficients[2])
        b[[imp]]<-  cbind(lm(perm_fit[[imp]]$variates$Y[,1] ~ perm_fit[[imp]]$variates$X[,1]+perm_fit[[imp]]$variates$M[,1])$coefficients[3], lm(perm_fit[[imp]]$variates$Y[,2] ~ perm_fit[[imp]]$variates$X[,2]+perm_fit[[imp]]$variates$M[,2])$coefficients[3])
        perm_ab[[imp]]<- data.frame(a[[imp]]*b[[imp]])
      }
      else{
        permuted_cov_LV[[imp]] <- rbind(permuted_cov_LV[[imp]],diag(cov(perm_fit[[imp]]$variates$X, perm_fit[[imp]]$variates$Y)))
        permuted_cor_LV[[imp]] <- rbind(permuted_cor_LV[[imp]],diag(cor(perm_fit[[imp]]$variates$X, perm_fit[[imp]]$variates$Y)))
        perm_fit[[imp]]$variates$X <- scale(perm_fit[[imp]]$variates$X)
        perm_fit[[imp]]$variates$M <- scale(perm_fit[[imp]]$variates$M)
        perm_fit[[imp]]$variates$Y <- scale(perm_fit[[imp]]$variates$Y)
        a[[imp]]<-  cbind(lm(perm_fit[[imp]]$variates$M[,1] ~ perm_fit[[imp]]$variates$X[,1])$coefficients[2], lm(perm_fit[[imp]]$variates$M[,2] ~ perm_fit[[imp]]$variates$X[,2])$coefficients[2])
        b[[imp]]<-  cbind(lm(perm_fit[[imp]]$variates$Y[,1] ~ perm_fit[[imp]]$variates$X[,1]+perm_fit[[imp]]$variates$M[,1])$coefficients[3], lm(perm_fit[[imp]]$variates$Y[,2] ~ perm_fit[[imp]]$variates$X[,2]+perm_fit[[imp]]$variates$M[,2])$coefficients[3])
        perm_ab[[imp]]<- rbind(perm_ab[[imp]], data.frame(a[[imp]]*b[[imp]]))
      }
    }
    
    for(comp in 1:ncomp){
      if(comp ==1){
        pval[[imp]] <- (sum(perm_ab[[imp]][,comp]>orig_ab[[imp]][comp])+1)/(nperm + 1)
        CI_upper[[imp]] <- data.frame(sort(perm_ab[[imp]][,comp], decreasing = TRUE))[n_CI,1]
      }
      else {
        pval[[imp]] <- rbind(pval[[imp]],(sum(perm_ab[[imp]][,comp]>orig_ab[[imp]][comp])+1)/(nperm + 1))
        CI_upper[[imp]] <- rbind(CI_upper[[imp]],data.frame(sort(perm_ab[[imp]][,comp], decreasing = TRUE))[n_CI,1])
      }
    }
    pval[[imp]] <- round(pval[[imp]], 5)
    pval[[imp]] <- data.frame(cbind(1:ncomp, pval[[imp]]))
  }  
  #Get pvals and CI for permutation
  pvalues <- NULL
  G <- NULL
  H <- NULL
  for(imp in 1:5){
    if(imp==1){
      pvalues <- pval[[imp]][,2]
      G <- CI_upper[[imp]]
      H <- t(orig_ab[[imp]])
      
    }
    if(imp>1){
      pvalues <-cbind(pvalues,pval[[imp]][,2])
      G <- G + CI_upper[[imp]]
      H <- H + t(orig_ab[[imp]])
    }
  }
  #Use Rubin's Z tranform method to pool p values
  p_pooled <- NULL 
  for(comp in 1:ncomp){
    z <- qnorm(pvalues[comp,])  # transform to z-scale
    zmean <- mean(z)
    imp_var <- sum((zmean-z)^2)/(5-1)
    total_var <- 1 + (1 + (1/5))*imp_var
    p_pooled[comp] <- 1-pnorm( zmean / sqrt(total_var)) # average and transform back
  }
  av_pval <- cbind(c(1,2), data.frame(p_pooled))
  av_CI_upper <- G/5
  av_orig_ab <- H/5
  colnames(av_orig_ab) <- 'av_orig_ab'
  
  fname <- paste(wd, '/' , A, '_', B,'_', C,'_pvals_dont_use.csv', sep='' )
  df <- cbind(av_pval,av_CI_upper, av_orig_ab)
  colnames(df) <- c('comp', 'p','Av.CI', 'Av.ab.path')
  write.csv(df, file = fname)
  
  fname <- paste(wd, '/' , A, '_', B,'_', C,'_pvals_each_imputation_dont_use.csv', sep='' )
  df <- cbind(c(1:ncomp), pvalues)
  colnames(df) <- c('comp', 'imp 1', 'imp 2', 'imp 3', 'imp 4', 'imp 5')
  write.csv(df, file = fname)

  #Get rule of thumb reliablility for loadings
  #This is based on the t-distribution ( t stat >2 is significant for two-tailed alpha=0.05 for this number of participants)
  av_pval_group <- data.frame(rep(1,nrow(av_pval))) 
  av_pval_group[which(av_pval[,2] > 0.05),] <- 0
  av_pval <- cbind(av_pval, av_pval_group)
  colnames(av_pval) <- c('X1', 'X2', 'significant')
  av_pval$significant <- as.factor(av_pval$significant)
  
  ggplot(av_pval, aes(X1,X2, fill=significant)) + 
    geom_bar(stat="identity") + 
    xlab("Component") + 
    ylab("P Value") +
    scale_x_continuous(breaks = seq(1, ncomp, by = 1))+
    scale_fill_manual(values=c("#999999", "#ec6c20")) +
    theme_minimal(base_size=15)+
    theme(axis.text.x = element_text(size=10, hjust=1), axis.text.y = element_text(size=10, hjust=1))+
    guides(fill=FALSE)
  fname <- paste(wd, '/' , A, '_', B,'_', C,'_pvals.png', sep='' )
  ggsave(fname, width = 15, height = 10, units = "cm")
  #Plot components with CI
  CI_lower = rep(0, ncomp)
  ggplot(data.frame(av_orig_ab), aes(av_orig_ab)) + 
    geom_line(aes(x=1:ncomp,y=av_orig_ab), colour="#ec6c20") + 
    geom_ribbon(aes(x=1:ncomp,ymin=CI_lower, ymax=av_CI_upper), alpha=0.2)  + 
    xlab("Component") + 
    ylab("Covariance between latent variables")+
    scale_x_continuous(breaks = seq(1, ncomp, by = 1)) +
    theme_minimal(base_size=15)+
    theme(axis.text.x = element_text(size=10, hjust=1), axis.text.y = element_text(size=10, hjust=1))
  fname <- paste(wd, '/' , A, '_', B,'_', C,'_SVs.png', sep='' )
  ggsave(fname, width = 15, height = 10, units = "cm")
}


###########Bootstrap #######################
ncomp = 1
nboot = 500
orig_fit <-perm_fit <-  NULL
X_boot_result <- M_boot_result <- Y_boot_result <-X_boot_star_result <- M_boot_star_result <- Y_boot_star_result <- ab_boot_result <- NULL
X_SE <- M_SE <- Y_SE <- X_star_SE <- M_star_SE <- Y_star_SE <- ab_SE <- NULL
x_orig <- y_orig <- m_orig <- data <- m_orig_all <- NULL
cols <- c('0' = '#999999','1'='#ec6c20') #set colours so that sig are orange
do_star <- FALSE
for(contrast in c(22:24)){
  for(imp in 1:5){
    if(contrast==1){
      x_orig[[imp]] <- questions[[imp]][,1:3]
      A <- 'SES'
      m_orig[[imp]] <- questions[[imp]][,c(4:6,8:22)]
      B <- 'Questions_excl_hours'
      y_orig[[imp]] <- cognitive[[imp]]
      C <- 'Cognitive'
      brain_vis <- FALSE
      
      
    } 
    if(contrast==2){  #PLS questionnaires- behaviour
      x_orig[[imp]] <- questions[[imp]][,1:3]
      A <- 'SES'
      m_orig[[imp]] <- questions[[imp]][,c(4:6,8:22)]
      B <- 'Questions_excl_hours'
      y_orig[[imp]] <- behaviour[[imp]]
      C <- 'Behaviour'
      brain_vis <- FALSE
      
    }
    if(contrast==3){  #PLS questionnaires- academic
      x_orig[[imp]] <- questions[[imp]]
      A <- 'Questions'
      m_orig[[imp]] <- cognitive[[imp]]
      B <- 'Cognitive'
      y_orig[[imp]] <- academic[[imp]]
      C <- 'Academic'
      brain_vis <- FALSE
      
    } 
    if(contrast==4){  #PLS questionnaires- academic
      x_orig[[imp]] <- questions[[imp]][,c(12,20,17,2,11)]
      A <- 'Questions'
      m_orig[[imp]] <- cognitive[[imp]]
      B <- 'Cognitive'
      y_orig[[imp]] <- behaviour[[imp]]
      C <- 'Behaviour'
      brain_vis <- FALSE
      
    } 
    if(contrast==5){  #PLS questionnaires- academic
      x_orig[[imp]] <- cognitive[[imp]]
      A <- 'Cognitive'
      m_orig[[imp]] <- questions[[imp]]
      B <- 'Questions'
      y_orig[[imp]] <- behaviour[[imp]]
      C <- 'Behaviour'
      brain_vis <- FALSE
      
    } 
    
    
    ##### MRI #####
    if(contrast==10){#PLS questions- global- academic
      x_orig[[imp]] <- questions[[imp]]
      A <- 'Questions'
      m_orig[[imp]] <- graph_measures
      B <- 'Global.connectome'
      y_orig[[imp]] <- academic[[imp]]
      C <- 'Academic'
      brain_vis <- FALSE
      
    } 
    if(contrast==11){#PLS questions- node degree- academic
      x_orig[[imp]] <- questions[[imp]]
      A <- 'Questions'
      m_orig[[imp]] <- node_degree
      B <- 'Node.degree'
      y_orig[[imp]] <- academic[[imp]]
      C <- 'Academic'
      brain_vis <- FALSE
      
    }
    if(contrast==12){#PLS questions- node_strength- academic
      x_orig[[imp]] <- questions[[imp]]
      A <- 'Questions'
      m_orig[[imp]] <- node_strength
      B <- 'Node.strength'
      y_orig[[imp]] <- academic[[imp]]
      C <- 'Academic'
      brain_vis <- FALSE
      
    }
    if(contrast==13){#PLS questions- connection strength- academic
      x_orig[[imp]] <- questions[[imp]]
      A <- 'Questions'
      m_orig[[imp]] <- half_nets
      B <- 'Connection.strength'
      nnodes <- 85
      labels_brain <- labels_85
      brain_vis <- TRUE
      y_orig[[imp]] <- academic[[imp]]
      C <- 'Academic'

    }
    if(contrast==14){#PLS questions- global- behaviour
      x_orig[[imp]] <- questions[[imp]]
      A <- 'Questions'
      m_orig[[imp]] <- graph_measures
      B <- 'Global.connectome'
      y_orig[[imp]] <- behaviour[[imp]]
      C <- 'Behaviour'
      brain_vis <- FALSE
      
    } 
    if(contrast==15){#PLS questions- node degree- behaviour
      x_orig[[imp]] <- questions[[imp]]
      A <- 'Questions'
      m_orig[[imp]] <- node_degree
      B <- 'Node.degree'
      y_orig[[imp]] <- behaviour[[imp]]
      C <- 'Behaviour'
      brain_vis <- FALSE
      
    }
    if(contrast==16){#PLS questions- node_strength- behaviour
      x_orig[[imp]] <- questions[[imp]]
      A <- 'Questions'
      m_orig[[imp]] <- node_strength
      B <- 'Node.strength'
      y_orig[[imp]] <- behaviour[[imp]]
      C <- 'Behaviour'
      brain_vis <- FALSE
      
    }
    if(contrast==17){#PLS questions- connection strength- Behaviour
      x_orig[[imp]] <- questions[[imp]]
      A <- 'Questions'
      m_orig[[imp]] <- half_nets
      B <- 'Connection.strength'
      nnodes <- 85
      labels_brain <- labels_85
      brain_vis <- TRUE
      y_orig[[imp]] <- behaviour[[imp]]
      C <- 'Behaviour'
    }
    if(contrast==18){
      x_orig[[imp]] <- questions[[imp]][,1:3]
      A <- 'SES'
      m_orig[[imp]] <- questions[[imp]][,c(4:6,8:22)]
      B <- 'Questions_excl_hours'
      y_orig[[imp]] <- graph_measures
      C <- 'Global.connectome'
      brain_vis <- FALSE
    } 
    if(contrast==19){
      x_orig[[imp]] <- questions[[imp]][,1:3]
      A <- 'SES'
      m_orig[[imp]] <- questions[[imp]][,c(4:6,8:22)]
      B <- 'Questions_excl_hours'
      y_orig[[imp]] <- node_degree
      C <- 'Node.degree'
      brain_vis <- FALSE
    } 
    if(contrast==20){
      x_orig[[imp]] <- questions[[imp]][,1:3]
      A <- 'SES'
      m_orig[[imp]] <- questions[[imp]][,c(4:6,8:22)]
      B <- 'Questions_excl_hours'
      y_orig[[imp]] <- node_strength
      C <- 'Node.strength'
      brain_vis <- FALSE
    } 
    if(contrast==21){
      x_orig[[imp]] <- questions[[imp]][,1:3]
      A <- 'SES'
      m_orig[[imp]] <- questions[[imp]][,c(4:6,8:22)]
      B <- 'Questions_excl_hours'
      y_orig[[imp]] <- half_nets
      C <- 'Connection.strength'
      nnodes <- 85
      labels_brain <- labels_85
      brain_vis <- TRUE
    } 
    
    ### MEG analysis ###
    if(contrast==22){#PLS connection strength - academic
      x_orig[[imp]] <- questions[[imp]][,1:3]
      A <- 'SES'
      m_orig[[imp]] <- questions[[imp]][,c(4:6,8:22)]
      B <- 'Questions_excl_hours'
      y_orig[[imp]] <- half_nets
      C <- paste(freq_range, '_', corr_type, sep='')
      nnodes <- 68
      labels_brain <- labels_68
      brain_vis <- TRUE
    }
    if(contrast==23){#PLS connection strength - academic
      x_orig[[imp]] <- questions[[imp]]
      A <- 'Questions'
      m_orig[[imp]] <- half_nets
      B <- paste(freq_range, '_', corr_type, sep='')
      y_orig[[imp]] <- academic[[imp]]
      C <- 'Academic'
      nnodes <- 68
      labels_brain <- labels_68
      brain_vis <- TRUE
    }
    if(contrast==24){#PLS connection strength - academic
      x_orig[[imp]] <- questions[[imp]]
      A <- 'Questions'
      m_orig[[imp]] <- half_nets
      B <- paste(freq_range, '_', corr_type, sep='')
      y_orig[[imp]] <- behaviour[[imp]]
      C <- 'Behaviour'
      nnodes <- 68
      labels_brain <- labels_68
      brain_vis <- TRUE
    }
    
    data[[imp]] = list(X=x_orig[[imp]] , M=m_orig[[imp]], Y = as.matrix(y_orig[[imp]]))
    
  }
      
  #### Run Block PLS ####  
  ncomp = min(ncol(x_orig[[imp]]),ncol(y_orig[[imp]])) #Number components to extract
  for(comp in 1:ncomp){
    for(imp in 1:5){
      print(paste(A,B,C))
      
      print(imp)
      boot_data <- cbind(x_orig[[imp]], m_orig[[imp]],y_orig[[imp]])
      orig_fit[[imp]] <- block.pls(data[[imp]], indY = 3,mode = "canonical",scheme='horst', ncomp =ncomp, scale=TRUE, max.iter = 1000)
      X_boot <- function(d,j){
        d2 <- d[j,]
        rownames(d2)<- NULL
        data[[imp]] = list(X=d2[,1:ncol(x_orig[[imp]])] , M=d2[,(ncol(x_orig[[imp]])+1):(ncol(x_orig[[imp]])+ncol(m_orig[[imp]]))], Y = d2[,((ncol(x_orig[[imp]])+ncol(m_orig[[imp]]))+1):ncol(d2)])
        perm_fit <- block.pls(data[[imp]], indY = 3,mode = "canonical",scheme='horst', ncomp =ncomp, scale=TRUE, max.iter = 1000)
        proc <- NULL
        proc$X <- procrustes(as.matrix(perm_fit$loadings$X[,comp]),as.matrix(orig_fit[[1]]$loadings$X[,comp]))
        return(proc$X$X.new)
      }
      X_boot_star <- function(d,j){
        d2 <- d[j,]
        rownames(d2)<- NULL
        data[[imp]] = list(X=d2[,1:ncol(x_orig[[imp]])] , M=d2[,(ncol(x_orig[[imp]])+1):(ncol(x_orig[[imp]])+ncol(m_orig[[imp]]))], Y = d2[,((ncol(x_orig[[imp]])+ncol(m_orig[[imp]]))+1):ncol(d2)])
        perm_fit <- block.pls(data[[imp]], indY = 3,mode = "canonical",scheme='horst', ncomp =ncomp, scale=TRUE, max.iter = 1000)
        sign_flip <- procrustes(as.matrix(perm_fit$loadings$X[,comp]),as.matrix(orig_fit[[1]]$loadings$X[,comp]))$R
        loadings.star <- ginv(perm_fit$X$X)%*%(sign_flip*perm_fit$variates$X[,comp])
        return(loadings.star)
      }
      M_boot <- function(d,j){
        d2 <- d[j,]
        rownames(d2)<- NULL
        data[[imp]] = list(X=d2[,1:ncol(x_orig[[imp]])] , M=d2[,(ncol(x_orig[[imp]])+1):(ncol(x_orig[[imp]])+ncol(m_orig[[imp]]))], Y = d2[,((ncol(x_orig[[imp]])+ncol(m_orig[[imp]]))+1):ncol(d2)])
        perm_fit <- block.pls(data[[imp]], indY = 3,mode = "canonical",scheme='horst', ncomp =ncomp, scale=TRUE, max.iter = 1000)
        proc <- NULL
        proc$M <- procrustes(as.matrix(perm_fit$loadings$M[,comp]),as.matrix(orig_fit[[1]]$loadings$M[,comp]))
        return(proc$M$X.new)
      }
      M_boot_star <- function(d,j){
        d2 <- d[j,]
        rownames(d2)<- NULL
        data[[imp]] = list(X=d2[,1:ncol(x_orig[[imp]])] , M=d2[,(ncol(x_orig[[imp]])+1):(ncol(x_orig[[imp]])+ncol(m_orig[[imp]]))], Y = d2[,((ncol(x_orig[[imp]])+ncol(m_orig[[imp]]))+1):ncol(d2)])
        perm_fit <- block.pls(data[[imp]], indY = 3,mode = "canonical",scheme='horst', ncomp =ncomp, scale=TRUE, max.iter = 1000)
        sign_flip <- procrustes(as.matrix(perm_fit$loadings$M[,comp]),as.matrix(orig_fit[[1]]$loadings$M[,comp]))$R
        loadings.star <- ginv(perm_fit$X$M)%*%(sign_flip*perm_fit$variates$M[,comp])
        return(loadings.star)
      }
      Y_boot <- function(d,j){
        d2 <- d[j,]
        rownames(d2)<- NULL
        data[[imp]] = list(X=d2[,1:ncol(x_orig[[imp]])] , M=d2[,(ncol(x_orig[[imp]])+1):(ncol(x_orig[[imp]])+ncol(m_orig[[imp]]))], Y = d2[,((ncol(x_orig[[imp]])+ncol(m_orig[[imp]]))+1):ncol(d2)])
        perm_fit <- block.pls(data[[imp]], indY = 3,mode = "canonical",scheme='horst', ncomp =ncomp, scale=TRUE, max.iter = 1000)
        proc <- NULL
        proc$Y <- procrustes(as.matrix(perm_fit$loadings$Y[,comp]),as.matrix(orig_fit[[1]]$loadings$Y[,comp]))
        return(proc$Y$X.new)
      }
      Y_boot_star <- function(d,j){
        d2 <- d[j,]
        rownames(d2)<- NULL
        data[[imp]] = list(X=d2[,1:ncol(x_orig[[imp]])] , M=d2[,(ncol(x_orig[[imp]])+1):(ncol(x_orig[[imp]])+ncol(m_orig[[imp]]))], Y = d2[,((ncol(x_orig[[imp]])+ncol(m_orig[[imp]]))+1):ncol(d2)])
        perm_fit <- block.pls(data[[imp]], indY = 3,mode = "canonical",scheme='horst', ncomp =ncomp, scale=TRUE, max.iter = 1000)
        sign_flip <- procrustes(as.matrix(perm_fit$loadings$Y[,comp]),as.matrix(orig_fit[[1]]$loadings$Y[,comp]))$R
        loadings.star <- ginv(perm_fit$X$Y)%*%(sign_flip*perm_fit$variates$Y[,comp])
        return(loadings.star)
      }
      ab_boot <- function(d,j){
        d2 <- d[j,]
        rownames(d2)<- NULL
        data[[imp]] = list(X=d2[,1:ncol(x_orig[[imp]])] , M=d2[,(ncol(x_orig[[imp]])+1):(ncol(x_orig[[imp]])+ncol(m_orig[[imp]]))], Y = d2[,((ncol(x_orig[[imp]])+ncol(m_orig[[imp]]))+1):ncol(d2)])
        perm_fit <- block.pls(data[[imp]], indY = 3,mode = "canonical",scheme='horst', ncomp =ncomp, scale=TRUE, max.iter = 1000)
        perm_fit$variates$X <- scale(perm_fit$variates$X)
        perm_fit$variates$M <- scale(perm_fit$variates$M)
        perm_fit$variates$Y <- scale(perm_fit$variates$Y)
        a<-  lm(perm_fit$variates$M[,comp] ~ perm_fit$variates$X[,comp])$coefficients[2]
        b<-  lm(perm_fit$variates$Y[,comp] ~ perm_fit$variates$X[,comp]+perm_fit$variates$M[,comp])$coefficients[3]
        c<-  lm(perm_fit$variates$Y[,comp] ~ perm_fit$variates$X[,comp])$coefficients[2]
        c_med<-  lm(perm_fit$variates$Y[,comp] ~ perm_fit$variates$X[,comp]+perm_fit$variates$M[,comp])$coefficients[2]
        b_before <- lm(perm_fit$variates$Y[,comp] ~ perm_fit$variates$M[,comp])$coefficients[2]
        ab <- a*b
        perm_ab<- rbind(c, c_med,a,b,ab, b_before)
        return(perm_ab)
      }
      X_boot_result[[imp]] <- boot(boot_data, X_boot, R=nboot)
      M_boot_result[[imp]] <- boot(boot_data, M_boot, R=nboot)
      Y_boot_result[[imp]] <- boot(boot_data, Y_boot, R=nboot)
      ab_boot_result[[imp]] <- boot(boot_data, ab_boot, R=nboot)
      X_SE[[imp]] <- diag(sqrt(var(X_boot_result[[imp]]$t)))
      M_SE[[imp]] <- diag(sqrt(var(M_boot_result[[imp]]$t)))
      Y_SE[[imp]] <- diag(sqrt(var(Y_boot_result[[imp]]$t)))
      ab_SE[[imp]] <- diag(sqrt(var(ab_boot_result[[imp]]$t)))
      
      if(do_star == TRUE){
        X_star_SE[[imp]] <- sqrt(diag(var(X_boot_star_result[[imp]]$t)))
        M_star_SE[[imp]] <- sqrt(diag(var(M_boot_star_result[[imp]]$t)))
        Y_star_SE[[imp]] <- sqrt(diag(var(Y_boot_star_result[[imp]]$t)))
        X_boot_star_result[[imp]] <- boot(boot_data, X_boot_star, R=nboot)
        M_boot_star_result[[imp]] <- boot(boot_data, M_boot_star, R=nboot)
        Y_boot_star_result[[imp]] <- boot(boot_data, Y_boot_star, R=nboot)
      }
      
    }
    G <- H <- I <- J <- K <- L<-M <- N <- G2 <- H2 <- I2 <- J2 <- K2 <- M2 <- ab_pool <-  NULL
    for(imp in 1:5){
      if(imp==1){
        G <- X_boot_result[[imp]]$t0
        H <- Y_boot_result[[imp]]$t0
        I <- (X_SE[[imp]])^2
        J <- (Y_SE[[imp]])^2
        K <- M_boot_result[[imp]]$t0
        L <- ab_boot_result[[imp]]$t0
        ab_pool <- ab_boot_result[[imp]]$t
        M <- (M_SE[[imp]])^2
        N <- (ab_SE[[imp]])^2
        if(do_star == TRUE){
          G2 <- X_boot_star_result[[imp]]$t0
          H2 <- Y_boot_star_result[[imp]]$t0
          I2 <- (X_star_SE[[imp]])^2
          J2 <- (Y_star_SE[[imp]])^2
          K2 <- M_boot_star_result[[imp]]$t0
          M2 <- (M_star_SE[[imp]])^2
        }
        
      }
      if(imp>1){
        G <- G  + X_boot_result[[imp]]$t0
        H <- H + Y_boot_result[[imp]]$t0
        I <- I + (X_SE[[imp]])^2
        J <- J + (Y_SE[[imp]])^2 
        K <- K  + M_boot_result[[imp]]$t0
        L <- L + ab_boot_result[[imp]]$t0
        ab_pool <- rbind(ab_pool,ab_boot_result[[imp]]$t)
        M <- M + (M_SE[[imp]])^2
        N <- N + (ab_SE[[imp]])^2 
        if(do_star == TRUE){
          G2 <- G2  + X_boot_star_result[[imp]]$t0
          H2 <- H2 + Y_boot_star_result[[imp]]$t0
          I2 <- I2 + (X_star_SE[[imp]])^2
          J2 <- J2 + (Y_star_SE[[imp]])^2 
          K2 <- K2  + M_boot_star_result[[imp]]$t0
          M2 <- M2 + (M_star_SE[[imp]])^2
        }
        
      }
    }
    av_X_boot_result <- G/5
    av_Y_boot_result <- H/5
    av_X_SE <- sqrt(I/5) #Note this is only the within imputation standard error
    av_Y_SE <- sqrt(J/5)
    av_M_boot_result <- K/5
    av_ab_boot_result <- L/5
    av_M_SE <- sqrt(M/5) 
    av_ab_SE <- sqrt(N/5)
    if(do_star == TRUE){
      av_X_boot_star_result <- G2/5
      av_Y_boot_star_result <- H2/5
      av_X_star_SE <- sqrt(I2/5) #Note this is only the within imputation standard error
      av_Y_star_SE <- sqrt(J2/5)
      av_M_boot_star_result <- K2/5
      av_M_star_SE <- sqrt(M2/5) 
    }
    ab_ci <- data.frame(lower_CI=rep(NA, length(av_ab_SE)), upper_CI=rep(NA, length(av_ab_SE)))
    for(path in 1:length(av_ab_SE)){
      ab_ci[path,1] <- sort(ab_pool[,path])[nboot*5*(0.5-(0.95/2))]
      ab_ci[path,2] <- sort(ab_pool[,path])[nboot*5*(0.5+(0.95/2))]
    }
    
    #update the standard error
    F <- G <- H<- I<- F2 <- G2 <- H2 <-  NULL
    for(imp in 1:5){
      if(imp==1){
        F <- (X_boot_result[[imp]]$t0-av_X_boot_result)^2
        G <- (Y_boot_result[[imp]]$t0-av_Y_boot_result)^2
        H <- (M_boot_result[[imp]]$t0-av_M_boot_result)^2
        I <- (ab_boot_result[[imp]]$t0-av_ab_boot_result)^2
        if(do_star == TRUE){
          F2 <- (X_boot_star_result[[imp]]$t0-av_X_boot_star_result)^2
          G2 <- (Y_boot_star_result[[imp]]$t0-av_Y_boot_star_result)^2
          H2 <- (M_boot_star_result[[imp]]$t0-av_M_boot_star_result)^2
        }

        
      }
      if(imp>1){
        F <- F + (X_boot_result[[imp]]$t0-av_X_boot_result)^2 
        G <- G + (Y_boot_result[[imp]]$t0-av_Y_boot_result)^2
        H <- H + (M_boot_result[[imp]]$t0-av_M_boot_result)^2 
        I <- I + (ab_boot_result[[imp]]$t0-av_ab_boot_result)^2
        if(do_star == TRUE){
          F2 <- F2 + (X_boot_star_result[[imp]]$t0-av_X_boot_star_result)^2
          G2 <- G2 + (Y_boot_star_result[[imp]]$t0-av_Y_boot_star_result)^2
          H2 <- H2 + (M_boot_star_result[[imp]]$t0-av_M_boot_star_result)^2
        }
      }
    }
    pooled_X_SE <- sqrt((av_X_SE)^2 + ((1+ (1/m))*(1/(m-1))*F)) #Remember var = SE^2
    pooled_Y_SE <- sqrt((av_Y_SE)^2 + ((1+ (1/m))*(1/(m-1))*G)) 
    pooled_M_SE <- sqrt((av_M_SE)^2 + ((1+ (1/m))*(1/(m-1))*H)) 
    pooled_ab_SE <- sqrt((av_ab_SE)^2 + ((1+ (1/m))*(1/(m-1))*I)) 
    if(do_star == TRUE){
      pooled_X_star_SE <- sqrt((av_X_star_SE)^2 + ((1+ (1/m))*(1/(m-1))*F2)) 
      pooled_Y_star_SE <- sqrt((av_Y_star_SE)^2 + ((1+ (1/m))*(1/(m-1))*G2)) 
      pooled_M_star_SE <- sqrt((av_M_star_SE)^2 + ((1+ (1/m))*(1/(m-1))*H2)) 
    }
    
    if(do_star == TRUE){
      if(comp==1){
        load_list <- c('load') #As load star= load for the first component
      }
      if(comp>1){
        load_list <- c('load', 'load_star')
      }
    }
    if(do_star == FALSE){
      load_list <- c('load')
    }

    for(load_type in load_list){
      if(load_type=='load_star'){
        av_X_boot_result <- av_X_boot_star_result
        pooled_X_SE <- pooled_X_star_SE
        av_M_boot_result <- av_M_boot_star_result
        pooled_M_SE <- pooled_M_star_SE
        av_Y_boot_result <- av_Y_boot_star_result
        pooled_Y_SE <- pooled_Y_star_SE
      }
      
      X_ratio <- data.frame(abs(av_X_boot_result)/pooled_X_SE)
      X_ratio_group <- data.frame(rep(1,nrow(X_ratio))) #To get rule of thumb reliable loadings using t-distribution
      X_ratio_group[which(X_ratio < qt(1-(0.05/2),(nrow(x_orig[[imp]])-1))),] <- 0
      Y_ratio <- data.frame(abs(av_Y_boot_result)/pooled_Y_SE)
      Y_ratio_group <- data.frame(rep(1,nrow(Y_ratio))) #To get rule of thumb reliable loadings
      Y_ratio_group[which(Y_ratio < qt(1-(0.05/2),(nrow(x_orig[[imp]])-1))),] <- 0
      M_ratio <- data.frame(abs(av_M_boot_result)/pooled_M_SE)
      M_ratio_group <- data.frame(rep(1,nrow(M_ratio))) #To get rule of thumb reliable loadings
      M_ratio_group[which(M_ratio < qt(1-(0.05/2),(nrow(x_orig[[imp]])-1))),] <- 0
      ab_ratio <- data.frame(abs(av_ab_boot_result)/pooled_ab_SE)
      ab_ratio_group <- data.frame(rep(1,nrow(ab_ratio))) #To get rule of thumb reliable loadings
      ab_ratio_group[which(ab_ratio < qt(1-(0.05/2),(nrow(x_orig[[imp]])-1))),] <- 0
      for_brain_plot <- av_M_boot_result*M_ratio_group
      
      flip <- sign(orig_fit[[1]]$loadings$Y[1,comp])
      av_X_boot_result <- flip*av_X_boot_result
      av_Y_boot_result <- flip*av_Y_boot_result
      av_M_boot_result <- flip*av_M_boot_result
      
      X_load <- NULL
      X_load <- data.frame(round(av_X_boot_result, 6),av_X_boot_result-pooled_X_SE, av_X_boot_result+pooled_X_SE, X_ratio_group,colnames(x_orig[[imp]]))
      colnames(X_load) <- c('loading', 'SE_low', 'SE_high','reliable', 'name')
      fname <- paste(wd,'/', A, '_', B, '_',  C, '_',comp, '_boot_X_loadings_', load_type, '.csv', sep='' )
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
        theme_minimal(base_size=10)+
        ylab('Weights')+
        theme(axis.text.x = element_text(angle = 90, size=10, hjust=1,vjust=0.4), axis.title.x = element_blank()) +
        guides(fill=FALSE)
      
      fname <- paste(wd,'/', A, '_', B, '_',  C, '_',comp,'_X_loadings_', load_type, '.png', sep='' )
      ggsave(fname, width = 15, height = 10, units = "cm")
      
      Y_load <- data.frame(round(av_Y_boot_result, 6),av_Y_boot_result-pooled_Y_SE, av_Y_boot_result+pooled_Y_SE, Y_ratio_group,colnames(y_orig[[imp]]))
      colnames(Y_load) <- c('loading', 'SE_low', 'SE_high','reliable', 'name')
      fname <- paste(wd,'/', A, '_', B, '_',  C, '_',comp,'_boot_Y_loadings_', load_type, '.csv', sep='' )
      write.csv(Y_load, file = fname)
      if(brain_vis == TRUE){
        adj <- matrix(data=NA,nrow=nnodes,ncol=nnodes)
        brain_load <- Y_load
        brain_load$reliable[brain_load$loading==0] <- 0
        brain_load$reliable[is.na(brain_load$loading)] <- 0
        
        var <- data.frame(brain_load$loading*brain_load$reliable)
        adj[lower.tri(adj)] <- as.numeric(var[,1]) #Note, has to be lower so that it gets filled in in the right order
        adj <- data.frame(adj)
        adj[is.na(adj)] <- 0
        adj = adj + t(adj)
        rownames(adj) <- t(labels_brain)
        colnames(adj) <- t(labels_brain)
        fname <- paste(wd,'/', A, '_', B, '_',C, '_', comp, '_adj_reliable_weights_', load_type, '.csv', sep='' )
        write.csv(adj, file = fname, row.names=FALSE)
        adj_bin <- sign(adj)
        adj_bin[adj_bin==1] <- 2
        adj_bin[adj_bin==-1] <- 1
        fname <- paste(wd,'/', A, '_', B, '_',C, '_', comp, '_adj_reliable_weights_bin_', load_type, '.csv', sep='' )
        write.csv(adj_bin, file = fname, row.names=FALSE)
        adj_pos <- adj
        adj_pos[adj_pos<0] <- 0
        fname <- paste(wd,'/', A, '_', B,'_',C,  '_', comp, '_adj_reliable_weights_positive_', load_type, '.csv', sep='' )
        write.csv(adj_pos, file = fname, row.names=FALSE)
        adj_neg <- adj
        adj_neg[adj_neg>0] <- 0
        adj_neg <- abs(adj_neg)
        
        fname <- paste(wd,'/', A, '_', B,'_',C,  '_', comp, '_adj_reliable_weights_negative_', load_type, '.csv', sep='' )
        write.csv(adj_neg, file = fname, row.names=FALSE)
      }
      Y_load <- Y_load[order(-abs(Y_load$loading)),]
      Y_load$name <- factor(Y_load$name, levels = Y_load$name[order(-abs(Y_load$loading))])
      Y_load$reliable <- factor(Y_load$reliable)
      
      ggplot(data = Y_load,
             aes(x = name, y = loading,  fill=reliable)) +
        geom_bar(stat = 'identity', position = 'identity', width = 0.9) +
        #geom_errorbar(aes(ymin=-1*CI_low,ymax=-1*CI_high)) +
        geom_errorbar(aes(ymin=SE_low,ymax=SE_high)) +
        scale_fill_manual(values=cols) +
        theme_minimal(base_size=10)+
        ylab('Weights')+
        theme(axis.text.x = element_text(angle = 90, size=10, hjust=1,vjust=0.4), axis.title.x = element_blank())+
        guides(fill=FALSE)
      
      fname <- paste(wd,'/', A, '_', B, '_',  C, '_',comp,'_Y_loadings_', load_type, '.png', sep='' )
      ggsave(fname, width = 15, height = 10, units = "cm")
      
      M_load <- data.frame(round(av_M_boot_result, 6),av_M_boot_result-pooled_M_SE, av_M_boot_result+pooled_M_SE, M_ratio_group,colnames(m_orig[[imp]]))
      colnames(M_load) <- c('loading', 'SE_low', 'SE_high','reliable', 'name')
      # Add in the variables that were removed #Not needed if full data (all variables, including NA) are used for bootstrap
      # na_loads <- data.frame(matrix(data=0, nrow=length(is_na), ncol=ncol(M_load)))
      # colnames(na_loads) <- colnames(M_load)
      # rownames(na_loads) <- names(is_na)
      # na_loads$name <- names(is_na)
      # full_loads <- rbind(M_load, na_loads)
      #M_load <- full_loads[colnames(m_orig_all[[1]]),]
      fname <- paste(wd,'/', A, '_', B, '_', C, '_',comp,'_boot_M_loadings_', load_type, '.csv', sep='' )
      write.csv(M_load, file = fname)
      
      ##Put Brain vis back in here
      
      M_load <- M_load[order(-abs(M_load$loading)),]
      M_load$name <- factor(M_load$name, levels = M_load$name[order(-abs(M_load$loading))])
      M_load$reliable <- factor(M_load$reliable)
      
      ggplot(data = M_load,
             aes(x = name, y = loading,  fill=reliable)) +
        geom_bar(stat = 'identity', position = 'identity', width = 0.9) +
        #geom_errorbar(aes(ymin=-1*CI_low,ymax=-1*CI_high)) +
        geom_errorbar(aes(ymin=SE_low,ymax=SE_high)) +
        scale_fill_manual(values=cols) +
        theme_minimal(base_size=10)+
        ylab('Weights')+
        theme(axis.text.x = element_text(angle = 90, size=10, hjust=1,vjust=0.4), axis.title.x = element_blank())+
        guides(fill=FALSE)
      
      fname <- paste(wd,'/', A, '_', B, '_', C, '_', comp,'_M_loadings_', load_type, '.png', sep='' )
      ggsave(fname, width = 15, height = 10, units = "cm")
      
      
      ab_load <- cbind(round(av_ab_boot_result, 6),pooled_ab_SE,av_ab_boot_result-pooled_ab_SE, av_ab_boot_result+pooled_ab_SE, ab_ratio_group,ab_ci )
      colnames(ab_load) <- c('loading', 'SE', 'SE_low', 'SE_high', 'reliable', 'CI_lower', 'CI_upper')
      fname <- paste(wd,'/', A, '_', B, '_', C, '_',comp,'_boot_ab_loadings_', '.csv', sep='' )
      write.csv(ab_load, file = fname)
      
      X_load$table <- rep(A, nrow(X_load))
      Y_load$table <- rep(C, nrow(Y_load))
      M_load$table <- rep(B, nrow(M_load))
      
      if(brain_vis==TRUE){
        all_load <- rbind(X_load,M_load)
      }
      if(brain_vis==FALSE){
        all_load <- rbind(X_load,Y_load, M_load)
        
      }
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
        theme_bw(base_size=10)+
        ylab('Weights')+
        theme(axis.text.x = element_text(angle = 90, size=10, hjust=1,vjust=0.4), axis.title.x = element_blank()) +
        guides(fill=FALSE)+
        theme(strip.background = element_blank(),
              strip.text.x = element_blank())
      
      fname <- paste(wd,'/', A, '_', B, '_', C, '_', comp,'_X_and_Y_loadings_', load_type, '.png', sep='' )
      ggsave(fname, width = 15, height = 10, units = "cm")
    }
  }
}

  



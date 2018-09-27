### Sparse PLS with Hyperparameter optimisation ###

#Setup
library(mixOmics)
library(reshape2)
library(psych)
library(Matrix)
library(gtools)
library(ggplot2)
library(MCMCpack)
library(caret)

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

wd <- "U:/Documents/ACE_Data/Thesis_Analysis/MI_PLS_MRI/results_sparse/"
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

wd <- ("U:/Documents/ACE_Data/Thesis_Analysis/MI_PLS_MEG/results_sparse/")
setwd(wd)

##### End MEG Connectome setup (MEG sample) ######


#Hyperparameter tuning
for(contrast in c(4)){
  test_cov_LV <- NULL
  test_cor_LV <- NULL
  train_cov_LV <- NULL
  train_cor_LV <- NULL
  test_dist_cov_LV <- NULL
  test_dist_cor_LV <- NULL
  train_dist_cov_LV <- NULL
  train_dist_cor_LV <- NULL
  x_orig <-y_orig<- NULL
  ncomp <- ncol(academic[[1]])
  comp <- 1
  varselect_perm <- 10
  nperm = 1000
  for(imp in 1:5){
    ### MRI ###
    if(contrast==1){#PLS Questions - connection strength 
      x_orig[[imp]] <- half_nets
      A <- 'Connection.strength'
      is_na <- which(colSums(!is.na(x_orig[[imp]])) == 0) #If needed to remove zero columns
      x_orig[[imp]] <- x_orig[[imp]][,colSums(!is.na(x_orig[[imp]])) > 0] #If needed to remove zero columns
      y_orig[[imp]] <- questions[[imp]]
      B <- 'Questions'

    }
    if(contrast==2){#PLS Questions - Node strength 
      x_orig[[imp]] <- node_strength
      A <- 'Node.strength'
      is_na <- which(colSums(!is.na(x_orig[[imp]])) == 0) #If needed to remove zero columns
      x_orig[[imp]] <- x_orig[[imp]][,colSums(!is.na(x_orig[[imp]])) > 0] #If needed to remove zero columns
      y_orig[[imp]] <- questions[[imp]]
      B <- 'Questions'
    }
    if(contrast==3){#PLS Questions - Node strength 
      x_orig[[imp]] <- node_degree
      A <- 'Node.degree'
      is_na <- which(colSums(!is.na(x_orig[[imp]])) == 0) #If needed to remove zero columns
      x_orig[[imp]] <- x_orig[[imp]][,colSums(!is.na(x_orig[[imp]])) > 0] #If needed to remove zero columns
      y_orig[[imp]] <- questions[[imp]]
      B <- 'Questions'
    }
 
    ### MEG ###
    if(contrast==4){#PLS connection strength- questions
      x_orig[[imp]] <- half_nets
      A <- paste(freq_range, '_', corr_type, sep='')
      is_na <- which(colSums(!is.na(x_orig[[imp]])) == 0) #If needed to remove zero columns
      x_orig[[imp]] <- x_orig[[imp]][,colSums(!is.na(x_orig[[imp]])) > 0] #If needed to remove zero columns
      y_orig[[imp]] <- questions[[imp]]
      B <- 'Questions'

    }

  }
  
  #fit sPLS with different selected variables
  ncomp=min(cbind(ncol(x_orig[[imp]]),ncol(y_orig[[imp]])))
  if(ncomp >5){
    ncomp <- 5
  }
  print(comp)
  for(j in 1:varselect_perm){ #repeat several times (e.g. 100)
    print(j)
    x_train <- NULL
    x_test <- NULL
    y_train <- NULL
    y_test <- NULL
    #Randomly split data 80:20
    inTrain <- createDataPartition(y =y_orig[[1]][,1], p=.8, list=FALSE) #Make sure use the same for each imp and param
    for(imp in 1:5){
      preProc  <- preProcess(x_orig[[imp]][ inTrain,], method=c("center", "scale")) #This ensures that the test data is scaled using the tranformation used for the training set
      x_train[[imp]] <- predict(preProc, x_orig[[imp]][ inTrain,])
      x_test[[imp]]  <- predict(preProc,  x_orig[[imp]][-inTrain,])
      preProc  <- preProcess(y_orig[[imp]][ inTrain,], method=c("center", "scale")) #This ensures that the test data is scaled using the tranformation used for the training set
      y_train[[imp]] <- predict(preProc, y_orig[[imp]][ inTrain,])
      y_test[[imp]]     <- predict(preProc,  y_orig[[imp]][-inTrain,])
      
      #Run sPLS with a range of keep X
      test_fit <- NULL
      pval <- list()
      CI_upper <- list()
      k <- 1
      keepX <- round(seq(2,ncol(x_orig[[imp]]),(ncol(x_orig[[imp]])-2)/20),0) #Sequence of tuning parameter
      
      for(i in keepX){
        train_fit <- spls(x_train[[imp]], y_train[[imp]],  mode = "canonical", ncomp =comp, scale=TRUE,  keepX=i)
        test_fit$variates$X <- as.matrix(x_test[[imp]])%*%as.matrix(train_fit$loadings.star[[1]])
        test_fit$variates$Y <- as.matrix(y_test[[imp]])%*%as.matrix(train_fit$loadings.star[[2]])
        
        test_cov_LV[[k]] <- diag(cov(test_fit$variates$X, test_fit$variates$Y)) #Covariance between LV = Singular Value
        test_cor_LV[[k]] <- diag(cor(test_fit$variates$X, test_fit$variates$Y)) #Correlation between LV 
        
        train_cov_LV[[k]] <- diag(cov(train_fit$variates$X, train_fit$variates$Y)) #Covariance between LV = Singular Value
        train_cor_LV[[k]] <- diag(cor(train_fit$variates$X, train_fit$variates$Y))
        k <- k+1
      }
      if(j==1){
        test_dist_cov_LV[[imp]] <- test_cov_LV
        test_dist_cor_LV[[imp]] <- test_cor_LV
        train_dist_cov_LV[[imp]] <- train_cov_LV
        train_dist_cor_LV[[imp]] <- train_cor_LV
        
      }
      if(j>1){
        test_dist_cov_LV[[imp]] <- cbind(test_dist_cov_LV[[imp]],test_cov_LV)
        test_dist_cor_LV[[imp]] <- cbind(test_dist_cor_LV[[imp]],test_cor_LV)
        train_dist_cov_LV[[imp]] <- cbind(train_dist_cov_LV[[imp]],train_cov_LV)
        train_dist_cor_LV[[imp]] <- cbind(train_dist_cor_LV[[imp]],train_cor_LV)
      }
    }
  }

  #Pool MI sets and select highest mean
  test_mean_cov <- rowMeans(do.call(cbind, test_dist_cov_LV))
  test_mean_cor <- rowMeans(do.call(cbind, test_dist_cor_LV))
  train_mean_cov <- rowMeans(do.call(cbind, train_dist_cov_LV))
  train_mean_cor <- rowMeans(do.call(cbind, train_dist_cor_LV))
  
  #Plot...
  df2 <- data.frame(cbind(test_mean_cov,test_mean_cor,train_mean_cov,train_mean_cor, keepX))
  ggplot(df2)+
    geom_line(size = 2, aes(x=keepX, y=train_mean_cov,  colour = 'Train'))+
    geom_line(size = 2, aes(x=keepX, y=test_mean_cov,  colour = 'Test'))+
    theme_minimal()+
    ylab('Mean covariance')
  fname <- paste(wd,'/', A, '_', B,'_component_', comp, '_select_x_train_and_test.png', sep='' )
  ggsave(fname, width = 10, height = 10, units = "cm")
  fname <- paste(wd,'/', A, '_', B,'_component_', comp, '_select_x_train_and_test.csv', sep='' )
  write.csv(df2, file = fname, row.names=FALSE)
  
  #Plot just test
  ggplot(df2)+
    geom_line(size = 2, aes(x=keepX, y=test_mean_cov,  colour = 'Test'))+
    theme_minimal()+
    ylab('Mean covariance')
  fname <- paste(wd,'/', A, '_', B,'_component_', comp, '_select_x_test.png', sep='' )
  ggsave(fname, width = 10, height = 10, units = "cm")
  
  #box plots
  df3 <- t(do.call(cbind, test_dist_cov_LV)) #Put what you want to look at here
  colnames(df3) <- keepX
  data_melt <- melt(df3, value.name='covariance')
  data_melt$Var2 <- as.factor(data_melt$Var2)
  ggplot(data_melt, aes(x=Var2, y=covariance, fill=Var2)) + 
    geom_boxplot()+
    theme_minimal()+
    guides(fill=FALSE)
  fname <- paste(wd,'/', A, '_', B,'_', comp, 'select_x_test_box_plot.png', sep='' )
  ggsave(fname, width = 15, height = 10, units = "cm")
  
  #Select keepX with highest mean covariance
  select_x <- keepX[which(test_mean_cov==max(test_mean_cov))]
  
  #Run sPLS with this
  orig_fit <- NULL
  orig_cov_LV <- NULL
  orig_cor_LV <- NULL
  for(imp in 1:5){
    #Fit PLS model with chosen keepX
    orig_fit[[imp]] <- spls(x_orig[[imp]], y_orig[[imp]],  keepX= select_x, mode = "canonical", ncomp =comp,  scale=TRUE)
    orig_cov_LV[[imp]] <- diag(cov(orig_fit[[imp]]$variates$X, orig_fit[[imp]]$variates$Y)) #Covariance between LV = Singular Value
    orig_cor_LV[[imp]] <- diag(cor(orig_fit[[imp]]$variates$X, orig_fit[[imp]]$variates$Y)) #Correlation between LV 
  }
  
  #Which are picked for each MI?
  select_loads <- cbind(orig_fit[[1]]$loadings$X[,comp], orig_fit[[2]]$loadings$X[,comp], orig_fit[[3]]$loadings$X[,comp], orig_fit[[4]]$loadings$X[,comp], orig_fit[[5]]$loadings$X[,comp])
  colnames(select_loads)<-c('Imp1', 'Imp2','Imp3','Imp4','Imp5')
  
  #constrain the keepX.constraint to the variables that are selected by all imputations
  #Note, this could be less stringent
  select_loads.bin <- select_loads
  select_loads.bin[abs(select_loads)>0]<-1
  keepX.constraint = list(colnames(x_orig[[1]])[which(rowSums(select_loads.bin)==5)])
  
  #Run for additional components, constraining the variables kept on previous components each time
  for(comp in 2:ncomp){
    test_dist_cov_LV<- NULL
    test_dist_cor_LV<- NULL
    train_dist_cov_LV<- NULL
    train_dist_cor_LV<- NULL
    print(comp)
    for(j in 1:varselect_perm){ #Permute 100 times
      print(j)
      x_train <- NULL
      x_test <- NULL
      y_train <- NULL
      y_test <- NULL
      #Randomly split data 80:20
      inTrain <- createDataPartition(y =y_orig[[1]][,1], p=.8, list=FALSE) #Make sure use the same for each imp and param
      for(imp in 1:5){
        preProc  <- preProcess(x_orig[[imp]][ inTrain,], method=c("center", "scale")) #This ensures that the test data is scaled using the tranfpormation used for the training set
        x_train[[imp]] <- predict(preProc, x_orig[[imp]][ inTrain,])
        x_test[[imp]]  <- predict(preProc,  x_orig[[imp]][-inTrain,])
        preProc  <- preProcess(y_orig[[imp]][ inTrain,], method=c("center", "scale")) #This ensures that the test data is scaled using the tranfpormation used for the training set
        y_train[[imp]] <- predict(preProc, y_orig[[imp]][ inTrain,])
        y_test[[imp]]     <- predict(preProc,  y_orig[[imp]][-inTrain,])
        
        #Run sPLS with a range of keep X
        test_fit <- NULL
        pval <- list()
        CI_upper <- list()
        k <- 1
        keepX <- round(seq(2,ncol(x_orig[[imp]]),(ncol(x_orig[[imp]])-2)/20),0) #Sequence of tuning parameter for given component
        
        for(i in keepX){
          train_fit <- spls(x_train[[imp]], y_train[[imp]],  mode = "canonical", ncomp =comp, keepX=i, keepX.constraint=keepX.constraint,scale=TRUE)
          test_fit$variates$X <- as.matrix(x_test[[imp]])%*%as.matrix(train_fit$loadings.star[[1]])
          test_fit$variates$Y <- as.matrix(y_test[[imp]])%*%as.matrix(train_fit$loadings.star[[2]])
          
          test_cov_LV[[k]] <- diag(cov(test_fit$variates$X, test_fit$variates$Y))[comp] #Covariance between LV = Singular Value
          test_cor_LV[[k]] <- diag(cor(test_fit$variates$X, test_fit$variates$Y))[comp] #Correlation between LV 
          
          train_cov_LV[[k]] <- diag(cov(train_fit$variates$X, train_fit$variates$Y))[comp] #Covariance between LV = Singular Value
          train_cor_LV[[k]] <- diag(cor(train_fit$variates$X, train_fit$variates$Y))[comp]
          k <- k+1
        }
        if(j==1){
          test_dist_cov_LV[[imp]] <- test_cov_LV
          test_dist_cor_LV[[imp]] <- test_cor_LV
          train_dist_cov_LV[[imp]] <- train_cov_LV
          train_dist_cor_LV[[imp]] <- train_cor_LV
          
        }
        if(j>1){
          test_dist_cov_LV[[imp]] <- cbind(test_dist_cov_LV[[imp]],test_cov_LV)
          test_dist_cor_LV[[imp]] <- cbind(test_dist_cor_LV[[imp]],test_cor_LV)
          train_dist_cov_LV[[imp]] <- cbind(train_dist_cov_LV[[imp]],train_cov_LV)
          
        }
      }
    }
    
    #Pool MI sets and select highest mean
    test_mean_cov <- rowMeans(do.call(cbind, test_dist_cov_LV))
    test_mean_cor <- rowMeans(do.call(cbind, test_dist_cor_LV))
    train_mean_cov <- rowMeans(do.call(cbind, train_dist_cov_LV))
    train_mean_cor <- rowMeans(do.call(cbind, train_dist_cor_LV))
    
    #Plot...
    df2 <- data.frame(cbind(test_mean_cov,test_mean_cor,train_mean_cov,train_mean_cor, keepX))
    ggplot(df2)+
      geom_line(size = 2, aes(x=keepX, y=train_mean_cov,  colour = 'Train'))+
      geom_line(size = 2, aes(x=keepX, y=test_mean_cov,  colour = 'Test'))+
      theme_minimal()+
      ylab('Mean covariance')
    fname <- paste(wd,'/', A, '_', B,'_component_', comp, '_select_x_train_and_test.png', sep='' )
    ggsave(fname, width = 15, height = 10, units = "cm")
    fname <- paste(wd,'/', A, '_', B,'_component_', comp, '_select_x_train_and_test.csv', sep='' )
    write.csv(df2, file = fname, row.names=FALSE)
    
    #Plot just test
    ggplot(df2)+
      geom_line(size = 2, aes(x=keepX, y=test_mean_cov,  colour = 'Test'))+
      theme_minimal()+
      ylab('Mean covariance')
    fname <- paste(wd,'/', A, '_', B,'_component_', comp, '_select_x_test.png', sep='' )
    ggsave(fname, width = 15, height = 10, units = "cm")
    
    #box plots
    df3 <- t(do.call(cbind, test_dist_cov_LV)) #Put what you want to look at here
    colnames(df3) <- keepX
    data_melt <- melt(df3, value.name='covariance')
    data_melt$Var2 <- as.factor(data_melt$Var2)
    ggplot(data_melt, aes(x=Var2, y=covariance, fill=Var2)) + 
      geom_boxplot()+
      theme_minimal()+
      guides(fill=FALSE)
    fname <- paste(wd,'/', A, '_', B,'_', comp, 'select_x_test_box_plot.png', sep='' )
    ggsave(fname, width = 15, height = 10, units = "cm")
    
    #Select keepX with highest mean covariance
    select_x <- keepX[which(test_mean_cov==max(test_mean_cov))]

    #Run sPLS with this
    orig_fit <- NULL
    orig_cov_LV <- NULL
    orig_cor_LV <- NULL
    for(imp in 1:5){
      #Fit PLS model with chosen keepX
      orig_fit[[imp]] <- spls(x_orig[[imp]], y_orig[[imp]],  keepX= select_x, mode = "canonical", ncomp =comp, keepX.constraint=keepX.constraint, scale=TRUE)
      orig_cov_LV[[imp]] <- diag(cov(orig_fit[[imp]]$variates$X, orig_fit[[imp]]$variates$Y)) #Covariance between LV = Singular Value
      orig_cor_LV[[imp]] <- diag(cor(orig_fit[[imp]]$variates$X, orig_fit[[imp]]$variates$Y)) #Correlation between LV 
    }
    #Which are picked for each MI?
    select_loads <- cbind(orig_fit[[1]]$loadings$X[,comp], orig_fit[[2]]$loadings$X[,comp], orig_fit[[3]]$loadings$X[,comp], orig_fit[[4]]$loadings$X[,comp], orig_fit[[5]]$loadings$X[,comp])
    colnames(select_loads)<-c('Imp1', 'Imp2','Imp3','Imp4','Imp5')
    #constrain the keepX.constraint to the variables that are selected by all imputations
    #Note, this could be less stringent
    select_loads.bin <- select_loads
    select_loads.bin[abs(select_loads)>0]<-1
    keepX.constraint = c(keepX.constraint,list(colnames(x_orig[[1]])[which(rowSums(select_loads.bin)==5)])) #Update for additional components each time
  }
  
  #Run sPLS with the final selected variables for each component and permute data to test significance of each component
  orig_fit <- NULL
  orig_cov_LV <- NULL
  orig_cor_LV <- NULL
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
  for(imp in 1:5){
    
    #Fit PLS model with just these selected variables
    orig_fit[[imp]] <- spls(x_orig[[imp]], y_orig[[imp]], keepX.constraint=keepX.constraint,mode = "canonical", ncomp =ncomp, scale=TRUE)
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
      perm_fit[[imp]] <- spls(x_perm, y_orig[[imp]], keepX.constraint=keepX.constraint,mode = "canonical", ncomp = ncomp, scale=TRUE)
      if(i ==1){
        permuted_cov_LV[[imp]] <- diag(cov(perm_fit[[imp]]$variates$X, perm_fit[[imp]]$variates$Y))
        permuted_cor_LV[[imp]] <- diag(cor(perm_fit[[imp]]$variates$X, perm_fit[[imp]]$variates$Y))
      }
      else{
        x_perm_orig <- x_orig[[1]][sample_for_perm[[i]],] #Permuted data
        perm_fit_orig <- spls(x_perm_orig, y_orig[[1]],  keepX.constraint=keepX.constraint,mode = "canonical", ncomp = ncomp, scale=TRUE)
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
    z <- qnorm(1-pvalues[comp,])  # transform to z-scale
    zmean <- mean(z)
    imp_var <- sum((zmean-z)^2)/(5-1)
    total_var <- 1 + (1 + (1/5))*imp_var
    p_pooled[comp] <- 1-pnorm( zmean / sqrt(total_var)) # average and transform back
  }
  av_pval <- cbind(c(1:ncomp), data.frame(p_pooled))
  av_CI_upper <- G/5
  av_orig_cov_LV <- H/5
  av_orig_cor_LV <- J/5
  
  #Save results 
  fname <- paste(wd, '/' , A, '_', B,'_pvals_CI_and_av_cov_LVs.csv', sep='' )
  df <- cbind(av_pval,av_CI_upper, av_orig_cov_LV, av_orig_cor_LV)
  colnames(df) <- c('comp', 'p','Av.CI', 'Av.cov.LVs', 'Av.cor.LVs')
  write.csv(df, file = fname)
  
  fname <- paste(wd, '/' , A, '_', B,'_pvals_each_imputation.csv', sep='' )
  df <- cbind(c(1:ncomp), pvalues)
  colnames(df) <- c('comp', 'imp 1', 'imp 2', 'imp 3', 'imp 4', 'imp 5')
  write.csv(df, file = fname)
  
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
  for(comp in 1:ncomp){
    fname <- paste(wd, '/' , A, '_', B, '_component_', comp, '_selected_variables.csv', sep='' )
    write.csv(keepX.constraint[[comp]], file = fname)
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
nboot = 500
orig_fit <- NULL
perm_fit <- NULL
X_boot_result <- NULL
Y_boot_result <- NULL
X_SE <- NULL
Y_SE <- NULL
x_orig <- NULL
y_orig <- NULL
ncomp <- 2
for(contrast in c(9)){
  #Get data 
  for(imp in 1:5){
    ### MRI ###
    if(contrast==1){#PLS Questions - connection strength 
      x_orig[[imp]] <- half_nets
      A <- 'Connection.strength'
      y_orig[[imp]] <- questions[[imp]]
      B <- 'Questions'
      brain_vis <- TRUE
      nnodes <- 85
      labels_brain <- labels_85
      
    }
    if(contrast==2){#PLS Questions - Node strength 
      x_orig[[imp]] <- node_strength
      A <- 'Node.strength'
      y_orig[[imp]] <- questions[[imp]]
      B <- 'Questions'
      brain_vis <- FALSE
    }
    if(contrast==3){#PLS Questions - Node strength 
      x_orig[[imp]] <- node_strength
      A <- 'Node.strength'
      y_orig[[imp]] <- questions[[imp]]
      B <- 'Questions'
      brain_vis <- FALSE
    }
    ### MEG ###
    if(contrast==4){#PLS connection strength- questions
      x_orig[[imp]] <- half_nets
      A <- paste(freq_range, '_', corr_type, sep='')
      y_orig[[imp]] <- questions[[imp]]
      B <- 'Questions'
      brain_vis <- TRUE
      nnodes <- 68
      labels_brain <- labels_68
      
    }
    
  }
  #Get constraints
  for(comp in 1:ncomp){ #Get list of variables to keep on x
    fname <- paste(wd, '/' , A, '_', B, '_component_', comp, '_selected_variables.csv', sep='' )
    if(comp ==1){
      keepX.constraint <- list(as.character(read.csv(file = fname)[,2]))
    }
    if(comp >1){
      keepX.constraint <- c(keepX.constraint,list(as.character(read.csv(file = fname)[,2])))
    }
  }

  for(comp in 1:ncomp){
    for(imp in 1:5){

      print(A)
      print(B)
      print(imp)
      data <- cbind(x_orig[[imp]], y_orig[[imp]])
      orig_fit[[imp]] <- spls(x_orig[[imp]], y_orig[[imp]], keepX.constraint=keepX.constraint,mode = "canonical", ncomp =ncomp, scale=TRUE)
      X_boot <- function(d,j){
        d2 <- d[j,]
        rownames(d2)<- NULL
        perm_fit<- spls(d2[,1:ncol(x_orig[[imp]])], d2[,(ncol(x_orig[[imp]])+1):ncol(d2)],  keepX.constraint=keepX.constraint,mode = "canonical", ncomp =ncomp, scale=TRUE)
        proc <- NULL
        #proc$X$X.new <- perm_fit$loadings$X
        proc$X <- procrustes(as.matrix(perm_fit$loadings$X[,comp]),as.matrix(orig_fit[[1]]$loadings$X[,comp]))
        return(proc$X$X.new)

      }
      Y_boot <- function(d,j){
        d2 <- d[j,]
        rownames(d2)<- NULL
        perm_fit <- spls(d2[,1:ncol(x_orig[[imp]])], d2[,(ncol(x_orig[[imp]])+1):ncol(d2)], keepX.constraint=keepX.constraint, mode = "canonical", ncomp =ncomp,  scale=TRUE)
        proc <- NULL
        #proc$Y$X.new <- perm_fit$loadings$Y
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
    pooled_X_SE <- sqrt((av_X_SE)^2 + ((1+ (1/m))*(1/(m-1))*F)) #Remember var = SE^2
    pooled_Y_SE <- sqrt((av_Y_SE)^2 + ((1+ (1/m))*(1/(m-1))*G)) #Remember var = SE^2
    
    X_ratio <- data.frame(abs(av_X_boot_result)/pooled_X_SE)
    X_ratio_group <- data.frame(rep(1,nrow(X_ratio))) #To get rule of thumb reliable loadings
    X_ratio_group[which(X_ratio < 2),] <- 0
    Y_ratio <- data.frame(abs(av_Y_boot_result)/pooled_Y_SE)
    Y_ratio_group <- data.frame(rep(1,nrow(Y_ratio))) #To get rule of thumb reliable loadings
    Y_ratio_group[which(Y_ratio < 2),] <- 0
    for_brain_plot <- av_X_boot_result*X_ratio_group
    
    #flip <- sign(orig_fit[[1]]$loadings$Y[1,1])
    av_X_boot_result <- flip*av_X_boot_result
    av_Y_boot_result <- flip*av_Y_boot_result
    X_load <- NULL
    X_load <- data.frame(round(av_X_boot_result, 6),av_X_boot_result-pooled_X_SE, av_X_boot_result+pooled_X_SE, X_ratio_group,colnames(x_orig[[imp]]))
    colnames(X_load) <- c('loading', 'SE_low', 'SE_high','reliable', 'name')
    fname <- paste(wd,'/s_', A, '_', B, '_', comp, '_boot_X_loadings.csv', sep='' )
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
      fname <- paste(wd,'/s_', A, '_', B, '_', comp, '_adj_reliable_loading.csv', sep='' )
      write.csv(adj, file = fname, row.names=FALSE)
      adj_bin <- sign(adj)
      adj_bin[adj_bin==1] <- 2
      adj_bin[adj_bin==-1] <- 1
      fname <- paste(wd,'/s_', A, '_', B, '_', comp, '_adj_reliable_loading_bin.csv', sep='' )
      write.csv(adj_bin, file = fname, row.names=FALSE)
      adj_pos <- adj
      adj_pos[adj_pos<0] <- 0
      fname <- paste(wd,'/s_', A, '_', B, '_', comp, '_adj_reliable_loading_positive.csv', sep='' )
      write.csv(adj_pos, file = fname, row.names=FALSE)
      adj_neg <- adj
      adj_neg[adj_neg>0] <- 0
      adj_neg <- abs(adj_neg)
      
      fname <- paste(wd,'/s_', A, '_', B, '_', comp, '_adj_reliable_loading_negative.csv', sep='' )
      write.csv(adj_neg, file = fname, row.names=FALSE)
    }
    
    #Sort for plotting
    X_load <- X_load[order(-abs(X_load$loading)),]
    X_load$name <- factor(X_load$name, levels = X_load$name[order(-abs(X_load$loading))])
    X_load$reliable <- factor(X_load$reliable)
    
    ggplot(data = X_load,
           aes(x = name, y = loading,  fill=reliable)) +
      geom_bar(stat = 'identity', position = 'identity', width = 0.9) +
      geom_errorbar(aes(ymin=SE_low,ymax=SE_high)) +
      scale_fill_manual(values=cols) +
      theme_minimal(base_size=10)+
      theme(axis.text.x = element_text(angle = 90, size=10, hjust=1,vjust=0.4), axis.title.x = element_blank()) +
      guides(fill=FALSE)
    
    fname <- paste(wd,'/s_', A, '_', B,'_', comp,'_X_loadings.png', sep='' )
    ggsave(fname, width = 15, height = 10, units = "cm")
    
    Y_load <- data.frame(round(av_Y_boot_result, 6),av_Y_boot_result-pooled_Y_SE, av_Y_boot_result+pooled_Y_SE, Y_ratio_group,colnames(y_orig[[imp]]))
    colnames(Y_load) <- c('loading', 'SE_low', 'SE_high','reliable', 'name')
    fname <- paste(wd,'/s_', A, '_', B, '_', comp,'_boot_Y_loadings.csv', sep='' )
    write.csv(Y_load, file = fname)
    Y_load <- Y_load[order(-abs(Y_load$loading)),]
    Y_load$name <- factor(Y_load$name, levels = Y_load$name[order(-abs(Y_load$loading))])
    Y_load$reliable <- factor(Y_load$reliable)
    
    ggplot(data = Y_load,
           aes(x = name, y = loading,  fill=reliable)) +
      geom_bar(stat = 'identity', position = 'identity', width = 0.9) +
      geom_errorbar(aes(ymin=SE_low,ymax=SE_high)) +
      scale_fill_manual(values=cols) +
      theme_minimal(base_size=10)+
      theme(axis.text.x = element_text(angle = 90, size=10, hjust=1,vjust=0.4), axis.title.x = element_blank())+
      guides(fill=FALSE)
    
    fname <- paste(wd,'/s_', A, '_', B, '_', comp,'_Y_loadings.png', sep='' )
    ggsave(fname, width = 15, height = 10, units = "cm")
    
    X_load$table <- rep(A, nrow(X_load))
    Y_load$table <- rep(B, nrow(Y_load))

    
    all_load <- rbind(X_load[1:length(keepX.constraint[[comp]]),], Y_load)
    ymin <- min(all_load$SE_low)
    if(ymin>0){
      ymin <- 0
    }
    ggplot(data = all_load,
           aes(x = name, y = loading,  fill=reliable)) +
      geom_bar(stat = 'identity', position = 'dodge', width = 0.9) +
      facet_grid(~table, space='free', scale="free")+
      geom_errorbar(aes(ymin=SE_low,ymax=SE_high)) +
      scale_fill_manual(values=cols) +
      ylim(ymin, max(all_load$SE_high)) +
      theme_bw(base_size=10)+
      theme(axis.text.x = element_text(angle = 90, size=10, hjust=1,vjust=0.4), axis.title.x = element_blank()) +
      guides(fill=FALSE)+
      theme(strip.background = element_blank(),
            strip.text.x = element_blank())
    
    fname <- paste(wd,'/s_', A, '_', B, '_', comp,'_X_and_Y_loadings.png', sep='' )
    ggsave(fname, width = 15, height = 10, units = "cm")
    

  }
}




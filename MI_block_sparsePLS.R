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
#Getting MEG connectome from adjacency matrices
wd <- "U:/Documents/ACE_Data/Thesis_Analysis/MI_PLS_MEG/"
setwd(wd)

#Load connectome data
#Select frequency range and correlation matrix type
freq_range <- 'beta_nets_LASSO' #'alpha_nets', 'beta_nets', 'beta_nets_LASSO'
corr_type <- 'envPartialCorrelationRegularized'
data_dir <- paste(wd, freq_range, '/', corr_type, '/', sep='' )
setwd(wd)
ID_brain <- data.frame(read.csv('IDs_for_MEG_over_5min.csv', header = TRUE)[,1])
colnames(ID_brain) <- 'ID'

k <- 1
for(i in ID_brain$ID){
  fname <- paste(data_dir, 'nets_bAdspm12_', i, '_mc_rest.csv', sep='')
  nets <- data.frame(read.csv(fname, header=FALSE))
  if(k ==1){
    all_nets <- nets[lower.tri(nets)]
  }
  if(k>1){
    all_nets <- rbind(all_nets,nets[lower.tri(nets)])
  }
  k <- k+1
}
half_nets <- scale(all_nets,  center=TRUE, scale=TRUE)

# #Or if you already have a csv file with connectome unravelled
# half_nets <- data.frame(read.csv('half_nets_beta_partial_z.csv', header=TRUE))
# half_nets <- half_nets[,-1]

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
half_nets <- cbind(half_nets, ID_brain)
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


# wd <- setwd("U:/Documents/ACE_Data/Thesis_Analysis/MI_PLS_MEG/sparse/over_5min")
# wd <- ("U:/Documents/ACE_Data/Thesis_Analysis/MI_PLS_MEG/sparse/over_5min")
wd <- setwd("U:/Documents/ACE_Data/Thesis_Analysis/MI_connectome/sparse/Connection_strength")
wd <- "U:/Documents/ACE_Data/Thesis_Analysis/MI_connectome/sparse/Connection_strength"

#Hyperparameter tuning
for(contrast in c(3:4)){
  test_cov_LV <- NULL
  test_cor_LV <- NULL
  train_cov_LV <- NULL
  train_cor_LV <- NULL
  test_dist_cov_LV <- NULL
  test_dist_cor_LV <- NULL
  train_dist_cov_LV <- NULL
  train_dist_cor_LV <- NULL
  train_a <- train_b <- train_ab <- test_a <- test_b <- test_ab <- NULL
  a <- b <- c<- c_med <- orig_ab <- NULL
  
  train_dist_ab <- test_dist_ab<- NULL
  data_full <- NULL
  ncomp <- ncol(academic[[1]])
  comp <- 1
  varselect_perm <- 100
  nperm = 1000
  x_orig<- m_orig<- y_orig <- NULL
  for(imp in 1:5){
    if(contrast==1){#PLS questions - MEG connection strength - academic
      x_orig[[imp]] <- questions[[imp]]
      A <- 'Questions'
      m_orig[[imp]] <- half_nets
      B <- paste(freq_range, '_', corr_type, sep='')
      is_na <- which(colSums(!is.na(m_orig[[imp]])) == 0) #If needed to remove zero columns
      m_orig[[imp]] <- m_orig[[imp]][,colSums(!is.na(m_orig[[imp]])) > 0] #If needed to remove zero columns
      y_orig[[imp]] <- academic[[imp]]
      C <- 'Academic'
      data_full[[imp]] = list(X=x_orig[[imp]] , M=m_orig[[imp]], Y=as.matrix(y_orig[[imp]]))
      
    }
    if(contrast==2){#PLS questions- MEG connection strength - behaviour
      x_orig[[imp]] <- questions[[imp]]
      A <- 'Questions'
      m_orig[[imp]] <- half_nets
      B <- paste(freq_range, '_', corr_type, sep='')
      is_na <- which(colSums(!is.na(m_orig[[imp]])) == 0) #If needed to remove zero columns
      m_orig[[imp]] <- m_orig[[imp]][,colSums(!is.na(m_orig[[imp]])) > 0] #If needed to remove zero columns
      y_orig[[imp]] <- behaviour[[imp]]
      C <- 'Behaviour'
      data_full[[imp]] = list(X=x_orig[[imp]] , M=m_orig[[imp]], Y=as.matrix(y_orig[[imp]]))
      
    }
    if(contrast==3){#PLS questions - MRI connection strength - academic
      x_orig[[imp]] <- questions[[imp]]
      A <- 'Questions'
      m_orig[[imp]] <- half_nets
      B <- 'Connection.strength'
      is_na <- which(colSums(!is.na(m_orig[[imp]])) == 0) #If needed to remove zero columns
      m_orig[[imp]] <- m_orig[[imp]][,colSums(!is.na(m_orig[[imp]])) > 0] #If needed to remove zero columns
      y_orig[[imp]] <- academic[[imp]]
      C <- 'Academic'
      data_full[[imp]] = list(X=x_orig[[imp]] , M=m_orig[[imp]], Y=as.matrix(y_orig[[imp]]))
      
    }
    if(contrast==4){#PLS questions- MRI connection strength - behaviour
      x_orig[[imp]] <- questions[[imp]]
      A <- 'Questions'
      m_orig[[imp]] <- half_nets
      B <- 'Connection.strength'
      is_na <- which(colSums(!is.na(m_orig[[imp]])) == 0) #If needed to remove zero columns
      m_orig[[imp]] <- m_orig[[imp]][,colSums(!is.na(m_orig[[imp]])) > 0] #If needed to remove zero columns
      y_orig[[imp]] <- behaviour[[imp]]
      C <- 'Behaviour'
      data_full[[imp]] = list(X=x_orig[[imp]] , M=m_orig[[imp]], Y=as.matrix(y_orig[[imp]]))
    }
  }
  print(comp)
  for(j in 1:varselect_perm){ #repeat several times (e.g. 100)
    print(j)
    x_train <- m_train <-y_train<- NULL
    x_test <- m_test <-y_test<- NULL
    train_data <- NULL
    #Randomly split data 80:20
    inTrain <- createDataPartition(y =y_orig[[1]][,1], p=.8, list=FALSE) #Make sure use the same for each imp and param
    for(imp in 1:5){
      preProc  <- preProcess(x_orig[[imp]][ inTrain,], method=c("center", "scale")) #This ensures that the test data is scaled using the tranfpormation used for the training set
      x_train[[imp]] <- predict(preProc, x_orig[[imp]][ inTrain,])
      x_test[[imp]]  <- predict(preProc,  x_orig[[imp]][-inTrain,])
      preProc  <- preProcess(m_orig[[imp]][ inTrain,], method=c("center", "scale")) #This ensures that the test data is scaled using the tranfpormation used for the training set
      m_train[[imp]] <- predict(preProc, m_orig[[imp]][ inTrain,])
      m_test[[imp]]  <- predict(preProc,  m_orig[[imp]][-inTrain,])
      preProc  <- preProcess(y_orig[[imp]][ inTrain,], method=c("center", "scale")) #This ensures that the test data is scaled using the tranfpormation used for the training set
      y_train[[imp]] <- predict(preProc, y_orig[[imp]][ inTrain,])
      y_test[[imp]]     <- predict(preProc,  y_orig[[imp]][-inTrain,])
      train_data[[imp]] = list(X=x_train[[imp]] , M=m_train[[imp]], Y=as.matrix(y_train[[imp]]))
      
      
      #Run sPLS with a range of keep X
      test_fit <- NULL
      pval <- list()
      CI_upper <- list()
      k <- 1
      keepM <- round(seq(2,ncol(m_orig[[imp]]),(ncol(m_orig[[imp]])-2)/20),0) #Sequence of tuning parameter
      
      for(i in keepM){
        train_fit <- block.spls(train_data[[imp]], indY = 3,mode = "canonical",scheme='horst', ncomp =ncomp, scale=TRUE, keepX=list(X=ncol(x_orig[[1]]),M=i,Y=ncol(y_orig[[1]])))
        test_fit$variates$X <- as.matrix(x_test[[imp]])%*%as.matrix(ginv(train_fit$X$X)%*%train_fit$variates$X)
        test_fit$variates$M <- as.matrix(m_test[[imp]])%*%as.matrix(ginv(train_fit$X$M)%*%train_fit$variates$M)
        test_fit$variates$Y <- as.matrix(y_test[[imp]])%*%as.matrix(ginv(train_fit$X$Y)%*%train_fit$variates$Y)
        
        #Standardise the LVs and run linear models to calculate the indirect path coefficient
        train_fit$variates$X <- scale(train_fit$variates$X)
        train_fit$variates$M <- scale(train_fit$variates$M)
        train_fit$variates$Y <- scale(train_fit$variates$Y)
        test_fit$variates$X <- scale(test_fit$variates$X)
        test_fit$variates$M <- scale(test_fit$variates$M)
        test_fit$variates$Y <- scale(test_fit$variates$Y)
        
        train_a[[k]]<-  cbind(lm(train_fit$variates$M[,1] ~ train_fit$variates$X[,1])$coefficients[2], lm(train_fit$variates$M[,2] ~ train_fit$variates$X[,2])$coefficients[2])
        train_b[[k]]<-  cbind(lm(train_fit$variates$Y[,1] ~ train_fit$variates$X[,1]+train_fit$variates$M[,1])$coefficients[3], lm(train_fit$variates$Y[,2] ~ train_fit$variates$X[,2]+train_fit$variates$M[,2])$coefficients[3])
        train_ab[[k]]<- (train_a[[k]]*train_b[[k]])[comp]
        test_a[[k]]<-  cbind(lm(test_fit$variates$M[,1] ~ test_fit$variates$X[,1])$coefficients[2], lm(test_fit$variates$M[,2] ~ test_fit$variates$X[,2])$coefficients[2])
        test_b[[k]]<-  cbind(lm(test_fit$variates$Y[,1] ~ test_fit$variates$X[,1]+test_fit$variates$M[,1])$coefficients[3], lm(test_fit$variates$Y[,2] ~ test_fit$variates$X[,2]+test_fit$variates$M[,2])$coefficients[3])
        test_ab[[k]]<- (test_a[[k]]*test_b[[k]])[comp]
        k <- k+1
        
      }
      if(j==1){
        train_dist_ab[[imp]] <- train_ab
        test_dist_ab[[imp]] <- test_ab
      }
      if(j>1){
        train_dist_ab[[imp]] <- cbind(train_dist_ab[[imp]],train_ab)
        test_dist_ab[[imp]] <- cbind(test_dist_ab[[imp]],test_ab)

      }
    }
  }
  
  #Pool MI sets and select highest mean
  train_mean_ab <- rowMeans(do.call(cbind, train_dist_ab))
  test_mean_ab <- rowMeans(do.call(cbind, test_dist_ab))

  
  #Plot...
  df2 <- data.frame(cbind(train_mean_ab,test_mean_ab, keepM))
  ggplot(df2)+
    geom_line(size = 2, aes(x=keepM, y=train_mean_ab,  colour = 'Train'))+
    geom_line(size = 2, aes(x=keepM, y=test_mean_ab,  colour = 'Test'))+
    theme_minimal()+
    ylab('Mean covariance')
  fname <- paste(wd,'/',A, '_', B,'_', C,'_', C,'_component_', comp, '_select_x_train_and_test.png', sep='' )
  ggsave(fname, width = 20, height = 15, units = "cm")
  fname <- paste(wd,'/',A, '_', B,'_', C,'_component_', comp, '_select_x_train_and_test.csv', sep='' )
  write.csv(df2, file = fname, row.names=FALSE)
  
  #Plot just test
  ggplot(df2)+
    geom_line(size = 2, aes(x=keepM, y=test_mean_ab,  colour = 'Test'))+
    theme_minimal()+
    ylab('Mean ab')
  fname <- paste(wd,'/',A, '_', B,'_', C,'_component_', comp, '_select_x_test.png', sep='' )
  ggsave(fname, width = 20, height = 15, units = "cm")
  
  #box plots
  df3 <- t(do.call(cbind, test_dist_ab)) #Put what you want to look at here
  colnames(df3) <- keepM
  data_melt <- melt(df3, value.name='ab')
  data_melt$Var2 <- as.factor(data_melt$Var2)
  ggplot(data_melt, aes(x=Var2, y=ab, fill=Var2)) + 
    geom_boxplot()+
    theme_minimal()+
    guides(fill=FALSE)
  fname <- paste(wd,'/',A, '_', B,'_', C,'_', comp, 'select_m_test_box_plot.png', sep='' )
  ggsave(fname, width = 20, height = 15, units = "cm")
  
  #Select keepX with highest mean ab
  select_m <- keepM[which(test_mean_ab==max(test_mean_ab))]
  
  #Run sPLS with this
  orig_fit <- NULL
  orig_cov_LV <- NULL
  orig_cor_LV <- NULL
  for(imp in 1:5){
    #Fit PLS model with chosen keepX
    orig_fit[[imp]] <- block.spls(data_full[[imp]], indY = 3, keepX=list(X=ncol(x_orig[[1]]),M=select_m,Y=ncol(y_orig[[1]])),mode = "canonical",scheme='horst', ncomp =ncomp, scale=TRUE)
    }
  
  #Which are picked for each MI?
  select_loads <- cbind(orig_fit[[1]]$loadings$M[,comp], orig_fit[[2]]$loadings$M[,comp], orig_fit[[3]]$loadings$M[,comp], orig_fit[[4]]$loadings$M[,comp], orig_fit[[5]]$loadings$M[,comp])
  colnames(select_loads)<-c('Imp1', 'Imp2','Imp3','Imp4','Imp5')
  
  #constrain the keepX.constraint to the variables that are selected by all imputations
  #Note, this could be less stringent
  select_loads.bin <- select_loads
  select_loads.bin[abs(select_loads)>0]<-1
  keepX.constraint=list(X=list(colnames(x_orig[[imp]])) , M=list(colnames(m_orig[[1]])[which(rowSums(select_loads.bin)==5)]), Y=list(colnames(y_orig[[imp]])))
  
  #Run for additional components, constraining the variables kept on previous components each time
  for(comp in 2:ncomp){
    test_dist_cov_LV<- NULL
    test_dist_cor_LV<- NULL
    train_dist_cov_LV<- NULL
    train_dist_cor_LV<- NULL
    print(comp)
    for(j in 1:varselect_perm){ #repeat several times (e.g. 100)
      print(j)
      x_train <- m_train <-y_train<- NULL
      x_test <- m_test <-y_test<- NULL
      train_data <- NULL
      #Randomly split data 80:20
      inTrain <- createDataPartition(y =y_orig[[1]][,1], p=.8, list=FALSE) #Make sure use the same for each imp and param
      for(imp in 1:5){
        preProc  <- preProcess(x_orig[[imp]][ inTrain,], method=c("center", "scale")) #This ensures that the test data is scaled using the tranfpormation used for the training set
        x_train[[imp]] <- predict(preProc, x_orig[[imp]][ inTrain,])
        x_test[[imp]]  <- predict(preProc,  x_orig[[imp]][-inTrain,])
        preProc  <- preProcess(m_orig[[imp]][ inTrain,], method=c("center", "scale")) #This ensures that the test data is scaled using the tranfpormation used for the training set
        m_train[[imp]] <- predict(preProc, m_orig[[imp]][ inTrain,])
        m_test[[imp]]  <- predict(preProc,  m_orig[[imp]][-inTrain,])
        preProc  <- preProcess(y_orig[[imp]][ inTrain,], method=c("center", "scale")) #This ensures that the test data is scaled using the tranfpormation used for the training set
        y_train[[imp]] <- predict(preProc, y_orig[[imp]][ inTrain,])
        y_test[[imp]]     <- predict(preProc,  y_orig[[imp]][-inTrain,])
        train_data[[imp]] = list(X=x_train[[imp]] , M=m_train[[imp]], Y=as.matrix(y_train[[imp]]))
        
        
        #Run sPLS with a range of keep X
        test_fit <- NULL
        pval <- list()
        CI_upper <- list()
        k <- 1
        keepM <- round(seq(2,ncol(m_orig[[imp]]),(ncol(m_orig[[imp]])-2)/20),0) #Sequence of tuning parameter
        
        for(i in keepM){
          train_fit <- block.spls(train_data[[imp]], indY = 3,mode = "canonical",scheme='horst', ncomp =comp, keepX.constraint=keepX.constraint,scale=TRUE, keepX=list(X=ncol(x_orig[[1]]),M=i,Y=ncol(y_orig[[1]])))
          test_fit$variates$X <- as.matrix(x_test[[imp]])%*%as.matrix(ginv(train_fit$X$X)%*%train_fit$variates$X) #The ginv part calculates loading star
          test_fit$variates$M <- as.matrix(m_test[[imp]])%*%as.matrix(ginv(train_fit$X$M)%*%train_fit$variates$M)
          test_fit$variates$Y <- as.matrix(y_test[[imp]])%*%as.matrix(ginv(train_fit$X$Y)%*%train_fit$variates$Y)
          
          #Standardise the LVs and run linear models to calculate the indirect path coefficient
          train_fit$variates$X <- scale(train_fit$variates$X)
          train_fit$variates$M <- scale(train_fit$variates$M)
          train_fit$variates$Y <- scale(train_fit$variates$Y)
          test_fit$variates$X <- scale(test_fit$variates$X)
          test_fit$variates$M <- scale(test_fit$variates$M)
          test_fit$variates$Y <- scale(test_fit$variates$Y)
          
          train_a[[k]]<-  cbind(lm(train_fit$variates$M[,1] ~ train_fit$variates$X[,1])$coefficients[2], lm(train_fit$variates$M[,2] ~ train_fit$variates$X[,2])$coefficients[2])
          train_b[[k]]<-  cbind(lm(train_fit$variates$Y[,1] ~ train_fit$variates$X[,1]+train_fit$variates$M[,1])$coefficients[3], lm(train_fit$variates$Y[,2] ~ train_fit$variates$X[,2]+train_fit$variates$M[,2])$coefficients[3])
          train_ab[[k]]<- (train_a[[k]]*train_b[[k]])[comp]
          test_a[[k]]<-  cbind(lm(test_fit$variates$M[,1] ~ test_fit$variates$X[,1])$coefficients[2], lm(test_fit$variates$M[,2] ~ test_fit$variates$X[,2])$coefficients[2])
          test_b[[k]]<-  cbind(lm(test_fit$variates$Y[,1] ~ test_fit$variates$X[,1]+test_fit$variates$M[,1])$coefficients[3], lm(test_fit$variates$Y[,2] ~ test_fit$variates$X[,2]+test_fit$variates$M[,2])$coefficients[3])
          test_ab[[k]]<- (test_a[[k]]*test_b[[k]])[comp]
          k <- k+1
          
        }
        if(j==1){
          train_dist_ab[[imp]] <- train_ab
          test_dist_ab[[imp]] <- test_ab
        }
        if(j>1){
          train_dist_ab[[imp]] <- cbind(train_dist_ab[[imp]],train_ab)
          test_dist_ab[[imp]] <- cbind(test_dist_ab[[imp]],test_ab)
          
        }
      }
    }
    
    
    #Pool MI sets and select highest mean
    train_mean_ab <- rowMeans(do.call(cbind, train_dist_ab))
    test_mean_ab <- rowMeans(do.call(cbind, test_dist_ab))
    
    
    #Plot...
    df2 <- data.frame(cbind(train_mean_ab,test_mean_ab, keepM))
    ggplot(df2)+
      geom_line(size = 2, aes(x=keepM, y=train_mean_ab,  colour = 'Train'))+
      geom_line(size = 2, aes(x=keepM, y=test_mean_ab,  colour = 'Test'))+
      theme_minimal()+
      ylab('Mean covariance')
    fname <- paste(wd,'/',A, '_', B,'_', C,'_', C,'_component_', comp, '_select_x_train_and_test.png', sep='' )
    ggsave(fname, width = 20, height = 15, units = "cm")
    fname <- paste(wd,'/',A, '_', B,'_', C,'_component_', comp, '_select_x_train_and_test.csv', sep='' )
    write.csv(df2, file = fname, row.names=FALSE)
    
    #Plot just test
    ggplot(df2)+
      geom_line(size = 2, aes(x=keepM, y=test_mean_ab,  colour = 'Test'))+
      theme_minimal()+
      ylab('Mean ab')
    fname <- paste(wd,'/',A, '_', B,'_', C,'_component_', comp, '_select_x_test.png', sep='' )
    ggsave(fname, width = 20, height = 15, units = "cm")
    
    #box plots
    df3 <- t(do.call(cbind, test_dist_ab)) #Put what you want to look at here
    colnames(df3) <- keepM
    data_melt <- melt(df3, value.name='ab')
    data_melt$Var2 <- as.factor(data_melt$Var2)
    ggplot(data_melt, aes(x=Var2, y=ab, fill=Var2)) + 
      geom_boxplot()+
      theme_minimal()+
      guides(fill=FALSE)
    fname <- paste(wd,'/',A, '_', B,'_', C,'_', comp, 'select_m_test_box_plot.png', sep='' )
    ggsave(fname, width = 20, height = 15, units = "cm")
    
    #Select keepX with highest mean ab
    select_m <- keepM[which(test_mean_ab==max(test_mean_ab))]
    
    
    #Run sPLS with this
    orig_fit <- NULL
    orig_cov_LV <- NULL
    orig_cor_LV <- NULL
    for(imp in 1:5){
      #Fit PLS model with chosen keepX
      orig_fit[[imp]] <- block.spls(data_full[[imp]], indY = 3,  keepX.constraint=keepX.constraint,keepX=list(X=ncol(x_orig[[1]]),M=select_m,Y=ncol(y_orig[[1]])),mode = "canonical",scheme='horst', ncomp =ncomp, scale=TRUE)
    }
    
    #Which are picked for each MI?
    select_loads <- cbind(orig_fit[[1]]$loadings$M[,comp], orig_fit[[2]]$loadings$M[,comp], orig_fit[[3]]$loadings$M[,comp], orig_fit[[4]]$loadings$M[,comp], orig_fit[[5]]$loadings$M[,comp])
    colnames(select_loads)<-c('Imp1', 'Imp2','Imp3','Imp4','Imp5')
    
    #constrain the keepX.constraint to the variables that are selected by all imputations
    #Note, this could be less stringent
    select_loads.bin <- select_loads
    select_loads.bin[abs(select_loads)>0]<-1
    keepX.constraint=list(X=c(keepX.constraint$X, list(colnames(x_orig[[imp]]))) , M=c(keepX.constraint$M,list(colnames(m_orig[[1]])[which(rowSums(select_loads.bin)==5)])), Y=c(keepX.constraint$Y,list(colnames(y_orig[[imp]]))))
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
    orig_fit[[imp]] <- block.spls(data_full[[imp]], indY = 3,  keepX.constraint=keepX.constraint,mode = "canonical",scheme='horst', ncomp =ncomp, scale=TRUE)
    orig_cov_LV[[imp]] <- diag(cov(orig_fit[[imp]]$variates$X, orig_fit[[imp]]$variates$Y)) #Covariance between LV = Singular Value
    orig_cor_LV[[imp]] <- diag(cor(orig_fit[[imp]]$variates$X, orig_fit[[imp]]$variates$Y)) #Correlation between LV 
    orig_fit[[imp]]$variates$X <- scale(orig_fit[[imp]]$variates$X)
    orig_fit[[imp]]$variates$M <- scale(orig_fit[[imp]]$variates$M)
    orig_fit[[imp]]$variates$Y <- scale(orig_fit[[imp]]$variates$Y)
    
    a[[imp]]<-  cbind(lm(orig_fit[[imp]]$variates$M[,1] ~ orig_fit[[imp]]$variates$X[,1])$coefficients[2], lm(orig_fit[[imp]]$variates$M[,2] ~ orig_fit[[imp]]$variates$X[,2])$coefficients[2])
    b[[imp]]<-  cbind(lm(orig_fit[[imp]]$variates$Y[,1] ~ orig_fit[[imp]]$variates$X[,1]+orig_fit[[imp]]$variates$M[,1])$coefficients[3], lm(orig_fit[[imp]]$variates$Y[,2] ~ orig_fit[[imp]]$variates$X[,2]+orig_fit[[imp]]$variates$M[,2])$coefficients[3])
    c[[imp]]<-  cbind(lm(orig_fit[[imp]]$variates$Y[,1] ~ orig_fit[[imp]]$variates$X[,1])$coefficients[2], lm(orig_fit[[imp]]$variates$Y[,2] ~ orig_fit[[imp]]$variates$X[,2])$coefficients[2])
    c_med[[imp]]<-  cbind(lm(orig_fit[[imp]]$variates$Y[,1] ~ orig_fit[[imp]]$variates$X[,1]+orig_fit[[imp]]$variates$M[,1])$coefficients[2], lm(orig_fit[[imp]]$variates$Y[,2] ~ orig_fit[[imp]]$variates$X[,2]+orig_fit[[imp]]$variates$M[,2])$coefficients[2])
    orig_ab[[imp]]<- (a[[imp]]*b[[imp]])  
  }
  
  #Save results 
  exp_var <- cbind(data.frame(round(orig_fit[[1]]$explained_variance$X,4)),data.frame(round(orig_fit[[1]]$explained_variance$Y,4)))
  First_LV <- cbind(data.frame(orig_fit[[1]]$variates$X[,1]),data.frame(orig_fit[[1]]$variates$M[,1]), data.frame(orig_fit[[1]]$variates$Y[,1]))
  X_1st_loadings <- data.frame(orig_fit[[1]]$loadings$X[,1])
  M_1st_loadings <- data.frame(orig_fit[[1]]$loadings$M[,1])
  Y_1st_loadings <- data.frame(orig_fit[[1]]$loadings$Y[,1])
  path_model <- rbind(c[[1]], c_med[[1]], a[[1]], b[[1]], orig_ab[[1]])
  rownames(path_model) <- c(paste(A, '->', C ),paste(A, '->', C, 'with_mediator'),paste(A, '->', B ),paste(B, '->', C ),paste(A, '->', B,'->', C  ))
  
  for(imp in 2:5){
    exp_var <- cbind(exp_var, data.frame(round(orig_fit[[imp]]$explained_variance$X,4)),data.frame(round(orig_fit[[imp]]$explained_variance$Y,4)))
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
  for(comp in 1:ncomp){
    fname <- paste(wd, '/' ,A, '_', B,'_', C, '_component_', comp, '_selected_variables.csv', sep='' )
    write.csv(keepX.constraint$M[[comp]], file = fname)
  }
  
  #Permuted PLS to get significance of each component
  permuted_cor_LV <- NULL
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
      perm_fit[[imp]] <- block.spls(data_perm, indY = 3,keepX.constraint=keepX.constraint,mode = "canonical",scheme='horst', ncomp =ncomp, scale=TRUE, max.iter = 1000)
      
      if(i ==1){
        permuted_cor_LV[[imp]] <- diag(cor(perm_fit[[imp]]$variates$X, perm_fit[[imp]]$variates$Y))
        perm_fit[[imp]]$variates$X <- scale(perm_fit[[imp]]$variates$X)
        perm_fit[[imp]]$variates$M <- scale(perm_fit[[imp]]$variates$M)
        perm_fit[[imp]]$variates$Y <- scale(perm_fit[[imp]]$variates$Y)
        a[[imp]]<-  cbind(lm(perm_fit[[imp]]$variates$M[,1] ~ perm_fit[[imp]]$variates$X[,1])$coefficients[2], lm(perm_fit[[imp]]$variates$M[,2] ~ perm_fit[[imp]]$variates$X[,2])$coefficients[2])
        b[[imp]]<-  cbind(lm(perm_fit[[imp]]$variates$Y[,1] ~ perm_fit[[imp]]$variates$X[,1]+perm_fit[[imp]]$variates$M[,1])$coefficients[3], lm(perm_fit[[imp]]$variates$Y[,2] ~ perm_fit[[imp]]$variates$X[,2]+perm_fit[[imp]]$variates$M[,2])$coefficients[3])
        perm_ab[[imp]]<- data.frame(a[[imp]]*b[[imp]])
      }
      else{
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
    z <- qnorm(1-pvalues[comp,])  # transform to z-scale
    zmean <- mean(z)
    imp_var <- sum((zmean-z)^2)/(5-1)
    total_var <- 1 + (1 + (1/5))*imp_var
    p_pooled[comp] <- 1-pnorm( zmean / sqrt(total_var)) # average and transform back
  }
  av_pval <- cbind(c(1:ncomp), data.frame(p_pooled))
  av_CI_upper <- G/5
  av_orig_ab <- H/5
  colnames(av_orig_ab) <- 'av_orig_ab'
  fname <- paste(wd, '/' , A, '_', B,'_', C,'_pvals.csv', sep='' )
  df <- cbind(av_pval,av_CI_upper)
  colnames(df) <- c('p','SV_CI')
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
  ggsave(fname, width = 10, height = 10, units = "cm")
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
  ggsave(fname, width = 20, height = 10, units = "cm")
  
  
}



###########Bootstrap #######################
library(boot)
ncomp = 2
nboot = 500
orig_fit <-perm_fit <-  NULL
X_boot_result <- M_boot_result <- Y_boot_result <- ab_boot_result <- NULL
X_SE <- M_SE <- Y_SE <- ab_SE <- NULL
x_orig <- y_orig <- NULL
cols <- c('0' = '#999999','1'='#ec6c20') #set colours so that sig are orange

for(contrast in c(3:4)){
  for(comp in 1:ncomp){
    for(imp in 1:5){
      if(contrast==1){#PLS questions - MEG connection strength - academic
        x_orig[[imp]] <- questions[[imp]]
        A <- 'Questions'
        m_orig[[imp]] <- half_nets
        B <- paste(freq_range, '_', corr_type, sep='')
        is_na <- which(colSums(!is.na(m_orig[[imp]])) == 0) #If needed to remove zero columns
        m_orig[[imp]] <- m_orig[[imp]][,colSums(!is.na(m_orig[[imp]])) > 0] #If needed to remove zero columns
        y_orig[[imp]] <- academic[[imp]]
        C <- 'Academic'
        data_full[[imp]] = list(X=x_orig[[imp]] , M=m_orig[[imp]], Y=as.matrix(y_orig[[imp]]))
        brain_vis <- TRUE
        nnodes <- 68
        labels_brain <- labels_68
        order_brain_vis <-read.csv(file = 'U:/Documents/ACE_Data/Thesis_Analysis/brain_visualisation/reordered_for_brain_vis/attributes_68_reordered_brain_vis.csv')[,1]
        
        
      }
      if(contrast==2){#PLS questions- MEG connection strength - behaviour
        x_orig[[imp]] <- questions[[imp]]
        A <- 'Questions'
        m_orig[[imp]] <- half_nets
        B <- paste(freq_range, '_', corr_type, sep='')
        is_na <- which(colSums(!is.na(m_orig[[imp]])) == 0) #If needed to remove zero columns
        m_orig[[imp]] <- m_orig[[imp]][,colSums(!is.na(m_orig[[imp]])) > 0] #If needed to remove zero columns
        y_orig[[imp]] <- behaviour[[imp]]
        C <- 'Behaviour'
        data_full[[imp]] = list(X=x_orig[[imp]] , M=m_orig[[imp]], Y=as.matrix(y_orig[[imp]]))
        brain_vis <- TRUE
        nnodes <- 68
        labels_brain <- labels_68
        order_brain_vis <-read.csv(file = 'U:/Documents/ACE_Data/Thesis_Analysis/brain_visualisation/reordered_for_brain_vis/attributes_68_reordered_brain_vis.csv')[,1]
        
      }
      if(contrast==3){#PLS questions - MRI connection strength - academic
        x_orig[[imp]] <- questions[[imp]]
        A <- 'Questions'
        m_orig[[imp]] <- half_nets
        B <- 'Connection.strength'
        is_na <- which(colSums(!is.na(m_orig[[imp]])) == 0) #If needed to remove zero columns
        m_orig[[imp]] <- m_orig[[imp]][,colSums(!is.na(m_orig[[imp]])) > 0] #If needed to remove zero columns
        y_orig[[imp]] <- academic[[imp]]
        C <- 'Academic'
        data_full[[imp]] = list(X=x_orig[[imp]] , M=m_orig[[imp]], Y=as.matrix(y_orig[[imp]]))
        brain_vis <- TRUE
        nnodes <- 85
        labels_brain <- labels_85
        order_brain_vis <-read.csv(file = 'U:/Documents/ACE_Data/Thesis_Analysis/brain_visualisation/reordered_for_brain_vis/attributes_85_reordered_brain_vis.csv')[,1]
        
      }
      if(contrast==4){#PLS questions- MRI connection strength - behaviour
        x_orig[[imp]] <- questions[[imp]]
        A <- 'Questions'
        m_orig[[imp]] <- half_nets
        B <- 'Connection.strength'
        is_na <- which(colSums(!is.na(m_orig[[imp]])) == 0) #If needed to remove zero columns
        m_orig[[imp]] <- m_orig[[imp]][,colSums(!is.na(m_orig[[imp]])) > 0] #If needed to remove zero columns
        y_orig[[imp]] <- behaviour[[imp]]
        C <- 'Behaviour'
        data_full[[imp]] = list(X=x_orig[[imp]] , M=m_orig[[imp]], Y=as.matrix(y_orig[[imp]]))
        brain_vis <- TRUE
        nnodes <- 85
        labels_brain <- labels_85
        order_brain_vis <-read.csv(file = 'U:/Documents/ACE_Data/Thesis_Analysis/brain_visualisation/reordered_for_brain_vis/attributes_85_reordered_brain_vis.csv')[,1]
        
      }
      #Get constraints
      for(comp_vars in 1:ncomp){ #Get list of variables to keep on x
        fname <- paste(wd, '/' ,A, '_', B,'_', C, '_component_', comp_vars, '_selected_variables.csv', sep='' )
        if(comp_vars ==1){
          keepX.constraint=list(X=list(colnames(x_orig[[imp]])) , M=list(as.character(read.csv(file = fname)[,2])), Y=list(colnames(y_orig[[imp]])))
        }
        if(comp_vars >1){
          keepX.constraint=list(X=c(keepX.constraint$X, list(colnames(x_orig[[imp]]))) , M=c(keepX.constraint$M,list(as.character(read.csv(file = fname)[,2]))), Y=c(keepX.constraint$Y,list(colnames(y_orig[[imp]]))))
        }
      }
      print(paste(A,B,C))
      
      print(imp)
      boot_data <- cbind(x_orig[[imp]], m_orig[[imp]],y_orig[[imp]])
      orig_fit[[imp]] <- block.spls(data_full[[imp]], indY = 3, keepX.constraint=keepX.constraint,mode = "canonical",scheme='horst', ncomp =ncomp, scale=TRUE)
      X_boot <- function(d,j){
        d2 <- d[j,]
        rownames(d2)<- NULL
        data_full[[imp]] = list(X=d2[,1:ncol(x_orig[[imp]])] , M=d2[,(ncol(x_orig[[imp]])+1):(ncol(x_orig[[imp]])+ncol(m_orig[[imp]]))], Y = d2[,((ncol(x_orig[[imp]])+ncol(m_orig[[imp]]))+1):ncol(d2)])
        perm_fit <- block.spls(data_full[[imp]], indY = 3, keepX.constraint=keepX.constraint,mode = "canonical",scheme='horst', ncomp =ncomp, scale=TRUE)
        proc <- NULL
        proc$X <- procrustes(as.matrix(perm_fit$loadings$X[,comp]),as.matrix(orig_fit[[1]]$loadings$X[,comp]))
        return(proc$X$X.new)
      }
      M_boot <- function(d,j){
        d2 <- d[j,]
        rownames(d2)<- NULL
        data_full[[imp]] = list(X=d2[,1:ncol(x_orig[[imp]])] ,  M=d2[,(ncol(x_orig[[imp]])+1):(ncol(x_orig[[imp]])+ncol(m_orig[[imp]]))], Y = d2[,((ncol(x_orig[[imp]])+ncol(m_orig[[imp]]))+1):ncol(d2)])
        perm_fit <- block.spls(data_full[[imp]], indY = 3, keepX.constraint=keepX.constraint,mode = "canonical",scheme='horst', ncomp =ncomp, scale=TRUE)
        proc <- NULL
        proc$M <- procrustes(as.matrix(perm_fit$loadings$M[,comp]),as.matrix(orig_fit[[1]]$loadings$M[,comp]))
        return(proc$M$X.new)
      }
      Y_boot <- function(d,j){
        d2 <- d[j,]
        rownames(d2)<- NULL
        data_full[[imp]] = list(X=d2[,1:ncol(x_orig[[imp]])] , M=d2[,(ncol(x_orig[[imp]])+1):(ncol(x_orig[[imp]])+ncol(m_orig[[imp]]))], Y = d2[,((ncol(x_orig[[imp]])+ncol(m_orig[[imp]]))+1):ncol(d2)])
        perm_fit <- block.spls(data_full[[imp]], indY = 3, keepX.constraint=keepX.constraint,mode = "canonical",scheme='horst', ncomp =ncomp, scale=TRUE)
        proc <- NULL
        proc$Y <- procrustes(as.matrix(perm_fit$loadings$Y[,comp]),as.matrix(orig_fit[[1]]$loadings$Y[,comp]))
        return(proc$Y$X.new)
      }
      ab_boot <- function(d,j){
        d2 <- d[j,]
        rownames(d2)<- NULL
        data_full[[imp]] = list(X=d2[,1:ncol(x_orig[[imp]])] , M=d2[,(ncol(x_orig[[imp]])+1):(ncol(x_orig[[imp]])+ncol(m_orig[[imp]]))], Y = d2[,((ncol(x_orig[[imp]])+ncol(m_orig[[imp]]))+1):ncol(d2)])
        perm_fit <- block.spls(data_full[[imp]], indY = 3, keepX.constraint=keepX.constraint,mode = "canonical",scheme='horst', ncomp =ncomp, scale=TRUE)
        proc <- NULL
        proc$X <- procrustes(as.matrix(perm_fit$loadings$X[,comp]),as.matrix(orig_fit[[1]]$loadings$X[,comp]))
        proc$M <- procrustes(as.matrix(perm_fit$loadings$M[,comp]),as.matrix(orig_fit[[1]]$loadings$M[,comp]))
        proc$Y <- procrustes(as.matrix(perm_fit$loadings$Y[,comp]),as.matrix(orig_fit[[1]]$loadings$Y[,comp]))
        perm_fit$variates$X <- as.matrix(d2[,1:ncol(x_orig[[imp]])])%*%as.matrix(ginv(perm_fit$X$X)%*%perm_fit$variates$X)
        perm_fit$variates$M <- as.matrix(d2[,(ncol(x_orig[[imp]])+1):(ncol(x_orig[[imp]])+ncol(m_orig[[imp]]))])%*%as.matrix(ginv(perm_fit$X$M)%*%perm_fit$variates$M)
        perm_fit$variates$Y <- as.matrix(d2[,((ncol(x_orig[[imp]])+ncol(m_orig[[imp]]))+1):ncol(d2)])%*%as.matrix(ginv(perm_fit$X$Y)%*%perm_fit$variates$Y)
        perm_fit$variates$X <- scale(perm_fit$variates$X)
        perm_fit$variates$M <- scale(perm_fit$variates$M)
        perm_fit$variates$Y <- scale(perm_fit$variates$Y)
        a<-  lm(perm_fit$variates$M[,comp] ~ perm_fit$variates$X[,comp])$coefficients[2]
        b<-  lm(perm_fit$variates$Y[,comp] ~ perm_fit$variates$X[,comp]+perm_fit$variates$M[,comp])$coefficients[3]
        c<-  lm(perm_fit$variates$Y[,comp] ~ perm_fit$variates$X[,comp])$coefficients[2]
        c_med<-  lm(perm_fit$variates$Y[,comp] ~ perm_fit$variates$X[,comp]+perm_fit$variates$M[,comp])$coefficients[2]
        ab <- a*b
        perm_ab<- rbind(c, c_med,a,b,ab)
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
    }
    G <- H <- I <- J <- K <- L<-NULL
    for(imp in 1:5){
      if(imp==1){
        G <- X_boot_result[[imp]]$t0
        H <- Y_boot_result[[imp]]$t0
        I <- (X_SE[[imp]])^2
        J <- (Y_SE[[imp]])^2
        K <- M_boot_result[[imp]]$t0
        L <- ab_boot_result[[imp]]$t0
        M <- (M_SE[[imp]])^2
        N <- (ab_SE[[imp]])^2
        
      }
      if(imp>1){
        G <- G  + X_boot_result[[imp]]$t0
        H <- H + Y_boot_result[[imp]]$t0
        I <- I + (X_SE[[imp]])^2
        J <- J + (Y_SE[[imp]])^2 
        K <- K  + M_boot_result[[imp]]$t0
        L <- L + ab_boot_result[[imp]]$t0
        M <- M + (M_SE[[imp]])^2
        N <- N + (ab_SE[[imp]])^2 
      }
    }
    av_X_boot_result <- G/5
    av_Y_boot_result <- H/5
    av_X_SE <- sqrt(I/5) #Note this is only the within imputation standard error
    av_Y_SE <- sqrt(J/5)
    av_M_boot_result <- K/5
    av_ab_boot_result <- L/5
    av_M_SE <- sqrt(M/5) #Note this is only the within imputation standard error
    av_ab_SE <- sqrt(N/5)
    
    #update the standard error
    F <- G <- H<- I<- NULL
    for(imp in 1:5){
      if(imp==1){
        F <- (X_boot_result[[imp]]$t0-av_X_boot_result)^2
        G <- (Y_boot_result[[imp]]$t0-av_Y_boot_result)^2
        H <- (M_boot_result[[imp]]$t0-av_M_boot_result)^2
        I <- (ab_boot_result[[imp]]$t0-av_ab_boot_result)^2
        
      }
      if(imp>1){
        F <- F + (X_boot_result[[imp]]$t0-av_X_boot_result)^2 
        G <- G + (Y_boot_result[[imp]]$t0-av_Y_boot_result)^2
        H <- H + (M_boot_result[[imp]]$t0-av_M_boot_result)^2 
        I <- I + (ab_boot_result[[imp]]$t0-av_ab_boot_result)^2
        
      }
    }
    pooled_X_SE <- sqrt((av_X_SE)^2 + ((1+ (1/m))*(1/(m-1))*F)) #Remember var = SE^2
    pooled_Y_SE <- sqrt((av_Y_SE)^2 + ((1+ (1/m))*(1/(m-1))*G)) #Remember var = SE^2
    pooled_M_SE <- sqrt((av_M_SE)^2 + ((1+ (1/m))*(1/(m-1))*H)) #Remember var = SE^2
    pooled_ab_SE <- sqrt((av_ab_SE)^2 + ((1+ (1/m))*(1/(m-1))*I)) #Remember var = SE^2
    
    X_ratio <- data.frame(abs(av_X_boot_result)/pooled_X_SE)
    X_ratio_group <- data.frame(rep(1,nrow(X_ratio))) #To get rule of thumb reliable loadings
    X_ratio_group[which(X_ratio < 2),] <- 0
    Y_ratio <- data.frame(abs(av_Y_boot_result)/pooled_Y_SE)
    Y_ratio_group <- data.frame(rep(1,nrow(Y_ratio))) #To get rule of thumb reliable loadings
    Y_ratio_group[which(Y_ratio < 2),] <- 0
    M_ratio <- data.frame(abs(av_M_boot_result)/pooled_M_SE)
    M_ratio_group <- data.frame(rep(1,nrow(M_ratio))) #To get rule of thumb reliable loadings
    M_ratio_group[which(M_ratio < 2),] <- 0
    ab_ratio <- data.frame(abs(av_ab_boot_result)/pooled_ab_SE)
    ab_ratio_group <- data.frame(rep(1,nrow(ab_ratio))) #To get rule of thumb reliable loadings
    ab_ratio_group[which(ab_ratio < 2),] <- 0
    for_brain_plot <- av_M_boot_result*M_ratio_group
    
    # flip <- sign(orig_fit[[1]]$loadings$Y[1,comp])
    # av_X_boot_result <- flip*av_X_boot_result
    # av_Y_boot_result <- flip*av_Y_boot_result
    # av_M_boot_result <- flip*av_M_boot_result
    # av_ab_boot_result <- flip*av_ab_boot_result
    
    X_load <- NULL
    X_load <- data.frame(round(av_X_boot_result, 6),av_X_boot_result-pooled_X_SE, av_X_boot_result+pooled_X_SE, X_ratio_group,colnames(x_orig[[imp]]))
    colnames(X_load) <- c('loading', 'SE_low', 'SE_high','reliable', 'name')
    fname <- paste(wd,'/', A, '_', B, '_',  C, '_',comp, '_boot_X_loadings.csv', sep='' )
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
      theme(axis.text.x = element_text(angle = 90, size=10, hjust=1,vjust=0.4), axis.title.x = element_blank()) +
      guides(fill=FALSE)
    
    fname <- paste(wd,'/', A, '_', B, '_',  C, '_',comp,'_X_loadings.png', sep='' )
    ggsave(fname, width = 15, height = 10, units = "cm")
    
    Y_load <- data.frame(round(av_Y_boot_result, 6),av_Y_boot_result-pooled_Y_SE, av_Y_boot_result+pooled_Y_SE, Y_ratio_group,colnames(y_orig[[imp]]))
    colnames(Y_load) <- c('loading', 'SE_low', 'SE_high','reliable', 'name')
    fname <- paste(wd,'/', A, '_', B, '_',  C, '_',comp,'_boot_Y_loadings.csv', sep='' )
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
      theme_minimal(base_size=10)+
      theme(axis.text.x = element_text(angle = 90, size=10, hjust=1,vjust=0.4), axis.title.x = element_blank())+
      guides(fill=FALSE)
    
    fname <- paste(wd,'/', A, '_', B, '_',  C, '_',comp,'_Y_loadings.png', sep='' )
    ggsave(fname, width = 15, height = 10, units = "cm")
    
    M_load <- data.frame(round(av_M_boot_result, 6),av_M_boot_result-pooled_M_SE, av_M_boot_result+pooled_M_SE, M_ratio_group,colnames(m_orig[[imp]]))
    colnames(M_load) <- c('loading', 'SE_low', 'SE_high','reliable', 'name')
    # Add in the variables that were removed #Not needed if full data (all variables, including NA) are used for bootstrap
    na_loads <- data.frame(matrix(data=0, nrow=length(is_na), ncol=ncol(M_load)))
    colnames(na_loads) <- colnames(M_load)
    full_loads <- rbind(M_load, na_loads)
    M_load <- full_loads[colnames(m_orig[[1]]),]
    fname <- paste(wd,'/', A, '_', B, '_', C, '_',comp,'_boot_M_loadings.csv', sep='' )
    write.csv(M_load, file = fname)
    if(brain_vis == TRUE){
      adj <- matrix(data=NA,nrow=nnodes,ncol=nnodes)
      brain_load <- M_load
      brain_load$reliable[brain_load$loading==0] <- 0
      brain_load$reliable[is.na(brain_load$loading)] <- 0
      
      var <- data.frame(brain_load$loading*brain_load$reliable)
      adj[lower.tri(adj)] <- as.numeric(var[,1]) #Note, has to be lower so that it gets filled in in the right order
      adj <- data.frame(adj)
      adj[is.na(adj)] <- 0
      adj = adj + t(adj)
      rownames(adj) <- t(labels_brain)
      colnames(adj) <- t(labels_brain)
      adj <- adj[order_brain_vis,order_brain_vis]
      fname <- paste(wd,'/s_', A, '_', B, '_',C, '_', comp, '_adj_reliable_loading_reordered.csv', sep='' )
      write.csv(adj, file = fname, row.names=FALSE)
      adj_bin <- sign(adj)
      adj_bin[adj_bin==1] <- 2
      adj_bin[adj_bin==-1] <- 1
      fname <- paste(wd,'/s_', A, '_', B, '_', C, '_', comp, '_adj_reliable_loading_bin_reordered.csv', sep='' )
      write.csv(adj_bin, file = fname, row.names=FALSE)
      adj_pos <- adj
      adj_pos[adj_pos<0] <- 0
      fname <- paste(wd,'/s_', A, '_', B, '_', C, '_', comp, '_adj_reliable_loading_positive_reordered.csv', sep='' )
      write.csv(adj_pos, file = fname, row.names=FALSE)
      adj_neg <- adj
      adj_neg[adj_neg>0] <- 0
      adj_neg <- abs(adj_neg)
      
      fname <- paste(wd,'/s_', A, '_', B, '_', C, '_', comp, '_adj_reliable_loading_negative_reordered.csv', sep='' )
      write.csv(adj_neg, file = fname, row.names=FALSE)
    }
    M_load <- M_load[order(-abs(M_load$loading)),]
    M_load$name <- factor(M_load$name, levels = M_load$name[order(-abs(M_load$loading))])
    M_load$reliable <- factor(M_load$reliable)
    
    ggplot(data = M_load[1:300,],
           aes(x = name, y = loading,  fill=reliable)) +
      geom_bar(stat = 'identity', position = 'identity', width = 0.9) +
      #geom_errorbar(aes(ymin=-1*CI_low,ymax=-1*CI_high)) +
      geom_errorbar(aes(ymin=SE_low,ymax=SE_high)) +
      scale_fill_manual(values=cols) +
      theme_minimal(base_size=10)+
      theme(axis.text.x = element_text(angle = 90, size=10, hjust=1,vjust=0.4), axis.title.x = element_blank())+
      guides(fill=FALSE)
    
    fname <- paste(wd,'/', A, '_', B, '_', C, '_', comp,'_M_loadings.png', sep='' )
    ggsave(fname, width = 15, height = 10, units = "cm")
    
    ab_load <- cbind(round(av_ab_boot_result, 6),av_ab_boot_result-pooled_ab_SE, av_ab_boot_result+pooled_ab_SE, ab_ratio_group)
    colnames(ab_load) <- c('loading', 'SE_low', 'SE_high', 'reliable')
    fname <- paste(wd,'/', A, '_', B, '_', C, '_',comp,'_boot_ab_loadings.csv', sep='' )
    write.csv(ab_load, file = fname)
    
    X_load$table <- rep(A, nrow(X_load))
    Y_load$table <- rep(C, nrow(Y_load))
    M_load$table <- rep(B, nrow(M_load))
    
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
      theme_bw(base_size=10)+
      theme(axis.text.x = element_text(angle = 90, size=10, hjust=1,vjust=0.4), axis.title.x = element_blank()) +
      guides(fill=FALSE)+
      theme(strip.background = element_blank(),
            strip.text.x = element_blank())
    
    fname <- paste(wd,'/', A, '_', B, '_', C, '_', comp,'_X_and_Y_loadings.png', sep='' )
    ggsave(fname, width = 15, height = 10, units = "cm")
    
  }
}


#Extras

#Plot all lines for path of choosing keepX
df2 <- t(do.call(cbind, test_dist_cov_LV)) #Put what you want to look at here
colnames(df2) <- keepX
df2 <- data.frame(cbind(df2, seq(1,nrow(df2))))
colnames(df2)[ncol(df2)] <- 'V1'
df2$V1 <- factor(df2$V1) #maintains order
data_melt <- melt(df2, id="V1", variable.name="Var2", value.name="Covariance")
data_melt$Var2 <- as.numeric(data_melt$Var2)
ggplot(data_melt)+
  geom_line(size = 0.5, aes(x=Var2, y=Covariance,  colour = V1)) +
  #scale_color_distiller(palette = 'Spectral')+
  guides(colour=FALSE)+
  theme_minimal()

#select_loads[select_loads>0]<-1
select_loads.m <- melt(select_loads[1:100,], variable.name="Imputation", value.name="Loading")
ggplot(select_loads.m, aes(Var2, Var1)) +
  geom_tile(aes(fill = Loading),colour = "white") +
  scale_fill_gradient(low = "white",high = "steelblue")+
  theme_minimal()+
  theme(axis.text.y=element_blank())

#Compare euclidean distance between nodes with loading strength
coords <- read.csv("coords_85.csv", header = FALSE)
euclid <- as.matrix(dist(coords, diag = TRUE, upper = TRUE))
half_euclid <- euclid[lower.tri(euclid)]
data <-cbind(half_euclid,X_load$reliable*X_load$loading, unravelled_labels_85)
data2 <- data[which(data[,2]!=0),]
plot(data2[,1], data2[,2])




L50_sta <- function(sample1, sample2){
  z1 <- abs(sample1-mean(sample1))
  z2 <- abs(sample2-mean(sample2))
  n1 <- length(sample1)
  n2 <- length(sample2)
  miu1 <- mean(z1)
  miu2 <- mean(z2)
  var1 <- sd(z1)^2
  var2 <- sd(z2)^2
  L50 <- abs(log(miu1)-log(miu2))/sqrt(((var1/(miu1^2))/n1)+((var2/(miu2^2))/n2))
  return(L50)
}


PL50_var_judge <- function(dataset,yname="species",y=c(1,0),alpha = 0.05,iter = 1000){
  y1name <- y[1]
  y0name <- y[2]
  data1 <- dataset[which(dataset[yname]==y1name),]
  data0 <- dataset[which(dataset[yname]==y0name),]
  sameindex <- c()
  diffindex <- c()
  for(i in 1 : (ncol(data1)-1)){
    L50_0 <- L50_sta(data1[,i],data0[,i])
    pool <- c(data1[,i],data0[,i])
    L50_set <- c()
    for(j in 1 : iter){
      index <- sample(1:nrow(dataset),nrow(data1),replace = F)
      sample1 <- pool[index]
      sample2 <- pool[-index]
      sub_L50 <- L50_sta(sample1,sample2)
      L50_set <- c(L50_set,sub_L50)
    }
    pvalue <- sum(L50_set >= L50_0)/iter
    if(pvalue < alpha){
      diffindex <- c(diffindex,i)
    }
    else{
      sameindex <- c(sameindex,i)
    }
  }
  return(list(same = sameindex,diff = diffindex))
}

lr_rate <- function(dataset,iter = 200,subtrain=150,subtest=300){
  rate <- data.frame()
  pro_label <- c()
  time <- c()
  for(i in 1:iter){
    t1 <- proc.time()
    datanum <- nrow(dataset)
    if((subtrain+subtest)>datanum){
      print("subtrain and subtest more than the datanum")
    }
    else{
      sub_dataset <- dataset[sample(nrow(dataset),subtrain+subtest,replace = F),]
      index <- sample(nrow(sub_dataset),subtrain,replace=F)
      traindata <- sub_dataset[index,]
      testdata <- sub_dataset[-index,]
    }
    
    test_X <- testdata[,-which(names(testdata)=="species")]
    test_y <- testdata["species"]
    
    model <- glm(species~.,family=binomial(link='logit'),data=traindata)
    prob <- predict(model,test_X,type='response')
    prob <- round(prob,4)
    
    t2 <- proc.time()
    t <- t2-t1
    sub_pro_label <- cbind(prob,test_y)
    names(sub_pro_label) <- c("predictions","labels")
    
    #AUC
    dfpred <-prediction(sub_pro_label$predictions,sub_pro_label$labels)
    auc <- performance(dfpred,'auc')
    auc <- unlist(slot(auc,"y.values"))
    #pre
    pre_y <- round(prob,0)
    #accuracy
    accuracy <- round(sum(pre_y == test_y)/nrow(test_y),5)
    #FP
    count_0_1 = 0
    for(j in 1:nrow(test_y)){
      if(pre_y[j] == 1 && test_y[j,] == 0){
        count_0_1=count_0_1+1
      }
    }
    FP <- round(count_0_1/nrow(test_y),5)
    #FN
    count_1_0 = 0
    for(k in 1:nrow(test_y)){
      if(pre_y[k] == 0 && test_y[k,] == 1){
        count_1_0=count_1_0+1
      }
    }
    FN <- round(count_1_0/nrow(test_y),5)
    
    subrate <- c(accuracy,FP,FN,auc)
    
    rate <- rbind(rate,subrate)
    time <- c(time,t[3][[1]])
    pro_label <- rbind(pro_label,sub_pro_label)
  }
  names(rate) <- c("Accuracy","FP","FN","AUC")
  return(list(rate=rate,pro_label=pro_label,time=time))
}

fans_lr_rate <- function(dataset,iter = 200,subtrain=150,subtest=300,L=10){
  rate <- data.frame()
  pro_label <- c()
  time <- c()
  for(i in 1 : iter){
    t1 <- proc.time()
    sub_dataset <- dataset[sample(nrow(dataset),2*subtrain+subtest,replace = F),]
    index <- sample(nrow(sub_dataset),2*subtrain,replace=F)
    train <- sub_dataset[index,]
    test <- sub_dataset[-index,]
    test_y <- test["species"]
    
    sub_pro <- rep(0,length(test_y))
    cat(paste("第",i,"次迭代\n"))
    for(p in 1 : L){
      fans <- LMDRT(train,test,'species', c(1, 0), T)
      traindata <- fans$traindata
      testdata <- fans$testdata
      test_X <- testdata[,-which(names(testdata)=="species")]
      
      model <- glm(species~.,family=binomial(link='logit'),data=traindata)
      prob <- predict(model,test_X,type='response')
      prob <- round(prob,4)
      
      sub_pro <- sub_pro+prob
    }
    t2 <- proc.time()
    t <- t2-t1
    
    sub_pro <- sub_pro/L
    sub_pro_label <- cbind(sub_pro,test_y)
    names(sub_pro_label) <- c("predictions","labels")
    #pre
    pre_y <- round(sub_pro,0)
    #AUC
    dfpred <-prediction(sub_pro_label$predictions,sub_pro_label$labels)
    auc <- performance(dfpred,'auc')
    auc <- unlist(slot(auc,"y.values"))
    #accuracy
    accuracy <- round(sum(pre_y == test_y)/nrow(test_y),5)
    #FP
    count_0_1 = 0
    for(j in 1:nrow(test_y)){
      if(pre_y[j] == 1 && test_y[j,] == 0){
        count_0_1=count_0_1+1
      }
    }
    FP <- round(count_0_1/nrow(test_y),5)
    #FN
    count_1_0 = 0
    for(k in 1:nrow(test_y)){
      if(pre_y[k] == 0 && test_y[k,] == 1){
        count_1_0=count_1_0+1
      }
    }
    FN <- round(count_1_0/nrow(test_y),5)
    
    subrate <- c(accuracy,FP,FN,auc)
    rate <- rbind(rate,subrate)
    time <- c(time,t[3][[1]])
    pro_label <- rbind(pro_label,sub_pro_label)
  }
  names(rate) <- c("Accuracy","FP","FN","AUC")
  return(list(rate=rate,pro_label=pro_label,time=time))
}

afans_lr_rate <- function(dataset,iter = 200,subtrain=150,subtest=300,L=10,alpha_var=0.05,var_iter=1000){
  rate <- data.frame()
  pro_label <- c()
  time <- c()
  #FANS index
  fans_index <- PL50_var_judge(dataset = dataset,iter = var_iter,alpha = alpha_var)[[2]]
  ori_index <- setdiff(1:(ncol(dataset)-1),fans_index)
  for(i in 1 : iter){
    t1 <- proc.time()
    sub_dataset <- dataset[sample(nrow(dataset),2*subtrain+subtest,replace = F),]
    index <- sample(nrow(sub_dataset),2*subtrain,replace=F)
    train <- sub_dataset[index,]
    test <- sub_dataset[-index,]
    test_y <- test["species"]
    
    sub_pro <- rep(0,length(test_y))
    #FANS index
    #fans_index <- PL50_var_judge(dataset = train,iter = var_iter)[[2]]
    #ori_index <- setdiff(1:(ncol(dataset)-1),fans_index)
    cat(paste("第",i,"次迭代\n"))
    for(p in 1 : L){
      if(length(fans_index)>0){
        #cat(paste("第",p,"次循环\n"))
        fans <- LMDRT(cbind(train[,fans_index],train["species"]),cbind(test[,fans_index],test["species"]),'species', c(1, 0), T)
        fans_traindata <- fans$traindata
        fans_testdata <- fans$testdata
        sub_row_train_index <- which(rownames(train)%in%rownames(fans_traindata))
        traindata <- cbind(list(train[sub_row_train_index,ori_index]),fans_traindata)
        testdata <- cbind(test[,ori_index],fans_testdata)
        colnames(traindata) <- c(1:(ncol(dataset)-1),"species")
        colnames(testdata) <- c(1:(ncol(dataset)-1),"species")
      }
      else{
        #cat(paste("第",p,"次循环\n"))
        traindata <- train[sample(nrow(train),nrow(train)/2,replace = F),]
        testdata <- test
        colnames(traindata) <- c(1:(ncol(dataset)-1),"species")
        colnames(testdata) <- c(1:(ncol(dataset)-1),"species")
      }
      
      test_X <- testdata[,-which(names(testdata)=="species")]
      model <- glm(species~.,family=binomial(link='logit'),data=traindata)
      #cat(paste("第",p,"次循环模型成功建立\n"))
      prob <- predict(model,test_X,type='response')
      #cat(paste("第",p,"次循环模型成功预测\n"))
      prob <- round(prob,4)
      
      sub_pro <- sub_pro+prob
    }
    
    t2 <- proc.time()
    t <- t2-t1
    sub_pro <- sub_pro/L
    sub_pro_label <- cbind(sub_pro,test_y)
    names(sub_pro_label) <- c("predictions","labels")
    #pre
    pre_y <- round(sub_pro,0)
    #AUC
    dfpred <-prediction(sub_pro_label$predictions,sub_pro_label$labels)
    auc <- performance(dfpred,'auc')
    auc <- unlist(slot(auc,"y.values"))
    #accuracy
    accuracy <- round(sum(pre_y == test_y)/nrow(test_y),5)
    #FP
    count_0_1 = 0
    for(j in 1:nrow(test_y)){
      if(pre_y[j] == 1 && test_y[j,] == 0){
        count_0_1=count_0_1+1
      }
    }
    FP <- round(count_0_1/nrow(test_y),5)
    #FN
    count_1_0 = 0
    for(k in 1:nrow(test_y)){
      if(pre_y[k] == 0 && test_y[k,] == 1){
        count_1_0=count_1_0+1
      }
    }
    FN <- round(count_1_0/nrow(test_y),5)
    
    subrate <- c(accuracy,FP,FN,auc)
    rate <- rbind(rate,subrate)
    time <- c(time,t[3][[1]])
    
    pro_label <- rbind(pro_label,sub_pro_label)
  }
  names(rate) <- c("Accuracy","FP","FN","AUC")
  return(list(rate=rate,pro_label=pro_label,fans_index=fans_index,time=time))
}

fans2_lr_rate <- function(dataset,iter = 200,subtrain=150,subtest=300,L=10){
  rate <- data.frame()
  pro_label <- c()
  time <- c()
  for(i in 1 : iter){
    t1 <- proc.time()
    sub_dataset <- dataset[sample(nrow(dataset),2*subtrain+subtest,replace = F),]
    index <- sample(nrow(sub_dataset),2*subtrain,replace=F)
    train <- sub_dataset[index,]
    test <- sub_dataset[-index,]
    test_y <- test["species"]
    
    sub_pro <- rep(0,length(test_y))
    cat(paste("第",i,"次迭代\n"))
    for(p in 1 : L){
      fans <- LMDRT(train,test,'species', c(1, 0), T)
      fans_traindata <- fans$traindata
      fans_testdata <- fans$testdata
      sub_row_train_index <- which(rownames(train)%in%rownames(fans_traindata))
      traindata <- cbind(fans_traindata[,1:(ncol(dataset)-1)],list(train[sub_row_train_index,1:(ncol(dataset)-1)]),fans_traindata["species"])
      testdata <- cbind(fans_testdata[,1:(ncol(dataset)-1)],test[,1:(ncol(dataset)-1)],fans_testdata["species"])
      colnames(traindata) <- c(1:(2*(ncol(dataset)-1)),"species")
      colnames(testdata) <- c(1:(2*(ncol(dataset)-1)),"species")
      
      test_X <- testdata[,-which(names(testdata)=="species")]
      
      model <- glm(species~.,family=binomial(link='logit'),data=traindata)
      prob <- predict(model,test_X,type='response')
      prob <- round(prob,4)
      
      sub_pro <- sub_pro+prob
    }
    t2 <- proc.time()
    t <- t2-t1
    sub_pro <- sub_pro/L
    sub_pro_label <- cbind(sub_pro,test_y)
    names(sub_pro_label) <- c("predictions","labels")
    #pre
    pre_y <- round(sub_pro,0)
    #AUC
    dfpred <-prediction(sub_pro_label$predictions,sub_pro_label$labels)
    auc <- performance(dfpred,'auc')
    auc <- unlist(slot(auc,"y.values"))
    #accuracy
    accuracy <- round(sum(pre_y == test_y)/nrow(test_y),5)
    #FP
    count_0_1 = 0
    for(j in 1:nrow(test_y)){
      if(pre_y[j] == 1 && test_y[j,] == 0){
        count_0_1=count_0_1+1
      }
    }
    FP <- round(count_0_1/nrow(test_y),5)
    #FN
    count_1_0 = 0
    for(k in 1:nrow(test_y)){
      if(pre_y[k] == 0 && test_y[k,] == 1){
        count_1_0=count_1_0+1
      }
    }
    FN <- round(count_1_0/nrow(test_y),5)
    
    subrate <- c(accuracy,FP,FN,auc)
    rate <- rbind(rate,subrate)
    
    pro_label <- rbind(pro_label,sub_pro_label)
    time <- c(time,t[3][[1]])
  }
  names(rate) <- c("Accuracy","FP","FN","AUC")
  return(list(rate=rate,pro_label=pro_label,time=time))
}

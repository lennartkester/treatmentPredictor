
progress <- function(n) cat(sprintf("task %d is complete\n", n))

featureSelection <- function(inputdata,
                             allTestVars,
                             innerCV=F,
                             nFoldsInnerCV=10,
                             predictSurvival=F,
                             maxFeatures=10,
                             minFeatures=3,
                             fileName=NULL,
                             printProgress="screen",
                             bestN=3){
  if (predictSurvival){
    if ( is.null(inputdata$time) || is.null(inputdata$status)){
      stop("time, status variables have to be in the inputdata data frame")
    }
  }else if ( is.null(inputdata$time) || is.null(inputdata$status) || is.null(inputdata$treatment)){
    stop("time, status and treatment variables have to be in the inputdata data frame")
  }
  if (!is.data.frame(inputdata)){
    stop("inputdata has to be a data frame")
  }
  if (length(allTestVars) < 3){
    stop("at least 3 possible features have to be provided")
  }
  if ( sum(allTestVars %in% colnames(inputdata)) != length(allTestVars)){
    stop("all possible features have to be in the inputdata data frame")
  }
  if ( length(allTestVars) < maxFeatures){
    maxFeatures <- length(allTestVars)
    message("maxFeatures > number of possible features, maxFeatures reset to number of possible features")
  }
  if (!(printProgress %in% c("screen","file","no"))){
    message("printProgress has to be one of the following: screen, file, no")
  }
  strt <- Sys.time()
  if (!predictSurvival){
    #treatmentCode <- as.numeric(unique(inputdata$treatment))
    if ( sum(unique(inputdata$treatment) %in% c(0,1)) == 2 ) {
      treatmentCode <- c(0,1)
    }else if (is.factor(inputdata$treatment)){
      treatmentCode <- c(0,1)
      inputdata$treatment <- as.numeric(inputdata$treatment)-1
    }
  }
  res <- list()
  res$params <- list()
  res$params$allTestVars <- allTestVars
  res$params$predictSurvival <- predictSurvival
  res$params$maxFeatures <- maxFeatures
  res$params$minFeatures <- minFeatures
  res$params$nFoldsInnerCV <- nFoldsInnerCV
  res$params$innerCV <- innerCV
  res$params$fileName <- fileName
  if (!predictSurvival){
    res$params$treatmentCode <- treatmentCode
  }
  res$inputdata <- inputdata
  
  ss <- Surv(inputdata$time,inputdata$status)
  coxList <- cmListFinal <- coefficientList <- pvalueList <- aicList <- selectedFeatures <- list()
  bestAIC <- bestPv <- c()
  for ( i in 1:length(allTestVars)){
    if ( predictSurvival){
      coxList[[i]] <- coxph(Surv(time,status) ~ inputdata[,allTestVars[i]] ,data=inputdata)
      bestPv[i] <- summary(coxList[[i]])$logtest[3]
    }else{
      coxList[[i]] <- coxph(as.formula(paste(c("Surv(time,status) ~  (",paste(allTestVars[i],collapse=" + "), ") * treatment"))), data=inputdata)    
      #coxList[[i]] <- coxph(as.formula(paste(c("Surv(time,status) ~  (",paste(allTestVars[i],collapse=" + "), ") + treatment"))), data=inputdata)    
      bestPv[i] <- summary(coxList[[i]])$coefficients[3,5]
      #bestPv[i] <- summary(coxList[[i]])$logtest[3]
    }
    bestAIC[i] <- (2*3)-(2*coxList[[i]]$loglik[2])
  }
  cmListFinal[[1]] <- coxList
  names(bestAIC) <- names(bestPv) <- allTestVars
  best <- which.min(bestPv)
  selectedFeatures[[1]] <- names(bestPv)[best]
  pvalueList[[1]] <- bestPv
  aicList[[1]] <- bestAIC
  coefficientList[[1]] <- lapply(coxList,function(x) x$coefficients)
  if ( printProgress == "screen"){
    cat("Best combination of features:",selectedFeatures[[1]],"  p-value:",bestPv[best],"AIC:",aicList[[1]][best],"\n")
  }else if ( printProgress == "file"){
    cat("Best combination of features:",selectedFeatures[[1]],"  p-value:",bestPv[best],"AIC:",aicList[[1]][best],"\n",file=fileName,append=TRUE)
  }
  for ( i in 2:maxFeatures){
    if ( i == 2){
      bestThree <- list()
      for ( b in 1:bestN){
        if (innerCV){
          bestThree[[b]]  <- names(bestPv[order(bestPv,decreasing = F)][b])
        }else{
          bestThree[[b]]  <- names(bestPv[order(bestAIC)][b])  
        }
        
      }
    }
    tsFeature <- coxCoefficients <- combinationList <- list()
    pv <- aic <- c()
    if ( i == length(allTestVars)){
      possibleCombinations <- rbind(allTestVars,allTestVars)
    }else{
      for ( m in 1:bestN){
        startSet <- bestThree[[m]]
        restSet <- allTestVars[!(allTestVars %in% startSet)]
        for (n in 1:length(restSet)){
          if ( m == 1 & n == 1){
            possibleCombinations <- c(startSet,restSet[n])
          }else{
            possibleCombinations <- rbind(possibleCombinations,c(startSet,restSet[n]))
          }
        }
      }
    }
    possibleCombinations <- t(apply(possibleCombinations,1, function(x) x[order(x)]))
    possibleCombinations <- possibleCombinations[!duplicated(possibleCombinations),]
    if ( i == length(allTestVars)){
      possibleCombinations <- as.matrix(t(possibleCombinations))
    }
    for ( k in 1:nrow(possibleCombinations)){
      testVars <- possibleCombinations[k,]
      if ( predictSurvival){
        tsFeature[[k]] <- cph(as.formula(paste(c("Surv(time,status) ~ (",paste(testVars,collapse=" + "), ") "))), data=inputdata,x=T,y=T)
        if ( innerCV){
          testList <- cm <- list()
          set.seed(1)
          folds <- sample(c(1:nFoldsInnerCV),nrow(inputdata),replace = T)
          for ( n in unique(folds)){
            trainSet <- inputdata[folds != n,]
            testSet <- inputdata[folds == n,]
            cm[[n]] <- coxph(as.formula(paste(c("Surv(time,status) ~  (",paste(testVars,collapse=" + "), ") "))), data=trainSet)
            coxCoefficientsTemp <- cm[[n]]$coefficients
            testMat <- model.matrix(Surv(time,status) ~ . ,data=testSet[,c(testVars,"time","status")])[,-1]
            testSet$predIndex <- testMat[,testVars] %*% coxCoefficientsTemp
            testList[[n]] <- testSet
          }
          testData <- do.call(rbind,testList)
          coxCoefficients[[k]] <- apply(do.call(rbind,lapply(cm,function(x) x$coefficients)),2,mean)
          tsFeature[[k]] <- coxph(Surv(time,status) ~  predIndex, data=testData)
        }
        else{
          tsFeature[[k]] <- coxph(as.formula(paste(c("Surv(time,status) ~  (",paste(testVars,collapse=" + "), ") "))), data=inputdata)
          coxCoefficients[[k]] <- tsFeature[[k]]$coefficients
        }
      }else{
        if ( innerCV ){
          testList <- cm <- list()
          set.seed(1)
          folds <- sample(c(1:nFoldsInnerCV),nrow(inputdata),replace = T)
          for ( n in unique(folds)){
            trainSet <- inputdata[folds != n,]
            testSet <- inputdata[folds == n,]
            cm[[n]] <- coxph(as.formula(paste(c("Surv(time,status) ~  (",paste(testVars,collapse=" + "), ") * treatment "))), data=trainSet)
            coxCoefficientsTemp <- cm[[n]]$coefficients
            testMat <- model.matrix(Surv(time,status) ~ . ,data=testSet[,c(testVars,"time","status","treatment")])[,-1]
            predT1 <- treatmentCode[1] * coxCoefficientsTemp[(((length(coxCoefficientsTemp)-1)/2)+1)] + treatmentCode[1] * as.matrix(testMat[,-ncol(testMat)]) %*% as.matrix(coxCoefficientsTemp)[(((length(coxCoefficientsTemp)-1)/2)+2):length(coxCoefficientsTemp)]
            predT2 <- treatmentCode[2] * coxCoefficientsTemp[(((length(coxCoefficientsTemp)-1)/2)+1)] + treatmentCode[2] * as.matrix(testMat[,-ncol(testMat)]) %*% as.matrix(coxCoefficientsTemp)[(((length(coxCoefficientsTemp)-1)/2)+2):length(coxCoefficientsTemp)]
            testSet$predIndexDelta <- predT1 - predT2
            testList[[n]] <- testSet
          }
          testData <- do.call(rbind,testList)
          tsFeature[[k]] <- coxph(Surv(time,status) ~  predIndexDelta * treatment, data=testData)
          coxCoefficients[[k]] <- apply(do.call(rbind,lapply(cm,function(x) x$coefficients)),2,median)
        }else{
          tsFeature[[k]] <- coxph(as.formula(paste(c("Surv(time,status) ~  (",paste(testVars,collapse=" + "), ") : treatment"))), data=inputdata, x=T, y=T)  
          #tsFeature[[k]] <- cph(as.formula(paste(c("Surv(time,status) ~  (",paste(testVars,collapse=" + "), ") + treatment"))), data=inputdata, x=T,y=T)  
          coxCoefficients[[k]] <- tsFeature[[k]]$coefficients
        }
      }
      if ( predictSurvival){
        pv[k] <- summary(tsFeature[[k]])$logtest[3]    
      }else{
        pv[k] <- summary(tsFeature[[k]])$logtest[3]    
        #pv[k] <- summary(tsFeature[[k]])$coefficients[3,5]      
      }
      aic[k] <- (2*((length(testVars))+1))-(2*tsFeature[[k]]$loglik[2])
      timeElapsed <- Sys.time()-strt
    }
    pvalueList[[i]] <- pv
    combinationList[[i]] <- possibleCombinations
    cmListFinal[[i]] <- tsFeature 
    coefficientList[[i]] <- coxCoefficients
    if ( innerCV ){
      best <- which.min(pv)
    }else{
      best <- which.min(pv)  
    }
    selectedFeatures[[i]] <- possibleCombinations[best,]
    if ( i < length(allTestVars)){
      bestThree <- list()
      for ( b in 1:bestN){
        if ( innerCV ) {
          bestThree[[b]] <- possibleCombinations[order(pv,decreasing = F)[b],]  
        }else{
          bestThree[[b]] <- possibleCombinations[order(aic,decreasing = F)[b],]  
        }
      }
    }
    aicList[[i]] <- aic
    
    if ( printProgress == "screen"){
      cat("Best combination of",i,"features:",selectedFeatures[[i]],"  p-value:",pv[best]," AIC:",aic[best],"    time elapsed",timeElapsed,"\n")
    }else if ( printProgress == "file"){
      cat("Best combination of",i,"features:",selectedFeatures[[i]],"  p-value:",pv[best]," AIC:",aic[best],"    time elapsed",timeElapsed,"\n",file=fileName,append=TRUE)
    }
  }
  res$pvalueListForwardFeatureSelection <- pvalueList
  res$pvalueListForwardFeatureSelectionBest <- unlist(lapply(pvalueList,min))
  res$aicListForwardFeatureSelection <- aicList
  res$selectedForwardFeatures <- selectedFeatures
  res$models <- cmListFinal
  res$coefficients <- coefficientList
  return(res)
}


testModel <- function(res,testData){
  if (sum(res$params$allTestVars %in% colnames(testData)) != length(res$params$allTestVars)){
    stop("testVars used in model construction have to be present in testData")
  }
  predictSurvival <- res$params$predictSurvival
  if (predictSurvival){
    if ( is.null(testData$time) || is.null(testData$status)){
      stop("time, status variables have to be in the inputdata data frame")
    }
  }else if ( is.null(testData$time) || is.null(testData$status) || is.null(testData$treatment)){
    stop("time, status and treatment variables have to be in the inputdata data frame")
  }
  treatmentCode <- res$params$treatmentCode
  trainList <- testList <- coefList <- pList <- list()
  for ( i in 1:length(res$pvalueListForwardFeatureSelection)){ 
    for ( j in 1:length(res$pvalueListForwardFeatureSelection[[i]])){
      if(res$pvalueListForwardFeatureSelection[[i]][j] == min(unlist(res$pvalueListForwardFeatureSelection))){
        coxCoefficients <- res$coefficients[[i]][[j]]
        if (predictSurvival){
          testVars <- names(coxCoefficients)  
          testMat <- model.matrix(Surv(time,status) ~ . ,data=testData[,c(testVars,"time","status")])[,-1]
          trainMat <- model.matrix(Surv(time,status) ~ . ,data=res$inputdata[,c(testVars,"time","status")])[,-1]
          cm <- coxph(as.formula(paste(c("Surv(time,status) ~  (",paste(testVars,collapse=" + "), ") "))), data=res$inputdata)  
          coxCoefficients <- cm$coefficients
          testData$predIndex <- testMat[,testVars] %*% coxCoefficients
          res$inputdata$predIndex <- trainMat[,testVars] %*% coxCoefficients
        }else{
          testVars <- names(coxCoefficients)[1:((length(coxCoefficients)-1)/2)] 
          testMat <- model.matrix(Surv(time,status) ~ . ,data=testData[,c(testVars,"time","status","treatment")])[,-1]
          trainMat <- model.matrix(Surv(time,status) ~ . ,data=res$inputdata[,c(testVars,"time","status","treatment")])[,-1]
          cm <- coxph(as.formula(paste(c("Surv(time,status) ~  (",paste(testVars,collapse=" + "), ") * treatment "))), data=res$inputdata)  
          coxCoefficients <- cm$coefficients
          
          #predT1 <- treatmentCode[1] * coxCoefficients[(length(testVars)+1)] + treatmentCode[1] * as.matrix(testMat[,-ncol(testMat)]) %*% as.matrix(coxCoefficients)[grep(":treatment",names(coxCoefficients))]
          #predT2 <- treatmentCode[2] * coxCoefficients[(length(testVars)+1)] + treatmentCode[2] * as.matrix(testMat[,-ncol(testMat)]) %*% as.matrix(coxCoefficients)[grep(":treatment",names(coxCoefficients))]
          testData$predIndexDelta <- -(coxCoefficients["treatment"] +as.matrix(testMat[,testVars]) %*% as.matrix(coxCoefficients)[grep(":treatment",names(coxCoefficients))])
          #testData$predIndexDelta <- predT1 - predT2
          #predT1Train <- treatmentCode[1] * coxCoefficients[(length(testVars)+1)] + treatmentCode[1] * as.matrix(trainMat[,-ncol(trainMat)]) %*% as.matrix(coxCoefficients)[grep(":treatment",names(coxCoefficients))]
          #predT2Train <- treatmentCode[2] * coxCoefficients[(length(testVars)+1)] + treatmentCode[2] * as.matrix(trainMat[,-ncol(trainMat)]) %*% as.matrix(coxCoefficients)[grep(":treatment",names(coxCoefficients))]
          #res$inputdata$predIndexDelta <- predT1Train - predT2Train
          res$inputdata$predIndexDelta <- -(coxCoefficients["treatment"] + as.matrix(trainMat[,-ncol(trainMat)]) %*% as.matrix(coxCoefficients)[grep(":treatment",names(coxCoefficients))])
          
        }
      }
    }
  }
  res$testData <- testData
  res$bestCoefficients <- coxCoefficients
  return(res)
}


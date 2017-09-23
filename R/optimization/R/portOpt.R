library(plyr)
library(Rmosek)
library(icsUtil)
library(abind)

optWeight <- function(alpha, covMatrix, risk, riskLambda, lastWeight, turnoverThreshold=Inf,
                      buyThreshold=turnoverThreshold,sellThreshold=turnoverThreshold,
                      singleMaxWeight=Inf,upWgtLimit=singleMaxWeight, dnWgtLimit=0-singleMaxWeight,
                      longOnly=TRUE, longSideWgt = NULL, shortSideWgt = NULL, optMode= 'riskLambda',
                      quantileNum = 10, totalTurnover = NULL, leverageRate=NULL, minWgtChg=0,
                      tc =0, verbose=TRUE){
  
  # parse parameters
  if(is.null(totalTurnover))totalTurnover=Inf
  if(is.null(tc)) tc = 0
  if(is.null(longSideWgt)) longSideWgt = ifelse(longOnly,1,0.5)
  if(is.null(shortSideWgt)) shortSideWgt = ifelse(longOnly,1,-longSideWgt)
  if(is.null(leverageRate)) leverageRate = longSideWgt - shortSideWgt
  
  if(optMode == "riskLambda"){
    which.z = diag(covMatrix) < 1e-12
    if(any(which.z)){
      print('Trying to handle 0 values in covMatrix')
      diag(covMatrix)[which.z] = median(diag(covMatrix)[!which.z])*10000
    }
    deCompMatrix <- tryCatch(Matrix::chol(covMatrix), error = function(e){print(e);NULL})
    if(is.null(deCompMatrix) && riskLambda == 0){
      deCompMatrix <- array(0,dim = dim(covMatrix))
    }
    
    if(is.null(deCompMatrix)){
      print("Trying to handle non positive definite")
      addon = diag(diag(covMatrix)*0.01)
      i = 1
      while(is.null(deCompMatrix)){
        deCompMatrix = tryCatch(Matrix::chol(covMatrix + addon*i*i),error=function(e){print(e);NULL})
        i = i + 1
        stopifnot(i<100)
      }
    }
    
    alpha = matrix(alpha)
    tarNum = dim(alpha)[1]
    N = 3*tarNum
    M = tarNum+1
    controlLeverage = shortSideWgt <-1e-9
    if(controlLeverage){
      idxL = N
      N = N + 2*tarNum
      M = M + tarNum
    }
    controlRiskBound = !is.null(risk) && is.list(risk)
    totalWeight <- longSideWgt + shortSideWgt
    infWeight <- matrix(pmax(lastWeight - sellThreshold, dnWgtLimit))
    if(shortSideWgt >= -1e-9){
      infWeight = pmax(0,infWeight)
    }
    supWeight <- matrix(pmin(lastWeight+buyThreshold,upWgtLimit))
    infWeight[infWeight>supWeight]<- dnWgtLimit[infWeight > supWeight]
    EnlargeMatrix = matrix(0,nrow = N, ncol = N)
    EnlargeAlpha = matrix(0,N,1)
    EnlargeEquaMatrix = matrix(0,M,N)
    EnlargeEquaArray = matrix(0,M,1)
    
    EnlargeUniMatrix = matrix(0,nrow=1,ncol=N)
    EnlargeUniMatrix[,(tarNum+1):(3*tarNum)] = 1
    EnlargeUniBound = matrix(totalTurnover)
    
    EnlargedMatrix[1:tarNum,1:tarNum] = deCompMatrix
    EnlargeAlpha[1:tarNum] = alpha
    
    if(controlLeverage){
      totalLevArray = matrix(1,1,N)
      totalLevArray[1:(3*tarNum)] = 0
      EnlargeUniMatrix = rbind(EnlargeUniMatrix, totalLevArray)
      EnlargeUniBound = rbind(EnlargeUniBound, leverageRate)
    }
    
    if(controlRiskBound){
      Nrisk = length(risk)
      uniMatrix = matrix(0,nrow = Nrisk*2, ncol = N)
      uniBound = matrix(0,nrow = Nrisk*2, ncol = 1)
      for(i in 1:Nrisk){
        uniMatrix[i,1:tarNum] = risk[[i]]$exposure
        uniMatrix[i+Nrisk,1:tarNum] = 0 - risk[[i]]$exposure
        uniBound[i,1] = risk[[i]]$supThres
        uniBound[i+Nrisk,1] = 0 - risk[[i]]$infThres
      }
      EnlargeUniMatrix = rbind(EnlargeUniMatrix,uniMatrix)
      EnlargeUniBound = rbind(EnlargeUniBound, uniBound)
    }
    
    EnlargeEquaMatrix[1,1:tarNum] = 1
    EnlargeEquaArray[1,1] = totalWeight
    
    EnlargeUniMatrix[1:tarNum+1,1:tarNum] = diag(tarNum)
    EnlargeUniMatrix[1:tarNum+1,(tarNum+1):(2*tarNum)] = -diag(tarNum)
    EnlargeUniMatrix[1:tarNum+1,(2*tarNum+1):(3*tarNum)] = diag(tarNum)
    EnlargeEquaArray[1:tarNum+1] = lastWeight
    
    if(controlLeverage){
      EnlargeUniMatrix[(tarNum+2):(2*tarNum+1),1:tarNum] = diag(tarNum)
      EnlargeUniMatrix[(tarNum+2):(2*tarNum+1),(3*tarNum+1):(4*tarNum)] = -diag(tarNum)
      EnlargeUniMatrix[(tarNum+2):(2*tarNum+1),(4*tarNum+1):(5*tarNum)] = diag(tarNum)
    }
    
    EnlargeInfWeight <- matrix(0, nrow = N, ncol = 1)
    EnlargeInfWeight[1:tarNum] = infWeight
    EnlargeSupWeight = matrix(Inf, nrow = N, ncol = 1)
    EnlargeSupWeight[1:tarNum] = supWeight
    
    if(riskLambda != 0){
      prob = mosek_qptoprob(EnlargeMatrix, -EnlargeAlpha/riskLambda, EnlargeUniMatrix, EnlargeUniBound, 
                            EnlargeEquaMatrix, EnlargeEquaArray, as.numeric(EnlargeInfWeight), as.numeric(EnlargeSupWeight))
    }else{
      prob = mosek_lptoprob(-EnlargeAlpha,EnlargeUniMatrix,EnlargeUniBound,EnlargeEquaMatrix,EnlargeEquaArray,
                            as.numeric(EnlargeInfWeight), as.numeric(EnlargeSupWeight))
    }
    
    result <- mosek(prob, list(verbose= ifelse(verbose,4,1)))
    weight <- result$sol$itr$xx[1:tarNum]
    return(list(weight = weight, status = result$sol$itr$solsta))
  }else if(optMode = "quantile"){
    alpha = matrix(alpha)
    alphaOrder = order(alpha)
    tarNum = dim(alpha)[1]
    holdNum = floor(tarNum/quantileNum)
    shortPos = head(alphaOrder,holdNum)
    longPos = tail(alphaOrder, holdNum)
    weight = rep(0, tarNum)
    
    if(longOnly){
      v = leverageRate/holdNum
      weight[longPos] = v
    }else{
      v = 0.5*leverageRate/holdNum
      weight[shortPos] = -v
      weight[longPos] = v
    }
    return(list(weigth = weight, status = "SOLVED"))
  }
}
 
optProcess <- function(cfg, verbose = TRUE){
  if(is.null(cfg$optMode)) cfg$optMode = "riskLambda"
  if(is.null(cfg$upWgtLimit)) cfg$upWgtLimit = cfg$singleMaxWeight
  if(is.null(cfg$dnWgtLimit)) cfg$dnWgtLimit = 0-cfg$singleMaxWeight
  if(is.null(cfg$totalTurnover)) cfg$totalTurnover = Inf
  
  tarWeight0 = cfg$lastWeight
  tarWeight = tarWeight0
  
  result <- optWeight(alpha = cfg$alpha, covMatrix = cfg$covMatrix, risk = cfg$risk, riskLambda = cfg$riskLambda,
                      lastWeight = cfg$lastWeight, dnWgtLimit = cfg$dnWgtLimit, longOnly = cfg$longOnly,
                      longSideWgt = cfg$longSideWgt, shortSideWgt = cfg$shortSideWgt, optMode = cfg$optMode,
                      quantileNum = cfg$quantileNum, totalTurnover = cfg$totalTurnover, 
                      leverageRate = cfg$leverageRate, tc = cfg$tc, verbose = verbose)
  if(grepl("INFEASIBLE", result$status)){
    print("Warning: No feasible solution")
    
    # try loose turnover
    if(is.finite(cfg$totalTurnover)){
      cfg$totalTurnover = cfg$totalTurnover*10
      print("Trying to loose turnover")
      result <- optWeight(alpha = cfg$alpha, covMatrix = cfg$covMatrix, risk = cfg$risk, riskLambda = cfg$riskLambda,
                          lastWeight = cfg$lastWeight, dnWgtLimit = cfg$dnWgtLimit, longOnly = cfg$longOnly,
                          longSideWgt = cfg$longSideWgt, shortSideWgt = cfg$shortSideWgt, optMode = cfg$optMode,
                          quantileNum = cfg$quantileNum, totalTurnover = cfg$totalTurnover, 
                          leverageRate = cfg$leverageRate, tc = cfg$tc, verbose = verbose)
      if(!grepl("INFEASIBLE", result$status)){
        tarWeight <- tarWeight0 + 0.1*(result$weight - tarWeight0)
      }else{
        print("Error: No feasible solution")
      }
    }
  }else{
    tarWeight = result$weight
  }
  
  print(paste("Turnover would be ", round(sum(abs(tarWeight-tarweight0))/sum(1e-10 + abs(tarWeight0)),2)))
  tarWeight
}
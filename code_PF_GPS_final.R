#############################################################################
###                                                                       ###
### Filtrage Particulaire : GPS et fusion de données (capteurs multiples) ###
###                                                                       ###
#############################################################################

install.packages("SDMTools")
install.packages("plyr")
install.packages("hydroGOF")
library(SDMTools)
library(MASS)
library(plyr)
library(hydroGOF)

########################################
# Dictionnaire de fonctions génériques #
########################################


# Resampling #
##############

# Méthodes de tirages (systématique, stratifiée, etc.) #
systematicRS<-function(N,W)
{ 
  sW<-length(W)
  vectNbParticle<-rep(0,sW)
  i<-1
  k<-1
  cs<-0
  u<-runif(1,0,1/N)
  while(k < sW + 1)
  {
    cs<-cs + W[k]
    while(u < cs || ((i == N) && (u <= cs)))
    {
      vectNbParticle[k] = vectNbParticle[k] + 1
      i<-i+1
      u<-u + (1/N)
    }
    k<-k + 1
  } 
  return(vectNbParticle)
}

multinomialRS<-function(N,W)
{
  vectNbParticle<-rmultinom(1,N,W)
  return(vectNbParticle)
}

stratifiedRS<-function(N,W)
{
  sW<-length(W)
  vectNbParticle<-rep(0,sW)
  i<-1
  k<-1
  cs<-0
  u<-runif(1,0,1/N)
  while(k < sW + 1)
  {
    cs<-cs + W[k]
    while(u < cs || ((i == N) && (u <= cs)))
    {
      vectNbParticle[k] = vectNbParticle[k] + 1
      i<-i+1
      u<-runif(1,((i-1)/N),(i/N)) 
    }
    k<-k + 1
  } 
  return(vectNbParticle)
}

residualRS<-function(N,W)
{
  nbTilde<-trunc(N*W)
  WBar<-W - (nbTilde/N)
  nbBar<-rmultinom(1,N,WBar)
  vectNbParticle<-nbTilde + nbBar
  return(vectNbParticle)
}


# Choix des ancètres #
particleResampling<-function(resampleMethod,N,W)
{
  nbParticle<-rep(0,N)
  nbParticle<-resampleMethod(N,W)
  
  vectIndice<-rep(0,0)
  for(i in 1:N)
  {
    ind<-rep(i,nbParticle[i])
    vectIndice<-c(vectIndice,ind)
  }
  
  return(vectIndice)
}

# normalisation des poids #
###########################

normWeightRoutine<-function(vectWeight,topLogWeight)
{
  N<-length(vectWeight)
  vectNormWeight<-rep(0,N)
  if (topLogWeight == 1)
  {
    m<-max(vectWeight)
    w<-rep(0,N)
    sum<-0
    for (i in 1:N)
    {
      w[i]<-exp(vectWeight[i]- m)
      sum<-sum + w[i]
    }
  }
  else
  {
    w<-vectWeight
    sum<-sum(w)
  }
  
  for (i in 1:N)
  {
    vectNormWeight[i]<- w[i]/sum
  }
  return(vectNormWeight)
}

# Algorithme de filtrage #
##########################

ParticleFilteringAsync<-function(dataInit,matData,initRoutine,samplingRoutine,weightRoutine,RSMethod,stateResampling,nbStateVar,nbSensorFaulty,N,listModelParameters,topLogWeight,ESSBoundary,
                                 nbIt = 0, itDeb = 1, topInitDone = 0, matSampleOld = 0, vectWeightOld = 0, meanUKF = 0, varUKF = 0)
{
  library(SDMTools)
  if(nbIt == 0 || nbIt > nrow(matData))
  {
    nbIt<-nrow(matData)
  }
  
  progress_barPF<- create_progress_bar("text")
  progress_barPF$init(nbIt - itDeb + 1)
  
  # setting
  situation<-matData[,ncol(matData)]
  data<-matData[,-ncol(matData)]
  data<-data[,-1]
  deltaT<-matData[,1]
  PFResultMean<-matrix(0,(nbIt- itDeb + 1),nbStateVar+(nbSensorFaulty*2))
  PFResultMeanPrev<-matrix(0,(nbIt- itDeb + 1),nbStateVar+(nbSensorFaulty*2))
  PFResultVar<-matrix(0,(nbIt- itDeb + 1),nbStateVar+(nbSensorFaulty*2))
  PFResultVarPrev<-matrix(0,(nbIt- itDeb + 1),nbStateVar+(nbSensorFaulty*2))

  if (topInitDone == 0)
  {
  # init
    result<-initRoutine(dataInit,nbStateVar,nbSensorFaulty,listModelParameters,N)
    matSample<-result$initSample
    vectWeight<-rep((topLogWeight*log(1/N))+((1-topLogWeight)*(1/N)),N)
  }
  else
  {
    matSample<-matSampleOld
    vectWeight<-vectWeightOld
    result = list(meanUKF = meanUKF,varUKF = varUKF)
  }
  
  vectWeightNorm<-normWeightRoutine(vectWeight,topLogWeight)

  # iteration
  for(t in itDeb:nbIt)
  {
    # sampling
    matSamplePrev<-matSample
    result<-samplingRoutine(N,matSamplePrev,data[t,],deltaT[t],nbStateVar,nbSensorFaulty,situation[t],listModelParameters,result$meanUKF,result$varUKF,PFResultMean[t,])
    matSample<-result$matSample
    
    # weight
    vectWeight<-weightRoutine(N,vectWeight,matSample,matSamplePrev,data[t,],situation[t],listModelParameters,nbSensorFaulty,nbStateVar,deltaT[t])

    vectWeightNorm<-normWeightRoutine(vectWeight,topLogWeight)
    
    # output result
    for(j in 1:nbStateVar)
    {
      PFResultMean[(t-itDeb + 1),j]<-wt.mean(matSample[((nbSensorFaulty*3)+j),],vectWeightNorm)
      PFResultMeanPrev[(t-itDeb + 1),j]<-mean(matSample[((nbSensorFaulty*3)+j),])
      PFResultVar[(t-itDeb + 1),j]<-wt.var(matSample[((nbSensorFaulty*3)+j),],vectWeightNorm)
      PFResultVarPrev[(t-itDeb + 1),j]<-var(matSample[((nbSensorFaulty*3)+j),])
    }
    if (nbSensorFaulty != 0 && listModelParameters$topCensorFailDetection == 0)
    {
      for(j in 0:(nbSensorFaulty-1))
      {
        PFResultMean[(t-itDeb + 1),nbStateVar+j+1]<-wt.mean(matSample[(3*j)+1,],vectWeightNorm)
        PFResultMeanPrev[(t-itDeb + 1),nbStateVar+j+1]<-mean(matSample[(3*j)+1,])
        PFResultVar[(t-itDeb + 1),nbStateVar+j+1]<-wt.var(matSample[(3*j)+1,],vectWeightNorm)
        PFResultVarPrev[(t-itDeb + 1),nbStateVar+j+1]<-var(matSample[(3*j)+1,])
        PFResultMean[(t-itDeb + 1),nbStateVar+j+2]<-wt.mean(matSample[(3*j)+2,],vectWeightNorm)
        PFResultMeanPrev[(t-itDeb + 1),nbStateVar+j+2]<-mean(matSample[(3*j)+2,])
        PFResultVar[(t-itDeb + 1),nbStateVar+j+2]<-wt.var(matSample[(3*j)+2,],vectWeightNorm)
        PFResultVarPrev[(t-itDeb + 1),nbStateVar+j+2]<-var(matSample[(3*j)+2,])
      }
    }
    
    # resampling
    ESS<-1/sum(vectWeightNorm^2)
    if (ESS < ESSBoundary)
    {
      vectNewInd<-particleResampling(RSMethod,N,vectWeightNorm)
      result<-stateResampling(vectNewInd,matSample,result$meanUKF,result$varUKF,N)
      matSample<-result$newMatSample
      vectWeight<-rep(1/N,N)
    }
      
    progress_barPF$step()
    
  }
  
  PFResult<-list(mean = PFResultMean, varMC = PFResultVar, meanPrev = PFResultMeanPrev, varMCPrev = PFResultVarPrev, 
                 matSample = matSample, vectWeight = vectWeight, meanUKF = result$meanUKF, varUKF = result$varUKF)
  return(PFResult)
}


# fonctions outils divers #
###########################

normDensityLog<-function(x,mu,matVar)
{
  N<-length(x)
  mx<-as.matrix(x)
  mmu<-as.matrix(mu)
  p<-t(mx-mmu)%*%solve(matVar)%*%(mx-mmu)
  ln<-((-0.5*N)*log(2*pi))-(0.5*log(det(matVar)))-(0.5*p[1,1])
  return(ln)
}

computationMSELatLong<-function(dataGPSClean,MCmean,dataFinal)
{
  estimX<-rep(0,0)
  estimY<-rep(0,0)
  
  for (i in 2:nrow(MCmean))
  {
    if (dataFinal[i-1,6] == 2 || dataFinal[i-1,6] == 3)
    {
      estimX<-c(estimX,MCmean[i,1])
      estimY<-c(estimY,MCmean[i,2])
    }
  }
  
  result<-rep(0,2)
  
  result[1]<-mse(estimX,dataGPSClean[1:length(estimX),2])
  result[2]<-mse(estimY,dataGPSClean[1:length(estimY),3])
  
  return(result)
}


##########################################
# Application guidage GPS multi-capteurs #
##########################################

# on suppose que les erreurs dans le modèles des observations sont de variance fixées et non corrélées #

# Dictionnaire : importation et mise en forme des données #
###########################################################

### alteration des données pour test UKF ###

dataAlteration<-function(dataGPS,matFailWindow)
{
  alteredDataGPS<-dataGPS
  if (max(matFailWindow[,2]) >= nrow(alteredDataGPS))
  {
    return ("error")
  }
  else
  {
      a<-0.35
      b<-0.55
      for (i in matFailWindow[1,1]:matFailWindow[1,2])
      {
        alteredDataGPS[i,3]<-alteredDataGPS[i,3] + a*((2*rbinom(1,1,0.2))-1)
      }
      for (i in matFailWindow[2,1]:matFailWindow[2,2])
      {
        alteredDataGPS[i,2]<-alteredDataGPS[i,2] + b*((2*rbinom(1,1,0.5))-1)
      }
  }
  return(alteredDataGPS)
}

### Préparation table de données finale ###

dataPreparation<-function(dataGPS,dataSpeedSensor,dataManager)
{
  dataManagerWithoutLaser<-matrix(0,1,4)
  for (i in 1:nrow(dataManager))
  {
    if ((dataManager[i,2] != 3) && (dataManager[i,1] != dataManagerWithoutLaser[nrow(dataManagerWithoutLaser),1]))
    {
      if (dataManager[i,2] == 1)
      {
        v<-c(dataManager[i,1],dataManager[i,2],0,dataManager[i,3])
      }
      else if (dataManager[i,2] == 2)
      {
        v<-c(dataManager[i,1],dataManager[i,2],dataManager[i,3],0)
      }
      dataManagerWithoutLaser<-rbind(dataManagerWithoutLaser,v)
      dataManagerWithoutLaser[nrow(dataManagerWithoutLaser),2]<- 3 - dataManagerWithoutLaser[nrow(dataManagerWithoutLaser),2]
    }
    else if ((dataManager[i,2] != 3) && (dataManager[i,1] == dataManagerWithoutLaser[nrow(dataManagerWithoutLaser),1]))
    {
      dataManagerWithoutLaser[nrow(dataManagerWithoutLaser),2]<- 3
      if (dataManager[i,2] == 1)
      {
        dataManagerWithoutLaser[nrow(dataManagerWithoutLaser),4]<-dataManager[i,3]
      }
      else if (dataManager[i,2] == 2)
      {
        dataManagerWithoutLaser[nrow(dataManagerWithoutLaser),3]<-dataManager[i,3]
      }
    }
  }
  dataFinal<-matrix(0,nrow(dataManagerWithoutLaser),6)
  
  for (i in 2:nrow(dataManagerWithoutLaser))
  {
    dataFinal[i,1]<-dataManagerWithoutLaser[i,1]-dataManagerWithoutLaser[(i-1),1]
    dataFinal[i,6]<-dataManagerWithoutLaser[i,2]
    if (dataManagerWithoutLaser[i,4] == 0)
    {
      dataFinal[i,2]<-NA
      dataFinal[i,3]<-NA
    }
    else
    {
      dataFinal[i,2]<-dataGPS[dataManagerWithoutLaser[i,4],2]
      dataFinal[i,3]<-dataGPS[dataManagerWithoutLaser[i,4],3]
    }
    if (dataManagerWithoutLaser[i,3] == 0)
    {
      dataFinal[i,4]<-NA
      dataFinal[i,5]<-NA
    }
    else
    {
      dataFinal[i,4]<-dataSpeedSteer[dataManagerWithoutLaser[i,3],2]
      dataFinal[i,5]<- - dataSpeedSteer[dataManagerWithoutLaser[i,3],3] # because angle are opposite in our model with respect to the data
    }
  }
  dataFinal<-dataFinal[-1,]
  
  dataFinal[1,1]<-dataFinal[2,1]
  
  return(dataFinal)
}


# Dictionnaire : dynamique des varibles d'états #
#################################################

stateTransition<-function(stateVariablePrev,deltaT,a,b,L)
{
  stateVariable<-rep(0,length(stateVariablePrev))
  
  stateVariable[1]<-stateVariablePrev[1] + (deltaT*((stateVariablePrev[3]*cos(stateVariablePrev[4]))+(((b*cos(stateVariablePrev[4]))-(a*sin(stateVariablePrev[4])))*
                                                                                                      (stateVariablePrev[3]/L)*tan(stateVariablePrev[5]))))
  stateVariable[2]<-stateVariablePrev[2] + (deltaT*((stateVariablePrev[3]*sin(stateVariablePrev[4]))+(((b*sin(stateVariablePrev[4]))+(a*cos(stateVariablePrev[4])))*
                                                                                                      (stateVariablePrev[3]/L)*tan(stateVariablePrev[5]))))
  stateVariable[3]<-stateVariablePrev[3]
  stateVariable[4]<-stateVariablePrev[4] + (deltaT*(stateVariablePrev[3]/L)*tan(stateVariablePrev[5]))
  stateVariable[5]<-stateVariablePrev[5] + (deltaT*stateVariablePrev[6])
  stateVariable[6]<-stateVariablePrev[6]
  
  return(stateVariable)
}

# Dictionnaire : filtre Bootstrap #
###################################

### initialisation ###

initRoutineBootstrap<-function(dataInit,nbStateVar,nbSensorFaulty,listModelParameters,N)
{
  initSample<-matrix(0,nbStateVar,N)
  
  matQ<-listModelParameters$matVarStateModel
  
  initSample[1,]<-rnorm(N,dataInit[1],sqrt(0.2)) # X
  initSample[2,]<-rnorm(N,dataInit[2],sqrt(0.2)) # Y
  initSample[3,]<-rlnorm(N,dataInit[3]/(1+(tan(dataInit[4])*(H/L))),sqrt(0.2)) # V
  initSample[4,]<-rnorm(N,pi,sqrt(0.2)) # heading psi
  initSample[5,]<-rnorm(N,dataInit[4],sqrt(0.2)) # steering
  initSample[6,]<-rnorm(N,0,sqrt(0.2)) # deriv steering
  
  meanUKF<-0
  varUKF<-0
  
  result<-list(initSample = initSample, meanUKF = meanUKF, varUKF = varUKF)
  return(result)
}

### modèle d'état ###

# (situation = 1 : Speed/Stearing seul) ;  (situation = 2 : GPS seul) ; (situation = 3 : tout) #

# sampling routine #
samplingRoutineBootstrap<-function(N,matSamplePrev,vectObs,deltaT,nbStateVar,nbSensorFaulty,situation,listModelParameters,meanUKF,varUKF,PFResultMean)
{ 
  matQ<-listModelParameters$matVarStateModel
  L<-listModelParameters$L
  a<-listModelParameters$a
  b<-listModelParameters$b
  
  matSample<-matrix(0,nbStateVar,N)
  
  for (i in 1:N)
  {
    mu<-mvrnorm(1,c(0,0),matQ)
    matSample[,i]<-stateTransition(matSamplePrev[,i],deltaT,a,b,L)
    matSample[3,i]<-matSample[3,i] + (deltaT*mu[1])
    matSample[6,i]<-matSample[6,i] + (deltaT*mu[2])
  }
  result<-list(matSample = matSample, meanUKF = meanUKF, varUKF = varUKF)
  return(result)
}

# weight routine #
weightRoutineBootstrap<-function(N,oldWeight,matSample,matSamplePrev,vectObs,situation,listModelParameters,nbSensorFaulty,nbStateVar,deltaT)
{
  matR<-listModelParameters$matVarObs
  H<-listModelParameters$H
  L<-listModelParameters$L
  
  weight<-rep(0,N)
  
  if (situation == 1)
  {
    mu<-rep(0,2)
    matR<-matR[c(3,4),c(3,4)]
    for (i in 1:N)
    {
      mu[1]<-(matSample[3,i]*(1+(tan(matSample[5,i])*(H/L))))
      mu[2]<-matSample[5,i]
      
      weight[i]<-normDensityLog(vectObs[c(3,4)],mu,matR)
    }
  }
  else if (situation == 2)
  {
    mu<-rep(0,2)
    matR<-matR[c(1,2),c(1,2)]
    for (i in 1:N)
    {
      mu[1]<-matSample[1,i]
      mu[2]<-matSample[2,i]
      
      weight[i]<-normDensityLog(vectObs[c(1,2)],mu,matR)
    }
  }
  else
  {
    mu<-rep(0,4)
    for (i in 1:N)
    {
      mu[1]<-matSample[1,i]
      mu[2]<-matSample[2,i]
      mu[3]<-(matSample[3,i]*(1+(tan(matSample[5,i])*(H/L))))
      mu[4]<-matSample[5,i]
      weight[i]<-normDensityLog(vectObs,mu,matR)
    }
  }
    
  return(weight)
}

# state resampling #
stateResamplingBootstrap<-function(vectNewInd,matSample,meanUKF,varUKF,N)
{
  newMatSample<-matrix(0,nrow(matSample),N)
  for (i in 1:N)
  {
    newMatSample[,i]<-matSample[,vectNewInd[i]]
  }
  result = list(newMatSample = newMatSample, meanUKF = meanUKF, varUKF = varUKF)
}


# Dictionnaire : filtre UKF Rao Blackwellisé #
##############################################

### initialisation ###

initRoutineRB<-function(dataInit,nbStateVar,nbSensorFaulty,listModelParameters,N)
{
  initSample<-matrix(0,((3*nbSensorFaulty)+nbStateVar),N)
  
  lambda<-listModelParameters$hyperParaVar
  matQ<-listModelParameters$matVarStateModel
  H<-listModelParameters$H
  L<-listModelParameters$L
  
  initSample[1,]<-rep(1,N)
  initSample[3,]<-rlnorm(N,0,sqrt(lambda))
  u<-runif(N,0.98,0.99)
  for (i in 1:N)
  {
    initSample[2,i]<-u[i]
  }
  initSample[4,]<-rnorm(N,dataInit[1],sqrt(0.02)) # X
  initSample[5,]<-rnorm(N,dataInit[2],sqrt(0.02)) # Y
  initSample[6,]<-rlnorm(N,dataInit[3]/(1+(tan(dataInit[4])*(H/L))),sqrt(0.002)) # V
  initSample[7,]<-rnorm(N,pi,sqrt(0.002)) # heading psi
  initSample[8,]<-rnorm(N,dataInit[4],sqrt(0.002))
  initSample[9,]<-rnorm(N,0,sqrt(0.02)) # deriv steering 
  
  meanUKF<-initSample[4:9,]
  varUKF<-matrix(0,(nbStateVar^2),N)
  for (i in 0:(nbStateVar-1))
  {
    varUKF[(nbStateVar*i)+(1+i),]<-rlnorm(N,0.02,0.002)
  }
  result<-list(initSample = initSample, meanUKF = meanUKF, varUKF = varUKF)
  return(result)
}

initRoutineUKF<-function(dataInit,nbStateVar,nbSensorFaulty,listModelParameters,N)
{
  initSample<-matrix(0,((3*nbSensorFaulty)+nbStateVar),N)
  
  lambda<-listModelParameters$hyperParaVar
  matQ<-listModelParameters$matVarStateModel
  H<-listModelParameters$H
  L<-listModelParameters$L
  
  initSample[1,]<-rep(0,N)
  initSample[2,]<-rep(0,N)
  initSample[3,]<-rep(0,N)
  initSample[4,]<-rnorm(N,dataInit[1],sqrt(0.02)) # X
  initSample[5,]<-rnorm(N,dataInit[2],sqrt(0.02)) # Y
  initSample[6,]<-rlnorm(N,dataInit[3]/(1+(tan(dataInit[4])*(H/L))),sqrt(0.002)) # V
  initSample[7,]<-rnorm(N,pi,sqrt(0.002)) # heading psi
  initSample[8,]<-rnorm(N,dataInit[4],sqrt(0.002))
  initSample[9,]<-rnorm(N,0,sqrt(0.002)) # deriv steering 
  
  meanUKF<-matrix(0,nbStateVar,N)
  meanUKF[1:nbStateVar,]<-initSample[4:9,]
  varUKF<-matrix(0,(nbStateVar^2),N)
  for (i in 0:(nbStateVar-1))
  {
    varUKF[(nbStateVar*i)+(1+i),]<-rlnorm(N,0.02,0.002)
  }
  result<-list(initSample = initSample, meanUKF = meanUKF, varUKF = varUKF)
  return(result)
}

### modèle d'états ###

# (situation = 1 : Speed/Stearing seul) ;  (situation = 2 : GPS seul) ; (situation = 3 : tout) #

topNonFailSamplingRB<-function(vectSamplePrev,vectObs,spreadFail,matVarObs,deltaT,a,b,L,nbSensorFaulty,nbStateVar)
{
  # bernoulli
  alphaPrev<-vectSamplePrev[2]
  condExp<-stateTransition(vectSamplePrev[((3*nbSensorFaulty) +1):((3*nbSensorFaulty) + nbStateVar)],deltaT,a,b,L)
  p1<-(alphaPrev*dnorm(vectObs[1],condExp[1],sqrt(matVarObs[1,1]))*dnorm(vectObs[2],condExp[2],sqrt(matVarObs[2,2])))
  p0<-(1-alphaPrev)*dunif(vectObs[1],vectSamplePrev[((3*nbSensorFaulty) +1)]-spreadFail,vectSamplePrev[((3*nbSensorFaulty) +1)]+spreadFail)*
    dunif(vectObs[2],vectSamplePrev[((3*nbSensorFaulty) +2)]-spreadFail,vectSamplePrev[((3*nbSensorFaulty) +2)]+spreadFail)
  if (p0 == 0 && p1 == 0)
  {
    sample<-0
  }
  else
  {  
    sample<-rbinom(1,1,(p1/(p0+p1)))
  } 
  return(sample)
}

alphaSamplingRB<-function(alphaPrev,topNonFail,hyperParaPrev)
{
  # beta
  hyperParaPrevModified<-hyperParaPrev + 1
  alphaPrevModified<-(hyperParaPrev*alphaPrev/hyperParaPrevModified) + (topNonFail/hyperParaPrevModified)
  sample<-rbeta(1,(hyperParaPrevModified*(1-alphaPrevModified)),(hyperParaPrevModified*alphaPrevModified))
  if (sample > (1 - 0.0001))
  {
    sample<-1 - 0.0001
  }
  if (sample < 0.0001)
  {
    sample<-0.0001
  }
  return(sample)
}

hyperParaSamplingRB<-function(alpha,hyperParaAlphaPrev,alphaPrev,varHyperPara)
{
  # Independant Metropolis Hasting
  sample<-hyperParaAlphaPrev
  for (mht in 1:1000)
  {
    sigmaLogNorm<-rlnorm(1,log(hyperParaAlphaPrev),sqrt(varHyperPara))
    r<-min(1,(dbeta(alpha,(sigmaLogNorm*(1-alphaPrev)),(sigmaLogNorm*alphaPrev))/dbeta(alpha,(sample*(1-alphaPrev)),(sample*alphaPrev))))
    u<-rbinom(1,1,r)
    sample<-(u*sigmaLogNorm)+((1-u)*sample)
  }
  return(sample)
}

UKF<-function(meanUKF,varUKF,vectSamplePrev,vectObs,situation,matQ,matR,nbStateVar,a,b,L,H,deltaT,nbSensorFaulty,topObsFail,topCensorFailDetection)
{
  sizeQ<-nrow(matQ)
  matVarUKF<-matrix(varUKF,nbStateVar,nbStateVar,byrow = F)
  
  topNonFail<- (topCensorFailDetection*(1 - topObsFail)) + ((1-topCensorFailDetection)*vectSamplePrev[1])
  
  if (situation == 1 || (situation == 3 && topNonFail == 0))
  {
    matR<-matR[3:4,3:4]
    sizeR<-nrow(matR)
    obs<-vectObs[3:4]
    oldMean<-rep(0,(nbStateVar + sizeQ + sizeR))
    oldMean[1:nbStateVar]<-meanUKF
  }
  else if (situation == 3 && topNonFail == 1)
  {
    sizeR<-nrow(matR)
    obs<-vectObs
    oldMean<-rep(0,(nbStateVar + sizeQ + sizeR))
    oldMean[1:nbStateVar]<-meanUKF
  }
  else if (situation == 2 && topNonFail == 1)
  {
    matR<-matR[1:2,1:2]
    sizeR<-nrow(matR)
    obs<-vectObs[1:2]
    oldMean<-rep(0,(nbStateVar + sizeQ + sizeR))
    oldMean[1:nbStateVar]<-meanUKF
  }
  else
  {
    sizeR<-0
    oldMean<-rep(0,(nbStateVar + sizeQ + sizeR))
    oldMean[1:nbStateVar]<-meanUKF
  }
  
  varOld<-matrix(0,(nbStateVar + sizeQ + sizeR),(nbStateVar + sizeQ + sizeR))
  varOld[1:nbStateVar,1:nbStateVar]<-matVarUKF
  varOld[(nbStateVar + 1):(nbStateVar + sizeQ) ,(nbStateVar + 1):(nbStateVar + sizeQ)]<-matQ
  if (situation != 2 || (topNonFail != 0))
  {
    varOld[(nbStateVar + sizeQ + 1):(nbStateVar + sizeQ + sizeR) ,(nbStateVar + sizeQ + 1):(nbStateVar + sizeQ + sizeR)]<-matR
  }
  cholVar<-t(chol(varOld))
  
  # weight UKF definition
  alpha<-0.001
  beta<-2
  M<-nbStateVar + sizeQ + sizeR
  lambda<-(alpha*alpha*M) - M
  weightMeanUKF<-rep((0.5/(M + lambda)),(2*M + 1))
  weightMeanUKF[(2*M) + 1]<-lambda/(lambda + M)
  weightCovUKF<-weightMeanUKF
  weightCovUKF[(2*M) + 1]<-weightCovUKF[(2*M) + 1] + (1 - (alpha*alpha) + beta)
  eta<-sqrt(lambda + M)
  
  sigmaPoint<-matrix(0,nbStateVar + sizeQ + sizeR,((2*M) + 1))
  Fsigma<-matrix(0,nbStateVar,((2*M)+1))
  sigmaPoint[,((2*M)+1)]<-oldMean
  Fsigma[,((2*M)+1)]<-stateTransition(sigmaPoint[(1:nbStateVar),((2*M)+1)],deltaT,a,b,L)
  Fsigma[3,((2*M)+1)]<-Fsigma[3,((2*M)+1)] + (deltaT*sigmaPoint[(nbStateVar + 1),((2*M)+1)])
  Fsigma[6,((2*M)+1)]<-Fsigma[6,((2*M)+1)] + (deltaT*sigmaPoint[(nbStateVar + 2),((2*M)+1)])
  for (i in 1:M)
  {
    sigmaPoint[,i]<-oldMean + (eta*cholVar[,i])
    Fsigma[,i]<-stateTransition(sigmaPoint[(1:nbStateVar),i],deltaT,a,b,L)
    Fsigma[3,i]<-Fsigma[3,i] + (deltaT*sigmaPoint[(nbStateVar + 1),i])
    Fsigma[6,i]<-Fsigma[6,i] + (deltaT*sigmaPoint[(nbStateVar + 2),i])
    sigmaPoint[,(i + M)]<-oldMean - (eta*cholVar[,i])
    Fsigma[,(i + M)]<-stateTransition(sigmaPoint[(1:nbStateVar),(i + M)],deltaT,a,b,L)
    Fsigma[3,(i + M)]<-Fsigma[3,(i + M)] + (deltaT*sigmaPoint[(nbStateVar + 1),(i + M)])
    Fsigma[6,(i + M)]<-Fsigma[6,(i + M)] + (deltaT*sigmaPoint[(nbStateVar + 2),(i + M)])
  }
  
  mux<-rep(0,nbStateVar)
  for (i in 1:((2*M)+1))
  {
    mux <- mux + weightMeanUKF[i]*Fsigma[,i]
  }
  
  pxx<-matrix(0,nbStateVar,nbStateVar)
  for (i in 1:(2*M))
  {
    pxx <- pxx + (weightCovUKF[i]*((as.matrix((Fsigma[,i] - Fsigma[,((2*M)+1)])))%*%(t(as.matrix(Fsigma[,i] - Fsigma[,((2*M)+1)]))))) # version modifiée pour rendre matrice cov definie positive
  }
  
  if (situation == 3 && topNonFail == 1)
  {
    HFsigma<-matrix(0,sizeR,((2*M)+1))
    for (i in 1:((2*M)+1))
    {
      HFsigma[,i]<-c(Fsigma[1,i]+ sigmaPoint[(nbStateVar + sizeQ + 1),i],Fsigma[2,i] + sigmaPoint[(nbStateVar + sizeQ + 2),i], ((1 + (tan(Fsigma[5,i])*(H/L)))*Fsigma[3,i]) + 
                       sigmaPoint[(nbStateVar + sizeQ + 3),i],Fsigma[5,i] + sigmaPoint[(nbStateVar + sizeQ + 4),i])
    }
  }
  else if (situation == 1 || (situation == 3 && topNonFail == 0))
  {
    HFsigma<-matrix(0,sizeR,((2*M)+1))
    for (i in 1:((2*M)+1))
    {
      HFsigma[,i]<-c(((1 + (tan(Fsigma[5,i])*(H/L)))*Fsigma[3,i]) + sigmaPoint[(nbStateVar + sizeQ + 1),i],Fsigma[5,i] + sigmaPoint[(nbStateVar + sizeQ + 2),i])
    }
  }
  else if ((situation == 2) && (topNonFail == 1))
  {
    HFsigma<-matrix(0,sizeR,((2*M)+1))
    for (i in 1:((2*M)+1))
    {
      HFsigma[,i]<-c(Fsigma[1,i]+ sigmaPoint[(nbStateVar + sizeQ + 1),i],Fsigma[2,i] + sigmaPoint[(nbStateVar + sizeQ + 2),i])
    }
  }
  
  if ((situation != 2) || (topNonFail != 0))
  {
    muy<-rep(0,sizeR)
    for (i in 1:((2*M)+1))
    {
      muy <- muy + weightMeanUKF[i]*HFsigma[,i]
    }
    pxy<-matrix(0,nbStateVar,sizeR)
    pyy<-matrix(0,sizeR,sizeR)
    for (i in 1:(2*M))
    {
      pxy <- pxy + (weightCovUKF[i]*((as.matrix((Fsigma[,i] - Fsigma[,((2*M)+1)])))%*%(t(as.matrix(HFsigma[,i] - HFsigma[,((2*M)+1)]))))) # version modifiée pour rendre matrice cov definie positive
      pyy <- pyy + (weightCovUKF[i]*((as.matrix((HFsigma[,i] - HFsigma[,((2*M)+1)])))%*%(t(as.matrix(HFsigma[,i] - HFsigma[,((2*M)+1)])))))
    }
  }
  
  if ((situation == 2) && (topNonFail == 0))
  {
    meanUKF[1:nbStateVar]<-mux
    matVarUKF<-pxx
  }
  else
  {
    mux<-as.matrix(mux)
    muy<-as.matrix(muy)
    meanUKF<- mux + (pxy%*% (solve(pyy))%*%(as.matrix(obs) - muy))
    matVarUKF<- pxx - (pxy%*%(solve(pyy))%*%(t(pxy)))
  }
  varUKF<-rep(0,0)
  for (i in 1:ncol(matVarUKF))
  {
    varUKF<-c(varUKF,matVarUKF[,i])
  }
  
  resultUKF<-list(meanUKF = meanUKF, varUKF = varUKF)
  return(resultUKF)
}


# sampling routine #
samplingRoutineRB<-function(N,matSamplePrev,vectObs,deltaT,nbStateVar,nbSensorFaulty,situation,listModelParameters,meanUKF,varUKF,PFResultMean)
{
  matR<-listModelParameters$matVarObs
  matQ<-listModelParameters$matVarStateModel
  lambda<-listModelParameters$hyperParaVar
  spread<-listModelParameters$spreadFail
  a<-listModelParameters$a
  b<-listModelParameters$b
  H<-listModelParameters$H
  L<-listModelParameters$L
  topCensorFailDetection<-listModelParameters$topCensorFailDetection
  
  matSample<-matrix(0,nrow(matSamplePrev),ncol(matSamplePrev))
  
  if(situation == 1)
  {
    for (i in 1:N)
    {
      matSample[1,i]<-matSamplePrev[1,i]
      matSample[2,i]<-matSamplePrev[2,i]
      matSample[3,i]<-matSamplePrev[3,i]
      UKFResult<-UKF(meanUKF[,i],varUKF[,i],matSamplePrev[,i],vectObs,situation,matQ,matR,nbStateVar,a,b,L,H,deltaT,nbSensorFaulty,0,topCensorFailDetection)
      meanUKF[,i]<-UKFResult$meanUKF
      varUKF[,i]<-UKFResult$varUKF
      matSample[((3*nbSensorFaulty)+1):((3*nbSensorFaulty) + nbStateVar),i]<-meanUKF[,i]
    }
  }
  else if (situation == 2)
  { 
    for (i in 1:N)
    {
      matSample[1,i]<-topNonFailSamplingRB(matSamplePrev[,i],vectObs,spread,matR,deltaT,a,b,L,nbSensorFaulty,nbStateVar)
      matSample[2,i]<-alphaSamplingRB(matSamplePrev[2,i],matSample[1,i],matSamplePrev[3,i])
      matSample[3,i]<-hyperParaSamplingRB(matSample[2,i],matSamplePrev[3,i],matSamplePrev[2,i],lambda)
      UKFResult<-UKF(meanUKF[,i],varUKF[,i],matSamplePrev[,i],vectObs,situation,matQ,matR,nbStateVar,a,b,L,H,deltaT,nbSensorFaulty,0,topCensorFailDetection)
      meanUKF[,i]<-UKFResult$meanUKF
      varUKF[,i]<-UKFResult$varUKF
      matSample[((3*nbSensorFaulty)+1):((3*nbSensorFaulty) + nbStateVar),i]<-meanUKF[,i]
    }
  }
  else
  {
    for (i in 1:N)
    {
      matSample[1,i]<-topNonFailSamplingRB(matSamplePrev[,i],vectObs,spread,matR,deltaT,a,b,L,nbSensorFaulty,nbStateVar)
      matSample[2,i]<-alphaSamplingRB(matSamplePrev[2,i],matSample[1,i],matSamplePrev[3,i])
      matSample[3,i]<-hyperParaSamplingRB(matSample[2,i],matSamplePrev[3,i],matSamplePrev[2,i],lambda)
      UKFResult<-UKF(meanUKF[,i],varUKF[,i],matSamplePrev[,i],vectObs,situation,matQ,matR,nbStateVar,a,b,L,H,deltaT,nbSensorFaulty,0,topCensorFailDetection)
      meanUKF[,i]<-UKFResult$meanUKF
      varUKF[,i]<-UKFResult$varUKF
      matSample[((3*nbSensorFaulty)+1):((3*nbSensorFaulty) + nbStateVar),i]<-meanUKF[,i]
    }
  }
  
  result<-list(matSample = matSample, meanUKF = meanUKF, varUKF = varUKF)
  return(result)
}

# sampling routine UKF only with Fail detection #
samplingRoutineUKF<-function(N,matSamplePrev,vectObs,deltaT,nbStateVar,nbSensorFaulty,situation,listModelParameters,meanUKF,varUKF,PFResultMean)
{
  matR<-listModelParameters$matVarObs
  matQ<-listModelParameters$matVarStateModel
  lambda<-listModelParameters$hyperParaVar
  spread<-listModelParameters$spreadFail
  a<-listModelParameters$a
  b<-listModelParameters$b
  H<-listModelParameters$H
  L<-listModelParameters$L
  topCensorFailDetection<-listModelParameters$topCensorFailDetection
  riskLevelCensorFailDetection<-listModelParameters$riskLevelCensorFailDetection
  
  matSample<-matrix(0,nrow(matSamplePrev),ncol(matSamplePrev))
  
  topObsFail<-0
  
  if (topCensorFailDetection == 1)
  # optional failure detection based on quadratic innovation for GPS sensor 
  if (situation == 2 || situation == 3)
  {
    S<-matR[1:2,1:2]
    newPred<-stateTransition(meanUKF,deltaT,a,b,L)
    v<-as.matrix(vectObs[1:2] - newPred[1:2])
    chi<-t(v)%*%solve(S)%*%v
    if (chi > qchisq((1 - riskLevelCensorFailDetection), df=2))
    {
      topObsFail<-1
    }
    else
    {
      topObsFail<-0
    }
  }
  
  if(situation == 1)
  {
    for (i in 1:N)
    {
      UKFResult<-UKF(meanUKF[,i],varUKF[,i],matSamplePrev[,i],vectObs,situation,matQ,matR,nbStateVar,a,b,L,H,deltaT,nbSensorFaulty,topObsFail,topCensorFailDetection)
      meanUKF[,i]<-UKFResult$meanUKF
      varUKF[,i]<-UKFResult$varUKF
      matSample[((3*nbSensorFaulty)+1):((3*nbSensorFaulty) + nbStateVar),i]<-meanUKF[,i]
      matSample[1,i]<-1-topObsFail
    }
  }
  else if (situation == 2)
  { 
    for (i in 1:N)
    {
      UKFResult<-UKF(meanUKF[,i],varUKF[,i],matSamplePrev[,i],vectObs,situation,matQ,matR,nbStateVar,a,b,L,H,deltaT,nbSensorFaulty,topObsFail,topCensorFailDetection)
      meanUKF[,i]<-UKFResult$meanUKF
      varUKF[,i]<-UKFResult$varUKF
      matSample[((3*nbSensorFaulty)+1):((3*nbSensorFaulty) + nbStateVar),i]<-meanUKF[,i]
      matSample[1,i]<-1-topObsFail
    }
  }
  else
  {
    for (i in 1:N)
    {
      UKFResult<-UKF(meanUKF[,i],varUKF[,i],matSamplePrev[,i],vectObs,situation,matQ,matR,nbStateVar,a,b,L,H,deltaT,nbSensorFaulty,topObsFail,topCensorFailDetection)
      meanUKF[,i]<-UKFResult$meanUKF
      varUKF[,i]<-UKFResult$varUKF
      matSample[((3*nbSensorFaulty)+1):((3*nbSensorFaulty) + nbStateVar),i]<-meanUKF[,i]
      matSample[1,i]<-1-topObsFail
    }
  }
  
  result<-list(matSample = matSample, meanUKF = meanUKF, varUKF = varUKF)
  return(result)
}


### calcul des poids ###

weightRoutineRB<-function(N,oldWeight,matSample,matSamplePrev,vectObs,situation,listModelParameters,nbSensorFaulty,nbStateVar,deltaT)
{
  matR<-listModelParameters$matVarObs
  matQ<-listModelParameters$matVarStateModel
  lambda<-listModelParameters$hyperParaVar
  spread<-listModelParameters$spreadFail
  H<-listModelParameters$H
  L<-listModelParameters$L
  a<-listModelParameters$a
  b<-listModelParameters$b
  
  
  lweight<-oldWeight
  if (situation == 1 || situation == 3)
  {
    for (i in 1:N)
    {
      lweight[i]<-lweight[i] + normDensityLog(c(vectObs[3],vectObs[4]),c(((1+(tan(matSample[((nbSensorFaulty*3)+5),i])*(H/L)))*matSample[((nbSensorFaulty*3)+3),i]),matSample[((nbSensorFaulty*3)+5),i]),matR[3:4,3:4])
    }
  }
  
  if (situation == 2 || situation == 3)
  {
    for (i in 1:N)
    {
      condExp<-stateTransition(matSamplePrev[((3*nbSensorFaulty) +1):((3*nbSensorFaulty) + nbStateVar),i],deltaT,a,b,L)
      hyperParaPrevModified<-matSamplePrev[3,i] + 1
      alphaPrevModified<-(matSamplePrev[3,i]*matSamplePrev[2,i]/hyperParaPrevModified) + (matSample[1,i]/hyperParaPrevModified)
      
      lweight[i]<-lweight[i] + (matSample[1,i]*(log(matSample[2,i]) - log(matSamplePrev[2,i]) + 
                      normDensityLog(c(vectObs[1],vectObs[2]),c(matSample[((nbSensorFaulty*3)+1),i],matSample[((nbSensorFaulty*3)+2),i]),matR[1:2,1:2]) -
                      normDensityLog(c(vectObs[1],vectObs[2]),condExp[1:2],matR[1:2,1:2]))) + ((1-matSample[1,i])*(log(1-matSample[2,i]) - log(1-matSamplePrev[2,i])))
      
      lweight[i]<-lweight[i] + log(dbeta(matSample[2,i],matSamplePrev[3,i]*(1-matSamplePrev[2,i]),matSamplePrev[3,i]*matSamplePrev[2,i])) -
                      log(dbeta(matSample[2,i],(hyperParaPrevModified*(1-alphaPrevModified)),(hyperParaPrevModified*alphaPrevModified)))
    }
  }
  return(lweight)
}

weightRoutineUKF<-function(N,oldWeight,matSample,matSamplePrev,vectObs,situation,listModelParameters,nbSensorFaulty,nbStateVar,deltaT)
{ 
  matR<-listModelParameters$matVarObs
  matQ<-listModelParameters$matVarStateModel
  lambda<-listModelParameters$hyperParaVar
  spread<-listModelParameters$spreadFail
  H<-listModelParameters$H
  L<-listModelParameters$L
  a<-listModelParameters$a
  b<-listModelParameters$b
  
  lweight<-oldWeight
  if (situation == 1 || situation == 3)
  {
    for (i in 1:N)
    {
      lweight[i]<-lweight[i] + normDensityLog(c(vectObs[3],vectObs[4]),c(((1+(tan(matSample[((nbSensorFaulty*3)+5),i])*(H/L)))*matSample[((nbSensorFaulty*3)+3),i]),
                                                                         matSample[((nbSensorFaulty*3)+5),i]),matR[3:4,3:4])
    }
  }
  
  if (situation == 2 || situation == 3)
  {
    for (i in 1:N)
    { 
      lweight[i]<-lweight[i] + (matSample[1,i]*normDensityLog(c(vectObs[1],vectObs[2]),c(matSample[((nbSensorFaulty*3)+1),i],matSample[((nbSensorFaulty*3)+2),i]),matR[1:2,1:2])) + 
        ((1- matSample[1,i])*log(1/2500))
    }
  }
  return(lweight)
}

# state resampling #
stateResamplingRB<-function(vectNewInd,matSample,meanUKF,varUKF,N)
{
  newMatSample<-matrix(0,nrow(matSample),N)
  oldMean<-meanUKF
  oldVar<-varUKF
  for (i in 1:N)
  {
    newMatSample[,i]<- matSample[,vectNewInd[i]]
    meanUKF[,i]<-oldMean[,vectNewInd[i]]
    varUKF[,i]<-oldVar[,vectNewInd[i]]
    
  }
  result = list(newMatSample = newMatSample, meanUKF = meanUKF, varUKF = varUKF)
}


# Estimation du modèles d'états #
#################################

### importation des données ###

dataGPS<-read.table(file.choose()) # fichier GPS.txt
dataSpeedSteer<-read.table(file.choose()) # fichier DRS.txt
dataManager<-read.table(file.choose()) # fichier Sensors_manager.txt

matFailWindow<-matrix(c(35,40,62,67),2,2, byrow = T)
alteredDataGPS2<-dataAlteration(dataGPS,matFailWindow)

dataFinal<-dataPreparation(dataGPS,dataSpeedSteer,dataManager)
alteredDataFinal<-dataPreparation(alteredDataGPS,dataSpeedSteer,dataManager)


# paramètres du modèle
a<-3.78
b<-0.5
H<-0.76
L<-2.83
hyperParaVar<-0.02
matVarStateModel<-matrix(c(0.1,0,0,0.1),2,2)
matVarObs<-diag(c(0.002,0.002,0.1,0.1))
spreadFail<-25
topCensorFailDetection<-0
riskLevelCensorFailDetection<-0
lmodPara<-list(a = a, b = b, H = H, L = L, spreadFail = spreadFail, hyperParaVar = hyperParaVar, matVarStateModel = matVarStateModel, matVarObs = matVarObs,
               topCensorFailDetection = topCensorFailDetection, riskLevelCensorFailDetection = riskLevelCensorFailDetection)
dataInit<-c(4.2,5.1,2.12,0.27)


### filtre Bootstrap ###

#---------------------------------------------------------------------------------------------------------------------------------------------------------------#

# test over 80 GPS data (720 combined data) with GPS failure on [35,40] and on [62,67] GPS observation index with Bootstrap Filter
# 1000 monte Carlo simulation per step
# length execution: 5 min
stateVariableMCBootstrap<-ParticleFilteringAsync(dataInit,alteredDataFinal,initRoutineBootstrap,samplingRoutineBootstrap,weightRoutineBootstrap,systematicRS,stateResamplingBootstrap,
                                                    6,0,1000,lmodPara,1,1000,720)

plot(stateVariableMCBootstrap$mean[,2],stateVariableMCBootstrap$mean[,1], col = "dodgerblue", ylab = "latitude X (mètre)", xlab = "longitude Y (mètre)",pch = 3,
     ylim = c(min(min(stateVariableMCBootstrap$mean[,1]),min(alteredDataGPS[1:80,2])),max(max(stateVariableMCBootstrap$mean[,1]),max(alteredDataGPS[1:80,2]))),
     xlim = c(min(min(stateVariableMCBootstrap$mean[,2]),min(alteredDataGPS[1:80,3])),max(max(stateVariableMCBootstrap$mean[,2]),max(alteredDataGPS[1:80,3]))))
points(alteredDataGPS[1:80,3],alteredDataGPS[1:80,2], col = "red", pch = 5)
points(dataGPS[c(35:40,62:67),3],dataGPS[c(35:40,62:67),2], col = "forestgreen", pch = 15)
legend("topright", legend = c("Bootstrap Filter Estimation", "Data (with failure)","Data (without failure)"), col = c("dodgerblue", "red", "forestgreen"), pch = 15, bty = "n",
       pt.cex = 2, cex = 0.8, text.col = "black", horiz = FALSE, inset = c(0.01, 0.01))

# statbility test
stateVariableMCBootstrap2<-ParticleFilteringAsync(dataInit,alteredDataFinal,initRoutineBootstrap,samplingRoutineBootstrap,weightRoutineBootstrap,systematicRS,stateResamplingBootstrap,
                                                 6,0,1000,lmodPara,1,1000,720)
stateVariableMCBootstrap3<-ParticleFilteringAsync(dataInit,alteredDataFinal,initRoutineBootstrap,samplingRoutineBootstrap,weightRoutineBootstrap,systematicRS,stateResamplingBootstrap,
                                                 6,0,1000,lmodPara,1,1000,720)
stateVariableMCBootstrap4<-ParticleFilteringAsync(dataInit,alteredDataFinal,initRoutineBootstrap,samplingRoutineBootstrap,weightRoutineBootstrap,systematicRS,stateResamplingBootstrap,
                                                 6,0,1000,lmodPara,1,1000,720)
stateVariableMCBootstrap5<-ParticleFilteringAsync(dataInit,alteredDataFinal,initRoutineBootstrap,samplingRoutineBootstrap,weightRoutineBootstrap,systematicRS,stateResamplingBootstrap,
                                                 6,0,1000,lmodPara,1,1000,720)
stateVariableMCBootstrap6<-ParticleFilteringAsync(dataInit,alteredDataFinal,initRoutineBootstrap,samplingRoutineBootstrap,weightRoutineBootstrap,systematicRS,stateResamplingBootstrap,
                                                 6,0,1000,lmodPara,1,1000,720)
stateVariableMCBootstrap7<-ParticleFilteringAsync(dataInit,alteredDataFinal,initRoutineBootstrap,samplingRoutineBootstrap,weightRoutineBootstrap,systematicRS,stateResamplingBootstrap,
                                                 6,0,1000,lmodPara,1,1000,720)
stateVariableMCBootstrap8<-ParticleFilteringAsync(dataInit,alteredDataFinal,initRoutineBootstrap,samplingRoutineBootstrap,weightRoutineBootstrap,systematicRS,stateResamplingBootstrap,
                                                 6,0,1000,lmodPara,1,1000,720)
stateVariableMCBootstrap9<-ParticleFilteringAsync(dataInit,alteredDataFinal,initRoutineBootstrap,samplingRoutineBootstrap,weightRoutineBootstrap,systematicRS,stateResamplingBootstrap,
                                                 6,0,1000,lmodPara,1,1000,720)
stateVariableMCBootstrap10<-ParticleFilteringAsync(dataInit,alteredDataFinal,initRoutineBootstrap,samplingRoutineBootstrap,weightRoutineBootstrap,systematicRS,stateResamplingBootstrap,
                                                 6,0,1000,lmodPara,1,1000,720)

meanY<-rep(0,720)
meanX<-rep(0,720)
for (t in 1:720)
{
  meanY[t]<-mean(stateVariableMCBootstrap$mean[t,2],stateVariableMCBootstrap2$mean[t,2],stateVariableMCBootstrap3$mean[t,2],stateVariableMCBootstrap4$mean[t,2],
                 stateVariableMCBootstrap5$mean[t,2],stateVariableMCBootstrap$mean6[t,2],stateVariableMCBootstrap7$mean[t,2],stateVariableMCBootstrap8$mean[t,2],
                 stateVariableMCBootstrap9$mean[t,2],stateVariableMCBootstrap10$mean[t,2])
  meanX[t]<-mean(stateVariableMCBootstrap$mean[t,1],stateVariableMCBootstrap2$mean[t,1],stateVariableMCBootstrap3$mean[t,1],stateVariableMCBootstrap4$mean[t,1],
                 stateVariableMCBootstrap5$mean[t,1],stateVariableMCBootstrap$mean6[t,1],stateVariableMCBootstrap7$mean[t,1],stateVariableMCBootstrap8$mean[t,1],
                 stateVariableMCBootstrap9$mean[t,1],stateVariableMCBootstrap10$mean[t,1])
  
}

#setEPS()
#postscript("rapport_PF_SMC/figure/Bootstrap_it720_dataGPS80_alterationGPS35to40and62to67_MC1000_ESSBound1000_failRandom.eps")
plot(meanY,meanX, col = "dodgerblue", ylab = "latitude X (mètre)", xlab = "longitude Y (mètre)",pch = 3,
     ylim = c(min(min(meanX),min(alteredDataGPS[1:80,2])),max(max(meanX),max(alteredDataGPS[1:80,2]))),
     xlim = c(min(min(meanY),min(alteredDataGPS[1:80,3])),max(max(meanY),max(alteredDataGPS[1:80,3]))))
points(alteredDataGPS[1:80,3],alteredDataGPS[1:80,2], col = "red", pch = 5)
points(dataGPS[c(35:40,62:67),3],dataGPS[c(35:40,62:67),2], col = "forestgreen", pch = 15)
legend("topleft", legend = c("Bootstrap Filter mean Estimation", "Data (with failure)","Data (without failure)"), col = c("dodgerblue", "red", "forestgreen"), pch = 15, bty = "n",
       pt.cex = 2, cex = 0.8, text.col = "black", horiz = FALSE, inset = c(0.01, 0.01))
#dev.off()

#---------------------------------------------------------------------------------------------------------------------------------------------------------------#


### RB UKF (article algorithm) ###

#---------------------------------------------------------------------------------------------------------------------------------------------------------------#

# test over 500 combined data (56 GPS data) with GPS failure on [35,40] GPS observation index with bayesian RB UKF
# 3000 monte Carlo simulation per step
# length execution: 2h
stateVariableMCRBUKF<-ParticleFilteringAsync(dataInit,alteredDataFinal,initRoutineRB,samplingRoutineRB,weightRoutineRB,systematicRS,stateResamplingRB,
                                             6,1,1000,lmodPara,1,1000,500)

#setEPS()
#postscript("rapport_PF_SMC/figure/RBUKF_it500_dataGPS56_alterationGPS35to40_MC1000_ESSBound1000_failRandom.eps")
plot(stateVariableMCRBUKF$mean[,2],stateVariableMCRBUKF$mean[,1], col = "dodgerblue", ylab = "latitude X (mètre)", xlab = "longitude Y (mètre)",pch = 3,
     ylim = c(min(min(stateVariableMCRBUKF$mean[,1]),min(alteredDataGPS[,2])),max(max(stateVariableMCRBUKF$mean[,1]),max(alteredDataGPS[1:56,2]))),
     xlim = c(min(min(stateVariableMCRBUKF$mean[,2]),min(alteredDataGPS[,3])),max(max(stateVariableMCRBUKF$mean[,2]),max(alteredDataGPS[1:56,3]))))
points(alteredDataGPS2[1:56,3],alteredDataGPS2[1:56,2], col = "red", pch = 5)
points(dataGPS[35:40,3],dataGPS[35:40,2], col = "forestgreen", pch = 15)
legend("topleft", legend = c("RB UKF Estimation", "Data (with failure)","Data (without failure)"), col = c("dodgerblue", "red", "forestgreen"), pch = 15, bty = "n",
       pt.cex = 2, cex = 0.8, text.col = "black", horiz = FALSE, inset = c(0.01, 0.01))
#dev.off()

MSE_RBUKF<-computationMSELatLong(dataGPS,stateVariableMCRBUKF$mean,dataFinal)
MSE_RBUKF

cMeanResult<-stateVariableMCRBUKF$mean[,7]
heatMapC<-matrix(0,length(cMeanResult),2)
heatMapC[,1]<-1-cMeanResult
heatMapC[,2]<-cMeanResult
colnames(heatMapC) <- c("c=0", "c=1")

image(seq(nrow(heatMapC)), seq(ncol(heatMapC)), heatMapC, col = gray(seq(0,1, 0.01)),breaks = seq(0, 1.01, 0.01), axes = FALSE, xlab = "", ylab = "")
axis(1, at = seq(nrow(heatMapC)))
axis(2, at = seq(ncol(heatMapC)),labels = colnames(heatMapC), las = 2)

alphaMeanResult<-rep(0,0)
for (i in 1:length(stateVariableMCRBUKF$mean[,8]))
{
  if (dataFinal[i,6] == 2 || dataFinal[i,6] == 3)
  {
    alphaMeanResult<-c(alphaMeanResult,stateVariableMCRBUKF$mean[i,8])
  }
}

plot(alphaMeanResult,col = "forestgreen", type = "l", ylim = c(0,1), ylab = "alpha", lwd = 1.5)
abline(v =35, col = "red", lwd = 1.5)
abline(v =40, col = "red", lwd = 1.5)


#---------------------------------------------------------------------------------------------------------------------------------------------------------------#


### UKF only (without failure detection) ###

#---------------------------------------------------------------------------------------------------------------------------------------------------------------#

# test over 500 combined data (80 GPS data) with GPS failure on [35,40] GPS observation index with UKF only
stateVariableMCUKF<-ParticleFilteringAsync(dataInit,alteredDataFinal,initRoutineUKF,samplingRoutineUKF,weightRoutineUKF,systematicRS,stateResamplingRB,
                                             6,1,1,lmodPara,1,1,500)
#setEPS()
#postscript("rapport_PF_SMC/figure/UKF_it500_dataGPS56_alterationGPS35to40_noFailDetection_failRandom.eps")
plot(stateVariableMCUKF$mean[,2],stateVariableMCUKF$mean[,1], col = "dodgerblue", ylab = "latitude X (mètre)", xlab = "longitude Y (mètre)",pch = 3,
     ylim = c(min(min(stateVariableMCUKF$mean[,1]),min(alteredDataGPS[1:56,2])),max(max(stateVariableMCUKF$mean[,1]),max(alteredDataGPS[1:56,2]))),
     xlim = c(min(min(stateVariableMCUKF$mean[,2]),min(alteredDataGPS[1:56,3])),max(max(stateVariableMCUKF$mean[,2]),max(alteredDataGPS[1:56,3]))))
points(alteredDataGPS[1:56,3],alteredDataGPS[1:56,2], col = "red", pch = 5)
points(dataGPS[35:40,3],dataGPS[35:40,2], col = "forestgreen", pch = 3)
legend("topleft", legend = c("UKF Filter", "Data (with failure)","Data (without failure)"), col = c("dodgerblue", "red", "forestgreen"), pch = 15, bty = "n",
       pt.cex = 2, cex = 0.8, text.col = "black", horiz = FALSE, inset = c(0.01, 0.01))
#dev.off()

MSE_UKF<-computationMSELatLong(dataGPS,stateVariableMCUKF$mean,dataFinal)
MSE_UKF

#---------------------------------------------------------------------------------------------------------------------------------------------------------------#


### UKF only (with failure detection) ### 

topCensorFailDetection<-1
riskLevelCensorFailDetection<-0.05

#---------------------------------------------------------------------------------------------------------------------------------------------------------------#

# test over 500 combined data (56 GPS data) with GPS failure on [35,40] GPS observation index with UKF only and with Sensor Failure detector based on quadratice innovation

# NO SUCCESS
stateVariableMCUKF_FD<-ParticleFilteringAsync(dataInit,alteredDataFinal,initRoutineUKF,samplingRoutineUKF,weightRoutineUKF,systematicRS,stateResamplingRB,
                                           6,1,1,lmodPara,1,1,500)
#setEPS()
#postscript("rapport_PF_SMC/figure/UKF_it500_dataGPS56_alterationGPS35to40_FailDetection_failRandom.eps")
plot(stateVariableMCUKF_FD$mean[,2],stateVariableMCUKF_FD$mean[,1], col = "dodgerblue", ylab = "latitude X (mètre)", xlab = "longitude Y (mètre)",pch = 3,
     ylim = c(min(min(stateVariableMCUKF_FD$mean[,1]),min(alteredDataGPS[1:56,2])),max(max(stateVariableMCUKF_FD$mean[,1]),max(alteredDataGPS[1:56,2]))),
     xlim = c(min(min(stateVariableMCUKF_FD$mean[,2]),min(alteredDataGPS[1:56,3])),max(max(stateVariableMCUKF_FD$mean[,2]),max(alteredDataGPS[1:56,3]))))
points(alteredDataGPS[1:56,3],alteredDataGPS[1:56,2], col = "red", pch = 5)
points(dataGPS[35:40,3],dataGPS[35:40,2], col = "forestgreen", pch = 3)
legend("topleft", legend = c("UKF Filter Fail Detec", "Data (with failure)","Data (without failure)"), col = c("dodgerblue", "red", "forestgreen"), pch = 15, bty = "n",
       pt.cex = 2, cex = 0.8, text.col = "black", horiz = FALSE, inset = c(0.01, 0.01))
#dev.off()

MSE_UKF_FD<-computationMSELatLong(dataGPS,stateVariableMCUKF_FD$mean,dataFinal)
MSE_UKF_FD

#---------------------------------------------------------------------------------------------------------------------------------------------------------------#









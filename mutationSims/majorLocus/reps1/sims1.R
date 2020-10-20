library(quantPop)
reps <- 1:100                          # vector of repetition identifiers
changePhenOpt <- 10
mutRate <-  7.35e-09                  # polygenic mutation rate to get to Vg = 0.6 from minor mutations (from Bürger & Lande 1994)
startHardSelec <- 1201                # the default generation to start density dependent population growth
c <- 6                                # standard deviation (width) of the fitness function
runTime <- 1300                       # number of generations to run the simulations
phenoLimit <- rep(NA,length(reps))    # save the phenotypic evolutionary limit
fitFun <- function(phen,maxFit,fitPhen,fitSd){maxFit*exp(-(((phen-fitPhen)^2)/(2*fitSd^2)))}   # the Gaussian fitness function
alphas <- seq(3,8,0.001)              # vector of potential large allelic effects. We sample one value randomly from this vector for each simulation repetition
freqIndicator <- rep(c(TRUE,FALSE),length(reps))      #TRUE for high frequency (p>0.5), FALSE for low frequency (p<0.5) large effect allele
phenOptMat <- matrix(NA,nrow=length(reps),ncol=runTime)   # matrix with each row representing a vector of phenotype optima through time for each simulation repetition
for(i in 1:length(reps)){
  if(i %in% seq(1,1000,2)) phenOptMat[i,] <- phenOptVec <- c(rep(0,1100),rep(12,100),rep(20,100)) # for high frequency (p>0.5) major alleles
  if(i %in% seq(2,1000,2)) phenOptMat[i,] <- phenOptVec <- c(rep(0,1100),rep(6,100),rep(14,100))  # for low frequency (p<0.5) major alleles
}
###########################################################
# initialize objects to store the results
###########################################################
popSizeMat <- matrix(NA,nrow=length(reps),ncol=runTime)
h2Mat <- matrix(NA,nrow=length(reps),ncol=runTime)
addVarMat <- matrix(NA,nrow=length(reps),ncol=runTime)
phenVarMat <- matrix(NA,nrow=length(reps),ncol=runTime)
numSegMat <- matrix(NA,nrow=length(reps),ncol=runTime)
meanPhenMat <- matrix(NA,nrow=length(reps),ncol=runTime)
NeMat <- matrix(NA,nrow=length(reps),ncol=runTime)
effSizeMat <- NULL
majorFreqMat <- matrix(NA,nrow=length(reps),ncol=runTime)
majorLocusFreqMat <- matrix(NA,nrow=length(reps),ncol=runTime)
hardSelecGenVec <- rep(NA,length(reps))
phenOptVec <- rep(NA,length(reps))    # save the realized phenotypic optima from the simulations
freqList <- list()               # save all the allele frequencies
effectSize <- sample(alphas,1)   # effect size of the large-effect locus (randomly selected from 'alphas')

effSizeVec <- rep(NA,length(reps))
i <- 1   # index repetitions

while(i %in% 1:length(reps)){
  effectSize <- sample(alphas,1)   # effect size of the large-effect locus (randomly selected from 'alphas')
  effSizeVec[i] <- effectSize

  outFreqs <- NULL
  attempt <- 1
  thisHardSelecGen <- startHardSelec
  while( (is.null(outFreqs) && attempt < 10) | thisHardSelecGen >= startHardSelec) {
    try(
      logQuant_mut_h2Dep(burnin=1100,gens=runTime,genSize=1000000,phen0=0,phenOpt=phenOptMat[i,],phenShift=10,
            c=c,mu=mutRate,muOff=startHardSelec-1,minSize=-0.5,maxSize=0.5,N=500,Ve=4,numLgEff=1,lgSize=effectSize,
            hardSelecGen=startHardSelec,desiredH2 = 0.6,h2Diff=0.05,highMajorFreq=freqIndicator[i],K=1000,lambda=1.5,f=4)
    )
    thisHardSelecGen <- hardSelecGen
    attempt <- attempt + 1
  }
  phenOptVec[i] <- max(phenOpt)
  hardSelecGenVec[i] <- thisHardSelecGen
  freqList [[i]] <- outFreqs
  h2Mat[i,1:length(phenVar)] <- addVar/phenVar
  addVarMat[i,1:length(addVar)] <- addVar
  phenVarMat[i,1:length(phenVar)] <- phenVar
  numSegMat[i,1:length(numSegVec)] <- numSegVec
  meanPhenMat[i,] <- colMeans(phenMat,na.rm=TRUE)
  popSizeMat[i,] <- NVec
  NeMat[i,] <- Ne
  majorLocusFreqMat[i,] <- majorLocusFreqs
  hardSelecGenVec[i] <- hardSelecGen
  #########################################
  # get the phenotypic evolutionary limit
  #########################################
  segLoci <- which(outFreqs[,hardSelecGen] > 0)
  segFreqs <- outFreqs[segLoci,hardSelecGen]
  segEffs <- freqMat[segLoci,3]
  Ve <- 4                                         # environmental variance in the phenotype
  zBar_2 <- rep(NA,length(segLoci))               # mean phenotype conditional on having two beneficial alleles
  zBar_1 <- rep(NA,length(segLoci))               # mean phenotype conditional on having one beneficial allele
  zBar_0 <- rep(NA,length(segLoci))               # mean phenotype conditional on having no beneficial alleles

  wBar_2 <- rep(NA,length(segLoci))               # mean fitness of each of the three possible genotypes
  wBar_1 <- rep(NA,length(segLoci))
  wBar_0 <- rep(NA,length(segLoci))

  sVec <- rep(NA,length(segLoci))
  for(j in 1:length(segLoci)){
    zBar_2 [j] <- 2*segEffs[j] + sum((segFreqs[-j]^2)*2*segEffs[-j] + 2*segFreqs[-j]*(1-segFreqs[-j])*segEffs[-j])
    zBar_1 [j] <- segEffs[j] + sum((segFreqs[-j]^2)*2*segEffs[-j] + 2*segFreqs[-j]*(1-segFreqs[-j])*segEffs[-j])
    zBar_0 [j] <- 0 + sum((segFreqs[-j]^2)*2*segEffs[-j] + 2*segFreqs[-j]*(1-segFreqs[-j])*segEffs[-j])
    rawFit_2 <- fitFun(phen=c(zBar_2[j]),maxFit = 1,fitPhen=phenOptVec[i],fitSd=c)
    rawFit_1 <- fitFun(phen=c(zBar_1[j]),maxFit = 1,fitPhen=phenOptVec[i],fitSd=c)
    rawFit_0 <- fitFun(phen=c(zBar_0[j]),maxFit = 1,fitPhen=phenOptVec[i],fitSd=c)

    wBar_2 [j] <- rawFit_2/(max(c(rawFit_0,rawFit_1,rawFit_2)))
    wBar_1 [j] <- rawFit_1/(max(c(rawFit_0,rawFit_1,rawFit_2)))
    wBar_0 [j] <- rawFit_0/(max(c(rawFit_0,rawFit_1,rawFit_2)))
    sVec[j] <- wBar_2[j] - wBar_0[j]
  }
  fixProbs <- (1-exp(-4*Ne[hardSelecGenVec[i]]*sVec*segFreqs))/(1-exp(-4*Ne[hardSelecGenVec[i]]*sVec))
  phenoLimit [i] <- sum(fixProbs * 2*segEffs)       # the expected selection limit
  freqMat <- cbind(rep(i,nrow(freqMat)),freqMat)
  effSizeMat <- rbind(effSizeMat,freqMat)

  # save the allele frequencies if you want to
  # write.table(outFreqs,file = paste("alleleFrequencies_effectSize_",effectSize,"_optimPhen_",optimPhen,"_mutRate_",
  #                                  mutRate,"_fitFunSd_",c,"_rep",reps[i],sep=""))
  i <- i + 1   # iterator
}

phenLimOut <- NULL
phenLimOut <- cbind(phenLimOut,phenoLimit)

###############################################
# save the results
###############################################
# write files
write.table(majorLocusFreqMat,file=paste("majorAlleleFrequenciesAcrossSimulations_effectSize_",effectSize,"_mutRate_",mutRate,"_fitFunSd_",c,"_",min(reps),
                                    "_",max(reps),sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(h2Mat,file=paste("h2TimeSeriesAcrossSimulations_effectSize_",effectSize,"_mutRate_",mutRate,"_fitFunSd_",c,"_",min(reps),"_",max(reps),sep=""),
            quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(numSegMat,file=paste("numSegregatingQTLsAcrossSimulations_effectSize_",effectSize,"_mutRate_",mutRate,"_fitFunSd_",c,"_",min(reps),"_",max(reps),sep=""),
            quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(phenVarMat,file=paste("phenotypicVarianceSeriesAcrossSimulations_effectSize_",effectSize,"_mutRate_",mutRate,"_fitFunSd_",c,"_",min(reps),"_",max(reps),sep=""),
            quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(addVarMat,file=paste("additiveVarianceAcrossSimulations_effectSize_",effectSize,"_mutRate_",mutRate,"_fitFunSd_",c,"_",min(reps),"_",max(reps),sep=""),
            quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(popSizeMat,file=paste("populationSize_effectSize_",effectSize,"_mutRate_",mutRate,"_fitFunSd_",c,"_",min(reps),"_",max(reps),sep=""),
            quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(meanPhenMat,file=paste("meanPhenotype_effectSize_",effectSize,"_mutRate_",mutRate,"_fitFunSd_",c,"_",min(reps),"_",max(reps),sep=""),
            quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(phenLimOut,file=paste("phenoLimit_effectSize_",effectSize,"_mutRate_",mutRate,"_fitFunSd_",c,"_",min(reps),"_",max(reps),sep=""),
            quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(NeMat,file=paste("Ne_effectSize_",effectSize,"_mutRate_",mutRate,"_fitFunSd_",c,"_",min(reps),"_",max(reps),sep=""),
            quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(phenLimOut,file=paste("phenoLimit_effectSize_",effectSize,"_mutRate_",mutRate,"_fitFunSd_",c,"_",min(reps),"_",max(reps),sep=""),
            quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(phenOptVec,file=paste("phenoOptima_effectSize_",effectSize,"_mutRate_",mutRate,"_fitFunSd_",c,"_",min(reps),"_",max(reps),sep=""),
            quote=FALSE,row.names=FALSE,col.names=FALSE)

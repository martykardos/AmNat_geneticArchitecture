library(quantPop)
reps <- 601:700                         # vector of repetition identifiers
effectSize <- 0                      # effect size of the large-effect locus
optimPhen <- 10                      # optimum phenotype at the end of the simulation
lastGenMu <- 1200                    # last generation mutations are allowed
hardSelecGen <- 1201                 # the generation to start density dependent population growth
c <- 6                               # standard deviation of the fitness function
runTime <- 1300                      # number of generations to run the simulations for
phenoLimit <- rep(NA,length(reps))   # save the phenotypic evolutionary limit


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
fitFun <- function(phen,maxFit,fitPhen,fitSd){maxFit*exp(-(((phen-fitPhen)^2)/(2*fitSd^2)))}   # the Gaussian fitness function
mutRate <- 7.35e-08        # mutation rate
for(i in 1:length(reps)){

  Ne <- NULL           # initialize an Ne vector
  while( is.null(Ne)) {            
    try(
      logQuant_mut(burnin=1,gens=runTime,genSize=1000000,phen0=0,phenOpt=c(rep(0,1200),rep(optimPhen,100)),c=c,mu=mutRate,muOff=lastGenMu,
                   minSize=-0.5,maxSize=0.5,N=500,Ve=4,numLgEff=0,lgSize=effectSize,hardSelecGen=hardSelecGen,K=1000,lambda=1.5,f=4)
    )
  } 
  
  h2Mat[i,1:length(phenVar)] <- addVar/phenVar
  addVarMat[i,1:length(addVar)] <- addVar
  phenVarMat[i,1:length(phenVar)] <- phenVar
  numSegMat[i,1:length(numSegVec)] <- numSegVec
  meanPhenMat[i,] <- colMeans(phenMat,na.rm=TRUE)
  popSizeMat[i,] <- NVec
  NeMat[i,] <- Ne
  
  #########################################
  # get the phenotypic evolutionary limit
  #########################################
  segLoci <- which(outFreqs[,1200] > 0)
  segFreqs <- outFreqs[segLoci,1200]
  segEffs <- freqMat[segLoci,3]
  Ve <- 4                                          # environmental variance in the phenotype
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
    rawFit_2 <- fitFun(phen=c(zBar_2[j]),maxFit = 1,fitPhen=optimPhen,fitSd=c)
    rawFit_1 <- fitFun(phen=c(zBar_1[j]),maxFit = 1,fitPhen=optimPhen,fitSd=c)
    rawFit_0 <- fitFun(phen=c(zBar_0[j]),maxFit = 1,fitPhen=optimPhen,fitSd=c)
    
    wBar_2 [j] <- rawFit_2/(max(c(rawFit_0,rawFit_1,rawFit_2)))
    wBar_1 [j] <- rawFit_1/(max(c(rawFit_0,rawFit_1,rawFit_2)))
    wBar_0 [j] <- rawFit_0/(max(c(rawFit_0,rawFit_1,rawFit_2)))
    sVec[j] <- wBar_2[j] - wBar_0[j]
  }
  fixProbs <- (1-exp(-4*Ne[1201]*sVec*segFreqs))/(1-exp(-4*Ne[1201]*sVec))
  phenoLimit [i] <- sum(fixProbs * 2*segEffs)       # the expected selection limit
  freqMat <- cbind(rep(i,nrow(freqMat)),freqMat)
  effSizeMat <- rbind(effSizeMat,freqMat)
  
  # save the allele frequencies
  #write.table(outFreqs,file = paste("alleleFrequencies_effectSize_",effectSize,"_optimPhen_",optimPhen,"_mutRate_",
  #                                  mutRate,"_fitFunSd_",c,"_rep",reps[i],sep=""))
}



phenLimOut <- NULL
phenLimOut <- cbind(phenLimOut,phenoLimit)

###############################################
# save the results
###############################################


# write files
write.table(h2Mat,file=paste("h2TimeSeriesAcrossSimulations_effectSize_",effectSize,"_optimPhen_",optimPhen,"_mutRate_",mutRate,"_fitFunSd_",c,"_",min(reps),"_",max(reps),sep=""),
            quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(numSegMat,file=paste("numSegregatingQTLsAcrossSimulations_effectSize_",effectSize,"_optimPhen_",optimPhen,"_mutRate_",mutRate,"_fitFunSd_",c,"_",min(reps),"_",max(reps),sep=""),
            quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(phenVarMat,file=paste("phenotypicVarianceSeriesAcrossSimulations_effectSize_",effectSize,"_optimPhen_",optimPhen,"_mutRate_",mutRate,"_fitFunSd_",c,"_",min(reps),"_",max(reps),sep=""),
            quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(addVarMat,file=paste("additiveVarianceAcrossSimulations_effectSize_",effectSize,"_optimPhen_",optimPhen,"_mutRate_",mutRate,"_fitFunSd_",c,"_",min(reps),"_",max(reps),sep=""),
            quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(popSizeMat,file=paste("populationSize_effectSize_",effectSize,"_optimPhen_",optimPhen,"_mutRate_",mutRate,"_fitFunSd_",c,"_",min(reps),"_",max(reps),sep=""),
            quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(meanPhenMat,file=paste("meanPhenotype_effectSize_",effectSize,"_optimPhen_",optimPhen,"_mutRate_",mutRate,"_fitFunSd_",c,"_",min(reps),"_",max(reps),sep=""),
            quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(phenLimOut,file=paste("phenoLimit_effectSize_",effectSize,"_optimPhen_",optimPhen,"_mutRate_",mutRate,"_fitFunSd_",c,"_",min(reps),"_",max(reps),sep=""),
            quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(NeMat,file=paste("Ne_effectSize_",effectSize,"_optimPhen_",optimPhen,"_mutRate_",mutRate,"_fitFunSd_",c,"_",min(reps),"_",max(reps),sep=""),
            quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(phenLimOut,file=paste("phenoLimit_effectSize_",effectSize,"_optimPhen_",optimPhen,"_mutRate_",mutRate,"_fitFunSd_",c,"_",min(reps),"_",max(reps),sep=""),
            quote=FALSE,row.names=FALSE,col.names=FALSE)





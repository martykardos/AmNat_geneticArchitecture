


##################################################################
#-----------------------------------------------------------------
# calculate evolutionary potential versus time for large mammals
# and major loci
#-----------------------------------------------------------------
##################################################################

gensEvolPot <- 80      # number of generations from the start to quantify evolutionary potential (i.e., the potential to evolove over ten generations)

#fitness function
fitFun <- function(phen,phenMu,phenSd,fitPhen,fitSd,maxFit){(1/(phenSd*sqrt(2*pi)))*maxFit*exp(-(((phen-fitPhen)^2)/(2*(fitSd^2)))-(((phen-phenMu)^2)/(2*(phenSd^2))))}
fit_sd <- 6     # standard deviation of the fitness function
Ve <- 4         # environmental phenotypic variance component
lambda <- 1.5   # intrinsic population growth rate (i.e., perfect adaptation and density near zero)
nLoci <- 100    # number of loci affecting the selected phenotype
setwd("~/Documents/work/rapidAdaptation/AmNatRevision/simulations/originalSims/h2Pt6/fecundity4/major")

###### major loci first
qtlInfo <- read.table("qtlInfo_N0_500_K_1000_propMaj_0.9",row.names=NULL)[,2:4]
rowStarts <- seq(1,nrow(qtlInfo),100)
rowEnds <- seq(100,nrow(qtlInfo),100)
allFreqs <- read.table("allFreqs_N0_500_K_1000_propMaj_0.9",header=TRUE)

potMatMajor <- matrix(NA,nrow=500,ncol=80)
for(i in 1:500){
  #####################################################
  # get the simulated data for the ith simulation rep
  #####################################################
  thisInfo <- qtlInfo[rowStarts[i]:rowEnds[i],]
  theseFreqs <- allFreqs[which(allFreqs[,1] == i),2:(ncol(allFreqs)-1)]

  evPotVec <- rep(NA,gensEvolPot)     # save evolutionary potential each generation
  genValMat <- cbind(2*thisInfo[,2],thisInfo[,2],rep(0,nrow(thisInfo)))      # genetic values of the three possible genotypes

  ############################################################
  # calculate the short term evolutionary potential
  # of this popualtion. This is the expected evolutionary
  # response with a constant difference between the mean
  # phenotype and the optimum phenotype and no genetic drift
  ############################################################
  evPotVec <- rep(NA,gensEvolPot)                                              # vector to store evolutionary potentials
  for(lm in 1:gensEvolPot){                                                    # loop over the generations you want to get the ten generation evolutionary potential for
    if(is.na(theseFreqs[1,lm]) == FALSE){
      newFreqs <- matrix(NA,nrow=nrow(theseFreqs),ncol=ncol(theseFreqs))         # projection of allele frequencies to calculate short term evolutionary potential
      newFreqs[,1] <- theseFreqs[,lm]
      phenVec <- rep(NA,10)                    # mean phenotype each generation
      vpVec <- rep(NA,gensEvolPot)             # phenotype varianbce

      for(j in 1:10){                          # calculate the 10 generation evolutionary potential as the total deterministic/ideal change in the phenotype over 10 generations
        genoFreqMat <- cbind(newFreqs[,j]^2,2*newFreqs[,j]*(1-newFreqs[,j]),(1-newFreqs[,j])^2)     # expected genotype frequencies in generation j
        phenVec [j] <- sum(rowSums(genValMat*genoFreqMat))                   # mean expected phenotype in the jth generation
        vpVec [j] <- sum(genoFreqMat[,2]*genValMat[,2]^2) + Ve               # total expected phenotypic variance this generation
        phenRange <- c(phenVec[j]-4*sqrt(vpVec[j]),phenVec[j]+4*sqrt(vpVec[j]))  # range across which to integreate the fitness function just below here
        pNext <- rep(NA,nrow(genoFreqMat))      # store allele frequencies for the next generation
        for(k in 1:nrow(genoFreqMat)){          # loop over loci and get the expected allele frequency in the next generation
          #-----------------------------------------------------------------------------------------------------------------
          # mean phenotype for each possible genotype at the kth locus, considering genotype distributions at all other loci
          #------------------------------------------------------------------------------------------------------------------

          phen_A1A1  <-  genValMat[k,3] + sum(rowSums(genValMat[-k,]*genoFreqMat[-k,]))  # mean phenotype conditional on being homozygous for the less fit allele at one locus
          phen_A1A2  <-  genValMat[k,2] + sum(rowSums(genValMat[-k,]*genoFreqMat[-k,]))  # mean phenotype conditional on being heterozygous at one locus
          phen_A2A2  <-  genValMat[k,1] + sum(rowSums(genValMat[-k,]*genoFreqMat[-k,]))  # mean phenotype conditional on being homozygous for the most fit allele at one locus

          #-------------------------------------------------------------------------------------------------
          # mean fitness for each possible genotype, considering genotype distributions at all other loci
          #-------------------------------------------------------------------------------------------------
          cond_Vp <- sum(genoFreqMat[-k,2]*genValMat[-k,2]^2) + Ve # phenotypic variance conditional on holding genotype constant at the kth locus

          # A1A1
          phenRange_A1A1 <- c(phen_A1A1-4*sqrt(cond_Vp),phen_A1A1+4*sqrt(cond_Vp))  # range across which to integreate the fitness function just below here
          meanFit_A1A1 <- integrate(fitFun,lower=phenRange_A1A1[1],upper=phenRange_A1A1[2],fitPhen=phenVec[j] + 10,fitSd=fit_sd,maxFit=lambda,phenMu=phen_A1A1,phenSd=sqrt(cond_Vp))[1][[1]] # mean fitness for individuals with A1A1 genotypes

          # A1A2
          phenRange_A1A2 <- c(phen_A1A2-4*sqrt(cond_Vp),phen_A1A2+4*sqrt(cond_Vp))  # range across which to integreate the fitness function just below here
          meanFit_A1A2 <- integrate(fitFun,lower=phenRange_A1A2[1],upper=phenRange_A1A2[2],fitPhen=phenVec[j] + 10,fitSd=fit_sd,maxFit=lambda,phenMu=phen_A1A2,phenSd=sqrt(cond_Vp))[1][[1]] # mean fitness for individuals with A1A2 genotypes

          # A2A2
          phenRange_A2A2 <- c(phen_A2A2-4*sqrt(cond_Vp),phen_A2A2+4*sqrt(cond_Vp))  # range across which to integreate the fitness function just below here
          meanFit_A2A2 <- integrate(fitFun,lower=phenRange_A2A2[1],upper=phenRange_A2A2[2],fitPhen=phenVec[j] + 10,fitSd=fit_sd,maxFit=lambda,phenMu=phen_A2A2,phenSd=sqrt(cond_Vp))[1][[1]] # mean fitness for individuals with A2A2 genotypes
          maxFit <- max(c( meanFit_A1A1, meanFit_A1A2, meanFit_A2A2))
          meanFit_A1A1 <-  meanFit_A1A1/maxFit
          meanFit_A1A2 <-  meanFit_A1A2/maxFit
          meanFit_A2A2 <-  meanFit_A2A2/maxFit
          pNext[k] <- ((newFreqs[k,j]^2)*meanFit_A2A2 + (newFreqs[k,j]*(1-newFreqs[k,j])*meanFit_A1A2))/mean(c( meanFit_A1A1, meanFit_A1A2, meanFit_A2A2))
        }
        if(sum(pNext > 1) > 0) pNext[which(pNext > 1)] <- 1
        newFreqs[,j+1] <- pNext
      }
      evPotVec[lm] <- phenVec[length(phenVec)] - phenVec[1]
    }
  }
  potMatMajor[i,] <- evPotVec
  print(paste("done with simulation repetition ",i))
}

write.table(potMatMajor,file="evolutionarLimits_N0_500_K_1000_propMaj_0.9",quote=FALSE,row.names=FALSE)


plot(c(0,80),c(0,15),type="n")
for(i in 1:10){
  lines(1:80,potMatMajor[i,])
}

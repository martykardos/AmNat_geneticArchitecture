
library(quantPop)
loci <- 100
reps <- 500
propMaj <- 0.5
majFreq <- 0.3
N_0 <- 500
h20 <- 0.8
carryCap <- 1000
NMat <- NULL    # store population sizes
phenMat <- NULL
h2Mat <- NULL
p0Vec <- NULL

  for(j in 1:reps){
    logQuantIndiv_varP0(N0=N_0,K=carryCap,t=80,lambda=1.5,f=4,sexType="herm",nLoci=2,phen_0=100,phen_opt=110,
                  fit_sd=6,h2_0=h20,Vp_0=10,beta1=0.5,beta2=0.5,freqBounds=c(0.05,0.95),nMajor=1,propMajor=propMaj,freqMajor=majFreq,immRate=0)
    NMat <- rbind(NMat,NVec)
    phenMat <- rbind(phenMat,phenVec)
    h2Mat <- rbind(h2Mat,h2Vec)
    print(j)
  }


write.table(NMat,file=paste("popSize_N0_",N_0,"_K_",carryCap,"_propMaj_",propMaj,"_freqMaj_",majFreq,sep=""),quote=FALSE,row.names=FALSE)
write.table(phenMat,file=paste("phen_N0_",N_0,"_K_",carryCap,"_propMaj_",propMaj,"_freqMaj_",majFreq,sep=""),quote=FALSE,row.names=FALSE)
write.table(h2Mat,file=paste("h2_N0_",N_0,"_K_",carryCap,"_propMaj_",propMaj,"_freqMaj_",majFreq,sep=""),quote=FALSE,row.names=FALSE)


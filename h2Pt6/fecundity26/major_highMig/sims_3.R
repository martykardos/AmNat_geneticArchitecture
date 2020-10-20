
library(quantPop)
loci <- 100
reps <- 500
propMaj <- 0.7
N_0 <- 5000
h20 <- 0.6
carryCap <- 10000
imm <- 8
NMat <- NULL    # store population sizes
phenMat <- NULL
h2Mat <- NULL
p0Vec <- NULL

  for(j in 1:reps){
    logQuantIndiv_varP0(N0=N_0,K=carryCap,t=80,lambda=1.3,f=26,sexType="herm",nLoci=loci,phen_0=100,phen_opt=110,
                  fit_sd=6,h2_0=h20,Vp_0=10,beta1=0.5,beta2=0.5,freqBounds=c(0.05,0.95),nMajor=1,propMajor=propMaj,freqMajor=runif(1,min=0.05,max=0.95),immRate=imm)
    NMat <- rbind(NMat,NVec)
    phenMat <- rbind(phenMat,phenVec)
    h2Mat <- rbind(h2Mat,h2Vec)
    print(j)
  }


write.table(NMat,file=paste("popSize_N0_",N_0,"_K_",carryCap,"_propMaj_",propMaj,sep=""),quote=FALSE,row.names=FALSE)
write.table(phenMat,file=paste("phen_N0_",N_0,"_K_",carryCap,"_propMaj_",propMaj,sep=""),quote=FALSE,row.names=FALSE)
write.table(h2Mat,file=paste("h2_N0_",N_0,"_K_",carryCap,"_propMaj_",propMaj,sep=""),quote=FALSE,row.names=FALSE)


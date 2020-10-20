library(quantPop)
loci <- 100
reps <- 500
N_0 <- 5000
carryCap <- 10000
p_0 <- 0.1
NMat <- NULL    # store population sizes
phenMat <- NULL
h2Mat <- NULL
p0Vec <- NULL

  for(j in 1:reps){
    logQuantIndiv(N0=N_0,K=carryCap,t=80,lambda=1.3,f=26,sexType="herm",nLoci=loci,phen_0=100,phen_opt=110,
                  fit_sd=6,h2_0=0.4,Vp_0=10,p0=p_0)
    NMat <- rbind(NMat,NVec)
    phenMat <- rbind(phenMat,phenVec)
    h2Mat <- rbind(h2Mat,h2Vec)
    print(j)
  }


write.table(NMat,file=paste("popSize_p0_",p_0,"_N0_",N_0,"_K_",carryCap,sep=""),quote=FALSE,row.names=FALSE)
write.table(phenMat,file=paste("phen_p0_",p_0,"_N0_",N_0,"_K_",carryCap,sep=""),quote=FALSE,row.names=FALSE)
write.table(h2Mat,file=paste("h2_p0_",p_0,"_N0_",N_0,"_K_",carryCap,sep=""),quote=FALSE,row.names=FALSE)


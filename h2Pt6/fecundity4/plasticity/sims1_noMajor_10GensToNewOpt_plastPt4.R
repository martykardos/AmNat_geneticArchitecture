library(quantPop)
setwd("~/Documents/work/rapidAdaptation/AmNatRevision/simulations/originalSims/h2Pt6/fecundity4/plasticity")
loci <- 100
reps <- 500
propMaj <- 0.01
N_0 <- 500
h20 <- 0.6
plast <- 0.4
carryCap <- 1000
NMat <- NULL    # store population sizes
phenMat <- NULL
h2Mat <- NULL
optPhenMat <- NULL
p0Vec <- NULL
gensToOpt <- 10
NeSave <- rep(NA,reps)
qtlInfoMat <- NULL
freqOut <- NULL
for(j in 1:reps){
  logQuantIndiv_varP0_stocPhenOpt_plasticity(N0=N_0,K=carryCap,t=80,lambda=1.5,f=4,sexType="herm",nLoci=loci,phen_0=100,phen_opt=110,plastic=plast,timeToOpt=gensToOpt,optSd=2,  
                                             fit_sd=6,h2_0=h20,Vp_0=10,beta1=0.5,beta2=0.5,freqBounds=c(0.05,0.95),nMajor=1,propMajor=propMaj,freqMajor=runif(1,min=0.05,max=0.95),immRate=0)
  NMat <- rbind(NMat,NVec)
  freqOut <- rbind(freqOut,freqMat)
  phenMat <- rbind(phenMat,phenVec)
  h2Mat <- rbind(h2Mat,h2Vec)  
  optPhenMat <- rbind(optPhenMat,optPhenVec)
  NeSave[j] <- NeVec[1]
  qtlInfoMat <- rbind(qtlInfoMat,qtlInfo)
  print(j)
}

write.table(freqOut,file=paste("allFreqs_N0_",N_0,"_K_",carryCap,"_propMaj_",propMaj,"_gensToNewOpt_",gensToOpt,"_plast_",plast,sep=""),quote=FALSE,row.names=FALSE)
write.table(qtlInfoMat,file=paste("qtlInfo_N0_",N_0,"_K_",carryCap,"_propMaj_",propMaj,"_gensToNewOpt_",gensToOpt,"_plast_",plast,sep=""),quote=FALSE,row.names=FALSE)
write.table(NeSave,file=paste("Ne_N0_",N_0,"_K_",carryCap,"_propMaj_",propMaj,"_gensToNewOpt_",gensToOpt,"_plast_",plast,sep=""),quote=FALSE,row.names=FALSE)
write.table(NMat,file=paste("popSize_N0_",N_0,"_K_",carryCap,"_propMaj_",propMaj,"_gensToNewOpt_",gensToOpt,"_plast_",plast,sep=""),quote=FALSE,row.names=FALSE)
write.table(phenMat,file=paste("phen_N0_",N_0,"_K_",carryCap,"_propMaj_",propMaj,"_gensToNewOpt_",gensToOpt,"_plast_",plast,sep=""),quote=FALSE,row.names=FALSE)
write.table(h2Mat,file=paste("h2_N0_",N_0,"_K_",carryCap,"_propMaj_",propMaj,"_gensToNewOpt_",gensToOpt,"_plast_",plast,sep=""),quote=FALSE,row.names=FALSE)
write.table(optPhenMat,file=paste("optPhen_N0_",N_0,"_K_",carryCap,"_propMaj_",propMaj,"_gensToNewOpt_",gensToOpt,"_plast_",plast,sep=""),quote=FALSE,row.names=FALSE)

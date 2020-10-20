# This script will reproduce the figures in Kardos and Luikart (2020)

##################################################################
#-----------------------------------------------------------------
# fitness function
# Figure 1
#-----------------------------------------------------------------
##################################################################


sdZ <- sqrt(10)
vZ <- 10
zVec <- seq(90,112,0.01)
zBar <- 100
Pz <- rep(NA,length(zVec))
fitZ <- rep(NA,length(zVec))
wMax <- 1.5                    # maximum intrinsic fitness
theta <- 110
c <- 6    #sd fitness function
for(i in 1:length(Pz)){
  Pz [i] <- (1/(sdZ*sqrt(2*pi)))*exp(-((zVec[i]-zBar)^2)/(2*vZ))  # phenotype distribution
  fitZ [i] <- wMax*exp(-(((zVec[i] - theta)^2)/(2*(c^2))))
}

# plot the phenotype density
par(mar=c(5,5,1,5),xpd=TRUE)
plot(Pz~zVec,ylab="",xlab="Phenotype",type="n",cex.lab=1.3)
lines(Pz~zVec,lwd=3,col="darkred")

# add mean phenotype
lines(c(100,100),c(0,0.13),lty="dashed")

# add another vertical axis
axis(side=4,at=seq(0,0.12,0.12/5),labels=seq(0,1.5,1.5/5))
text(x=117,y=0.065,labels="Fitness",srt=90,cex=1.3,col="darkgray")
text(x=85,y=0.065,labels="Probability",srt=90,cex=1.3,col="darkred")

# add the fitness function
lines(0.12*(fitZ/1.5)  ~zVec,lty="dashed",lwd=3,col="darkgray")
lines(c(110,110),c(0,0.13),lty="dotted")


##################################################################
#-----------------------------------------------------------------
# fitness function
# Figure 1
#-----------------------------------------------------------------
##################################################################


sdZ <- sqrt(10)
vZ <- 10
zVec <- seq(90,112,0.01)
zBar <- 100
Pz <- rep(NA,length(zVec))
fitZ <- rep(NA,length(zVec))
fitZ1 <- rep(NA,length(zVec))
theta1 <- 100
wMax <- 1.5                    # maximum intrinsic fitness
theta <- 110
c <- 6    #sd fitness function
for(i in 1:length(Pz)){
  Pz [i] <- (1/(sdZ*sqrt(2*pi)))*exp(-((zVec[i]-zBar)^2)/(2*vZ))  # phenotype distribution
  fitZ [i] <- wMax*exp(-(((zVec[i] - theta)^2)/(2*(c^2))))
  fitZ1 [i] <- wMax*exp(-(((zVec[i] - theta1)^2)/(2*(c^2))))
}

# plot the phenotype density
par(mar=c(5,5,1,5),xpd=TRUE)
plot(Pz~zVec,ylab="",xlab="Phenotype",type="n",cex.lab=1.3)
lines(Pz~zVec,lwd=3,col="darkred")
lines(0.12*(fitZ1/1.5)  ~zVec,lty="dashed",lwd=3,col="darkgray")
# add mean phenotype
lines(c(100,100),c(0,0.13),lty="dashed")


# add another vertical axis
axis(side=4,at=seq(0,0.12,0.12/5),labels=seq(0,1.5,1.5/5))
text(x=117,y=0.065,labels="Fitness",srt=90,cex=1.3,col="darkgray")
text(x=85,y=0.065,labels="Probability",srt=90,cex=1.3,col="darkred")





#------------------------------------------------
# move the optimum phenotype
#------------------------------------------------
# plot the phenotype density
par(mar=c(5,5,1,5),xpd=TRUE)
plot(Pz~zVec,ylab="",xlab="Phenotype",type="n",cex.lab=1.3)
lines(Pz~zVec,lwd=3,col="darkred")
lines(0.12*(fitZ/1.5)  ~zVec,lty="dashed",lwd=3,col="darkgray")
# add mean phenotype
lines(c(100,100),c(0,0.13),lty="dashed")


# add another vertical axis
axis(side=4,at=seq(0,0.12,0.12/5),labels=seq(0,1.5,1.5/5))
text(x=117,y=0.065,labels="Fitness",srt=90,cex=1.3,col="darkgray")
text(x=85,y=0.065,labels="Probability",srt=90,cex=1.3,col="darkred")


# add the fitness function

lines(c(110,110),c(0,0.13),lty="dotted")






#####################################################
# plot deterministic dynamics with 
# initial heritability 0.6
# Figure 2
#####################################################

library(quantPop)
startN <- 500
carryCap <- 1000
runTime <- 80
intGrowth <- 1.5                 # lambda for a perfectly adapted population at N -> 0.
startPheno <- 100
optPheno <- 110
sdFitFun <- 6
h2Start <- 0.6
phenVar <- 10
startFreq <- c(0.1,0.25,0.5,0.75,0.9)
t <- runTime

######### major locus

majorNMat <- NULL
majorPhenMat <- NULL
majorH2Mat <- NULL
majorPMat <- NULL
for(i in 1:length(startFreq)){
  logQuant(N0=startN,K=carryCap,t=runTime,lambda=intGrowth,nLoci=1,
           phen_0=startPheno,phen_opt=optPheno,fit_sd=sdFitFun,h2_0=h2Start,Vp_0=phenVar,p0=startFreq[i])
  majorNMat <- rbind(majorNMat,NVec)
  majorPhenMat <- rbind(majorPhenMat,phenVec)
  majorH2Mat <- rbind(majorH2Mat,h2Vec)
  majorPMat <- rbind(majorPMat,freqVec)
  NVec <- NULL
  phenVec <- NULL
  h2Vec <- NULL
  freqVec <- NULL
}

######### 2 major loci

major2NMat <- NULL
major2PhenMat <- NULL
major2H2Mat <- NULL
major2PMat <- NULL
for(i in 1:length(startFreq)){
  logQuant(N0=startN,K=carryCap,t=runTime,lambda=intGrowth,nLoci=2,
           phen_0=startPheno,phen_opt=optPheno,fit_sd=sdFitFun,h2_0=h2Start,Vp_0=phenVar,p0=startFreq[i])
  major2NMat <- rbind(major2NMat,NVec)
  major2PhenMat <- rbind(major2PhenMat,phenVec)
  major2H2Mat <- rbind(major2H2Mat,h2Vec)
  major2PMat <- rbind(major2PMat,freqVec)
  NVec <- NULL
  phenVec <- NULL
  h2Vec <- NULL
  freqVec <- NULL
}


######## polygenic 
polyNMat <- NULL
polyPhenMat <- NULL
polyH2Mat <- NULL
polyPMat <- NULL
for(i in 1:length(startFreq)){
  logQuant(N0=startN,K=carryCap,t=runTime,lambda=intGrowth,nLoci=100,
           phen_0=startPheno,phen_opt=optPheno,fit_sd=sdFitFun,h2_0=h2Start,Vp_0=phenVar,p0=startFreq[i])
  polyNMat <- rbind(polyNMat,NVec)
  polyPhenMat <- rbind(polyPhenMat,phenVec)
  polyH2Mat <- rbind(polyH2Mat,h2Vec)
  polyPMat <- rbind(polyPMat,freqVec)
}

#########################################################################
# make heritabilities and phenotypes missing for generations with N < 2
#########################################################################

for(i in 1:nrow(majorNMat)){
  if(sum(majorNMat[i,] < 2) > 0){
    majorH2Mat[i,which(majorNMat[i,] < 2)-1] <- NA
    majorPhenMat[i,which(majorNMat[i,] < 2)-1] <- NA
    majorNMat[i,which(majorNMat[i,] < 2)-1] <- NA
  }
  if(sum(major2NMat[i,] < 2) > 0){
    major2H2Mat[i,which(major2NMat[i,] < 2)-1] <- NA
    major2PhenMat[i,which(major2NMat[i,] < 2)-1] <- NA
    major2NMat[i,which(major2NMat[i,] < 2)-1] <- NA
  }
  if(sum(polyNMat[i,] < 2) > 0){
    polyH2Mat[i,which(polyNMat[i,] < 2)-1] <- NA
    polyPhenMat[i,which(polyNMat[i,] < 2)-1] <- NA
    polyNMat[i,which(polyNMat[i,] < 2)-1] <- NA
  }
}


########################################
#################### Make the figure
########################################

# population size
par(mfrow=c(3,3),xpd=TRUE,mar=c(4,5,1,2))
plot(1:t,majorNMat[1:t],ylim=c(0,carryCap),type="n",cex.axis=0.8,
     cex.lab=1.5,ylab=expression(italic(""*N*"")),xlab="")
text(x=-22,y=1050,labels="A",cex=2)
lineVec <- c("solid","dashed","dotted","dotdash","twodash")
for(i in 1:nrow(majorNMat)){
  lines(1:t,majorNMat[i,1:t],lty=lineVec[i])
}

plot(1:t,major2NMat[1:t],ylim=c(0,carryCap),type="n",cex.axis=0.8,
     cex.lab=1.5,ylab="",xlab="")
text(x=-22,y=1050,labels="B",cex=2)
lineVec <- c("solid","dashed","dotted","dotdash","twodash")
for(i in 1:nrow(major2NMat)){
  lines(1:t,major2NMat[i,1:t],lty=lineVec[i])
}

plot(1:t,polyNMat[1:t],ylim=c(0,carryCap),type="n",cex.axis=0.8,
     cex.lab=1.5,ylab="",xlab="")
text(x=-22,y=1050,labels="C",cex=2)
for(i in 1:nrow(polyNMat)){
  lines(1:t,polyNMat[i,1:t],lty=lineVec[i])
}
legend(x=30,y=0,lty=lineVec,legend=c(expression(paste(italic(""*p*"")[0]," = 0.1",sep="")),
                                     expression(paste(italic(""*p*"")[0]," = 0.25",sep="")),
                                     expression(paste(italic(""*p*"")[0]," = 0.5",sep="")),
                                     expression(paste(italic(""*p*"")[0]," = 0.75",sep="")),
                                     expression(paste(italic(""*p*"")[0]," = 0.9",sep=""))),
       bty="n",xjust=FALSE,yjust=FALSE,cex=1)

# heritability
plot(1:t,majorNMat[1:t],ylim=c(0,0.8),type="n",cex.axis=0.8,
     cex.lab=1.5,ylab=expression(italic(""*h*"")^2),xlab="")
for(i in 1:nrow(majorH2Mat)){
  lines(1:t,majorH2Mat[i,1:t],lty=lineVec[i])
}
plot(1:t,major2NMat[1:t],ylim=c(0,0.8),type="n",cex.axis=0.8,
     cex.lab=1.5,ylab="",xlab="")
for(i in 1:nrow(major2H2Mat)){
  lines(1:t,major2H2Mat[i,1:t],lty=lineVec[i])
}
plot(1:t,polyH2Mat[1:t],ylim=c(0,0.8),type="n",cex.axis=0.8,
     cex.lab=1.5,ylab="",xlab="")
for(i in 1:nrow(polyH2Mat)){
  lines(1:t,polyH2Mat[i,1:t],lty=lineVec[i])
}


# phenotype
plot(1:t,majorNMat[1:t],ylim=c(100,111),type="n",cex.axis=0.8,
     cex.lab=1.5,ylab="Mean phenotype",xlab="Generation")
for(i in 1:nrow(majorPhenMat)){
  lines(1:t,majorPhenMat[i,1:t],lty=lineVec[i])
}
plot(1:t,majorNMat[1:t],ylim=c(100,111),type="n",cex.axis=0.8,
     cex.lab=1.5,ylab="",xlab="Generation")
for(i in 1:nrow(major2PhenMat)){
  lines(1:t,major2PhenMat[i,1:t],lty=lineVec[i])
}
plot(1:t,majorNMat[1:t],ylim=c(100,111),type="n",cex.axis=0.8,
     cex.lab=1.5,ylab="",xlab="Generation")
for(i in 1:nrow(polyPhenMat)){
  lines(1:t,polyPhenMat[i,1:t],lty=lineVec[i])
}






####################################################################
#-------------------------------------------------------------------
# fixed p0, single locus versus polygenic for mammals with h20 = 0.6,
# individual-based simulations
# FIGURE 3
#-------------------------------------------------------------------
####################################################################
# get polygenic data
setwd("~/Documents/work/rapidAdaptation/AmNatRevision/simulations/originalSims/constantP0/h2Pt6/noMajor")
p0Vec <- c(0.1,0.25,0.5,0.75,0.9)

popFileNames <- paste("popSize_p0_",p0Vec,"_N0_500_K_1000",sep="")
phenFileNames <- paste("phen_p0_",p0Vec,"_N0_500_K_1000",sep="")
h2FileNames <- paste("h2_p0_",p0Vec,"_N0_500_K_1000",sep="")
polyPopSize <- list()
polyPhen <- list()
polyH2 <- list()

for(i in 1:5){
  polyPopSize [[i]] <- read.table(popFileNames[i],header=TRUE)
  polyPhen [[i]] <- read.table(phenFileNames[i],header=TRUE)
  polyH2 [[i]] <- read.table(h2FileNames[i],header=TRUE)
}


# get single major locus data

setwd("~/Documents/work/rapidAdaptation/AmNatRevision/simulations/originalSims/constantP0/h2Pt6/major")
popFileNames <- paste("popSize_N0_500_K_1000_propMaj_0.9999999_freqMaj_",p0Vec,sep="")
phenFileNames <- paste("phen_N0_500_K_1000_propMaj_0.9999999_freqMaj_",p0Vec,sep="")
h2FileNames <- paste("h2_N0_500_K_1000_propMaj_0.9999999_freqMaj_",p0Vec,sep="")
majorPopSize <- list()
majorPhen <- list()
majorH2 <- list()

for(i in 1:5){
  majorPopSize [[i]] <- read.table(popFileNames[i],header=TRUE)
  majorPhen [[i]] <- read.table(phenFileNames[i],header=TRUE)
  majorH2 [[i]] <- read.table(h2FileNames[i],header=TRUE)
}



# get  two major loci data

setwd("~/Documents/work/rapidAdaptation/AmNatRevision/simulations/originalSims/constantP0/h2Pt6/twoMajor")
popFileNames <- paste("popSize_N0_500_K_1000_propMaj_0.5_freqMaj_",p0Vec,sep="")
phenFileNames <- paste("phen_N0_500_K_1000_propMaj_0.5_freqMaj_",p0Vec,sep="")
h2FileNames <- paste("h2_N0_500_K_1000_propMaj_0.5_freqMaj_",p0Vec,sep="")
twoMajorPopSize <- list()
twoMajorPhen <- list()
twoMajorH2 <- list()

for(i in 1:5){
  twoMajorPopSize [[i]] <- read.table(popFileNames[i],header=TRUE)
  twoMajorPhen [[i]] <- read.table(phenFileNames[i],header=TRUE)
  twoMajorH2 [[i]] <- read.table(h2FileNames[i],header=TRUE)
}

#-----------------------------------------------
# make NA population sizes zero
#-----------------------------------------------

for(i in 1:5){
  for(j in 1:ncol(polyPopSize[[i]])){
    polyPopSize[[i]][which(is.na(polyPopSize[[i]][,j])),j] <- 0
    majorPopSize[[i]][which(is.na(majorPopSize[[i]][,j])),j] <- 0
    twoMajorPopSize[[i]][which(is.na(twoMajorPopSize[[i]][,j])),j] <- 0
  }
}

#----------------------------------------------------------------------------
# calculate mean parmameters and extinction rate
# for major, two major an polygenic architectures
#----------------------------------------------------------------------------

twoMajorExtProp <- list()
majorExtProp <- list()
polyExtProp <- list()

twoMajorMeanPop <- list()
twoMajorMeanPhen <- list()
twoMajorMeanH2 <- list()

majorMeanPop <- list()
majorMeanPhen <- list()
majorMeanH2 <- list()

polyMeanPop <- list()
polyMeanPhen <- list()
polyMeanH2 <- list()

for(i in 1:5){
  polyMeanPop [[i]] <- colMeans(polyPopSize[[i]])
  polyMeanPhen [[i]] <- colMeans(polyPhen[[i]],na.rm=TRUE)
  polyMeanH2  [[i]] <- colMeans(polyH2[[i]],na.rm=TRUE)

  majorMeanPop [[i]] <- colMeans(majorPopSize[[i]])
  majorMeanPhen [[i]] <- colMeans(majorPhen[[i]],na.rm=TRUE)
  majorMeanH2  [[i]] <- colMeans(majorH2[[i]],na.rm=TRUE)

  twoMajorMeanPop [[i]] <- colMeans(twoMajorPopSize[[i]])
  twoMajorMeanPhen [[i]] <- colMeans(twoMajorPhen[[i]],na.rm=TRUE)
  twoMajorMeanH2  [[i]] <- colMeans(twoMajorH2[[i]],na.rm=TRUE)

  polyExtProp [[i]] <- colSums(polyPopSize[[i]] == 0)/500
  majorExtProp [[i]] <- colSums(majorPopSize[[i]] == 0)/500
  twoMajorExtProp [[i]] <- colSums(twoMajorPopSize[[i]] == 0)/500
}



########################################
# Make the figure
########################################

# population size
t <- 81
carryCap <- 1000
# one major locus
par(mfrow=c(4,3),xpd=TRUE,mar=c(2,5,1.5,1))
plot(1:t,majorMeanPop[[1]],ylim=c(0,carryCap),type="n",cex.axis=0.7, main ="1 major locus",cex.main=1.4,
     cex.lab=1.5,ylab=expression(italic(""*N*"")),xlab="")
text(x=-30,y=1100,labels="A",cex=2)
lineVec <- c("solid","dashed","dotted","dotdash","twodash")
for(i in 1:length(majorMeanPop)){
  lines(1:t,majorMeanPop[[i]],lty=lineVec[i],lwd=1.3)
}

# two major loci
plot(1:t,twoMajorMeanPop[[1]],ylim=c(0,carryCap),type="n",cex.axis=0.7, main ="2 major loci",cex.main=1.4,
     cex.lab=1.5,xlab="",ylab="")
text(x=-30,y=1100,labels="B",cex=2)
lineVec <- c("solid","dashed","dotted","dotdash","twodash")
for(i in 1:length(majorMeanPop)){
  lines(1:t,twoMajorMeanPop[[i]],lty=lineVec[i],lwd=1.3)
}

# polygenic
plot(1:t,polyMeanPop[[1]],ylim=c(0,carryCap),type="n",cex.axis=0.7, main ="100 small-effect loci",cex.main=1.4,
     cex.lab=1.5,xlab="",ylab="")
text(x=-30,y=1100,labels="C",cex=2)
lineVec <- c("solid","dashed","dotted","dotdash","twodash")
for(i in 1:length(majorMeanPop)){
  lines(1:t,polyMeanPop[[i]],lty=lineVec[i],lwd=1.3)
}

############# heritability

# one major locus
par(mar=c(2,5,1,1))
plot(1:t,majorMeanPop[[1]],ylim=c(0,0.8),type="n",cex.axis=0.7,
     cex.lab=1.5,ylab=expression(italic(""*h*"")^2),xlab="")
text(x=-35,y=1100,labels="A",cex=2)
lineVec <- c("solid","dashed","dotted","dotdash","twodash")
for(i in 1:length(majorMeanPop)){
  lines(1:(t-1),majorMeanH2[[i]],lty=lineVec[i],lwd=1.3)
}

# two major loci
plot(1:t,twoMajorMeanPop[[1]],ylim=c(0,0.8),type="n",cex.axis=0.7,
     cex.lab=1.5,xlab="",ylab="")
text(x=-35,y=1100,labels="A",cex=2)
lineVec <- c("solid","dashed","dotted","dotdash","twodash")
for(i in 1:length(majorMeanPop)){
  lines(1:(t-1),twoMajorMeanH2[[i]],lty=lineVec[i],lwd=1.3)
}

# polygenic
plot(1:t,polyMeanPop[[1]],ylim=c(0,0.8),type="n",cex.axis=0.7,
     cex.lab=1.5,xlab="",ylab="")
text(x=-35,y=1100,labels="A",cex=2)
lineVec <- c("solid","dashed","dotted","dotdash","twodash")
for(i in 1:length(majorMeanPop)){
  lines(1:(t-1),polyMeanH2[[i]],lty=lineVec[i],lwd=1.3)
}


############# phenotype

# one major locus

plot(1:t,majorMeanPop[[1]],ylim=c(100,111),type="n",cex.axis=0.7,
     cex.lab=1.5,ylab="Mean phenotype",xlab="")
text(x=-35,y=1100,labels="A",cex=2)
lineVec <- c("solid","dashed","dotted","dotdash","twodash")
for(i in 1:length(majorMeanPop)){
  lines(1:(t-1),majorMeanPhen[[i]],lty=lineVec[i],lwd=1.3)
}

# two major loci
plot(1:t,majorMeanPop[[1]],ylim=c(100,111),type="n",cex.axis=0.7,
     cex.lab=1.5,ylab="",xlab="")
text(x=-35,y=1100,labels="A",cex=2)
lineVec <- c("solid","dashed","dotted","dotdash","twodash")
for(i in 1:length(majorMeanPop)){
  lines(1:(t-1),twoMajorMeanPhen[[i]],lty=lineVec[i],lwd=1.3)
}

# polygenic
plot(1:t,majorMeanPop[[1]],ylim=c(100,111),type="n",cex.axis=0.7,
     cex.lab=1.5,ylab="",xlab="")
text(x=-35,y=1100,labels="A",cex=2)
lineVec <- c("solid","dashed","dotted","dotdash","twodash")
for(i in 1:length(majorMeanPop)){
  lines(1:(t-1),polyMeanPhen[[i]],lty=lineVec[i],lwd=1.3)
}

legend(x=30,y=100,lty=lineVec,legend=c(expression(paste(italic(""*p*"")[0]," = 0.1",sep="")),
                                       expression(paste(italic(""*p*"")[0]," = 0.25",sep="")),
                                       expression(paste(italic(""*p*"")[0]," = 0.5",sep="")),
                                       expression(paste(italic(""*p*"")[0]," = 0.75",sep="")),
                                       expression(paste(italic(""*p*"")[0]," = 0.9",sep=""))),
       bty="n",xjust=FALSE,yjust=FALSE)



#---------------------------------------------------
# plot extinction rates for the three architectures
#---------------------------------------------------

par(mar=c(4,5,1,1))
# single major locus
plot(c(0,80),y=c(0,1),type="n",xlab="Generation",ylab="Proportion extinct",
     cex.lab=1.5,main ="",cex.axis=0.7)
for(i in 1:5){
  lines(1:81,majorExtProp[[i]],lty=lineVec[i],lwd=1.3)
}


# two major loci
plot(c(0,80),y=c(0,1),type="n",xlab="Generation",
     ylab="",cex.lab=1.5,main ="",cex.axis=0.7)
for(i in 1:5){
  lines(1:81,twoMajorExtProp[[i]],lty=lineVec[i],lwd=1.3)
}

# 100 loci
plot(c(0,80),y=c(0,1),type="n",xlab="Generation",
     ylab="",cex.lab=1.5,main ="",cex.axis=0.7)
for(i in 1:5){
  lines(1:81,polyExtProp[[i]],lty=lineVec[i],lwd=1.3)
}



########################################################
########################################################
# analyze results where the loci have variable
# p0 for coral and large mammals and the large effect
# locus is responsible for 90% of the genetic variance
# FIGURE 4
########################################################
########################################################

propAnalyze <- 0.9
####----------------------
#### large mammals first
####----------------------
setwd("~/Documents/work/rapidAdaptation/AmNatRevision/simulations/originalSims/h2Pt6/fecundity4/major")

popFileName <- paste("popSize_N0_500_K_1000_propMaj_",propAnalyze,sep="")
phenFileName <-  paste("phen_N0_500_K_1000_propMaj_",propAnalyze,sep="")
h2FileName <-  paste("h2_N0_500_K_1000_propMaj_",propAnalyze,sep="")

mam_majorPopSize <- read.table(popFileName,header=TRUE)
mam_majorPhen <- read.table(phenFileName,header=TRUE)
mam_majorH2 <- read.table(h2FileName,header=TRUE)


setwd("~/Documents/work/rapidAdaptation/AmNatRevision/simulations/originalSims/h2Pt6/fecundity4/noMajor")

mam_noMajorPopSize <- read.table("popSize_N0_500_K_1000_propMaj_0.01",header=TRUE)
mam_noMajorPhen <- read.table("phen_N0_500_K_1000_propMaj_0.01",header=TRUE)
mam_noMajorH2 <- read.table("h2_N0_500_K_1000_propMaj_0.01",header=TRUE)

#-----------------------------------------------
# make NA population sizes zero
#-----------------------------------------------


for(j in 1:ncol(mam_noMajorPopSize)){
  mam_noMajorPopSize[which(is.na(mam_noMajorPopSize[,j])),j] <- 0
  mam_majorPopSize[which(is.na(mam_majorPopSize[,j])),j] <- 0
}


#######################
# check for extinction
#######################

mam_polyExt <- rep(NA,81)
mam_OligExt <- rep(NA,81)

for(i in 1:81){
  mam_polyExt[i] <- sum(mam_noMajorPopSize[,i] == 0)
  mam_OligExt[i] <- sum(mam_majorPopSize[,i] == 0)
}


#############################################
# get bootstrap CIs for proportion extinct
#############################################

bootReps <- 1000
lowPolyBoots <- matrix(NA,nrow=bootReps,ncol=81)
lowOligBoots <- matrix(NA,nrow=bootReps,ncol=81)

for(i in 1:bootReps){
  lowPolyDat <- mam_noMajorPopSize[sample(1:nrow(mam_noMajorPopSize),nrow(mam_noMajorPopSize),replace=TRUE),]
  lowOligDat <- mam_majorPopSize[sample(1:nrow(mam_majorPopSize),nrow(mam_majorPopSize),replace=TRUE),]
  lowPolyBoots [i,] <- colSums(lowPolyDat == 0)/nrow(lowPolyDat )
  lowOligBoots [i,] <- colSums(lowOligDat == 0)/nrow(lowOligDat )
  print(i)
}

mamlowPolyCI <- NULL
mamlowOligCI <- NULL

for(i in 1:81){
  mamlowPolyCI <- cbind(mamlowPolyCI,quantile(lowPolyBoots[,i],c(0.025,0.975)))
  mamlowOligCI <- cbind(mamlowOligCI,quantile(lowOligBoots[,i],c(0.025,0.975)))
}


#-----------------------------------------
# now corals
#-----------------------------------------

setwd("~/Documents/work/rapidAdaptation/AmNatRevision/simulations/originalSims/h2Pt6/fecundity26/N010000/major")

popFileName <- paste("popSize_N0_10000_K_20000_propMaj_0.9",sep="")
phenFileName <-  paste("phen_N0_10000_K_20000_propMaj_0.9",sep="")
h2FileName <-  paste("h2_N0_10000_K_20000_propMaj_0.9",sep="")

cor_majorPopSize <- read.table(popFileName,header=TRUE)
cor_majorPhen <- read.table(phenFileName,header=TRUE)
cor_majorH2 <- read.table(h2FileName,header=TRUE)


setwd("~/Documents/work/rapidAdaptation/AmNatRevision/simulations/originalSims/h2Pt6/fecundity26/N010000/noMajor")


cor_noMajorPopSize <- read.table("popSize_N0_10000_K_20000_propMaj_0.01",header=TRUE)
cor_noMajorPhen <- read.table("phen_N0_10000_K_20000_propMaj_0.01",header=TRUE)
cor_noMajorH2 <- read.table("h2_N0_10000_K_20000_propMaj_0.01",header=TRUE)

#-----------------------------------------------
# make NA population sizes zero
#-----------------------------------------------


for(j in 1:ncol(cor_noMajorPopSize)){
  cor_noMajorPopSize[which(is.na(cor_noMajorPopSize[,j])),j] <- 0
  cor_majorPopSize[which(is.na(cor_majorPopSize[,j])),j] <- 0
}


#######################
# check for extinction
#######################
cor_polyExt <- rep(NA,81)
cor_OligExt <- rep(NA,81)

for(i in 1:81){
  cor_polyExt[i] <- sum(cor_noMajorPopSize[,i] == 0)
  cor_OligExt[i] <- sum(cor_majorPopSize[,i] == 0)
}


#############################################
# get bootstrap CIs for proportion extinct
#############################################

bootReps <- 1000
cor_lowPolyBoots <- matrix(NA,nrow=bootReps,ncol=81)
cor_lowOligBoots <- matrix(NA,nrow=bootReps,ncol=81)

for(i in 1:bootReps){
  lowPolyDat <- cor_noMajorPopSize[sample(1:nrow(cor_noMajorPopSize),nrow(cor_noMajorPopSize),replace=TRUE),]
  lowOligDat <- cor_majorPopSize[sample(1:nrow(cor_majorPopSize),nrow(cor_majorPopSize),replace=TRUE),]
  cor_lowPolyBoots [i,] <- colSums(lowPolyDat == 0)/nrow(lowPolyDat )
  cor_lowOligBoots [i,] <- colSums(lowOligDat == 0)/nrow(lowOligDat )
  print(i)
}

corlowPolyCI <- NULL
corlowOligCI <- NULL

for(i in 1:81){
  corlowPolyCI <- cbind(corlowPolyCI,quantile(cor_lowPolyBoots[,i],c(0.025,0.975)))
  corlowOligCI <- cbind(corlowOligCI,quantile(cor_lowOligBoots[,i],c(0.025,0.975)))
}

############################################
# plot results
############################################
par(mfrow=c(3,2),mar=c(2,5,3,0.5),xpd=FALSE)
library(scales)
carryCap <- 1000
##### mammal population size
plot(c(0,81),c(0,carryCap + 200)/1000,type="n",ylab=expression(paste(italic(""*N*"")," (thousands)",sep="")),xlab="",cex.lab=1.5,
     cex.main=1.5,cex.axis=0.9)
for(i in 1:nrow(mam_noMajorPopSize)){
  lines(1:81,mam_noMajorPopSize[i,]/1000,col=alpha("orange",alpha=0.05))
  lines(1:81,mam_majorPopSize[i,]/1000,col=alpha("darkblue",alpha=0.05))
}
polySizeMean <- rep(NA,81)
oligSizeMean <- rep(NA,81)
for(i in 1:81){
  polySizeMean[i] <- mean(mam_noMajorPopSize[,i])/1000
  oligSizeMean[i] <- mean(mam_majorPopSize[,i])/1000
}
lines(1:81,polySizeMean,lwd=3,col="orange")
lines(1:81,oligSizeMean,lwd=3,col="darkblue")
lines(c(0,200),c(carryCap/1000,carryCap/1000),lty="dashed")
par(xpd=TRUE)
text(x=-25,y=1.4,labels="A",cex=2)
par(xpd=FALSE)
##### coral population size
plot(c(0,81),c(0,20000 + 4000)/1000,type="n",ylab="",xlab="",cex.lab=1.5,
     cex.main=1.5,cex.axis=0.9)
carryCap <- 20000
for(i in 1:nrow(cor_noMajorPopSize)){
  lines(1:81,cor_noMajorPopSize[i,]/1000,col=alpha("orange",alpha=0.05))
  lines(1:81,cor_majorPopSize[i,]/1000,col=alpha("darkblue",alpha=0.05))
}
polySizeMean <- rep(NA,81)
oligSizeMean <- rep(NA,81)
for(i in 1:81){
  polySizeMean[i] <- mean(cor_noMajorPopSize[,i]/1000)
  oligSizeMean[i] <- mean(cor_majorPopSize[,i]/1000)
}
lines(1:81,polySizeMean,lwd=3,col="orange")
lines(1:81,oligSizeMean,lwd=3,col="darkblue")
lines(c(0,200),c(carryCap/1000,carryCap/1000),lty="dashed")
par(xpd=TRUE)
text(x=-25,y=28,labels="B",cex=2)
par(xpd=FALSE)

#---------------------------------
##### plot the phenotypes
#---------------------------------
par(mar=c(3,5,2,0.5))
idealPhen <- 110

####### mammals
plot(c(0,80),c(100,112),type="n",ylab="Mean Phenotype",xlab="",cex.lab=1.5,cex.axis=0.8)
for(i in 1:500){
  lines(1:80,mam_noMajorPhen[i,],col=alpha("orange",alpha=0.05))
  lines(1:80,mam_majorPhen[i,],col=alpha("darkblue",alpha=0.05))
}

# get rolling mean
polyPhenMean <- rep(NA,80)
oligPhenMean <- rep(NA,80)
for(i in 1:80){
  polyPhenMean[i] <- mean(mam_noMajorPhen[,i],na.rm=TRUE)
  oligPhenMean[i] <- mean(mam_majorPhen[,i],na.rm=TRUE)
}
lines(1:80,polyPhenMean,lwd=3,col="orange")
lines(1:80,oligPhenMean,lwd=3,col="darkblue")
lines(c(0,80),c(110,110),lty="dashed")


par(mar=c(3,5,2,0.5))
idealPhen <- 110

####### corals
plot(c(0,80),c(100,112),type="n",ylab="",xlab="",cex.lab=1.5,cex.axis=0.8)
for(i in 1:500){
  lines(1:80,cor_noMajorPhen[i,],col=alpha("orange",alpha=0.05))
  lines(1:80,cor_majorPhen[i,],col=alpha("darkblue",alpha=0.05))
}

# get rolling mean
polyPhenMean <- rep(NA,80)
oligPhenMean <- rep(NA,80)
for(i in 1:80){
  polyPhenMean[i] <- mean(cor_noMajorPhen[,i],na.rm=TRUE)
  oligPhenMean[i] <- mean(cor_majorPhen[,i],na.rm=TRUE)
}
lines(1:80,polyPhenMean,lwd=3,col="orange")
lines(1:80,oligPhenMean,lwd=3,col="darkblue")
lines(c(0,80),c(110,110),lty="dashed")


#####################################################################
# plot the fraction of extinct population versus time
#####################################################################
par(mar=c(4,5,1,0.5))
# extinct proportion for mammals
plot(c(1,81),c(0,1),type="n",xlab="Generation",ylab="Proportion Extinct",cex.lab=1.5)
lines(1:81,mam_polyExt/nrow(lowPolyDat),col="orange",lwd=3)
lines(1:81,mam_OligExt/nrow(lowOligDat),col="darkblue",lwd=3)

polygon(x=c(1:81,rev(1:81)),y=c(mamlowPolyCI[1,1:81],rev(mamlowPolyCI[2,1:81])),col=adjustcolor("orange",alpha.f=0.5),border=NA)
polygon(x=c(1:81,rev(1:81)),y=c(mamlowOligCI[1,1:81],rev(mamlowOligCI[2,1:81])),col=adjustcolor("darkblue",alpha.f=0.5),border=NA)
# get pictures

setwd("~/Documents/work/rapidAdaptation/AmNatRevision/simulations/originalSims/h2Pt6")


library(jpeg)
img<-readJPEG("wildebeest.jpg")
rasterImage(img,xleft=-3, ybottom=0.8, xright=25, ytop=1.1)


# corals


par(mar=c(4,5,1,0.5))
# extinct proportion for mammals
plot(c(1,81),c(0,1),type="n",xlab="Generation",ylab="",cex.lab=1.5)
lines(1:81,cor_polyExt/nrow(lowPolyDat),col="orange",lwd=3)
lines(1:81,cor_OligExt/nrow(lowOligDat),col="darkblue",lwd=3)

legend(x=20,y=0.0,xjust=FALSE,yjust=FALSE,legend=c("No Major Locus","Major Locus"),
       lwd=3,col=c("orange","darkblue"),bty="n")

polygon(x=c(1:81,rev(1:81)),y=c(corlowPolyCI[1,1:81],rev(corlowPolyCI[2,1:81])),col=adjustcolor("orange",alpha.f=0.5),border=NA)
polygon(x=c(1:81,rev(1:81)),y=c(corlowOligCI[1,1:81],rev(corlowOligCI[2,1:81])),col=adjustcolor("darkblue",alpha.f=0.5),border=NA)

img<-readJPEG("coralPic.jpg")
rasterImage(img,xleft=-3, ybottom=0.8, xright=25, ytop=1.1)


##############################################################################
##############################################################################
# plot population size versus initial large-effect allele frequency
# FIGURE 5
##############################################################################
##############################################################################

############ mammals
setwd("~/Documents/work/rapidAdaptation/AmNatRevision/simulations/originalSims/h2Pt6/fecundity4/major")
mammPopFiles <- c("popSize_N0_500_K_1000_propMaj_0.9",paste("popSize_N0_500_K_1000_propMaj_0.9","_",2:4,sep=""))
mammQtlFiles <- c("qtlInfo_N0_500_K_1000_propMaj_0.9",paste("qtlInfo_N0_500_K_1000_propMaj_0.9","_",2:4,sep=""))

popDat <- NULL
qtlDat <- list()
qtlIter <- 1
for(i in 1:4){
  thisPopDat <- read.table(mammPopFiles[i],header=TRUE,row.names=NULL)
  popDat <- rbind(popDat,thisPopDat)
  thisQtlDat <- read.table(mammQtlFiles[i],header=TRUE,row.names=NULL)
  for(j in 1:500){
    qtlDat[[qtlIter]] <- thisQtlDat[which(thisQtlDat[,1] == j),]
    qtlIter <- qtlIter + 1
  }
}

# get the frequency of the major allele
mamm_majFreq0 <- rep(NA,length(qtlDat))
for (i in 1:length(qtlDat)){
  mamm_majFreq0[i] <-  qtlDat[[i]][which(qtlDat[[i]][,4] == max(qtlDat[[i]][,4])),2]
}

##### get rid of simulation reps with extreme allele frequencies
extremeFreqs <- which(mamm_majFreq0 < 0.05 | mamm_majFreq0 > 0.95)
mamm_majFreq0 <- mamm_majFreq0[-extremeFreqs]

# zero out NA population sizes
for(i in 1:81){
  if(sum(is.na(popDat[,i]) > 0)){popDat[which(is.na(popDat[,i])),i] <- 0}
}

mammFinalSizes <- popDat[,81]
mammFinalSizes <- mammFinalSizes[-extremeFreqs]
#------------------------------------------------------------------------------
# get bootstrap CIs for mean population size across initial allele frequencies
#------------------------------------------------------------------------------
boots <- 1000
stepSize <- 0.05
freqStarts <- seq(0.05,0.9,stepSize)
freqEnds <- freqStarts + stepSize

mamm_popMeans <- rep(NA,length(freqStarts))
mamm_bootLow <- rep(NA,length(freqStarts))
mamm_bootHigh  <- rep(NA,length(freqStarts))

for(i in 1:length(freqStarts)){
  thesePops <- mammFinalSizes[which(mamm_majFreq0 >freqStarts[i] & mamm_majFreq0 <= freqEnds[i])]
  mamm_popMeans[i] <- mean(thesePops)
  mamm_bootEsts <- rep(NA,boots)
  for(j in 1:boots){
    theseSamps <- sample(thesePops,length(thesePops),replace=TRUE)
    mamm_bootEsts[j] <- mean(theseSamps)
  }
  quants <- quantile(mamm_bootEsts,probs=c(0.05,0.95))
  mamm_bootHigh[i] <- quants[2]
  mamm_bootLow[i] <- quants[1]
}



############ corals
setwd("~/Documents/work/rapidAdaptation/AmNatRevision/simulations/originalSims/h2Pt6/fecundity26/major")
corPopFiles <- paste("popSize_N0_10000_K_20000_propMaj_0.9_",2:11,sep="")
corQtlFiles <- paste("qtlInfo_N0_10000_K_20000_propMaj_0.9_",2:11,sep="")

popDat <- NULL
qtlDat <- list()
qtlIter <- 1
for(i in 1:10){
  thisPopDat <- read.table(corPopFiles[i],header=TRUE,row.names=NULL)
  popDat <- rbind(popDat,thisPopDat)
  thisQtlDat <- read.table(corQtlFiles[i],header=TRUE,row.names=NULL)
  for(j in 1:200){
    qtlDat[[qtlIter]] <- thisQtlDat[which(thisQtlDat[,1] == j),]
    qtlIter <- qtlIter + 1
  }
}

# get the frency of the major allele
cor_majFreq0 <- rep(NA,length(qtlDat))
for (i in 1:length(qtlDat)){
  cor_majFreq0[i] <-  qtlDat[[i]][which(qtlDat[[i]][,4] == max(qtlDat[[i]][,4])),2]
}

# zero out NA population sizes
for(i in 1:80){
  if(sum(is.na(popDat[,i]) > 0)){popDat[which(is.na(popDat[,i])),i] <- 0}
}

corFinalSizes <- popDat[,80]



#------------------------------------------------------------------------------
# get bootstrap CIs for mean population size across initial allele frequencies
#------------------------------------------------------------------------------
boots <- 1000
stepSize <- 0.05
freqStarts <- seq(0.05,0.9,stepSize)
freqEnds <- freqStarts + stepSize

cor_popMeans <- rep(NA,length(freqStarts))
cor_bootLow <- rep(NA,length(freqStarts))
cor_bootHigh  <- rep(NA,length(freqStarts))
for(i in 1:length(freqStarts)){
  thesePops <- corFinalSizes[which(cor_majFreq0 >freqStarts[i] & cor_majFreq0 <= freqEnds[i])]
  cor_popMeans[i] <- mean(thesePops)
  cor_bootEsts <- rep(NA,boots)
  for(j in 1:boots){
    theseSamps <- sample(thesePops,length(thesePops),replace=TRUE)
    cor_bootEsts[j] <- mean(theseSamps)
  }
  quants <- quantile(cor_bootEsts,probs=c(0.05,0.95))
  cor_bootHigh[i] <- quants[2]
  cor_bootLow[i] <- quants[1]
}

#-----------------------------------
# make the plot
#-----------------------------------
library(scales)
par(mfrow=c(1,2),mar=c(5,5,3,1))

#mammals
plot(mamm_majFreq0,mammFinalSizes,pch=1,col=alpha("darkgray",alpha=0.3),cex.lab=1.3,cex.axis=0.7,
     ylab=expression(paste("Final population size (",italic(""*N*"")[80],")",sep="")),xlab=expression(italic(""*p*"")[0]))

lines(freqEnds - 0.5*stepSize,mamm_popMeans,lwd=2)
lines(freqEnds - 0.5*stepSize,mamm_bootLow,lty="dashed",lwd=2)
lines(freqEnds - 0.5*stepSize,mamm_bootHigh,lty="dashed",lwd=2)

setwd("~/Documents/work/rapidAdaptation/AmNatRevision/simulations/originalSims/h2Pt6")
library(jpeg)
img<-readJPEG("wildebeest.jpg")
rasterImage(img,xleft=.7, ybottom=990, xright=1, ytop=1300)

#corals
par(mar=c(5,3,3,3))
plot(cor_majFreq0,corFinalSizes,pch=1,col=alpha("darkgray",alpha=0.3),cex.lab=1.3,cex.axis=0.7,
     ylab="",xlab=expression(italic(""*p*"")[0]))

lines(freqEnds - 0.5*stepSize,cor_popMeans,lwd=2)
lines(freqEnds - 0.5*stepSize,cor_bootLow,lty="dashed",lwd=2)
lines(freqEnds - 0.5*stepSize,cor_bootHigh,lty="dashed",lwd=2)

setwd("~/Documents/work/rapidAdaptation/AmNatRevision/simulations/originalSims/h2Pt6")
library(jpeg)
img<-readJPEG("coralPic.jpg")
rasterImage(img,xleft=.7, ybottom=17000, xright=1, ytop=23000)


###########################################################################
#--------------------------------------------------------------------------
# plot the phenotype distributions for different numbers of loci
# to see how well they conform to the assumption of a normal distribiution
# Figure S1
#--------------------------------------------------------------------------
###########################################################################
phen_0 <- 100
h2_0<-0.6
Vp_0<-10
p0<-startFreq[i]
pseudoN<-50000
Ve <- Vp_0 - (h2_0*Vp_0)
par(mfrow=c(3,3))
startFreq <- c(0.1,0.25,0.5)
locusNums <- c(1,2,100)

i <- 1   # we're only considering the phenotype distribution in the first generation

for(j in 1:length(locusNums)){    # loop through the desired numbers of loci
  nLoci <- locusNums[j]
  for(k in 1:length(startFreq)){  # loop throught the desired starting frequencies of the positively selected allele(s)
    # the following genotpes and phenotype are directly from the function logQuantPseudoGenos in R package quantPop
    p0 <- startFreq[k]
    vQTL <- (h2_0*Vp_0)/nLoci   # genetic variance attributed to each QTL
    a <- sqrt((vQTL/(2*p0*(1-p0))))  # allelic effect
    geno1Mat <- matrix(sample(c(0,1),pseudoN*nLoci,replace=TRUE,prob=c(1-p0,p0)),nrow=nLoci,ncol=pseudoN)
    geno2Mat <- matrix(sample(c(0,1),pseudoN*nLoci,replace=TRUE,prob=c(1-p0,p0)),nrow=nLoci,ncol=pseudoN)
    modInt <- phen_0 - weighted.mean( c(0,a,2*a),c((1-p0)^2,2*p0*(1-p0),(p0)^2))*nLoci # control the mean phenotype in generation 1
    phens <-  modInt + colSums(geno1Mat + geno2Mat)*a + rnorm(n=ncol(geno1Mat),mean=0,sd=sqrt(Ve))   # individual phenotypes
    if(j == 1)hist(phens,breaks=50,xlab="Phenotype",main=paste("p0 = ",startFreq[k],"; ",locusNums[j]," locus",sep=""),border=NULL,col="darkred")
    if(j > 1)hist(phens,breaks=50,xlab="Phenotype",main=paste("p0 = ",startFreq[k],"; ",locusNums[j]," loci",sep=""),border=NULL,col="darkred")
    print(paste("mean pheno = ",mean(phens)))
    print(paste(" pheno variance = ",var(phens)))
  }
}


#########################################################
#--------------------------------------------------------
# plot population size, h2 and phenotype through time
# from individual-based simulations for populations with
# h20 = 0.4
# FIGURES S6
#--------------------------------------------------------
#########################################################

# get polygenic data
setwd("~/Documents/work/rapidAdaptation/AmNatRevision/simulations/originalSims/constantP0/h2Pt4/noMajor")
p0Vec <- c(0.1,0.25,0.5,0.75,0.9)

popFileNames <- paste("popSize_p0_",p0Vec,"_N0_500_K_1000",sep="")
phenFileNames <- paste("phen_p0_",p0Vec,"_N0_500_K_1000",sep="")
h2FileNames <- paste("h2_p0_",p0Vec,"_N0_500_K_1000",sep="")
polyPopSize <- list()
polyPhen <- list()
polyH2 <- list()

for(i in 1:5){
  polyPopSize [[i]] <- read.table(popFileNames[i],header=TRUE)
  polyPhen [[i]] <- read.table(phenFileNames[i],header=TRUE)
  polyH2 [[i]] <- read.table(h2FileNames[i],header=TRUE)
}


# get single major locus data

setwd("~/Documents/work/rapidAdaptation/AmNatRevision/simulations/originalSims/constantP0/h2Pt4/major")
popFileNames <- paste("popSize_N0_500_K_1000_propMaj_0.9999999_freqMaj_",p0Vec,sep="")
phenFileNames <- paste("phen_N0_500_K_1000_propMaj_0.9999999_freqMaj_",p0Vec,sep="")
h2FileNames <- paste("h2_N0_500_K_1000_propMaj_0.9999999_freqMaj_",p0Vec,sep="")
majorPopSize <- list()
majorPhen <- list()
majorH2 <- list()

for(i in 1:5){
  majorPopSize [[i]] <- read.table(popFileNames[i],header=TRUE)
  majorPhen [[i]] <- read.table(phenFileNames[i],header=TRUE)
  majorH2 [[i]] <- read.table(h2FileNames[i],header=TRUE)
}



# get  two major loci data

setwd("~/Documents/work/rapidAdaptation/AmNatRevision/simulations/originalSims/constantP0/h2Pt4/twoMajor")
popFileNames <- paste("popSize_N0_500_K_1000_propMaj_0.5_freqMaj_",p0Vec,sep="")
phenFileNames <- paste("phen_N0_500_K_1000_propMaj_0.5_freqMaj_",p0Vec,sep="")
h2FileNames <- paste("h2_N0_500_K_1000_propMaj_0.5_freqMaj_",p0Vec,sep="")
twoMajorPopSize <- list()
twoMajorPhen <- list()
twoMajorH2 <- list()

for(i in 1:5){
  twoMajorPopSize [[i]] <- read.table(popFileNames[i],header=TRUE)
  twoMajorPhen [[i]] <- read.table(phenFileNames[i],header=TRUE)
  twoMajorH2 [[i]] <- read.table(h2FileNames[i],header=TRUE)
}

#-----------------------------------------------
# make NA population sizes zero
#-----------------------------------------------

for(i in 1:5){
  for(j in 1:ncol(polyPopSize[[i]])){
    polyPopSize[[i]][which(is.na(polyPopSize[[i]][,j])),j] <- 0
    majorPopSize[[i]][which(is.na(majorPopSize[[i]][,j])),j] <- 0
    twoMajorPopSize[[i]][which(is.na(twoMajorPopSize[[i]][,j])),j] <- 0
  }
}


#------------------------------------------------------------------
# make missing heritabilities for populations N < 10
#------------------------------------------------------------------
for (i in 1:5){
  #### polygenic
  thisPop <- polyPopSize[[i]][,1:80]
  thisH2 <- polyH2[[i]]
  for(j in 1:nrow(thisPop)){
    if(sum(thisPop[j,] < 10,na.rm=TRUE) > 0){
      thisH2[j,which(thisPop[j,] < 10)] <- NA
    }
  }
  polyH2[[i]] <- thisH2

  #### single major locus
  thisPop <- majorPopSize[[i]][,1:80]
  thisH2 <- majorH2[[i]]
  for(j in 1:nrow(thisPop)){
    if(sum(thisPop[j,] < 10,na.rm=TRUE) > 0){
      thisH2[j,which(thisPop[j,] < 10)] <- NA
    }
  }
  majorH2[[i]] <- thisH2

  #### two major loci
  thisPop <- twoMajorPopSize[[i]][,1:80]
  thisH2 <- twoMajorH2[[i]]
  for(j in 1:nrow(thisPop)){
    if(sum(thisPop[j,] < 10,na.rm=TRUE) > 0){
      thisH2[j,which(thisPop[j,] < 10)] <- NA
    }
  }
  twoMajorH2[[i]] <- thisH2

}





#----------------------------------------------------------------------------
# calculate mean parmameters and extinction rate
# for major, two major an polygenic architectures
#----------------------------------------------------------------------------

twoMajorExtProp <- list()
majorExtProp <- list()
polyExtProp <- list()

twoMajorMeanPop <- list()
twoMajorMeanPhen <- list()
twoMajorMeanH2 <- list()

majorMeanPop <- list()
majorMeanPhen <- list()
majorMeanH2 <- list()

polyMeanPop <- list()
polyMeanPhen <- list()
polyMeanH2 <- list()

for(i in 1:5){
  polyMeanPop [[i]] <- colMeans(polyPopSize[[i]])
  polyMeanPhen [[i]] <- colMeans(polyPhen[[i]],na.rm=TRUE)
  polyMeanH2  [[i]] <- colMeans(polyH2[[i]],na.rm=TRUE)

  majorMeanPop [[i]] <- colMeans(majorPopSize[[i]])
  majorMeanPhen [[i]] <- colMeans(majorPhen[[i]],na.rm=TRUE)
  majorMeanH2  [[i]] <- colMeans(majorH2[[i]],na.rm=TRUE)

  twoMajorMeanPop [[i]] <- colMeans(twoMajorPopSize[[i]])
  twoMajorMeanPhen [[i]] <- colMeans(twoMajorPhen[[i]],na.rm=TRUE)
  twoMajorMeanH2  [[i]] <- colMeans(twoMajorH2[[i]],na.rm=TRUE)

  polyExtProp [[i]] <- colSums(polyPopSize[[i]] == 0)/500
  majorExtProp [[i]] <- colSums(majorPopSize[[i]] == 0)/500
  twoMajorExtProp [[i]] <- colSums(twoMajorPopSize[[i]] == 0)/500
}

########################################
# Make the figure
########################################

# population size
t <- 81
carryCap <- 1000
# one major locus
par(mfrow=c(3,3),xpd=TRUE,mar=c(4,5,1.5,1))
plot(1:t,majorMeanPop[[1]],ylim=c(0,carryCap),type="n",cex.axis=0.7, main ="1 major locus",cex.main=1.4,
     cex.lab=1.5,ylab=expression(italic(""*N*"")),xlab="")
text(x=-30,y=1100,labels="A",cex=2)
lineVec <- c("solid","dashed","dotted","dotdash","twodash")
for(i in 1:length(majorMeanPop)){
  lines(1:t,majorMeanPop[[i]],lty=lineVec[i],lwd=1.3)
}

# two major loci
plot(1:t,twoMajorMeanPop[[1]],ylim=c(0,carryCap),type="n",cex.axis=0.7, main ="2 major loci",cex.main=1.4,
     cex.lab=1.5,xlab="",ylab="")
text(x=-30,y=1100,labels="B",cex=2)
lineVec <- c("solid","dashed","dotted","dotdash","twodash")
for(i in 1:length(majorMeanPop)){
  lines(1:t,twoMajorMeanPop[[i]],lty=lineVec[i],lwd=1.3)
}

# polygenic
plot(1:t,polyMeanPop[[1]],ylim=c(0,carryCap),type="n",cex.axis=0.7, main ="100 small-effect loci",cex.main=1.4,
     cex.lab=1.5,xlab="",ylab="")
text(x=-30,y=1100,labels="C",cex=2)
lineVec <- c("solid","dashed","dotted","dotdash","twodash")
for(i in 1:length(majorMeanPop)){
  lines(1:t,polyMeanPop[[i]],lty=lineVec[i],lwd=1.3)
}

############# heritability

# one major locus
plot(1:t,majorMeanPop[[1]],ylim=c(0,0.8),type="n",cex.axis=0.7,
     cex.lab=1.5,ylab=expression(italic(""*h*"")^2),xlab="")
text(x=-35,y=1100,labels="A",cex=2)
lineVec <- c("solid","dashed","dotted","dotdash","twodash")
for(i in 1:length(majorMeanPop)){
  lines(1:(t-1),majorMeanH2[[i]],lty=lineVec[i],lwd=1.3)
}

# two major loci
plot(1:t,twoMajorMeanPop[[1]],ylim=c(0,0.8),type="n",cex.axis=0.7,
     cex.lab=1.5,xlab="",ylab="")
text(x=-35,y=1100,labels="A",cex=2)
lineVec <- c("solid","dashed","dotted","dotdash","twodash")
for(i in 1:length(majorMeanPop)){
  lines(1:(t-1),twoMajorMeanH2[[i]],lty=lineVec[i],lwd=1.3)
}

# polygenic
plot(1:t,polyMeanPop[[1]],ylim=c(0,0.8),type="n",cex.axis=0.7,
     cex.lab=1.5,xlab="",ylab="")
text(x=-35,y=1100,labels="A",cex=2)
lineVec <- c("solid","dashed","dotted","dotdash","twodash")
for(i in 1:length(majorMeanPop)){
  lines(1:(t-1),polyMeanH2[[i]],lty=lineVec[i],lwd=1.3)
}


############# phenotype

# one major locus

plot(1:t,majorMeanPop[[1]],ylim=c(100,111),type="n",cex.axis=0.7,
     cex.lab=1.5,ylab="Mean phenotype",xlab="Generation")
text(x=-35,y=1100,labels="A",cex=2)
lineVec <- c("solid","dashed","dotted","dotdash","twodash")
for(i in 1:length(majorMeanPop)){
  lines(1:(t-1),majorMeanPhen[[i]],lty=lineVec[i],lwd=1.3)
}

# two major loci
plot(1:t,majorMeanPop[[1]],ylim=c(100,111),type="n",cex.axis=0.7,
     cex.lab=1.5,ylab="",xlab="Generation")
text(x=-35,y=1100,labels="A",cex=2)
lineVec <- c("solid","dashed","dotted","dotdash","twodash")
for(i in 1:length(majorMeanPop)){
  lines(1:(t-1),twoMajorMeanPhen[[i]],lty=lineVec[i],lwd=1.3)
}

# polygenic
plot(1:t,majorMeanPop[[1]],ylim=c(100,111),type="n",cex.axis=0.7,
     cex.lab=1.5,ylab="",xlab="Generation")
text(x=-35,y=1100,labels="A",cex=2)
lineVec <- c("solid","dashed","dotted","dotdash","twodash")
for(i in 1:length(majorMeanPop)){
  lines(1:(t-1),polyMeanPhen[[i]],lty=lineVec[i],lwd=1.3)
}

legend(x=30,y=100,lty=lineVec,legend=c(expression(paste(italic(""*p*"")[0]," = 0.1",sep="")),
                                       expression(paste(italic(""*p*"")[0]," = 0.25",sep="")),
                                       expression(paste(italic(""*p*"")[0]," = 0.5",sep="")),
                                       expression(paste(italic(""*p*"")[0]," = 0.75",sep="")),
                                       expression(paste(italic(""*p*"")[0]," = 0.9",sep=""))),
       bty="n",xjust=FALSE,yjust=FALSE)



#---------------------------------------------------
# plot extinction rates for the three architectures
#---------------------------------------------------

par(mar=c(4,5,2,1),mfrow=c(1,3))
# single major locus
plot(c(0,80),y=c(0,1),type="n",xlab="Generation",ylab="Proportion extinct",
     cex.lab=1.5,main ="1 major locus",cex.axis=0.7)
for(i in 1:5){
  lines(1:81,majorExtProp[[i]],lty=lineVec[i],lwd=1.3)
}


# two major loci
plot(c(0,80),y=c(0,1),type="n",xlab="Generation",
     ylab="",cex.lab=1.5,main ="2 major loci",cex.axis=0.7)
for(i in 1:5){
  lines(1:81,twoMajorExtProp[[i]],lty=lineVec[i],lwd=1.3)
}

# 100 loci
plot(c(0,80),y=c(0,1),type="n",xlab="Generation",
     ylab="",cex.lab=1.5,main ="100 small-effect loci",cex.axis=0.7)
for(i in 1:5){
  lines(1:81,polyExtProp[[i]],lty=lineVec[i],lwd=1.3)
}
legend(x=30,y=0,lty=lineVec,legend=c(expression(paste(italic(""*p*"")[0]," = 0.1",sep="")),
                                       expression(paste(italic(""*p*"")[0]," = 0.25",sep="")),
                                       expression(paste(italic(""*p*"")[0]," = 0.5",sep="")),
                                       expression(paste(italic(""*p*"")[0]," = 0.75",sep="")),
                                       expression(paste(italic(""*p*"")[0]," = 0.9",sep=""))),
       bty="n",xjust=FALSE,yjust=FALSE)

#########################################################
# plot population size, h2 and phenotype through time
# from individual-based simulations for populations with
# h20 = 0.8
# FIGURE S8
#########################################################
# get polygenic data
setwd("~/Documents/work/rapidAdaptation/AmNatRevision/simulations/originalSims/constantP0/h2Pt8/noMajor")
p0Vec <- c(0.1,0.25,0.5,0.75,0.9)

popFileNames <- paste("popSize_p0_",p0Vec,"_N0_500_K_1000",sep="")
phenFileNames <- paste("phen_p0_",p0Vec,"_N0_500_K_1000",sep="")
h2FileNames <- paste("h2_p0_",p0Vec,"_N0_500_K_1000",sep="")
polyPopSize <- list()
polyPhen <- list()
polyH2 <- list()

for(i in 1:5){
  polyPopSize [[i]] <- read.table(popFileNames[i],header=TRUE)
  polyPhen [[i]] <- read.table(phenFileNames[i],header=TRUE)
  polyH2 [[i]] <- read.table(h2FileNames[i],header=TRUE)
}


# get single major locus data

setwd("~/Documents/work/rapidAdaptation/AmNatRevision/simulations/originalSims/constantP0/h2Pt8/major")
popFileNames <- paste("popSize_N0_500_K_1000_propMaj_0.9999999_freqMaj_",p0Vec,sep="")
phenFileNames <- paste("phen_N0_500_K_1000_propMaj_0.9999999_freqMaj_",p0Vec,sep="")
h2FileNames <- paste("h2_N0_500_K_1000_propMaj_0.9999999_freqMaj_",p0Vec,sep="")
majorPopSize <- list()
majorPhen <- list()
majorH2 <- list()

for(i in 1:5){
  majorPopSize [[i]] <- read.table(popFileNames[i],header=TRUE)
  majorPhen [[i]] <- read.table(phenFileNames[i],header=TRUE)
  majorH2 [[i]] <- read.table(h2FileNames[i],header=TRUE)
}



# get  two major loci data

setwd("~/Documents/work/rapidAdaptation/AmNatRevision/simulations/originalSims/constantP0/h2Pt8/twoMajor")
popFileNames <- paste("popSize_N0_500_K_1000_propMaj_0.5_freqMaj_",p0Vec,sep="")
phenFileNames <- paste("phen_N0_500_K_1000_propMaj_0.5_freqMaj_",p0Vec,sep="")
h2FileNames <- paste("h2_N0_500_K_1000_propMaj_0.5_freqMaj_",p0Vec,sep="")
twoMajorPopSize <- list()
twoMajorPhen <- list()
twoMajorH2 <- list()

for(i in 1:5){
  twoMajorPopSize [[i]] <- read.table(popFileNames[i],header=TRUE)
  twoMajorPhen [[i]] <- read.table(phenFileNames[i],header=TRUE)
  twoMajorH2 [[i]] <- read.table(h2FileNames[i],header=TRUE)
}

#-----------------------------------------------
# make NA population sizes zero
#-----------------------------------------------

for(i in 1:5){
  for(j in 1:ncol(polyPopSize[[i]])){
    polyPopSize[[i]][which(is.na(polyPopSize[[i]][,j])),j] <- 0
    majorPopSize[[i]][which(is.na(majorPopSize[[i]][,j])),j] <- 0
    twoMajorPopSize[[i]][which(is.na(twoMajorPopSize[[i]][,j])),j] <- 0
  }
}


#------------------------------------------------------------------
# make missing heritabilities for populations N < 10
#------------------------------------------------------------------
for (i in 1:5){
  #### polygenic
  thisPop <- polyPopSize[[i]][,1:80]
  thisH2 <- polyH2[[i]]
  for(j in 1:nrow(thisPop)){
    if(sum(thisPop[j,] < 10,na.rm=TRUE) > 0){
      thisH2[j,which(thisPop[j,] < 10)] <- NA
    }
  }
  polyH2[[i]] <- thisH2

  #### single major locus
  thisPop <- majorPopSize[[i]][,1:80]
  thisH2 <- majorH2[[i]]
  for(j in 1:nrow(thisPop)){
    if(sum(thisPop[j,] < 10,na.rm=TRUE) > 0){
      thisH2[j,which(thisPop[j,] < 10)] <- NA
    }
  }
  majorH2[[i]] <- thisH2

  #### two major loci
  thisPop <- twoMajorPopSize[[i]][,1:80]
  thisH2 <- twoMajorH2[[i]]
  for(j in 1:nrow(thisPop)){
    if(sum(thisPop[j,] < 10,na.rm=TRUE) > 0){
      thisH2[j,which(thisPop[j,] < 10)] <- NA
    }
  }
  twoMajorH2[[i]] <- thisH2

}

#----------------------------------------------------------------------------
# calculate mean parmameters and extinction rate
# for major, two major an polygenic architectures
#----------------------------------------------------------------------------

twoMajorExtProp <- list()
majorExtProp <- list()
polyExtProp <- list()

twoMajorMeanPop <- list()
twoMajorMeanPhen <- list()
twoMajorMeanH2 <- list()

majorMeanPop <- list()
majorMeanPhen <- list()
majorMeanH2 <- list()

polyMeanPop <- list()
polyMeanPhen <- list()
polyMeanH2 <- list()

for(i in 1:5){
  polyMeanPop [[i]] <- colMeans(polyPopSize[[i]])
  polyMeanPhen [[i]] <- colMeans(polyPhen[[i]],na.rm=TRUE)
  polyMeanH2  [[i]] <- colMeans(polyH2[[i]],na.rm=TRUE)

  majorMeanPop [[i]] <- colMeans(majorPopSize[[i]])
  majorMeanPhen [[i]] <- colMeans(majorPhen[[i]],na.rm=TRUE)
  majorMeanH2  [[i]] <- colMeans(majorH2[[i]],na.rm=TRUE)

  twoMajorMeanPop [[i]] <- colMeans(twoMajorPopSize[[i]])
  twoMajorMeanPhen [[i]] <- colMeans(twoMajorPhen[[i]],na.rm=TRUE)
  twoMajorMeanH2  [[i]] <- colMeans(twoMajorH2[[i]],na.rm=TRUE)

  polyExtProp [[i]] <- colSums(polyPopSize[[i]] == 0)/500
  majorExtProp [[i]] <- colSums(majorPopSize[[i]] == 0)/500
  twoMajorExtProp [[i]] <- colSums(twoMajorPopSize[[i]] == 0)/500
}

########################################
# Make the figure
########################################

# population size
t <- 81
carryCap <- 1000
# one major locus
par(mfrow=c(3,3),xpd=TRUE,mar=c(4,5,1.5,1))
plot(1:t,majorMeanPop[[1]],ylim=c(0,carryCap),type="n",cex.axis=0.7, main ="1 major locus",cex.main=1.4,
     cex.lab=1.5,ylab=expression(italic(""*N*"")),xlab="")

lineVec <- c("solid","dashed","dotted","dotdash","twodash")
for(i in 1:length(majorMeanPop)){
  lines(1:t,majorMeanPop[[i]],lty=lineVec[i],lwd=1.3)
}
text(x=-35,y=1100,labels="A",cex=2)
# two major loci
plot(1:t,twoMajorMeanPop[[1]],ylim=c(0,carryCap),type="n",cex.axis=0.7, main ="2 major loci",cex.main=1.4,
     cex.lab=1.5,xlab="",ylab="")
text(x=-30,y=1100,labels="B",cex=2)
lineVec <- c("solid","dashed","dotted","dotdash","twodash")
for(i in 1:length(majorMeanPop)){
  lines(1:t,twoMajorMeanPop[[i]],lty=lineVec[i],lwd=1.3)
}

# polygenic
plot(1:t,polyMeanPop[[1]],ylim=c(0,carryCap),type="n",cex.axis=0.7, main ="100 small-effect loci",cex.main=1.4,
     cex.lab=1.5,xlab="",ylab="")
text(x=-30,y=1100,labels="C",cex=2)
lineVec <- c("solid","dashed","dotted","dotdash","twodash")
for(i in 1:length(majorMeanPop)){
  lines(1:t,polyMeanPop[[i]],lty=lineVec[i],lwd=1.3)
}

par(xpd=FALSE)
############# heritability

# one major locus
plot(1:t,majorMeanPop[[1]],ylim=c(0,1),type="n",cex.axis=0.7,
     cex.lab=1.5,ylab=expression(italic(""*h*"")^2),xlab="")
text(x=-35,y=1100,labels="A",cex=2)
lineVec <- c("solid","dashed","dotted","dotdash","twodash")
for(i in 1:length(majorMeanPop)){
  lines(1:(t-1),majorMeanH2[[i]],lty=lineVec[i],lwd=1.3)
}

# two major loci
plot(1:t,twoMajorMeanPop[[1]],ylim=c(0,1),type="n",cex.axis=0.7,
     cex.lab=1.5,xlab="",ylab="")
text(x=-35,y=1100,labels="A",cex=2)
lineVec <- c("solid","dashed","dotted","dotdash","twodash")
for(i in 1:length(majorMeanPop)){
  lines(1:(t-1),twoMajorMeanH2[[i]],lty=lineVec[i],lwd=1.3)
}

# polygenic
plot(1:t,polyMeanPop[[1]],ylim=c(0,1),type="n",cex.axis=0.7,
     cex.lab=1.5,xlab="",ylab="")
text(x=-35,y=1100,labels="A",cex=2)
lineVec <- c("solid","dashed","dotted","dotdash","twodash")
for(i in 1:length(majorMeanPop)){
  lines(1:(t-1),polyMeanH2[[i]],lty=lineVec[i],lwd=1.3)
}


############# phenotype

# one major locus

plot(1:t,majorMeanPop[[1]],ylim=c(100,112),type="n",cex.axis=0.7,
     cex.lab=1.5,ylab="Mean phenotype",xlab="Generation")
text(x=-35,y=1100,labels="A",cex=2)
lineVec <- c("solid","dashed","dotted","dotdash","twodash")
for(i in 1:length(majorMeanPop)){
  lines(1:(t-1),majorMeanPhen[[i]],lty=lineVec[i],lwd=1.3)
}

# two major loci
plot(1:t,majorMeanPop[[1]],ylim=c(100,112),type="n",cex.axis=0.7,
     cex.lab=1.5,ylab="",xlab="Generation")
text(x=-35,y=1100,labels="A",cex=2)
lineVec <- c("solid","dashed","dotted","dotdash","twodash")
for(i in 1:length(majorMeanPop)){
  lines(1:(t-1),twoMajorMeanPhen[[i]],lty=lineVec[i],lwd=1.3)
}

# polygenic
plot(1:t,majorMeanPop[[1]],ylim=c(100,112),type="n",cex.axis=0.7,
     cex.lab=1.5,ylab="",xlab="Generation")
text(x=-35,y=1100,labels="A",cex=2)
lineVec <- c("solid","dashed","dotted","dotdash","twodash")
for(i in 1:length(majorMeanPop)){
  lines(1:(t-1),polyMeanPhen[[i]],lty=lineVec[i],lwd=1.3)
}

legend(x=30,y=100,lty=lineVec,legend=c(expression(paste(italic(""*p*"")[0]," = 0.1",sep="")),
                                       expression(paste(italic(""*p*"")[0]," = 0.25",sep="")),
                                       expression(paste(italic(""*p*"")[0]," = 0.5",sep="")),
                                       expression(paste(italic(""*p*"")[0]," = 0.75",sep="")),
                                       expression(paste(italic(""*p*"")[0]," = 0.9",sep=""))),
       bty="n",xjust=FALSE,yjust=FALSE)


####################################################
#---------------------------------------------------
# plot extinction rates for the three architectures
# Figure S9
#---------------------------------------------------
####################################################
par(mar=c(4,5,2,1))
# single major locus
par(mfrow=c(1,3))
plot(c(0,80),y=c(0,1),type="n",xlab="Generation",ylab="Proportion extinct",
     cex.lab=1.5,main ="1 major locus",cex.axis=0.7)
for(i in 1:5){
  lines(1:81,majorExtProp[[i]],lty=lineVec[i],lwd=1.3)
}


# two major loci
plot(c(0,80),y=c(0,1),type="n",xlab="Generation",
     ylab="",cex.lab=1.5,main ="2 major loci",cex.axis=0.7)
for(i in 1:5){
  lines(1:81,twoMajorExtProp[[i]],lty=lineVec[i],lwd=1.3)
}

# 100 loci
plot(c(0,80),y=c(0,1),type="n",xlab="Generation",
     ylab="",cex.lab=1.5,main ="100 small-effect loci",cex.axis=0.7)
for(i in 1:5){
  lines(1:81,polyExtProp[[i]],lty=lineVec[i],lwd=1.3)
}
legend(x=30,y=0.4,lty=lineVec,legend=c(expression(paste(italic(""*p*"")[0]," = 0.1",sep="")),
                                       expression(paste(italic(""*p*"")[0]," = 0.25",sep="")),
                                       expression(paste(italic(""*p*"")[0]," = 0.5",sep="")),
                                       expression(paste(italic(""*p*"")[0]," = 0.75",sep="")),
                                       expression(paste(italic(""*p*"")[0]," = 0.9",sep=""))),
       bty="n",xjust=FALSE,yjust=FALSE)


########################################################
########################################################
# analyze results where the loci have variable
# p0 for coral and large mammals and the large effect
# locus is responsible for 90% of the genetic variance
# for popultions with LOW IMMIGRATION from a population with
# a different phenotypic optimum
# FIGURE S10
########################################################
########################################################

propAnalyze <- 0.9

setwd("~/Documents/work/rapidAdaptation/AmNatRevision/simulations/originalSims/h2Pt6/fecundity4/major_lowMig")
popFileName <- paste("popSize_N0_500_K_1000_propMaj_",propAnalyze,sep="")
phenFileName <-  paste("phen_N0_500_K_1000_propMaj_",propAnalyze,sep="")
h2FileName <-  paste("h2_N0_500_K_1000_propMaj_",propAnalyze,sep="")

mam_majorPopSize <- read.table(popFileName,header=TRUE)
mam_majorPhen <- read.table(phenFileName,header=TRUE)
mam_majorH2 <- read.table(h2FileName,header=TRUE)


setwd("~/Documents/work/rapidAdaptation/AmNatRevision/simulations/originalSims/h2Pt6/fecundity4/noMajor_lowMig")

mam_noMajorPopSize <- read.table("popSize_N0_500_K_1000_propMaj_0.01",header=TRUE)
mam_noMajorPhen <- read.table("phen_N0_500_K_1000_propMaj_0.01",header=TRUE)
mam_noMajorH2 <- read.table("h2_N0_500_K_1000_propMaj_0.01",header=TRUE)

#-----------------------------------------------
# make NA population sizes zero
#-----------------------------------------------


for(j in 1:ncol(mam_noMajorPopSize)){
  mam_noMajorPopSize[which(is.na(mam_noMajorPopSize[,j])),j] <- 0
  mam_majorPopSize[which(is.na(mam_majorPopSize[,j])),j] <- 0
}


#######################
# check for extinction
#######################

mam_polyExt <- rep(NA,81)
mam_OligExt <- rep(NA,81)

for(i in 1:81){
  mam_polyExt[i] <- sum(mam_noMajorPopSize[,i] == 0)
  mam_OligExt[i] <- sum(mam_majorPopSize[,i] == 0)
}


#############################################
# get bootstrap CIs for proportion extinct
#############################################

bootReps <- 1000
lowPolyBoots <- matrix(NA,nrow=bootReps,ncol=81)
lowOligBoots <- matrix(NA,nrow=bootReps,ncol=81)

for(i in 1:bootReps){
  lowPolyDat <- mam_noMajorPopSize[sample(1:nrow(mam_noMajorPopSize),nrow(mam_noMajorPopSize),replace=TRUE),]
  lowOligDat <- mam_majorPopSize[sample(1:nrow(mam_majorPopSize),nrow(mam_majorPopSize),replace=TRUE),]
  lowPolyBoots [i,] <- colSums(lowPolyDat == 0)/nrow(lowPolyDat )
  lowOligBoots [i,] <- colSums(lowOligDat == 0)/nrow(lowOligDat )
  print(i)
}

mamlowPolyCI <- NULL
mamlowOligCI <- NULL

for(i in 1:81){
  mamlowPolyCI <- cbind(mamlowPolyCI,quantile(lowPolyBoots[,i],c(0.025,0.975)))
  mamlowOligCI <- cbind(mamlowOligCI,quantile(lowOligBoots[,i],c(0.025,0.975)))
}


############################################
# plot results
############################################
par(mfrow=c(3,1),mar=c(2,5,3,0.5),xpd=FALSE)
library(scales)
carryCap <- 1000
##### mammal population size
plot(c(0,81),c(0,carryCap + 200)/1000,type="n",ylab=expression(paste(italic(""*N*"")," (thousands)",sep="")),xlab="",cex.lab=1.5,
     cex.main=1.5,cex.axis=0.9)
for(i in 1:nrow(mam_noMajorPopSize)){
  lines(1:81,mam_noMajorPopSize[i,]/1000,col=alpha("orange",alpha=0.05))
  lines(1:81,mam_majorPopSize[i,]/1000,col=alpha("darkblue",alpha=0.05))
}
polySizeMean <- rep(NA,81)
oligSizeMean <- rep(NA,81)
for(i in 1:81){
  polySizeMean[i] <- mean(mam_noMajorPopSize[,i])/1000
  oligSizeMean[i] <- mean(mam_majorPopSize[,i])/1000
}
lines(1:81,polySizeMean,lwd=3,col="orange")
lines(1:81,oligSizeMean,lwd=3,col="darkblue")
lines(c(0,200),c(carryCap/1000,carryCap/1000),lty="dashed")

#---------------------------------
##### plot the phenotypes
#---------------------------------
par(mar=c(3,5,2,0.5))
idealPhen <- 110

####### mammals
plot(c(0,80),c(100,112),type="n",ylab="Mean Phenotype",xlab="",cex.lab=1.5,cex.axis=0.8)
for(i in 1:500){
  lines(1:80,mam_noMajorPhen[i,],col=alpha("orange",alpha=0.05))
  lines(1:80,mam_majorPhen[i,],col=alpha("darkblue",alpha=0.05))
}

# get rolling mean
polyPhenMean <- rep(NA,80)
oligPhenMean <- rep(NA,80)
for(i in 1:80){
  polyPhenMean[i] <- mean(mam_noMajorPhen[,i],na.rm=TRUE)
  oligPhenMean[i] <- mean(mam_majorPhen[,i],na.rm=TRUE)
}
lines(1:80,polyPhenMean,lwd=3,col="orange")
lines(1:80,oligPhenMean,lwd=3,col="darkblue")
lines(c(0,80),c(110,110),lty="dashed")


par(mar=c(3,5,2,0.5))
idealPhen <- 110

#####################################################################
# plot the fraction of extinct population versus time
#####################################################################

par(mar=c(4,5,1,0.5))
# extinct proportion for mammals
plot(c(1,81),c(0,1),type="n",xlab="Generation",ylab="Proportion Extinct",cex.lab=1.5)
lines(1:81,mam_polyExt/nrow(lowPolyDat),col="orange",lwd=3)
lines(1:81,mam_OligExt/nrow(lowOligDat),col="darkblue",lwd=3)

polygon(x=c(1:81,rev(1:81)),y=c(mamlowPolyCI[1,1:81],rev(mamlowPolyCI[2,1:81])),col=adjustcolor("orange",alpha.f=0.5),border=NA)
polygon(x=c(1:81,rev(1:81)),y=c(mamlowOligCI[1,1:81],rev(mamlowOligCI[2,1:81])),col=adjustcolor("darkblue",alpha.f=0.5),border=NA)
# get pictures

setwd("~/Documents/work/rapidAdaptation/AmNatRevision/simulations/originalSims/h2Pt6")
library(jpeg)
img<-readJPEG("wildebeest.jpg")
rasterImage(img,xleft=-3, ybottom=0.8, xright=25, ytop=1.1)


########################################################
########################################################
# analyze results where the loci have variable
# p0 for coral and large mammals and the large effect
# locus is responsible for 90% of the genetic variance
# for popultions with HIGH IMMIGRATION from a population with
# a different phenotypic optimum
# FIGURE S11
########################################################
########################################################

propAnalyze <- 0.9
####----------------------
#### large mammals first
####----------------------
setwd("~/Documents/work/rapidAdaptation/AmNatRevision/simulations/originalSims/h2Pt6/fecundity4/major_highMig")
popFileName <- paste("popSize_N0_500_K_1000_propMaj_",propAnalyze,sep="")
phenFileName <-  paste("phen_N0_500_K_1000_propMaj_",propAnalyze,sep="")
h2FileName <-  paste("h2_N0_500_K_1000_propMaj_",propAnalyze,sep="")

mam_majorPopSize <- read.table(popFileName,header=TRUE)
mam_majorPhen <- read.table(phenFileName,header=TRUE)
mam_majorH2 <- read.table(h2FileName,header=TRUE)


setwd("~/Documents/work/rapidAdaptation/AmNatRevision/simulations/originalSims/h2Pt6/fecundity4/noMajor_highMig")

mam_noMajorPopSize <- read.table("popSize_N0_500_K_1000_propMaj_0.01",header=TRUE)
mam_noMajorPhen <- read.table("phen_N0_500_K_1000_propMaj_0.01",header=TRUE)
mam_noMajorH2 <- read.table("h2_N0_500_K_1000_propMaj_0.01",header=TRUE)

#-----------------------------------------------
# make NA population sizes zero
#-----------------------------------------------


for(j in 1:ncol(mam_noMajorPopSize)){
  mam_noMajorPopSize[which(is.na(mam_noMajorPopSize[,j])),j] <- 0
  mam_majorPopSize[which(is.na(mam_majorPopSize[,j])),j] <- 0
}


#######################
# check for extinction
#######################

mam_polyExt <- rep(NA,81)
mam_OligExt <- rep(NA,81)

for(i in 1:81){
  mam_polyExt[i] <- sum(mam_noMajorPopSize[,i] == 0)
  mam_OligExt[i] <- sum(mam_majorPopSize[,i] == 0)
}


#############################################
# get bootstrap CIs for proportion extinct
#############################################

bootReps <- 1000
lowPolyBoots <- matrix(NA,nrow=bootReps,ncol=81)
lowOligBoots <- matrix(NA,nrow=bootReps,ncol=81)

for(i in 1:bootReps){
  lowPolyDat <- mam_noMajorPopSize[sample(1:nrow(mam_noMajorPopSize),nrow(mam_noMajorPopSize),replace=TRUE),]
  lowOligDat <- mam_majorPopSize[sample(1:nrow(mam_majorPopSize),nrow(mam_majorPopSize),replace=TRUE),]
  lowPolyBoots [i,] <- colSums(lowPolyDat == 0)/nrow(lowPolyDat )
  lowOligBoots [i,] <- colSums(lowOligDat == 0)/nrow(lowOligDat )
  print(i)
}

mamlowPolyCI <- NULL
mamlowOligCI <- NULL

for(i in 1:81){
  mamlowPolyCI <- cbind(mamlowPolyCI,quantile(lowPolyBoots[,i],c(0.025,0.975)))
  mamlowOligCI <- cbind(mamlowOligCI,quantile(lowOligBoots[,i],c(0.025,0.975)))
}


############################################
# plot results
############################################
par(mfrow=c(3,1),mar=c(2,5,3,0.5),xpd=FALSE)
library(scales)
carryCap <- 1000
##### mammal population size
plot(c(0,81),c(0,carryCap + 200)/1000,type="n",ylab=expression(paste(italic(""*N*"")," (thousands)",sep="")),xlab="",cex.lab=1.5,
     cex.main=1.5,cex.axis=0.9)
for(i in 1:nrow(mam_noMajorPopSize)){
  lines(1:81,mam_noMajorPopSize[i,]/1000,col=alpha("orange",alpha=0.05))
  lines(1:81,mam_majorPopSize[i,]/1000,col=alpha("darkblue",alpha=0.05))
}
polySizeMean <- rep(NA,81)
oligSizeMean <- rep(NA,81)
for(i in 1:81){
  polySizeMean[i] <- mean(mam_noMajorPopSize[,i])/1000
  oligSizeMean[i] <- mean(mam_majorPopSize[,i])/1000
}
lines(1:81,polySizeMean,lwd=3,col="orange")
lines(1:81,oligSizeMean,lwd=3,col="darkblue")
lines(c(0,200),c(carryCap/1000,carryCap/1000),lty="dashed")

#---------------------------------
##### plot the phenotypes
#---------------------------------
par(mar=c(3,5,2,0.5))
idealPhen <- 110

####### mammals
plot(c(0,80),c(100,112),type="n",ylab="Mean Phenotype",xlab="",cex.lab=1.5,cex.axis=0.8)
for(i in 1:500){
  lines(1:80,mam_noMajorPhen[i,],col=alpha("orange",alpha=0.05))
  lines(1:80,mam_majorPhen[i,],col=alpha("darkblue",alpha=0.05))
}

# get rolling mean
polyPhenMean <- rep(NA,80)
oligPhenMean <- rep(NA,80)
for(i in 1:80){
  polyPhenMean[i] <- mean(mam_noMajorPhen[,i],na.rm=TRUE)
  oligPhenMean[i] <- mean(mam_majorPhen[,i],na.rm=TRUE)
}
lines(1:80,polyPhenMean,lwd=3,col="orange")
lines(1:80,oligPhenMean,lwd=3,col="darkblue")
lines(c(0,80),c(110,110),lty="dashed")


par(mar=c(3,5,2,0.5))
idealPhen <- 110

#####################################################################
# plot the fraction of extinct population versus time
#####################################################################

par(mar=c(4,5,1,0.5))
# extinct proportion for mammals
plot(c(1,81),c(0,1),type="n",xlab="Generation",ylab="Proportion Extinct",cex.lab=1.5)
lines(1:81,mam_polyExt/nrow(lowPolyDat),col="orange",lwd=3)
lines(1:81,mam_OligExt/nrow(lowOligDat),col="darkblue",lwd=3)

polygon(x=c(1:81,rev(1:81)),y=c(mamlowPolyCI[1,1:81],rev(mamlowPolyCI[2,1:81])),col=adjustcolor("orange",alpha.f=0.5),border=NA)
polygon(x=c(1:81,rev(1:81)),y=c(mamlowOligCI[1,1:81],rev(mamlowOligCI[2,1:81])),col=adjustcolor("darkblue",alpha.f=0.5),border=NA)
# get pictures

setwd("~/Documents/work/rapidAdaptation/AmNatRevision/simulations/originalSims/h2Pt6")
library(jpeg)
img<-readJPEG("wildebeest.jpg")
rasterImage(img,xleft=-3, ybottom=0.8, xright=25, ytop=1.1)





##################################################################
##################################################################
##################################################################
# dynamics in populations with different major locus effect sizes
# FIGURE S12
##################################################################
##################################################################
##################################################################

####----------------------
#### large mammals first
####----------------------
setwd("~/Documents/work/rapidAdaptation/AmNatRevision/simulations/originalSims/h2Pt6/fecundity4/major")
mam_majorPopSize_effSize5 <- read.table("popSize_N0_500_K_1000_propMaj_0.5",header=TRUE)
mam_majorPopSize_effSize7 <- read.table("popSize_N0_500_K_1000_propMaj_0.7",header=TRUE)
mam_majorPopSize_effSize9 <- read.table("popSize_N0_500_K_1000_propMaj_0.9",header=TRUE)


mam_phen_effSize5 <-  read.table("phen_N0_500_K_1000_propMaj_0.5",header=TRUE)
mam_phen_effSize7 <-  read.table("phen_N0_500_K_1000_propMaj_0.7",header=TRUE)
mam_phen_effSize9 <-  read.table("phen_N0_500_K_1000_propMaj_0.9",header=TRUE)

setwd("~/Documents/work/rapidAdaptation/AmNatRevision/simulations/originalSims/h2Pt6/fecundity4/noMajor")
mam_noMajorPopSize <- read.table("popSize_N0_500_K_1000_propMaj_0.01",header=TRUE)
mam_noMajorPhen <-  read.table("phen_N0_500_K_1000_propMaj_0.01",header=TRUE)



#-----------------------------------------------
# make NA population sizes zero
#-----------------------------------------------


for(j in 1:ncol(mam_noMajorPopSize)){
  mam_noMajorPopSize[which(is.na(mam_noMajorPopSize[,j])),j] <- 0
  mam_majorPopSize_effSize5[which(is.na(mam_majorPopSize_effSize5[,j])),j] <- 0
  mam_majorPopSize_effSize7[which(is.na(mam_majorPopSize_effSize7[,j])),j] <- 0
  mam_majorPopSize_effSize9[which(is.na(mam_majorPopSize_effSize9[,j])),j] <- 0
}

#######################
# check for extinction
#######################

mam_polyExt <- rep(NA,81)
mam_OligExt_effSize5 <- rep(NA,81)
mam_OligExt_effSize7 <- rep(NA,81)
mam_OligExt_effSize9 <- rep(NA,81)

for(i in 1:81){
  mam_polyExt[i] <- sum(mam_noMajorPopSize[,i] == 0)
  mam_OligExt_effSize5[i] <- sum(mam_majorPopSize_effSize5[,i] == 0)
  mam_OligExt_effSize7[i] <- sum(mam_majorPopSize_effSize7[,i] == 0)
  mam_OligExt_effSize9[i] <- sum(mam_majorPopSize_effSize9[,i] == 0)
}


############################################
# plot results
############################################
par(mfrow=c(3,1),mar=c(2,5,3,0.5),xpd=FALSE)
library(scales)
carryCap <- 1000
##### mammal population size
plot(c(0,81),c(0,carryCap + 200),type="n",ylab=expression(paste(italic(""*N*"")," (thousands)",sep="")),xlab="",cex.lab=1.5,
     cex.main=1.5,cex.axis=0.9)



oligSizeMean_effSize5 <- rep(NA,81)
oligSizeMean_effSize7 <- rep(NA,81)
oligSizeMean_effSize9 <- rep(NA,81)
polySizeMean <- rep(NA,81)

for(i in 1:81){
  polySizeMean[i] <- mean(mam_noMajorPopSize[,i])
  oligSizeMean_effSize5[i] <- mean(mam_majorPopSize_effSize5[,i])
  oligSizeMean_effSize7[i] <- mean(mam_majorPopSize_effSize7[,i])
  oligSizeMean_effSize9[i] <- mean(mam_majorPopSize_effSize9[,i])
}
lines(1:81,polySizeMean,lwd=3,col="orange")
lines(1:81,oligSizeMean_effSize5,lwd=3,col="darkblue")
lines(1:81,oligSizeMean_effSize7,lwd=3,col="darkred")
lines(1:81,oligSizeMean_effSize9,lwd=3,col="darkgreen")

lines(c(0,200),c(carryCap,carryCap),lty="dashed")

#---------------------------------
##### plot the phenotypes
#---------------------------------
par(mar=c(3,5,2,0.5))
idealPhen <- 110



plot(c(0,80),c(100,112),type="n",ylab="Mean phenotype",xlab="",cex.lab=1.5,
     cex.main=1.5,cex.axis=0.9)



oligPhenMean_effSize5 <- rep(NA,80)
oligPhenMean_effSize7 <- rep(NA,80)
oligPhenMean_effSize9 <- rep(NA,80)
polyPhenMean <- rep(NA,80)

for(i in 1:80){
  polyPhenMean[i] <- mean(mam_noMajorPhen[,i],na.rm=TRUE)
  oligPhenMean_effSize5[i] <- mean(mam_phen_effSize5[,i],na.rm=TRUE)
  oligPhenMean_effSize7[i] <- mean(mam_phen_effSize7[,i],na.rm=TRUE)
  oligPhenMean_effSize9[i] <- mean(mam_phen_effSize9[,i],na.rm=TRUE)
}
lines(1:80,polyPhenMean,lwd=3,col="orange")
lines(1:80,oligPhenMean_effSize5,lwd=3,col="darkblue")
lines(1:80,oligPhenMean_effSize7,lwd=3,col="darkred")
lines(1:80,oligPhenMean_effSize9,lwd=3,col="darkgreen")
lines(c(0,80),c(110,110),lty="dashed")

legend(lty="solid",col=c("orange","darkblue","darkred","darkgreen"),x=20,y=100,xjust=FALSE,yjust=FALSE,
       legend = c("Polygenic",paste("Vqtl/Vg = ",c("0.5","0.7","0.9"))),lwd=2,bty="n")



#####################################################################
# plot the fraction of extinct population versus time
#####################################################################

par(mar=c(4,5,1,0.5))
# extinct proportion for mammals
plot(c(1,81),c(0,0.7),type="n",xlab="Generation",ylab="Proportion Extinct",cex.lab=1.5)
lines(1:81,mam_polyExt/nrow(mam_majorPopSize_effSize5),col="orange",lwd=3)
lines(1:81,mam_OligExt_effSize5/nrow(mam_majorPopSize_effSize5),col="darkblue",lwd=3)
lines(1:81,mam_OligExt_effSize7/nrow(mam_majorPopSize_effSize5),col="darkred",lwd=3)
lines(1:81,mam_OligExt_effSize9/nrow(mam_majorPopSize_effSize5),col="darkgreen",lwd=3)


###################################################################
###################################################################
# analyze results where the phenotypic optimum is a moving target
# with the to the new phenothpic optimum being either 10 or 20 generations
# FIGURE S13
###################################################################
###################################################################
propAnalyze <- 0.9
####--------------------------------------
#### time is 10 generations to new optimum
####--------------------------------------
setwd("~/Documents/work/rapidAdaptation/AmNatRevision/simulations/originalSims/h2Pt6/fecundity4/movingPhenoOpt")
popFileName <- "popSize_N0_500_K_1000_propMaj_0.9_gensToNewOpt_10"
phenFileName <-  "phen_N0_500_K_1000_propMaj_0.9_gensToNewOpt_10"
h2FileName <-  "h2_N0_500_K_1000_propMaj_0.9_gensToNewOpt_10"

mam_majorPopSize <- read.table(popFileName,header=TRUE)
mam_majorPhen <- read.table(phenFileName,header=TRUE)
mam_majorH2 <- read.table(h2FileName,header=TRUE)



mam_noMajorPopSize <- read.table("popSize_N0_500_K_1000_propMaj_0.01_gensToNewOpt_10",header=TRUE)
mam_noMajorPhen <- read.table("phen_N0_500_K_1000_propMaj_0.01_gensToNewOpt_10",header=TRUE)
mam_noMajorH2 <- read.table("h2_N0_500_K_1000_propMaj_0.01_gensToNewOpt_10",header=TRUE)

#-----------------------------------------------
# make NA population sizes zero
#-----------------------------------------------


for(j in 1:ncol(mam_noMajorPopSize)){
  mam_noMajorPopSize[which(is.na(mam_noMajorPopSize[,j])),j] <- 0
  mam_majorPopSize[which(is.na(mam_majorPopSize[,j])),j] <- 0
}


#######################
# check for extinction
#######################

mam_polyExt <- rep(NA,81)
mam_OligExt <- rep(NA,81)

for(i in 1:81){
  mam_polyExt[i] <- sum(mam_noMajorPopSize[,i] == 0)
  mam_OligExt[i] <- sum(mam_majorPopSize[,i] == 0)
}


#############################################
# get bootstrap CIs for proportion extinct
#############################################

bootReps <- 1000
lowPolyBoots <- matrix(NA,nrow=bootReps,ncol=81)
lowOligBoots <- matrix(NA,nrow=bootReps,ncol=81)

for(i in 1:bootReps){
  lowPolyDat <- mam_noMajorPopSize[sample(1:nrow(mam_noMajorPopSize),nrow(mam_noMajorPopSize),replace=TRUE),]
  lowOligDat <- mam_majorPopSize[sample(1:nrow(mam_majorPopSize),nrow(mam_majorPopSize),replace=TRUE),]
  lowPolyBoots [i,] <- colSums(lowPolyDat == 0)/nrow(lowPolyDat )
  lowOligBoots [i,] <- colSums(lowOligDat == 0)/nrow(lowOligDat )
  print(i)
}

mamlowPolyCI <- NULL
mamlowOligCI <- NULL

for(i in 1:81){
  mamlowPolyCI <- cbind(mamlowPolyCI,quantile(lowPolyBoots[,i],c(0.025,0.975)))
  mamlowOligCI <- cbind(mamlowOligCI,quantile(lowOligBoots[,i],c(0.025,0.975)))
}


####--------------------------------------
#### time is 20 generations to new optimum
####--------------------------------------
popFileName <- "popSize_N0_500_K_1000_propMaj_0.9_gensToNewOpt_20"
phenFileName <-  "phen_N0_500_K_1000_propMaj_0.9_gensToNewOpt_20"
h2FileName <-  "h2_N0_500_K_1000_propMaj_0.9_gensToNewOpt_20"

cor_majorPopSize <- read.table(popFileName,header=TRUE)
cor_majorPhen <- read.table(phenFileName,header=TRUE)
cor_majorH2 <- read.table(h2FileName,header=TRUE)



cor_noMajorPopSize <- read.table( "popSize_N0_500_K_1000_propMaj_0.01_gensToNewOpt_20",header=TRUE)
cor_noMajorPhen <- read.table("phen_N0_500_K_1000_propMaj_0.01_gensToNewOpt_20",header=TRUE)
cor_noMajorH2 <- read.table("h2_N0_500_K_1000_propMaj_0.01_gensToNewOpt_20",header=TRUE)

#-----------------------------------------------
# make NA population sizes zero
#-----------------------------------------------


for(j in 1:ncol(cor_noMajorPopSize)){
  cor_noMajorPopSize[which(is.na(cor_noMajorPopSize[,j])),j] <- 0
  cor_majorPopSize[which(is.na(cor_majorPopSize[,j])),j] <- 0
}


#######################
# check for extinction
#######################

cor_polyExt <- rep(NA,81)
cor_OligExt <- rep(NA,81)

for(i in 1:81){
  cor_polyExt[i] <- sum(cor_noMajorPopSize[,i] == 0)
  cor_OligExt[i] <- sum(cor_majorPopSize[,i] == 0)
}


#############################################
# get bootstrap CIs for proportion extinct
#############################################

bootReps <- 1000
cor_lowPolyBoots <- matrix(NA,nrow=bootReps,ncol=81)
cor_lowOligBoots <- matrix(NA,nrow=bootReps,ncol=81)

for(i in 1:bootReps){
  lowPolyDat <- cor_noMajorPopSize[sample(1:nrow(cor_noMajorPopSize),nrow(cor_noMajorPopSize),replace=TRUE),]
  lowOligDat <- cor_majorPopSize[sample(1:nrow(cor_majorPopSize),nrow(cor_majorPopSize),replace=TRUE),]
  cor_lowPolyBoots [i,] <- colSums(lowPolyDat == 0)/nrow(lowPolyDat )
  cor_lowOligBoots [i,] <- colSums(lowOligDat == 0)/nrow(lowOligDat )
  print(i)
}

corlowPolyCI <- NULL
corlowOligCI <- NULL

for(i in 1:81){
  corlowPolyCI <- cbind(corlowPolyCI,quantile(cor_lowPolyBoots[,i],c(0.025,0.975)))
  corlowOligCI <- cbind(corlowOligCI,quantile(cor_lowOligBoots[,i],c(0.025,0.975)))
}

############################################
# plot results
############################################
par(mfrow=c(3,2),mar=c(2,5,3,0.5),xpd=FALSE)
library(scales)
carryCap <- 1000
##### 10 generation to the new optimum
plot(c(0,81),c(0,carryCap + 200)/1000,type="n",ylab=expression(paste(italic(""*N*"")," (thousands)",sep="")),xlab="",cex.lab=1.5,
     cex.main=1.5,cex.axis=0.9)
for(i in 1:nrow(mam_noMajorPopSize)){
  lines(1:81,mam_noMajorPopSize[i,]/1000,col=alpha("orange",alpha=0.05))
  lines(1:81,mam_majorPopSize[i,]/1000,col=alpha("darkblue",alpha=0.05))
}
polySizeMean <- rep(NA,81)
oligSizeMean <- rep(NA,81)
for(i in 1:81){
  polySizeMean[i] <- mean(mam_noMajorPopSize[,i])/1000
  oligSizeMean[i] <- mean(mam_majorPopSize[,i])/1000
}
lines(1:81,polySizeMean,lwd=3,col="orange")
lines(1:81,oligSizeMean,lwd=3,col="darkblue")
lines(c(0,200),c(carryCap/1000,carryCap/1000),lty="dashed")
par(xpd=TRUE)
text(x=-25,y=1.4,labels="A",cex=2)
par(xpd=FALSE)
##### 20 generations to the new optimum
plot(c(0,81),c(0,carryCap + 200)/1000,type="n",ylab="",xlab="",cex.lab=1.5,
     cex.main=1.5,cex.axis=0.9)
carryCap <- 10000
for(i in 1:nrow(cor_noMajorPopSize)){
  lines(1:81,cor_noMajorPopSize[i,]/1000,col=alpha("orange",alpha=0.05))
  lines(1:81,cor_majorPopSize[i,]/1000,col=alpha("darkblue",alpha=0.05))
}
polySizeMean <- rep(NA,81)
oligSizeMean <- rep(NA,81)
for(i in 1:81){
  polySizeMean[i] <- mean(cor_noMajorPopSize[,i]/1000)
  oligSizeMean[i] <- mean(cor_majorPopSize[,i]/1000)
}
lines(1:81,polySizeMean,lwd=3,col="orange")
lines(1:81,oligSizeMean,lwd=3,col="darkblue")
lines(c(0,200),c(20000/1000,20000/1000),lty="dashed")
par(xpd=TRUE)
text(x=-25,y=1.4,labels="B",cex=2)
par(xpd=FALSE)




#---------------------------------
##### plot the phenotypes
#---------------------------------
par(mar=c(3,5,2,0.5))
idealPhen <- 110

plot(c(0,80),c(95,112),type="n",ylab="Mean Phenotype",xlab="",cex.lab=1.5,cex.axis=0.8)
for(i in 1:500){
  lines(1:80,mam_noMajorPhen[i,],col=alpha("orange",alpha=0.05))
  lines(1:80,mam_majorPhen[i,],col=alpha("darkblue",alpha=0.05))
}

# get rolling mean
polyPhenMean <- rep(NA,80)
oligPhenMean <- rep(NA,80)
for(i in 1:80){
  polyPhenMean[i] <- mean(mam_noMajorPhen[,i],na.rm=TRUE)
  oligPhenMean[i] <- mean(mam_majorPhen[,i],na.rm=TRUE)
}
lines(1:80,polyPhenMean,lwd=3,col="orange")
lines(1:80,oligPhenMean,lwd=3,col="darkblue")
lines(c(0,80),c(110,110),lty="dashed")


par(mar=c(3,5,2,0.5))
idealPhen <- 110


plot(c(0,80),c(95,112),type="n",ylab="",xlab="",cex.lab=1.5,cex.axis=0.8)
for(i in 1:500){
  lines(1:80,cor_noMajorPhen[i,],col=alpha("orange",alpha=0.05))
  lines(1:80,cor_majorPhen[i,],col=alpha("darkblue",alpha=0.05))
}

# get rolling mean
polyPhenMean <- rep(NA,80)
oligPhenMean <- rep(NA,80)
for(i in 1:80){
  polyPhenMean[i] <- mean(cor_noMajorPhen[,i],na.rm=TRUE)
  oligPhenMean[i] <- mean(cor_majorPhen[,i],na.rm=TRUE)
}
lines(1:80,polyPhenMean,lwd=3,col="orange")
lines(1:80,oligPhenMean,lwd=3,col="darkblue")
lines(c(0,80),c(110,110),lty="dashed")


#####################################################################
# plot the fraction of extinct population versus time
#####################################################################

par(mar=c(4,5,1,0.5))
# extinct proportion for mammals
plot(c(1,81),c(0,0.8),type="n",xlab="Generation",ylab="Proportion Extinct",cex.lab=1.5)
lines(1:81,mam_polyExt/nrow(lowPolyDat),col="orange",lwd=3)
lines(1:81,mam_OligExt/nrow(lowOligDat),col="darkblue",lwd=3)

polygon(x=c(1:81,rev(1:81)),y=c(mamlowPolyCI[1,1:81],rev(mamlowPolyCI[2,1:81])),col=adjustcolor("orange",alpha.f=0.5),border=NA)
polygon(x=c(1:81,rev(1:81)),y=c(mamlowOligCI[1,1:81],rev(mamlowOligCI[2,1:81])),col=adjustcolor("darkblue",alpha.f=0.5),border=NA)

par(mar=c(4,5,1,0.5))
plot(c(1,81),c(0,0.8),type="n",xlab="Generation",ylab="",cex.lab=1.5)
lines(1:81,cor_polyExt/nrow(lowPolyDat),col="orange",lwd=3)
lines(1:81,cor_OligExt/nrow(lowOligDat),col="darkblue",lwd=3)

legend(x=20,y=-0.05,xjust=FALSE,yjust=FALSE,legend=c("No Major Locus","Major Locus"),
       lwd=3,col=c("orange","darkblue"),bty="n")

polygon(x=c(1:81,rev(1:81)),y=c(corlowPolyCI[1,1:81],rev(corlowPolyCI[2,1:81])),col=adjustcolor("orange",alpha.f=0.5),border=NA)
polygon(x=c(1:81,rev(1:81)),y=c(corlowOligCI[1,1:81],rev(corlowOligCI[2,1:81])),col=adjustcolor("darkblue",alpha.f=0.5),border=NA)



###################################################################
###################################################################
# analyze results with mammal life history and linked loci on
# 10 chromosomes
# FIGURE S14
###################################################################
###################################################################
propAnalyze <- 0.9
####--------------------------------------
#### collect mammal results
####--------------------------------------

setwd("~/Documents/work/rapidAdaptation/AmNatRevision/simulations/originalSims/h2Pt6/fecundity4/linkage/major")

pedFiles <- paste("pedObject_propMajor_0.9_100Gens_",1:500,sep="")
NMat_major <- matrix(NA,nrow=500,ncol=81)
phenMat_major <- matrix(NA,nrow=500,ncol=81)
for(i in 1:length(pedFiles)){
  thisDat <- read.table(pedFiles[i],header=TRUE)
  theseNs <- rep(NA,81)
  meanPhen <- rep(NA,81)
  for(j in 1:max(thisDat[,4])){
    theseNs [j]<- sum(thisDat[,4] == j)
    meanPhen[j] <- mean(thisDat[thisDat[,4] == j,6])
  }
  NMat_major[i,] <- theseNs
  phenMat_major[i,] <- meanPhen
}


setwd("~/Documents/work/rapidAdaptation/AmNatRevision/simulations/originalSims/h2Pt6/fecundity4/linkage/noMajor")
pedFiles <- paste("pedObject_propMajor_0.01_100Gens_",1:500,sep="")
pedFiles <- paste("pedObject_propMajor_0.01_100Gens_",1:500,sep="")
NMat_noMajor <- matrix(NA,nrow=500,ncol=81)
phenMat_noMajor <- matrix(NA,nrow=500,ncol=81)
for(i in 1:length(pedFiles)){
  thisDat <- read.table(pedFiles[i],header=TRUE)
  theseNs <- rep(NA,81)
  meanPhen <- rep(NA,81)
  for(j in 1:max(thisDat[,4])){
    theseNs [j]<- sum(thisDat[,4] == j)
    meanPhen[j] <- mean(thisDat[thisDat[,4] == j,6])
  }
  NMat_noMajor[i,] <- theseNs
  phenMat_noMajor[i,] <- meanPhen
}


mam_majorPopSize <- NMat_major
mam_majorPhen <- phenMat_major

mam_noMajorPopSize <- NMat_noMajor
mam_noMajorPhen <- phenMat_noMajor


#-----------------------------------------------
# make NA population sizes zero
#-----------------------------------------------


for(j in 1:ncol(mam_noMajorPopSize)){
  mam_noMajorPopSize[which(is.na(mam_noMajorPopSize[,j])),j] <- 0
  mam_majorPopSize[which(is.na(mam_majorPopSize[,j])),j] <- 0
}


#######################
# check for extinction
#######################

mam_polyExt <- rep(NA,81)
mam_OligExt <- rep(NA,81)

for(i in 1:81){
  mam_polyExt[i] <- sum(mam_noMajorPopSize[,i] == 0)
  mam_OligExt[i] <- sum(mam_majorPopSize[,i] == 0)
}


#############################################
# get bootstrap CIs for proportion extinct
#############################################

bootReps <- 1000
lowPolyBoots <- matrix(NA,nrow=bootReps,ncol=81)
lowOligBoots <- matrix(NA,nrow=bootReps,ncol=81)

for(i in 1:bootReps){
  lowPolyDat <- mam_noMajorPopSize[sample(1:nrow(mam_noMajorPopSize),nrow(mam_noMajorPopSize),replace=TRUE),]
  lowOligDat <- mam_majorPopSize[sample(1:nrow(mam_majorPopSize),nrow(mam_majorPopSize),replace=TRUE),]
  lowPolyBoots [i,] <- colSums(lowPolyDat == 0)/nrow(lowPolyDat )
  lowOligBoots [i,] <- colSums(lowOligDat == 0)/nrow(lowOligDat )
  print(i)
}

mamlowPolyCI <- NULL
mamlowOligCI <- NULL

for(i in 1:81){
  mamlowPolyCI <- cbind(mamlowPolyCI,quantile(lowPolyBoots[,i],c(0.025,0.975)))
  mamlowOligCI <- cbind(mamlowOligCI,quantile(lowOligBoots[,i],c(0.025,0.975)))
}




############################################
# plot results
############################################
par(mfrow=c(3,1),mar=c(2,5,3,0.5),xpd=FALSE)
library(scales)
carryCap <- 1000
##### mammal population size
plot(c(0,81),c(0,carryCap + 200)/1000,type="n",ylab=expression(paste(italic(""*N*"")," (thousands)",sep="")),xlab="",cex.lab=1.5,
     cex.main=1.5,cex.axis=0.9)
for(i in 1:nrow(mam_noMajorPopSize)){
  lines(1:81,mam_noMajorPopSize[i,]/1000,col=alpha("orange",alpha=0.05))
  lines(1:81,mam_majorPopSize[i,]/1000,col=alpha("darkblue",alpha=0.05))
}
polySizeMean <- rep(NA,81)
oligSizeMean <- rep(NA,81)
for(i in 1:81){
  polySizeMean[i] <- mean(mam_noMajorPopSize[,i])/1000
  oligSizeMean[i] <- mean(mam_majorPopSize[,i])/1000
}
lines(1:81,polySizeMean,lwd=3,col="orange")
lines(1:81,oligSizeMean,lwd=3,col="darkblue")
lines(c(0,200),c(carryCap/1000,carryCap/1000),lty="dashed")

#---------------------------------
##### plot the phenotypes
#---------------------------------
par(mar=c(3,5,2,0.5))
idealPhen <- 110

####### mammals
plot(c(0,80),c(100,112),type="n",ylab="Mean Phenotype",xlab="",cex.lab=1.5,cex.axis=0.8)
for(i in 1:500){
  lines(1:81,mam_noMajorPhen[i,],col=alpha("orange",alpha=0.05))
  lines(1:81,mam_majorPhen[i,],col=alpha("darkblue",alpha=0.05))
}

# get rolling mean
polyPhenMean <- rep(NA,80)
oligPhenMean <- rep(NA,80)
for(i in 1:80){
  polyPhenMean[i] <- mean(mam_noMajorPhen[,i],na.rm=TRUE)
  oligPhenMean[i] <- mean(mam_majorPhen[,i],na.rm=TRUE)
}
lines(1:80,polyPhenMean,lwd=3,col="orange")
lines(1:80,oligPhenMean,lwd=3,col="darkblue")
lines(c(0,80),c(110,110),lty="dashed")


par(mar=c(3,5,2,0.5))
idealPhen <- 110

#####################################################################
# plot the fraction of extinct population versus time
#####################################################################

par(mar=c(4,5,1,0.5))
# extinct proportion for mammals
plot(c(1,81),c(0,0.8),type="n",xlab="Generation",ylab="Proportion Extinct",cex.lab=1.5)
lines(1:81,mam_polyExt/nrow(lowPolyDat),col="orange",lwd=3)
lines(1:81,mam_OligExt/nrow(lowOligDat),col="darkblue",lwd=3)

polygon(x=c(1:81,rev(1:81)),y=c(mamlowPolyCI[1,1:81],rev(mamlowPolyCI[2,1:81])),col=adjustcolor("orange",alpha.f=0.5),border=NA)
polygon(x=c(1:81,rev(1:81)),y=c(mamlowOligCI[1,1:81],rev(mamlowOligCI[2,1:81])),col=adjustcolor("darkblue",alpha.f=0.5),border=NA)
# get pictures

setwd("~/Documents/work/rapidAdaptation/AmNatRevision/simulations/originalSims/h2Pt6")
library(jpeg)
img<-readJPEG("wildebeest.jpg")
rasterImage(img,xleft=-3, ybottom=0.62, xright=25, ytop=0.84)



###################################################################
###################################################################
# analyze results where the phenotypic optimum is a moving target
# and there is phenotypic plasticity
# FIGURE S15
###################################################################
###################################################################
####--------------------------------------
#### plasticity parameter m = 0.1
####--------------------------------------
setwd("~/Documents/work/rapidAdaptation/AmNatRevision/simulations/originalSims/h2Pt6/fecundity4/plasticity")
mam_majorPopSize <- read.table("popSize_N0_500_K_1000_propMaj_0.9_gensToNewOpt_10_plast_0.1",header=TRUE)
mam_majorPhen <- read.table("phen_N0_500_K_1000_propMaj_0.9_gensToNewOpt_10_plast_0.1",header=TRUE)
mam_majorH2 <- read.table("h2_N0_500_K_1000_propMaj_0.9_gensToNewOpt_10_plast_0.1",header=TRUE)
mam_noMajorPopSize <- read.table("popSize_N0_500_K_1000_propMaj_0.01_gensToNewOpt_10_plast_0.1",header=TRUE)
mam_noMajorPhen <- read.table("phen_N0_500_K_1000_propMaj_0.01_gensToNewOpt_10_plast_0.1",header=TRUE)
mam_noMajorH2 <- read.table("h2_N0_500_K_1000_propMaj_0.01_gensToNewOpt_10_plast_0.1",header=TRUE)

#-----------------------------------------------
# make NA population sizes zero
#-----------------------------------------------
for(j in 1:ncol(mam_noMajorPopSize)){
  mam_noMajorPopSize[which(is.na(mam_noMajorPopSize[,j])),j] <- 0
  mam_majorPopSize[which(is.na(mam_majorPopSize[,j])),j] <- 0
}


#######################
# check for extinction
#######################
mam_polyExt <- rep(NA,81)
mam_OligExt <- rep(NA,81)
for(i in 1:81){
  mam_polyExt[i] <- sum(mam_noMajorPopSize[,i] == 0)
  mam_OligExt[i] <- sum(mam_majorPopSize[,i] == 0)
}


#############################################
# get bootstrap CIs for proportion extinct
#############################################
bootReps <- 1000
lowPolyBoots <- matrix(NA,nrow=bootReps,ncol=81)
lowOligBoots <- matrix(NA,nrow=bootReps,ncol=81)
for(i in 1:bootReps){
  lowPolyDat <- mam_noMajorPopSize[sample(1:nrow(mam_noMajorPopSize),nrow(mam_noMajorPopSize),replace=TRUE),]
  lowOligDat <- mam_majorPopSize[sample(1:nrow(mam_majorPopSize),nrow(mam_majorPopSize),replace=TRUE),]
  lowPolyBoots [i,] <- colSums(lowPolyDat == 0)/nrow(lowPolyDat )
  lowOligBoots [i,] <- colSums(lowOligDat == 0)/nrow(lowOligDat )
  print(i)
}

mamlowPolyCI <- NULL
mamlowOligCI <- NULL
for(i in 1:81){
  mamlowPolyCI <- cbind(mamlowPolyCI,quantile(lowPolyBoots[,i],c(0.025,0.975)))
  mamlowOligCI <- cbind(mamlowOligCI,quantile(lowOligBoots[,i],c(0.025,0.975)))
}




####--------------------------------------
#### plasticity with me = 0.2
####--------------------------------------
cor_majorPopSize <- read.table("popSize_N0_500_K_1000_propMaj_0.9_gensToNewOpt_10_plast_0.2",header=TRUE)
cor_majorPhen <- read.table("phen_N0_500_K_1000_propMaj_0.9_gensToNewOpt_10_plast_0.2",header=TRUE)
cor_majorH2 <- read.table("h2_N0_500_K_1000_propMaj_0.9_gensToNewOpt_10_plast_0.2",header=TRUE)

cor_noMajorPopSize <- read.table("popSize_N0_500_K_1000_propMaj_0.01_gensToNewOpt_10_plast_0.2",header=TRUE)
cor_noMajorPhen <- read.table("phen_N0_500_K_1000_propMaj_0.01_gensToNewOpt_10_plast_0.2",header=TRUE)
cor_noMajorH2 <- read.table("h2_N0_500_K_1000_propMaj_0.01_gensToNewOpt_10_plast_0.2",header=TRUE)

#-----------------------------------------------
# make NA population sizes zero
#-----------------------------------------------
for(j in 1:ncol(cor_noMajorPopSize)){
  cor_noMajorPopSize[which(is.na(cor_noMajorPopSize[,j])),j] <- 0
  cor_majorPopSize[which(is.na(cor_majorPopSize[,j])),j] <- 0
}

#######################
# check for extinction
#######################
cor_polyExt <- rep(NA,81)
cor_OligExt <- rep(NA,81)
for(i in 1:81){
  cor_polyExt[i] <- sum(cor_noMajorPopSize[,i] == 0)
  cor_OligExt[i] <- sum(cor_majorPopSize[,i] == 0)
}

#############################################
# get bootstrap CIs for proportion extinct
#############################################
bootReps <- 1000
cor_lowPolyBoots <- matrix(NA,nrow=bootReps,ncol=81)
cor_lowOligBoots <- matrix(NA,nrow=bootReps,ncol=81)
for(i in 1:bootReps){
  lowPolyDat <- cor_noMajorPopSize[sample(1:nrow(cor_noMajorPopSize),nrow(cor_noMajorPopSize),replace=TRUE),]
  lowOligDat <- cor_majorPopSize[sample(1:nrow(cor_majorPopSize),nrow(cor_majorPopSize),replace=TRUE),]
  cor_lowPolyBoots [i,] <- colSums(lowPolyDat == 0)/nrow(lowPolyDat )
  cor_lowOligBoots [i,] <- colSums(lowOligDat == 0)/nrow(lowOligDat )
  print(i)
}

corlowPolyCI <- NULL
corlowOligCI <- NULL
for(i in 1:81){
  corlowPolyCI <- cbind(corlowPolyCI,quantile(cor_lowPolyBoots[,i],c(0.025,0.975)))
  corlowOligCI <- cbind(corlowOligCI,quantile(cor_lowOligBoots[,i],c(0.025,0.975)))
}


####--------------------------------------
#### plasticity with me = 0.4
####--------------------------------------
cor_majorPopSize_4 <- read.table("popSize_N0_500_K_1000_propMaj_0.9_gensToNewOpt_10_plast_0.4",header=TRUE)
cor_majorPhen_4 <- read.table("phen_N0_500_K_1000_propMaj_0.9_gensToNewOpt_10_plast_0.4",header=TRUE)
cor_majorH2_4 <- read.table("h2_N0_500_K_1000_propMaj_0.9_gensToNewOpt_10_plast_0.4",header=TRUE)

cor_noMajorPopSize_4 <- read.table("popSize_N0_500_K_1000_propMaj_0.01_gensToNewOpt_10_plast_0.4",header=TRUE)
cor_noMajorPhen_4 <- read.table("phen_N0_500_K_1000_propMaj_0.01_gensToNewOpt_10_plast_0.4",header=TRUE)
cor_noMajorH2_4 <- read.table("h2_N0_500_K_1000_propMaj_0.01_gensToNewOpt_10_plast_0.4",header=TRUE)

#-----------------------------------------------
# make NA population sizes zero
#-----------------------------------------------
for(j in 1:ncol(cor_noMajorPopSize)){
  cor_noMajorPopSize_4[which(is.na(cor_noMajorPopSize_4[,j])),j] <- 0
  cor_majorPopSize_4[which(is.na(cor_majorPopSize_4[,j])),j] <- 0
}

#######################
# check for extinction
#######################
cor_polyExt_4 <- rep(NA,81)
cor_OligExt_4 <- rep(NA,81)
for(i in 1:81){
  cor_polyExt_4[i] <- sum(cor_noMajorPopSize_4[,i] == 0)
  cor_OligExt_4[i] <- sum(cor_majorPopSize_4[,i] == 0)
}


#############################################
# get bootstrap CIs for proportion extinct
#############################################

bootReps <- 1000
cor_lowPolyBoots_4 <- matrix(NA,nrow=bootReps,ncol=81)
cor_lowOligBoots_4 <- matrix(NA,nrow=bootReps,ncol=81)

for(i in 1:bootReps){
  lowPolyDat_4 <- cor_noMajorPopSize_4[sample(1:nrow(cor_noMajorPopSize_4),nrow(cor_noMajorPopSize_4),replace=TRUE),]
  lowOligDat_4 <- cor_majorPopSize_4[sample(1:nrow(cor_majorPopSize_4),nrow(cor_majorPopSize_4),replace=TRUE),]
  cor_lowPolyBoots_4 [i,] <- colSums(lowPolyDat_4 == 0)/nrow(lowPolyDat_4 )
  cor_lowOligBoots_4 [i,] <- colSums(lowOligDat_4 == 0)/nrow(lowOligDat_4 )
  print(i)
}

corlowPolyCI_4 <- NULL
corlowOligCI_4 <- NULL

for(i in 1:81){
  corlowPolyCI_4 <- cbind(corlowPolyCI_4,quantile(cor_lowPolyBoots_4[,i],c(0.025,0.975)))
  corlowOligCI_4 <- cbind(corlowOligCI_4,quantile(cor_lowOligBoots_4[,i],c(0.025,0.975)))
}






############################################
# plot results
############################################
par(mfrow=c(3,3),mar=c(2,5,3,0.5),xpd=FALSE)
library(scales)
carryCap <- 1000
##### mammal population size
plot(c(0,81),c(0,carryCap + 200)/1000,type="n",ylab=expression(paste(italic(""*N*"")," (thousands)",sep="")),xlab="",cex.lab=1.5,
     cex.main=1.5,cex.axis=0.9)
for(i in 1:nrow(mam_noMajorPopSize)){
  lines(1:81,mam_noMajorPopSize[i,]/1000,col=alpha("orange",alpha=0.02))
  lines(1:81,mam_majorPopSize[i,]/1000,col=alpha("darkblue",alpha=0.02))
}
polySizeMean <- rep(NA,81)
oligSizeMean <- rep(NA,81)
for(i in 1:81){
  polySizeMean[i] <- mean(mam_noMajorPopSize[,i])/1000
  oligSizeMean[i] <- mean(mam_majorPopSize[,i])/1000
}
lines(1:81,polySizeMean,lwd=3,col="orange")
lines(1:81,oligSizeMean,lwd=3,col="darkblue")
lines(c(0,200),c(carryCap/1000,carryCap/1000),lty="dashed")
par(xpd=TRUE)
text(x=-25,y=1.4,labels="A",cex=2)
par(xpd=FALSE)


##### coral population size
plot(c(0,81),c(0,carryCap + 200)/1000,type="n",ylab="",xlab="",cex.lab=1.5,
     cex.main=1.5,cex.axis=0.9)
carryCap <- 10000
for(i in 1:nrow(cor_noMajorPopSize)){
  lines(1:81,cor_noMajorPopSize[i,]/1000,col=alpha("orange",alpha=0.02))
  lines(1:81,cor_majorPopSize[i,]/1000,col=alpha("darkblue",alpha=0.02))
}
polySizeMean <- rep(NA,81)
oligSizeMean <- rep(NA,81)
for(i in 1:81){
  polySizeMean[i] <- mean(cor_noMajorPopSize[,i]/1000)
  oligSizeMean[i] <- mean(cor_majorPopSize[,i]/1000)
}
lines(1:81,polySizeMean,lwd=3,col="orange")
lines(1:81,oligSizeMean,lwd=3,col="darkblue")
lines(c(0,200),c(1,1),lty="dashed")
par(xpd=TRUE)
text(x=-25,y=1.4,labels="B",cex=2)
par(xpd=FALSE)



##### m = 0.4
plot(c(0,81),c(0,1200)/1000,type="n",ylab="",xlab="",cex.lab=1.5,
     cex.main=1.5,cex.axis=0.9)
carryCap <- 10000
for(i in 1:nrow(cor_noMajorPopSize_4)){
  lines(1:81,cor_noMajorPopSize_4[i,]/1000,col=alpha("orange",alpha=0.02))
  lines(1:81,cor_majorPopSize_4[i,]/1000,col=alpha("darkblue",alpha=0.02))
}
polySizeMean <- rep(NA,81)
oligSizeMean <- rep(NA,81)
for(i in 1:81){
  polySizeMean[i] <- mean(cor_noMajorPopSize_4[,i]/1000)
  oligSizeMean[i] <- mean(cor_majorPopSize_4[,i]/1000)
}
lines(1:81,polySizeMean,lwd=3,col="orange")
lines(1:81,oligSizeMean,lwd=3,col="darkblue")
lines(c(0,200),c(1,1),lty="dashed")
par(xpd=TRUE)
text(x=-25,y=1.4,labels="C",cex=2)
par(xpd=FALSE)





#---------------------------------
##### plot the phenotypes
#---------------------------------
par(mar=c(3,5,2,0.5))
idealPhen <- 110

####### mammals
plot(c(0,80),c(100,112),type="n",ylab="Mean Phenotype",xlab="",cex.lab=1.5,cex.axis=0.8)
for(i in 1:500){
  lines(1:80,mam_noMajorPhen[i,],col=alpha("orange",alpha=0.02))
  lines(1:80,mam_majorPhen[i,],col=alpha("darkblue",alpha=0.02))
}

# get rolling mean
polyPhenMean <- rep(NA,80)
oligPhenMean <- rep(NA,80)
for(i in 1:80){
  polyPhenMean[i] <- mean(mam_noMajorPhen[,i],na.rm=TRUE)
  oligPhenMean[i] <- mean(mam_majorPhen[,i],na.rm=TRUE)
}
lines(1:80,polyPhenMean,lwd=3,col="orange")
lines(1:80,oligPhenMean,lwd=3,col="darkblue")
lines(c(0,80),c(110,110),lty="dashed")




####### m = 0.2
plot(c(0,80),c(100,112),type="n",ylab="",xlab="",cex.lab=1.5,cex.axis=0.8)
for(i in 1:500){
  lines(1:80,cor_noMajorPhen[i,],col=alpha("orange",alpha=0.02))
  lines(1:80,cor_majorPhen[i,],col=alpha("darkblue",alpha=0.02))
}

# get rolling mean
polyPhenMean <- rep(NA,80)
oligPhenMean <- rep(NA,80)
for(i in 1:80){
  polyPhenMean[i] <- mean(cor_noMajorPhen[,i],na.rm=TRUE)
  oligPhenMean[i] <- mean(cor_majorPhen[,i],na.rm=TRUE)
}
lines(1:80,polyPhenMean,lwd=3,col="orange")
lines(1:80,oligPhenMean,lwd=3,col="darkblue")
lines(c(0,80),c(110,110),lty="dashed")



####### m = 0.4
plot(c(0,80),c(100,112),type="n",ylab="",xlab="",cex.lab=1.5,cex.axis=0.8)
for(i in 1:500){
  lines(1:80,cor_noMajorPhen_4[i,],col=alpha("orange",alpha=0.02))
  lines(1:80,cor_majorPhen_4[i,],col=alpha("darkblue",alpha=0.02))
}

# get rolling mean
polyPhenMean <- rep(NA,80)
oligPhenMean <- rep(NA,80)
for(i in 1:80){
  polyPhenMean[i] <- mean(cor_noMajorPhen_4[,i],na.rm=TRUE)
  oligPhenMean[i] <- mean(cor_majorPhen_4[,i],na.rm=TRUE)
}
lines(1:80,polyPhenMean,lwd=3,col="orange")
lines(1:80,oligPhenMean,lwd=3,col="darkblue")
lines(c(0,80),c(110,110),lty="dashed")

#####################################################################
# plot the fraction of extinct population versus time
#####################################################################
par(mar=c(4,5,1,0.5))
# extinct proportion for mammals
plot(c(1,81),c(0,0.55),type="n",xlab="Generation",ylab="Proportion Extinct",cex.lab=1.5)
lines(1:81,mam_polyExt/nrow(lowPolyDat),col="orange",lwd=3)
lines(1:81,mam_OligExt/nrow(lowOligDat),col="darkblue",lwd=3)

polygon(x=c(1:81,rev(1:81)),y=c(mamlowPolyCI[1,1:81],rev(mamlowPolyCI[2,1:81])),col=adjustcolor("orange",alpha.f=0.5),border=NA)
polygon(x=c(1:81,rev(1:81)),y=c(mamlowOligCI[1,1:81],rev(mamlowOligCI[2,1:81])),col=adjustcolor("darkblue",alpha.f=0.5),border=NA)


# m = 0.2
par(mar=c(4,5,1,0.5))
plot(c(1,81),c(0,0.55),type="n",xlab="Generation",ylab="",cex.lab=1.5)
lines(1:81,cor_polyExt/nrow(lowPolyDat),col="orange",lwd=3)
lines(1:81,cor_OligExt/nrow(lowOligDat),col="darkblue",lwd=3)
polygon(x=c(1:81,rev(1:81)),y=c(corlowPolyCI[1,1:81],rev(corlowPolyCI[2,1:81])),col=adjustcolor("orange",alpha.f=0.5),border=NA)
polygon(x=c(1:81,rev(1:81)),y=c(corlowOligCI[1,1:81],rev(corlowOligCI[2,1:81])),col=adjustcolor("darkblue",alpha.f=0.5),border=NA)

# m = 0.4
par(mar=c(4,5,1,0.5))
plot(c(1,81),c(0,0.55),type="n",xlab="Generation",ylab="",cex.lab=1.5)
lines(1:81,cor_polyExt_4/nrow(lowPolyDat_4),col="orange",lwd=3)
lines(1:81,cor_OligExt_4/nrow(lowOligDat_4),col="darkblue",lwd=3)
legend(x=0,y=0.2,xjust=FALSE,yjust=FALSE,legend=c("No Major Locus","Major Locus"),
       lwd=3,col=c("orange","darkblue"),bty="n")
polygon(x=c(1:81,rev(1:81)),y=c(corlowPolyCI_4[1,1:81],rev(corlowPolyCI_4[2,1:81])),col=adjustcolor("orange",alpha.f=0.5),border=NA)
polygon(x=c(1:81,rev(1:81)),y=c(corlowOligCI_4[1,1:81],rev(corlowOligCI_4[2,1:81])),col=adjustcolor("darkblue",alpha.f=0.5),border=NA)




##################################################################
#-----------------------------------------------------------------
# simulations with a smaller shift
# in the phenotypic optimum (theta), from 100 to 105 and 107.5.
# FIGURE S16
#-----------------------------------------------------------------
##################################################################

propAnalyze <- 0.9
####----------------------
#### theta shift to 105
####----------------------
setwd("~/Documents/work/rapidAdaptation/AmNatRevision/simulations/originalSims/h2Pt6/fecundity4/smallThetaShift/major")
popFileName <- paste("popSize_N0_500_K_1000_propMaj_",propAnalyze,sep="")
phenFileName <-  paste("phen_N0_500_K_1000_propMaj_",propAnalyze,sep="")
h2FileName <-  paste("h2_N0_500_K_1000_propMaj_",propAnalyze,sep="")

mam_majorPopSize <- read.table(popFileName,header=TRUE)
mam_majorPhen <- read.table(phenFileName,header=TRUE)
mam_majorH2 <- read.table(h2FileName,header=TRUE)


setwd("~/Documents/work/rapidAdaptation/AmNatRevision/simulations/originalSims/h2Pt6/fecundity4/smallThetaShift/noMajor")
popFileName <- paste("popSize_N0_500_K_1000_propMaj_0.01",sep="")
phenFileName <-  paste("phen_N0_500_K_1000_propMaj_0.01",sep="")
h2FileName <-  paste("h2_N0_500_K_1000_propMaj_0.01",sep="")

mam_noMajorPopSize <- read.table(popFileName,header=TRUE)
mam_noMajorPhen <- read.table(phenFileName,header=TRUE)
mam_noMajorH2 <- read.table(h2FileName,header=TRUE)
#-----------------------------------------------
# make NA population sizes zero
#-----------------------------------------------


for(j in 1:ncol(mam_noMajorPopSize)){
  mam_noMajorPopSize[which(is.na(mam_noMajorPopSize[,j])),j] <- 0
  mam_majorPopSize[which(is.na(mam_majorPopSize[,j])),j] <- 0
}


#######################
# check for extinction
#######################

mam_polyExt <- rep(NA,81)
mam_OligExt <- rep(NA,81)

for(i in 1:81){
  mam_polyExt[i] <- sum(mam_noMajorPopSize[,i] == 0)
  mam_OligExt[i] <- sum(mam_majorPopSize[,i] == 0)
}


#############################################
# get bootstrap CIs for proportion extinct
#############################################

bootReps <- 1000
lowPolyBoots <- matrix(NA,nrow=bootReps,ncol=81)
lowOligBoots <- matrix(NA,nrow=bootReps,ncol=81)

for(i in 1:bootReps){
  lowPolyDat <- mam_noMajorPopSize[sample(1:nrow(mam_noMajorPopSize),nrow(mam_noMajorPopSize),replace=TRUE),]
  lowOligDat <- mam_majorPopSize[sample(1:nrow(mam_majorPopSize),nrow(mam_majorPopSize),replace=TRUE),]
  lowPolyBoots [i,] <- colSums(lowPolyDat == 0)/nrow(lowPolyDat )
  lowOligBoots [i,] <- colSums(lowOligDat == 0)/nrow(lowOligDat )
  print(i)
}

mamlowPolyCI <- NULL
mamlowOligCI <- NULL

for(i in 1:81){
  mamlowPolyCI <- cbind(mamlowPolyCI,quantile(lowPolyBoots[,i],c(0.025,0.975)))
  mamlowOligCI <- cbind(mamlowOligCI,quantile(lowOligBoots[,i],c(0.025,0.975)))
}


####----------------------
#### theta shift to 107.5
####----------------------
setwd("~/Documents/work/rapidAdaptation/AmNatRevision/simulations/originalSims/h2Pt6/fecundity4/smallThetaShift/major_phenOpt7Pt5")
popFileName <- paste("popSize_N0_500_K_1000_propMaj_",propAnalyze,sep="")
phenFileName <-  paste("phen_N0_500_K_1000_propMaj_",propAnalyze,sep="")
h2FileName <-  paste("h2_N0_500_K_1000_propMaj_",propAnalyze,sep="")

mam_majorPopSize_biggerShift <- read.table(popFileName,header=TRUE)
mam_majorPhen_biggerShift <- read.table(phenFileName,header=TRUE)
mam_majorH2_biggerShift <- read.table(h2FileName,header=TRUE)


setwd("~/Documents/work/rapidAdaptation/AmNatRevision/simulations/originalSims/h2Pt6/fecundity4/smallThetaShift/noMajor_phenOpt7Pt5")

popFileName_biggerShift <- paste("popSize_N0_500_K_1000_propMaj_0.01",sep="")
phenFileName_biggerShift <-  paste("phen_N0_500_K_1000_propMaj_0.01",sep="")
h2FileName_biggerShift <-  paste("h2_N0_500_K_1000_propMaj_0.01",sep="")

mam_noMajorPopSize_biggerShift <- read.table(popFileName_biggerShift,header=TRUE)
mam_noMajorPhen_biggerShift <- read.table(phenFileName_biggerShift,header=TRUE)
mam_noMajorH2_biggerShift <- read.table(h2FileName_biggerShift,header=TRUE)
#-----------------------------------------------
# make NA population sizes zero
#-----------------------------------------------


for(j in 1:ncol(mam_noMajorPopSize_biggerShift)){
  mam_noMajorPopSize_biggerShift[which(is.na(mam_noMajorPopSize_biggerShift[,j])),j] <- 0
  mam_majorPopSize_biggerShift[which(is.na(mam_majorPopSize_biggerShift[,j])),j] <- 0
}


#######################
# check for extinction
#######################

mam_polyExt_biggerShift <- rep(NA,81)
mam_OligExt_biggerShift <- rep(NA,81)

for(i in 1:81){
  mam_polyExt_biggerShift[i] <- sum(mam_noMajorPopSize_biggerShift[,i] == 0)
  mam_OligExt_biggerShift[i] <- sum(mam_majorPopSize_biggerShift[,i] == 0)
}


#############################################
# get bootstrap CIs for proportion extinct
#############################################

bootReps <- 1000
lowPolyBoots_biggerShift <- matrix(NA,nrow=bootReps,ncol=81)
lowOligBoots_biggerShift <- matrix(NA,nrow=bootReps,ncol=81)

for(i in 1:bootReps){
  lowPolyDat_biggerShift <- mam_noMajorPopSize_biggerShift[sample(1:nrow(mam_noMajorPopSize_biggerShift),nrow(mam_noMajorPopSize_biggerShift),replace=TRUE),]
  lowOligDat_biggerShift <- mam_majorPopSize_biggerShift[sample(1:nrow(mam_majorPopSize_biggerShift),nrow(mam_majorPopSize_biggerShift),replace=TRUE),]
  lowPolyBoots_biggerShift [i,] <- colSums(lowPolyDat_biggerShift == 0)/nrow(lowPolyDat_biggerShift )
  lowOligBoots_biggerShift [i,] <- colSums(lowOligDat_biggerShift == 0)/nrow(lowOligDat_biggerShift )
  print(i)
}

mamlowPolyCI_biggerShift <- NULL
mamlowOligCI_biggerShift <- NULL

for(i in 1:81){
  mamlowPolyCI_biggerShift <- cbind(mamlowPolyCI_biggerShift,quantile(lowPolyBoots_biggerShift[,i],c(0.025,0.975)))
  mamlowOligCI_biggerShift <- cbind(mamlowOligCI_biggerShift,quantile(lowOligBoots_biggerShift[,i],c(0.025,0.975)))
}


############################################
# plot results
############################################
par(mfrow=c(3,2),mar=c(2,5,3,0.5),xpd=FALSE)
library(scales)
carryCap <- 1000
##### small shift population size
plot(c(0,81),c(0,carryCap + 200)/1000,type="n",ylab=expression(paste(italic(""*N*"")," (thousands)",sep="")),xlab="",cex.lab=1.5,
     cex.main=1.5,cex.axis=0.9)
for(i in 1:nrow(mam_noMajorPopSize)){
  lines(1:81,mam_noMajorPopSize[i,]/1000,col=alpha("orange",alpha=0.05))
  lines(1:81,mam_majorPopSize[i,]/1000,col=alpha("darkblue",alpha=0.05))
}
polySizeMean <- rep(NA,81)
oligSizeMean <- rep(NA,81)
for(i in 1:81){
  polySizeMean[i] <- mean(mam_noMajorPopSize[,i])/1000
  oligSizeMean[i] <- mean(mam_majorPopSize[,i])/1000
}
lines(1:81,polySizeMean,lwd=3,col="orange")
lines(1:81,oligSizeMean,lwd=3,col="darkblue")
lines(c(0,200),c(carryCap/1000,carryCap/1000),lty="dashed")
par(xpd=TRUE)
text(x=-25,y=1.4,labels="A",cex=2)
par(xpd=FALSE)

##### larger shift population size
plot(c(0,81),c(0,carryCap + 200)/1000,type="n",ylab=expression(paste(italic(""*N*"")," (thousands)",sep="")),xlab="",cex.lab=1.5,
     cex.main=1.5,cex.axis=0.9)
for(i in 1:nrow(mam_noMajorPopSize_biggerShift)){
  lines(1:81,mam_noMajorPopSize_biggerShift[i,]/1000,col=alpha("orange",alpha=0.05))
  lines(1:81,mam_majorPopSize_biggerShift[i,]/1000,col=alpha("darkblue",alpha=0.05))
}
polySizeMean <- rep(NA,81)
oligSizeMean <- rep(NA,81)
for(i in 1:81){
  polySizeMean[i] <- mean(mam_noMajorPopSize_biggerShift[,i])/1000
  oligSizeMean[i] <- mean(mam_majorPopSize_biggerShift[,i])/1000
}
lines(1:81,polySizeMean,lwd=3,col="orange")
lines(1:81,oligSizeMean,lwd=3,col="darkblue")
lines(c(0,200),c(carryCap/1000,carryCap/1000),lty="dashed")
par(xpd=TRUE)
text(x=-25,y=1.4,labels="B",cex=2)
par(xpd=FALSE)

#---------------------------------
##### plot the phenotypes
#---------------------------------
par(mar=c(3,5,2,0.5))
idealPhen <- 105

####### mammals
plot(c(0,80),c(100,110),type="n",ylab="Mean Phenotype",xlab="",cex.lab=1.5,cex.axis=0.8)
for(i in 1:500){
  lines(1:80,mam_noMajorPhen[i,],col=alpha("orange",alpha=0.05))
  lines(1:80,mam_majorPhen[i,],col=alpha("darkblue",alpha=0.05))
}

# get rolling mean
polyPhenMean <- rep(NA,80)
oligPhenMean <- rep(NA,80)
for(i in 1:80){
  polyPhenMean[i] <- mean(mam_noMajorPhen[,i],na.rm=TRUE)
  oligPhenMean[i] <- mean(mam_majorPhen[,i],na.rm=TRUE)
}
lines(1:80,polyPhenMean,lwd=3,col="orange")
lines(1:80,oligPhenMean,lwd=3,col="darkblue")
lines(c(0,80),c(105,105),lty="dashed")

######## bigger shift
par(mar=c(3,5,2,0.5))
idealPhen <- 107.5


plot(c(0,80),c(100,110),type="n",ylab="Mean Phenotype",xlab="",cex.lab=1.5,cex.axis=0.8)
for(i in 1:500){
  lines(1:80,mam_noMajorPhen_biggerShift[i,],col=alpha("orange",alpha=0.05))
  lines(1:80,mam_majorPhen_biggerShift[i,],col=alpha("darkblue",alpha=0.05))
}

# get rolling mean
polyPhenMean <- rep(NA,80)
oligPhenMean <- rep(NA,80)
for(i in 1:80){
  polyPhenMean[i] <- mean(mam_noMajorPhen_biggerShift[,i],na.rm=TRUE)
  oligPhenMean[i] <- mean(mam_majorPhen_biggerShift[,i],na.rm=TRUE)
}
lines(1:80,polyPhenMean,lwd=3,col="orange")
lines(1:80,oligPhenMean,lwd=3,col="darkblue")
lines(c(0,80),c(107.5,107.5),lty="dashed")

#####################################################################
# plot the fraction of extinct population versus time
#####################################################################
par(mar=c(4,5,1,0.5))
# extinct proportion for mammals
plot(c(1,81),c(0,0.12),type="n",xlab="Generation",ylab="Proportion Extinct",cex.lab=1.5)
lines(1:81,mam_polyExt/nrow(lowPolyDat),col="orange",lwd=3)
lines(1:81,mam_OligExt/nrow(lowOligDat),col="darkblue",lwd=3)
polygon(x=c(1:81,rev(1:81)),y=c(mamlowPolyCI[1,1:81],rev(mamlowPolyCI[2,1:81])),col=adjustcolor("orange",alpha.f=0.5),border=NA)
polygon(x=c(1:81,rev(1:81)),y=c(mamlowOligCI[1,1:81],rev(mamlowOligCI[2,1:81])),col=adjustcolor("darkblue",alpha.f=0.5),border=NA)

# bigger shift in theta
plot(c(1,81),c(0,0.12),type="n",xlab="Generation",ylab="Proportion Extinct",cex.lab=1.5)
lines(1:81,mam_polyExt_biggerShift/nrow(lowPolyDat_biggerShift),col="orange",lwd=3)
lines(1:81,mam_OligExt_biggerShift/nrow(lowOligDat_biggerShift),col="darkblue",lwd=3)
polygon(x=c(1:81,rev(1:81)),y=c(mamlowPolyCI_biggerShift[1,1:81],rev(mamlowPolyCI_biggerShift[2,1:81])),col=adjustcolor("orange",alpha.f=0.5),border=NA)
polygon(x=c(1:81,rev(1:81)),y=c(mamlowOligCI_biggerShift[1,1:81],rev(mamlowOligCI_biggerShift[2,1:81])),col=adjustcolor("darkblue",alpha.f=0.5),border=NA)


########################################################
#-------------------------------------------------------
# Major locus simulations with mutation
# Figure S17
#-------------------------------------------------------
########################################################
propAnalyze <- 0.9
####------------------------------
#### major locus simulations first
####------------------------------

popFileName <- paste("populationSize_mutRate_7.35e-09_fitFunSd_6_",c("1_100","101_200","201_300","301_400","401_500","501_520"),sep="")
phenFileName <- paste("meanPhenotype_mutRate_7.35e-09_fitFunSd_6_",c("1_100","101_200","201_300","301_400","401_500","501_520"),sep="")
majorFreqNames <- paste("majorAlleleFrequenciesAcrossSimulations_mutRate_7.35e-09_fitFunSd_6_",c("1_100","101_200","201_300","301_400","401_500","501_520"),sep="")
h2FileName <-  paste("h2TimeSeriesAcrossSimulations_mutRate_7.35e-09_fitFunSd_6_",c("1_100","101_200","201_300","301_400","401_500","501_520"),sep="")
mam_majorPopSize <- NULL
mam_majorPhen <- NULL
mam_majorH2 <- NULL
mam_majorFreqs <- NULL
for(i in 1:length(popFileName)){
  setwd(paste("~/Documents/work/rapidAdaptation/AmNatRevision/simulations/mutationSims/majorLocus/reps",i,sep=""))

  thisMajorPopSize <- as.matrix(read.table(popFileName[i],header=FALSE))
  thisMajorPhen <- as.matrix(read.table(phenFileName[i],header=FALSE))
  thisMajorH2 <- as.matrix(read.table(h2FileName[i],header=FALSE))
  thisMajorFreqs <- as.matrix(read.table(majorFreqNames[i],header=FALSE))
  print(nrow(thisMajorPopSize))

  mam_majorPopSize <- rbind(mam_majorPopSize,thisMajorPopSize)
  mam_majorPhen <- rbind(mam_majorPhen,thisMajorPhen)
  mam_majorH2 <- rbind(mam_majorH2,thisMajorH2)
  mam_majorFreqs <- rbind(mam_majorFreqs,thisMajorFreqs)
}


####-------------------------------------
#### now simulations with no major locus
####-------------------------------------

nameExts <- paste("_",c(1,101,201,301,401,501,601,701,801,901),"_",c(100,200,300,400,500,600,700,800,900,1000),sep="")
mam_noMajorPopSize <- NULL
mam_noMajorPhen <- NULL
mam_noMajorH2 <- NULL

# compile the five batches of simulations
for(i in 1:10){
  setwd(paste("/Users/marty_kardos/Documents/work/rapidAdaptation/AmNatRevision/simulations/mutationSims/polygenic/reps",i,sep=""))
  thisPop <- read.table(paste("populationSize_effectSize_0_optimPhen_10_mutRate_7.35e-08_fitFunSd_6",nameExts[i],sep=""),header=FALSE)
  thisPhen <- read.table(paste("meanPhenotype_effectSize_0_optimPhen_10_mutRate_7.35e-08_fitFunSd_6",nameExts[i],sep=""),header=FALSE)
  thisH2 <- read.table(paste("h2TimeSeriesAcrossSimulations_effectSize_0_optimPhen_10_mutRate_7.35e-08_fitFunSd_6",nameExts[i],sep=""),header=FALSE)
  print(nrow(thisPop))

  mam_noMajorPopSize <- rbind(mam_noMajorPopSize,as.matrix(thisPop))
  mam_noMajorPhen <- rbind(mam_noMajorPhen,as.matrix(thisPhen))
  mam_noMajorH2 <- rbind(mam_noMajorH2,as.matrix(thisH2))
}

#-----------------------------------------------
# make NA population sizes zero
#-----------------------------------------------


for(j in 1:ncol(mam_noMajorPopSize)){
  mam_noMajorPopSize[which(is.na(mam_noMajorPopSize[,j])),j] <- 0
  mam_majorPopSize[which(is.na(mam_majorPopSize[,j])),j] <- 0
}



#############################################
# keep only polygenic simulations
# with h2 very close to 0.6 at generation
# 1200
#############################################
keepers <- which(mam_noMajorH2[,1200] > 0.55 & mam_noMajorH2[,1200] < 0.65)[1:500]
mam_noMajorPopSize <- mam_noMajorPopSize[keepers,]
mam_noMajorPhen <- mam_noMajorPhen[keepers,]
mam_noMajorH2 <- mam_noMajorH2[keepers,]


##################################################
# identify the beginning of hard selection
# in the major locus simulations
##################################################
hardGens <- rep(NA,500)
for(i in 1:nrow(mam_majorPopSize)){
  hardGens [i] <- which(mam_majorPopSize[i,] != 500) [1] -1
}

#######################
# check for extinction
#######################


# make a new population size Matrix for major locus simulations

majorPopSize2 <- matrix(NA,nrow=500,ncol=81)
majorFreqs2 <- matrix(NA,nrow=500,ncol=81)
for(i in 1:500){
  majorPopSize2[i,] <- mam_majorPopSize[i,hardGens[i]:(hardGens[i] + 80)]
  majorFreqs2[i,] <- mam_majorFreqs[i,hardGens[i]:(hardGens[i] + 80)]
}

# get extinction rates per gereration after the onse of hard selection
mam_polyExt <- rep(NA,81)
mam_OligExt <- rep(NA,81)

gensAnalyze <- 1200:1280
for(i in 1:81){
  mam_polyExt[i] <- sum(mam_noMajorPopSize[,gensAnalyze[i]] == 0)
  mam_OligExt[i] <- sum(majorPopSize2[,i] == 0)
}
#############################################
# get bootstrap CIs for proportion extinct
#############################################

bootReps <- 1000
lowPolyBoots <- matrix(NA,nrow=bootReps,ncol=81)
lowOligBoots <- matrix(NA,nrow=bootReps,ncol=81)

for(i in 1:bootReps){
  lowPolyDat <- mam_noMajorPopSize[sample(1:nrow(mam_noMajorPopSize),nrow(mam_noMajorPopSize),replace=TRUE),gensAnalyze]
  lowOligDat <- majorPopSize2[sample(1:nrow(majorPopSize2),nrow(majorPopSize2),replace=TRUE),]
  lowPolyBoots [i,] <- colSums(lowPolyDat == 0)/nrow(lowPolyDat )
  lowOligBoots [i,] <- colSums(lowOligDat == 0)/nrow(lowOligDat )
  print(i)
}

mamlowPolyCI <- NULL
mamlowOligCI <- NULL

for(i in 1:81){
  mamlowPolyCI <- cbind(mamlowPolyCI,quantile(lowPolyBoots[,i],c(0.025,0.975)))
  mamlowOligCI <- cbind(mamlowOligCI,quantile(lowOligBoots[,i],c(0.025,0.975)))
}

############################################
# plot results
############################################
par(mfrow=c(3,1),mar=c(2,5,3,0.5),xpd=FALSE)
library(scales)
carryCap <- 1000
##### mammal population size
plot(c(0,81),c(0,carryCap + 200)/1000,type="n",ylab=expression(paste(italic(""*N*"")," (thousands)",sep="")),xlab="",cex.lab=1.5,
     cex.main=1.5,cex.axis=0.9)
for(i in 1:nrow(mam_noMajorPopSize)){
  lines(1:81,mam_noMajorPopSize[i,1200:1280]/1000,col=alpha("darkred",alpha=0.05))
  lines(1:81,majorPopSize2[i,]/1000,col=alpha("darkgray",alpha=0.2))
}
polySizeMean <- rep(NA,81)
oligSizeMean <- rep(NA,81)
for(i in 1:81){
  polySizeMean[i] <- mean(mam_noMajorPopSize[,i+1199])/1000
  oligSizeMean[i] <- mean(majorPopSize2[,i])/1000
}
lines(1:81,polySizeMean,lwd=3,col="darkred")
lines(1:81,oligSizeMean,lwd=3,col="darkgray",lty="dashed")
lines(c(0,200),c(carryCap/1000,carryCap/1000),lty="dashed")

#---------------------------------
##### plot the phenotypes
#---------------------------------
par(mar=c(3,5,2,0.5))
idealPhen <- 110

# organize the phenotypes from simulations with a major locus
mam_majorPhen2 <- matrix(NA,nrow=500,ncol=81)
for(i in 1:500){
  mam_majorPhen2[i,] <- mam_majorPhen[i,hardGens[i]:(hardGens[i]+80)]
}

####### mammals
plot(c(1,81),c(0,13),type="n",ylab="Mean Phenotype",xlab="",cex.lab=1.5,cex.axis=0.8)
for(i in 1:500){
  lines(1:81,mam_noMajorPhen[i,1200:1280],col=alpha("darkred",alpha=0.02))
  lines(1:81,mam_majorPhen2[i,]-mam_majorPhen2[i,1],col=alpha("darkgray",alpha=0.1))
}

# get rolling mean
polyPhenMean <- rep(NA,80)
oligPhenMean <- rep(NA,80)
for(i in 1:80){
  polyPhenMean[i] <- mean(mam_noMajorPhen[,i+1199],na.rm=TRUE)
  oligPhenMean[i] <- mean(mam_majorPhen2[,i]-mam_majorPhen2[,1],na.rm=TRUE)
}
lines(1:80,polyPhenMean,lwd=3,col=alpha("darkred"))
lines(1:80,oligPhenMean,lwd=3,col=alpha("darkgray",alpha=0.8),lty="dashed")

lines(c(0,80),c(10,10),lty="dashed")


par(mar=c(3,5,2,0.5))
idealPhen <- 112

#####################################################################
# plot the fraction of extinct population versus time
#####################################################################

par(mar=c(4,5,1,0.5))
# extinct proportion
plot(c(1,81),c(0,0.5),type="n",xlab="Generation",ylab="Proportion Extinct",cex.lab=1.5)
lines(1:81,mam_polyExt/nrow(lowPolyDat),col="darkred",lwd=3)
lines(1:81,mam_OligExt/nrow(lowOligDat),col="darkgray",lty="dashed",lwd=3)

polygon(x=c(1:81,rev(1:81)),y=c(mamlowPolyCI[1,],rev(mamlowPolyCI[2,])),col=adjustcolor("darkred",alpha.f=0.5),border=NA)
polygon(x=c(1:81,rev(1:81)),y=c(mamlowOligCI[1,],rev(mamlowOligCI[2,])),col=adjustcolor("darkgray",alpha.f=0.5),border=NA)
# get pictures

legend(x=30,y=0.1,xjust=FALSE,yjust=FALSE,legend=c("No Major Locus","Major Locus"),
       lwd=3,col=c("darkred","darkgray"),bty="n")




#################################################################################
#--------------------------------------------------------------------------------
#################################################################################
# Make figures of evolutionary potential versus time for major and polygenic
# simulations of the scenarios in Figure S11
#################################################################################
#--------------------------------------------------------------------------------
#################################################################################
majorMeanL0 <- NULL # save the mean L0 for major and polygenic scenarios
polyMeanL0 <- NULL


pValVec_major <- NULL       # save the p values and regression coefficients for persistence versus selection limit
coefVec_major <- NULL        
pValVec_poly <- NULL       # save the p values and regression coefficients for persistence versus selection limit
coefVec_poly <- NULL        

# save all of the major locus and polygenic evolutionary limits and populations sizes

allMajorLims <- NULL
allPolyLims <- NULL
allMajorNs <- NULL
allPolyNs <- NULL



keepPolyLims <- NULL # keep polygenic limits for plotting
keepMajorLims <- NULL # keep polygenic limits for plotting


#####################################
#------------------------------------
############plot with effect size 0.9
############ Figure S19
#------------------------------------
#####################################
setwd("/Users/marty_kardos/Documents/work/rapidAdaptation/AmNatRevision/simulations/originalSims/h2Pt6/fecundity4/major")

majorLims <- read.table("evolutionarLimits_N0_500_K_1000_propMaj_0.9",header=TRUE,colClasses="numeric")
keepMajorLims <- rbind(keepMajorLims,majorLims)
majorNs <- read.table("popSize_N0_500_K_1000_propMaj_0.9",header=TRUE,colClasses="numeric")

polyLims <- read.table("evolutionarLimits_N0_500_K_1000_propMaj_0.01",header=TRUE,colClasses="numeric")
keepPolyLims <- rbind(keepPolyLims,polyLims)
polyNs <- read.table("popSize_N0_500_K_1000_propMaj_0.01",header=TRUE,colClasses="numeric")


################ major locus

par(mfrow=c(2,2),mar=c(4,6,2,2),xpd=FALSE)
library(scales)
plot(c(1,30),c(0,15),type="n",ylab=expression(italic(""*L*"")),xlab="Generation",cex.lab=1.3)
for(i in 1:500){
  lines(1:80,majorLims[i,],col=alpha("blue",alpha=0.05))
}
lines(1:80,colMeans(majorLims,na.rm=TRUE),col="blue",lwd=3)
par(xpd=TRUE)
text(x=-8,y=16.5,labels="A",cex=2)
par(xpd=FALSE)

majorFinalN <- majorNs[,81]
majorFinalN[is.na(majorFinalN)] <- 0
majorPersist <- majorFinalN > 0
plot(majorLims[,1],majorPersist,col=alpha("blue",alpha=0.2),xlab=expression(italic(""*L*"")[0]),ylab="Persistence",cex.lab=1.3)

par(xpd=TRUE)
text(x=0.5,y=1.1,labels="B",cex=2)
par(xpd=FALSE)

majorMeanL0 <- c(majorMeanL0,mean(majorLims[,1]))
#logistic model
limits <- majorLims[,1]
logMod <- glm(formula = majorPersist ~ limits, family = binomial(logit))

pValVec_major <- c(pValVec_major,summary(logMod)$coefficients[2,4])
coefVec_major <- c(coefVec_major,summary(logMod)$coefficients[2,1])


yPreds <- predict(logMod,type="response")
plotX <- limits[order(limits)]
plotY <- yPreds[order(limits)]
lines(plotX,plotY, col = "black", lwd = 2)

allMajorLims <- c(allMajorLims,limits)
allMajorNs <- c(allMajorNs,majorPersist)

############### polygenic

plot(c(1,30),c(0,15),type="n",ylab=expression(italic(""*L*"")),xlab="Generation",cex.lab=1.3)
for(i in 1:500){
  lines(1:80,polyLims[i,],col=alpha("orange",alpha=0.05))
}
lines(1:80,colMeans(polyLims,na.rm=TRUE),col="orange",lwd=3)

polyFinalN <- polyNs[,81]
polyFinalN[is.na(polyFinalN)] <- 0
polyPersist <- polyFinalN > 0
plot(polyLims[,1],polyPersist,col=alpha("orange",alpha=0.2),xlab=expression(italic(""*L*"")[0]),ylab="Persistence",cex.lab=1.3)

pVal <- summary(logMod)$coefficients[2,4]
coef1 <- summary(logMod)$coefficients[2,1]
coef2 <- summary(logMod)$coefficients[1,1]

polyMeanL0 <- c(polyMeanL0,mean(polyLims[,1]))
#logistic model
limits <- polyLims[,1]
logMod <- glm(formula = polyPersist ~ limits, family = binomial(logit))
yPreds <- predict(logMod,type="response")
plotX <- limits[order(limits)]
plotY <- yPreds[order(limits)]
lines(plotX,plotY, col = "black", lwd = 2)


allPolyLims <- c(allPolyLims,limits)
allPolyNs <- c(allPolyNs,polyPersist)
pValVec_poly <- c(pValVec_poly,summary(logMod)$coefficients[2,4])
coefVec_poly <- c(coefVec_poly,summary(logMod)$coefficients[2,1])



#####################################
#------------------------------------
############plot with effect size 0.7
############ Figure S20
#------------------------------------
#####################################
setwd("/Users/marty_kardos/Documents/work/rapidAdaptation/AmNatRevision/simulations/originalSims/h2Pt6/fecundity4/major")

majorLims <- read.table("evolutionarLimits_N0_500_K_1000_propMaj_0.7",header=TRUE,colClasses="numeric")
keepMajorLims <- rbind(keepMajorLims,majorLims)
majorNs <- read.table("popSize_N0_500_K_1000_propMaj_0.7",header=TRUE,colClasses="numeric")

polyLims <- read.table("evolutionarLimits_N0_500_K_1000_propMaj_0.01",header=TRUE,colClasses="numeric")
keepPolyLims <- rbind(keepPolyLims,polyLims)
polyNs <- read.table("popSize_N0_500_K_1000_propMaj_0.01",header=TRUE,colClasses="numeric")


# major locus
par(mfrow=c(2,2),mar=c(4,6,2,2),xpd=FALSE)
library(scales)
plot(c(1,30),c(0,15),type="n",ylab=expression(italic(""*L*"")),xlab="Generation",cex.lab=1.3)
for(i in 1:500){
  lines(1:80,majorLims[i,],col=alpha("blue",alpha=0.05))
}
lines(1:80,colMeans(majorLims,na.rm=TRUE),col="blue",lwd=3)
par(xpd=TRUE)
text(x=-8,y=16.5,labels="A",cex=2)
par(xpd=FALSE)

majorFinalN <- majorNs[,81]
majorFinalN[is.na(majorFinalN)] <- 0
majorPersist <- majorFinalN > 0
plot(majorLims[,1],majorPersist,col=alpha("blue",alpha=0.2),xlab=expression(italic(""*L*"")[0]),ylab="Persistence",cex.lab=1.3)

par(xpd=TRUE)
text(x=0.5,y=1.1,labels="B",cex=2)
par(xpd=FALSE)

majorMeanL0 <- c(majorMeanL0,mean(majorLims[,1]))
#logistic model
limits <- majorLims[,1]
logMod <- glm(formula = majorPersist ~ limits, family = binomial(logit))
yPreds <- predict(logMod,type="response")
plotX <- limits[order(limits)]
plotY <- yPreds[order(limits)]
lines(plotX,plotY, col = "black", lwd = 2)
pValVec_major <- c(pValVec_major,summary(logMod)$coefficients[2,4])
coefVec_major <- c(coefVec_major,summary(logMod)$coefficients[2,1])


allMajorLims <- c(allMajorLims,limits)
allMajorNs <- c(allMajorNs,majorPersist)

# polygenic
plot(c(1,30),c(0,15),type="n",ylab=expression(italic(""*L*"")),xlab="Generation",cex.lab=1.3)
for(i in 1:500){
  lines(1:80,polyLims[i,],col=alpha("orange",alpha=0.05))
}
lines(1:80,colMeans(polyLims,na.rm=TRUE),col="orange",lwd=3)
polyFinalN <- polyNs[,81]
polyFinalN[is.na(polyFinalN)] <- 0
polyPersist <- polyFinalN > 0
plot(polyLims[,1],polyPersist,col=alpha("orange",alpha=0.2),xlab=expression(italic(""*L*"")[0]),ylab="Persistence",cex.lab=1.3)

pVal <- summary(logMod)$coefficients[2,4]
coef1 <- summary(logMod)$coefficients[2,1]
coef2 <- summary(logMod)$coefficients[1,1]

polyMeanL0 <- c(polyMeanL0,mean(polyLims[,1]))
#logistic model
limits <- polyLims[,1]
logMod <- glm(formula = polyPersist ~ limits, family = binomial(logit))
yPreds <- predict(logMod,type="response")
plotX <- limits[order(limits)]
plotY <- yPreds[order(limits)]
lines(plotX,plotY, col = "black", lwd = 2)


#allPolyLims <- c(allPolyLims,limits)
#allPolyNs <- c(allPolyNs,polyPersist)
pValVec_poly <- c(pValVec_poly,summary(logMod)$coefficients[2,4])
coefVec_poly <- c(coefVec_poly,summary(logMod)$coefficients[2,1])




#####################################
#------------------------------------
############plot with effect size 0.5
############ Figure S21
#------------------------------------
#####################################
setwd("/Users/marty_kardos/Documents/work/rapidAdaptation/AmNatRevision/simulations/originalSims/h2Pt6/fecundity4/major")

majorLims <- read.table("evolutionarLimits_N0_500_K_1000_propMaj_0.5",header=TRUE,colClasses="numeric")
keepMajorLims <- rbind(keepMajorLims,majorLims)
majorNs <- read.table("popSize_N0_500_K_1000_propMaj_0.5",header=TRUE,colClasses="numeric")

polyLims <- read.table("evolutionarLimits_N0_500_K_1000_propMaj_0.01",header=TRUE,colClasses="numeric")
keepPolyLims <- rbind(keepPolyLims,polyLims)
polyNs <- read.table("popSize_N0_500_K_1000_propMaj_0.01",header=TRUE,colClasses="numeric")


# major locus
par(mfrow=c(2,2),mar=c(4,6,2,2),xpd=FALSE)
library(scales)
plot(c(1,30),c(0,15),type="n",ylab=expression(italic(""*L*"")),xlab="Generation",cex.lab=1.3)
for(i in 1:500){
  lines(1:80,majorLims[i,],col=alpha("blue",alpha=0.05))
}
lines(1:80,colMeans(majorLims,na.rm=TRUE),col="blue",lwd=3)
par(xpd=TRUE)
text(x=-8,y=16.5,labels="A",cex=2)
par(xpd=FALSE)

majorFinalN <- majorNs[,81]
majorFinalN[is.na(majorFinalN)] <- 0
majorPersist <- majorFinalN > 0
plot(majorLims[,1],majorPersist,col=alpha("blue",alpha=0.2),xlab=expression(italic(""*L*"")[0]),ylab="Persistence",cex.lab=1.3)

par(xpd=TRUE)
text(x=0.5,y=1.1,labels="B",cex=2)
par(xpd=FALSE)

majorMeanL0 <- c(majorMeanL0,mean(majorLims[,1]))
#logistic model
limits <- majorLims[,1]
logMod <- glm(formula = majorPersist ~ limits, family = binomial(logit))
yPreds <- predict(logMod,type="response")
plotX <- limits[order(limits)]
plotY <- yPreds[order(limits)]
lines(plotX,plotY, col = "black", lwd = 2)


allMajorLims <- c(allMajorLims,limits)
allMajorNs <- c(allMajorNs,majorPersist)
pValVec_major <- c(pValVec_major,summary(logMod)$coefficients[2,4])
coefVec_major <- c(coefVec_major,summary(logMod)$coefficients[2,1])


# polygenic
plot(c(1,30),c(0,15),type="n",ylab=expression(italic(""*L*"")),xlab="Generation",cex.lab=1.3)
for(i in 1:500){
  lines(1:80,polyLims[i,],col=alpha("orange",alpha=0.05))
}
lines(1:80,colMeans(polyLims,na.rm=TRUE),col="orange",lwd=3)

polyFinalN <- polyNs[,81]
polyFinalN[is.na(polyFinalN)] <- 0
polyPersist <- polyFinalN > 0
plot(polyLims[,1],polyPersist,col=alpha("orange",alpha=0.2),xlab=expression(italic(""*L*"")[0]),ylab="Persistence",cex.lab=1.3)

pVal <- summary(logMod)$coefficients[2,4]
coef1 <- summary(logMod)$coefficients[2,1]
coef2 <- summary(logMod)$coefficients[1,1]

polyMeanL0 <- c(polyMeanL0,mean(polyLims[,1]))
#logistic model
limits <- polyLims[,1]
logMod <- glm(formula = polyPersist ~ limits, family = binomial(logit))
yPreds <- predict(logMod,type="response")
plotX <- limits[order(limits)]
plotY <- yPreds[order(limits)]
lines(plotX,plotY, col = "black", lwd = 2)


#allPolyLims <- c(allPolyLims,limits)
#allPolyNs <- c(allPolyNs,polyPersist)
pValVec_poly <- c(pValVec_poly,summary(logMod)$coefficients[2,4])
coefVec_poly <- c(coefVec_poly,summary(logMod)$coefficients[2,1])

#######################################################
#------------------------------------------------------
############plot with 20 generations to the new theta
############ Figure S22
#------------------------------------------------------
#######################################################
setwd("~/Documents/work/rapidAdaptation/AmNatRevision/simulations/originalSims/h2Pt6/fecundity4/movingPhenoOpt")

majorLims <- read.table("evolutionarLimits_N0_500_K_1000_propMaj_0.9_gensToNewOpt_20",header=TRUE,colClasses="numeric")
keepMajorLims <- rbind(keepMajorLims,majorLims)
majorNs <- read.table("popSize_N0_500_K_1000_propMaj_0.9_gensToNewOpt_20",header=TRUE,colClasses="numeric")

polyLims <- read.table("evolutionarLimits_N0_500_K_1000_propMaj_0.01_gensToNewOpt_20",header=TRUE,colClasses="numeric")
keepPolyLims <- rbind(keepPolyLims,polyLims)
polyNs <- read.table("popSize_N0_500_K_1000_propMaj_0.01_gensToNewOpt_20",header=TRUE,colClasses="numeric")


################ major locus

par(mfrow=c(2,2),mar=c(4,6,2,2),xpd=FALSE)
library(scales)
plot(c(1,30),c(0,15),type="n",ylab=expression(italic(""*L*"")),xlab="Generation",cex.lab=1.3)
for(i in 1:500){
  lines(1:80,majorLims[i,],col=alpha("blue",alpha=0.05))
}
lines(1:80,colMeans(majorLims,na.rm=TRUE),col="blue",lwd=3)
par(xpd=TRUE)
text(x=-8,y=16.5,labels="A",cex=2)
par(xpd=FALSE)

majorFinalN <- majorNs[,81]
majorFinalN[is.na(majorFinalN)] <- 0
majorPersist <- majorFinalN > 0
plot(majorLims[,1],majorPersist,col=alpha("blue",alpha=0.2),xlab=expression(italic(""*L*"")[0]),ylab="Persistence",cex.lab=1.3)

par(xpd=TRUE)
text(x=0.5,y=1.1,labels="B",cex=2)
par(xpd=FALSE)
majorMeanL0 <- c(majorMeanL0,mean(majorLims[,1]))

#logistic model
limits <- majorLims[,1]
logMod <- glm(formula = majorPersist ~ limits, family = binomial(logit))
yPreds <- predict(logMod,type="response")
plotX <- limits[order(limits)]
plotY <- yPreds[order(limits)]
lines(plotX,plotY, col = "black", lwd = 2)
allMajorLims <- c(allMajorLims,limits)
allMajorNs <- c(allMajorNs,majorPersist)
pValVec_major <- c(pValVec_major,summary(logMod)$coefficients[2,4])
coefVec_major <- c(coefVec_major,summary(logMod)$coefficients[2,1])


############### polygenic

plot(c(1,30),c(0,15),type="n",ylab=expression(italic(""*L*"")),xlab="Generation",cex.lab=1.3)
for(i in 1:500){
  lines(1:80,polyLims[i,],col=alpha("orange",alpha=0.05))
}
lines(1:80,colMeans(polyLims,na.rm=TRUE),col="orange",lwd=3)

polyFinalN <- polyNs[,81]
polyFinalN[is.na(polyFinalN)] <- 0
polyPersist <- polyFinalN > 0
plot(polyLims[,1],polyPersist,col=alpha("orange",alpha=0.2),xlab=expression(italic(""*L*"")[0]),ylab="Persistence",cex.lab=1.3)

pVal <- summary(logMod)$coefficients[2,4]
coef1 <- summary(logMod)$coefficients[2,1]
coef2 <- summary(logMod)$coefficients[1,1]

polyMeanL0 <- c(polyMeanL0,mean(polyLims[,1]))
#logistic model
limits <- polyLims[,1]
logMod <- glm(formula = polyPersist ~ limits, family = binomial(logit))
yPreds <- predict(logMod,type="response")
plotX <- limits[order(limits)]
plotY <- yPreds[order(limits)]
lines(plotX,plotY, col = "black", lwd = 2)
allPolyLims <- c(allPolyLims,limits)
allPolyNs <- c(allPolyNs,polyPersist)
pValVec_poly <- c(pValVec_poly,summary(logMod)$coefficients[2,4])
coefVec_poly <- c(coefVec_poly,summary(logMod)$coefficients[2,1])



#######################################################
#------------------------------------------------------
############plot with 10 generations to the new theta
############ Figure S23
#------------------------------------------------------
#######################################################
setwd("~/Documents/work/rapidAdaptation/AmNatRevision/simulations/originalSims/h2Pt6/fecundity4/movingPhenoOpt")

majorLims <- read.table("evolutionarLimits_N0_500_K_1000_propMaj_0.9_gensToNewOpt_10",header=TRUE,colClasses="numeric")
keepMajorLims <- rbind(keepMajorLims,majorLims)
majorNs <- read.table("popSize_N0_500_K_1000_propMaj_0.9_gensToNewOpt_10",header=TRUE,colClasses="numeric")

polyLims <- read.table("evolutionarLimits_N0_500_K_1000_propMaj_0.01_gensToNewOpt_10",header=TRUE,colClasses="numeric")
keepPolyLims <- rbind(keepPolyLims,polyLims)
polyNs <- read.table("popSize_N0_500_K_1000_propMaj_0.01_gensToNewOpt_10",header=TRUE,colClasses="numeric")


################ major locus

par(mfrow=c(2,2),mar=c(4,6,2,2),xpd=FALSE)
library(scales)
plot(c(1,30),c(0,15),type="n",ylab=expression(italic(""*L*"")),xlab="Generation",cex.lab=1.3)
for(i in 1:500){
  lines(1:80,majorLims[i,],col=alpha("blue",alpha=0.05))
}
lines(1:80,colMeans(majorLims,na.rm=TRUE),col="blue",lwd=3)
par(xpd=TRUE)
text(x=-8,y=16.5,labels="A",cex=2)
par(xpd=FALSE)

majorFinalN <- majorNs[,81]
majorFinalN[is.na(majorFinalN)] <- 0
majorPersist <- majorFinalN > 0
plot(majorLims[,1],majorPersist,col=alpha("blue",alpha=0.2),xlab=expression(italic(""*L*"")[0]),ylab="Persistence",cex.lab=1.3)

par(xpd=TRUE)
text(x=0.5,y=1.1,labels="B",cex=2)
par(xpd=FALSE)

majorMeanL0 <- c(majorMeanL0,mean(majorLims[,1]))
#logistic model
limits <- majorLims[,1]
logMod <- glm(formula = majorPersist ~ limits, family = binomial(logit))
yPreds <- predict(logMod,type="response")
plotX <- limits[order(limits)]
plotY <- yPreds[order(limits)]
lines(plotX,plotY, col = "black", lwd = 2)

allMajorLims <- c(allMajorLims,limits)
allMajorNs <- c(allMajorNs,majorPersist)
pValVec_major <- c(pValVec_major,summary(logMod)$coefficients[2,4])
coefVec_major <- c(coefVec_major,summary(logMod)$coefficients[2,1])


############### polygenic

plot(c(1,30),c(0,15),type="n",ylab=expression(italic(""*L*"")),xlab="Generation",cex.lab=1.3)
for(i in 1:500){
  lines(1:80,polyLims[i,],col=alpha("orange",alpha=0.05))
}
lines(1:80,colMeans(polyLims,na.rm=TRUE),col="orange",lwd=3)

polyFinalN <- polyNs[,81]
polyFinalN[is.na(polyFinalN)] <- 0
polyPersist <- polyFinalN > 0
plot(polyLims[,1],polyPersist,col=alpha("orange",alpha=0.2),xlab=expression(italic(""*L*"")[0]),ylab="Persistence",cex.lab=1.3)

pVal <- summary(logMod)$coefficients[2,4]
coef1 <- summary(logMod)$coefficients[2,1]
coef2 <- summary(logMod)$coefficients[1,1]

polyMeanL0 <- c(polyMeanL0,mean(polyLims[,1]))
#logistic model
limits <- polyLims[,1]
logMod <- glm(formula = polyPersist ~ limits, family = binomial(logit))
yPreds <- predict(logMod,type="response")
plotX <- limits[order(limits)]
plotY <- yPreds[order(limits)]
lines(plotX,plotY, col = "black", lwd = 2)


allPolyLims <- c(allPolyLims,limits)
allPolyNs <- c(allPolyNs,polyPersist)
pValVec_poly <- c(pValVec_poly,summary(logMod)$coefficients[2,4])
coefVec_poly <- c(coefVec_poly,summary(logMod)$coefficients[2,1])


#######################################################
#------------------------------------------------------
############ plasticity parameter m = 0.1
############ Figure S24
#------------------------------------------------------

setwd("~/Documents/work/rapidAdaptation/AmNatRevision/simulations/originalSims/h2Pt6/fecundity4/plasticity")

majorLims <- read.table("evolutionarLimits_N0_500_K_1000_propMaj_0.9_gensToNewOpt_10_plast_0.1",header=TRUE,colClasses="numeric")
keepMajorLims <- rbind(keepMajorLims,majorLims)
majorNs <- read.table("popSize_N0_500_K_1000_propMaj_0.9_gensToNewOpt_10_plast_0.1",header=TRUE,colClasses="numeric")

polyLims <- read.table("evolutionarLimits_N0_500_K_1000_propMaj_0.01_gensToNewOpt_10_plast_0.1",header=TRUE,colClasses="numeric")
keepPolyLims <- rbind(keepPolyLims,polyLims)
polyNs <- read.table("popSize_N0_500_K_1000_propMaj_0.01_gensToNewOpt_10_plast_0.1",header=TRUE,colClasses="numeric")


################ major locus

par(mfrow=c(2,2),mar=c(4,6,2,2),xpd=FALSE)
library(scales)
plot(c(1,30),c(0,15),type="n",ylab=expression(italic(""*L*"")),xlab="Generation",cex.lab=1.3)
for(i in 1:500){
  lines(1:80,majorLims[i,],col=alpha("blue",alpha=0.05))
}
lines(1:80,colMeans(majorLims,na.rm=TRUE),col="blue",lwd=3)
par(xpd=TRUE)
text(x=-8,y=16.5,labels="A",cex=2)
par(xpd=FALSE)

majorFinalN <- majorNs[,81]
majorFinalN[is.na(majorFinalN)] <- 0
majorPersist <- majorFinalN > 0
plot(majorLims[,1],majorPersist,col=alpha("blue",alpha=0.2),xlab=expression(italic(""*L*"")[0]),ylab="Persistence",cex.lab=1.3)

par(xpd=TRUE)
text(x=0.5,y=1.1,labels="B",cex=2)
par(xpd=FALSE)

majorMeanL0 <- c(majorMeanL0,mean(majorLims[,1]))
#logistic model
limits <- majorLims[,1]
logMod <- glm(formula = majorPersist ~ limits, family = binomial(logit))
yPreds <- predict(logMod,type="response")
plotX <- limits[order(limits)]
plotY <- yPreds[order(limits)]
lines(plotX,plotY, col = "black", lwd = 2)

allMajorLims <- c(allMajorLims,limits)
allMajorNs <- c(allMajorNs,majorPersist)
pValVec_major <- c(pValVec_major,summary(logMod)$coefficients[2,4])
coefVec_major <- c(coefVec_major,summary(logMod)$coefficients[2,1])


############### polygenic

plot(c(1,30),c(0,15),type="n",ylab=expression(italic(""*L*"")),xlab="Generation",cex.lab=1.3)
for(i in 1:500){
  lines(1:80,polyLims[i,],col=alpha("orange",alpha=0.05))
}
lines(1:80,colMeans(polyLims,na.rm=TRUE),col="orange",lwd=3)

polyFinalN <- polyNs[,81]
polyFinalN[is.na(polyFinalN)] <- 0
polyPersist <- polyFinalN > 0
plot(polyLims[,1],polyPersist,col=alpha("orange",alpha=0.2),xlab=expression(italic(""*L*"")[0]),ylab="Persistence",cex.lab=1.3)

pVal <- summary(logMod)$coefficients[2,4]
coef1 <- summary(logMod)$coefficients[2,1]
coef2 <- summary(logMod)$coefficients[1,1]

polyMeanL0 <- c(polyMeanL0,mean(polyLims[,1]))
#logistic model
limits <- polyLims[,1]
logMod <- glm(formula = polyPersist ~ limits, family = binomial(logit))
yPreds <- predict(logMod,type="response")
plotX <- limits[order(limits)]
plotY <- yPreds[order(limits)]
lines(plotX,plotY, col = "black", lwd = 2)

allPolyLims <- c(allPolyLims,limits)
allPolyNs <- c(allPolyNs,polyPersist)
pValVec_poly <- c(pValVec_poly,summary(logMod)$coefficients[2,4])
coefVec_poly <- c(coefVec_poly,summary(logMod)$coefficients[2,1])

#######################################################
#------------------------------------------------------
############ plasticity parameter m = 0.2
############ Figure S25
#------------------------------------------------------
#######################################################
setwd("~/Documents/work/rapidAdaptation/AmNatRevision/simulations/originalSims/h2Pt6/fecundity4/plasticity")

majorLims <- read.table("evolutionarLimits_N0_500_K_1000_propMaj_0.9_gensToNewOpt_10_plast_0.2",header=TRUE,colClasses="numeric")
keepMajorLims <- rbind(keepMajorLims,majorLims)
majorNs <- read.table("popSize_N0_500_K_1000_propMaj_0.9_gensToNewOpt_10_plast_0.2",header=TRUE,colClasses="numeric")

polyLims <- read.table("evolutionarLimits_N0_500_K_1000_propMaj_0.01_gensToNewOpt_10_plast_0.2",header=TRUE,colClasses="numeric")
keepPolyLims <- rbind(keepPolyLims,polyLims)
polyNs <- read.table("popSize_N0_500_K_1000_propMaj_0.01_gensToNewOpt_10_plast_0.2",header=TRUE,colClasses="numeric")


################ major locus

par(mfrow=c(2,2),mar=c(4,6,2,2),xpd=FALSE)
library(scales)
plot(c(1,30),c(0,15),type="n",ylab=expression(italic(""*L*"")),xlab="Generation",cex.lab=1.3)
for(i in 1:500){
  lines(1:80,majorLims[i,],col=alpha("blue",alpha=0.05))
}
lines(1:80,colMeans(majorLims,na.rm=TRUE),col="blue",lwd=3)
par(xpd=TRUE)
text(x=-8,y=16.5,labels="A",cex=2)
par(xpd=FALSE)

majorFinalN <- majorNs[,81]
majorFinalN[is.na(majorFinalN)] <- 0
majorPersist <- majorFinalN > 0
plot(majorLims[,1],majorPersist,col=alpha("blue",alpha=0.2),xlab=expression(italic(""*L*"")[0]),ylab="Persistence",cex.lab=1.3)

par(xpd=TRUE)
text(x=0.5,y=1.1,labels="B",cex=2)
par(xpd=FALSE)

majorMeanL0 <- c(majorMeanL0,mean(majorLims[,1]))
#logistic model
limits <- majorLims[,1]
logMod <- glm(formula = majorPersist ~ limits, family = binomial(logit))
yPreds <- predict(logMod,type="response")
plotX <- limits[order(limits)]
plotY <- yPreds[order(limits)]
lines(plotX,plotY, col = "black", lwd = 2)
pValVec_major <- c(pValVec_major,summary(logMod)$coefficients[2,4])
coefVec_major <- c(coefVec_major,summary(logMod)$coefficients[2,1])


allMajorLims <- c(allMajorLims,limits)
allMajorNs <- c(allMajorNs,majorPersist)


############### polygenic

plot(c(1,30),c(0,15),type="n",ylab=expression(italic(""*L*"")),xlab="Generation",cex.lab=1.3)
for(i in 1:500){
  lines(1:80,polyLims[i,],col=alpha("orange",alpha=0.05))
}
lines(1:80,colMeans(polyLims,na.rm=TRUE),col="orange",lwd=3)

polyFinalN <- polyNs[,81]
polyFinalN[is.na(polyFinalN)] <- 0
polyPersist <- polyFinalN > 0
plot(polyLims[,1],polyPersist,col=alpha("orange",alpha=0.2),xlab=expression(italic(""*L*"")[0]),ylab="Persistence",cex.lab=1.3)

pVal <- summary(logMod)$coefficients[2,4]
coef1 <- summary(logMod)$coefficients[2,1]
coef2 <- summary(logMod)$coefficients[1,1]

polyMeanL0 <- c(polyMeanL0,mean(polyLims[,1]))
#logistic model
limits <- polyLims[,1]
logMod <- glm(formula = polyPersist ~ limits, family = binomial(logit))
yPreds <- predict(logMod,type="response")
plotX <- limits[order(limits)]
plotY <- yPreds[order(limits)]
lines(plotX,plotY, col = "black", lwd = 2)

pValVec_poly <- c(pValVec_poly,summary(logMod)$coefficients[2,4])
coefVec_poly <- c(coefVec_poly,summary(logMod)$coefficients[2,1])

allPolyLims <- c(allPolyLims,limits)
allPolyNs <- c(allPolyNs,polyPersist)


#######################################################
#------------------------------------------------------
############ plasticity parameter m = 0.4
############ Figure S26
#------------------------------------------------------
#######################################################
setwd("~/Documents/work/rapidAdaptation/AmNatRevision/simulations/originalSims/h2Pt6/fecundity4/plasticity")

majorLims <- read.table("evolutionarLimits_N0_500_K_1000_propMaj_0.9_gensToNewOpt_10_plast_0.4",header=TRUE,colClasses="numeric")
keepMajorLims <- rbind(keepMajorLims,majorLims)
majorNs <- read.table("popSize_N0_500_K_1000_propMaj_0.9_gensToNewOpt_10_plast_0.4",header=TRUE,colClasses="numeric")

polyLims <- read.table("evolutionarLimits_N0_500_K_1000_propMaj_0.01_gensToNewOpt_10_plast_0.4",header=TRUE,colClasses="numeric")
keepPolyLims <- rbind(keepPolyLims,polyLims)
polyNs <- read.table("popSize_N0_500_K_1000_propMaj_0.01_gensToNewOpt_10_plast_0.4",header=TRUE,colClasses="numeric")


################ major locus

par(mfrow=c(2,2),mar=c(4,6,2,2),xpd=FALSE)
library(scales)
plot(c(1,30),c(0,15),type="n",ylab=expression(italic(""*L*"")),xlab="Generation",cex.lab=1.3)
for(i in 1:500){
  lines(1:80,majorLims[i,],col=alpha("blue",alpha=0.05))
}
lines(1:80,colMeans(majorLims,na.rm=TRUE),col="blue",lwd=3)
par(xpd=TRUE)
text(x=-8,y=16.5,labels="A",cex=2)
par(xpd=FALSE)

majorFinalN <- majorNs[,81]
majorFinalN[is.na(majorFinalN)] <- 0
majorPersist <- majorFinalN > 0
plot(majorLims[,1],majorPersist,col=alpha("blue",alpha=0.2),xlab=expression(italic(""*L*"")[0]),ylab="Persistence",cex.lab=1.3)

par(xpd=TRUE)
text(x=0.5,y=1.1,labels="B",cex=2)
par(xpd=FALSE)

majorMeanL0 <- c(majorMeanL0,mean(majorLims[,1]))
#logistic model
limits <- majorLims[,1]
logMod <- glm(formula = majorPersist ~ limits, family = binomial(logit))
yPreds <- predict(logMod,type="response")
plotX <- limits[order(limits)]
plotY <- yPreds[order(limits)]
lines(plotX,plotY, col = "black", lwd = 2)
allMajorLims <- c(allMajorLims,limits)
allMajorNs <- c(allMajorNs,majorPersist)
pValVec_major <- c(pValVec_major,summary(logMod)$coefficients[2,4])
coefVec_major <- c(coefVec_major,summary(logMod)$coefficients[2,1])


############### polygenic

plot(c(1,30),c(0,15),type="n",ylab=expression(italic(""*L*"")),xlab="Generation",cex.lab=1.3)
for(i in 1:500){
  lines(1:80,polyLims[i,],col=alpha("orange",alpha=0.05))
}
lines(1:80,colMeans(polyLims,na.rm=TRUE),col="orange",lwd=3)

polyFinalN <- polyNs[,81]
polyFinalN[is.na(polyFinalN)] <- 0
polyPersist <- polyFinalN > 0
plot(polyLims[,1],polyPersist,col=alpha("orange",alpha=0.2),xlab=expression(italic(""*L*"")[0]),ylab="Persistence",cex.lab=1.3)

pVal <- summary(logMod)$coefficients[2,4]
coef1 <- summary(logMod)$coefficients[2,1]
coef2 <- summary(logMod)$coefficients[1,1]

polyMeanL0 <- c(polyMeanL0,mean(polyLims[,1]))
#logistic model
limits <- polyLims[,1]
logMod <- glm(formula = polyPersist ~ limits, family = binomial(logit))
yPreds <- predict(logMod,type="response")
plotX <- limits[order(limits)]
plotY <- yPreds[order(limits)]
lines(plotX,plotY, col = "black", lwd = 2)
pValVec_poly <- c(pValVec_poly,summary(logMod)$coefficients[2,4])
coefVec_poly <- c(coefVec_poly,summary(logMod)$coefficients[2,1])


poly_oddsRat <- exp(coefVec_poly)
major_oddsRat <- exp(coefVec_major)

poly_logResMat <- cbind(round(pValVec_poly,digits=3),round(coefVec_poly,digits=3),round(poly_oddsRat,digits=3))
major_logResMat <- cbind(round(pValVec_major,digits=3),round(coefVec_major,digits=3),round(major_oddsRat,digits=3))
colnames(poly_logResMat) <- c("pVal","coef","oddsRat")
colnames(major_logResMat) <- c("pVal","coef","oddsRat")


##########################################################
#---------------------------------------------------------
# plot odds ratio from polygenic versus from major locus
# Figure S27
#---------------------------------------------------------
##########################################################

par(mfrow=c(1,1))
plot(major_logResMat[,3],poly_logResMat[,3],xlim=c(0.8,3),ylim=c(0.8,3),
     xlab="Odds ratio (major locus)",ylab="Odds ratio (polygenic)",pch=16,col="darkgray",cex=1.5)
lines(c(0,4),c(0,4),lty="dashed")


################################################################
#---------------------------------------------------------------
# logistic model across all simulations
# Fiogure S28
#---------------------------------------------------------------
################################################################


allLims <- c(allMajorLims,allPolyLims)
allNs <- c(allMajorNs,allPolyNs)

colVec <- rep(NA,length(allLims))
colVec[1:length(allMajorLims)] <- "blue"
colVec[(length(allMajorLims)+1):length(allLims)] <- "orange"
par(mfrow=c(1,1))
plot(allLims,allNs,xlab=expression(italic(""*L*"")[0]),ylab="Persistence",col=alpha("darkgray",alpha=0.02))
#logistic model

logMod <- glm(formula = allNs ~ allLims, family = binomial(logit))
yPreds <- predict(logMod,type="response")
plotX <- allLims[order(allLims)]
plotY <- yPreds[order(allLims)]
lines(plotX,plotY, col = "black", lwd = 2)


setwd("~/Documents/work/rapidAdaptation/AmNatRevision/simulations/newFigs")
write.table(poly_logResMat,file="logisticRegResults_polygenic_26May2020.txt",quote=FALSE)
write.table(major_logResMat,file="logisticRegResults_major_26May2020.txt",quote=FALSE)






##################################################################
#-----------------------------------------------------------------
############ all polygenic and major locus sims plotted separately
############ Figure 6
#-----------------------------------------------------------------
##################################################################




################ major locus

par(mfrow=c(2,2),mar=c(4,6,2,2),xpd=FALSE)
library(scales)
plot(c(1,30),c(0,15),type="n",ylab=expression(italic(""*L*"")),xlab="Generation",cex.lab=1.3)
for(i in 1:nrow(keepMajorLims)){
  lines(1:80,keepMajorLims[i,],col=alpha("blue",alpha=0.015))
}
lines(1:80,colMeans(keepMajorLims,na.rm=TRUE),col="blue",lwd=3)
par(xpd=TRUE)
text(x=-8,y=16.5,labels="A",cex=2)
par(xpd=FALSE)


plot(allMajorLims,allMajorNs,col=alpha("blue",alpha=0.05),xlab=expression(italic(""*L*"")[0]),ylab="Persistence",cex.lab=1.3)

par(xpd=TRUE)
text(x=-3,y=1.1,labels="B",cex=2)
par(xpd=FALSE)


#logistic model
logMod <- glm(formula = allMajorNs ~ allMajorLims, family = binomial(logit))
summary(logMod)
yPreds <- predict(logMod,type="response")
plotX <- allMajorLims[order(allMajorLims)]
plotY <- yPreds[order(allMajorLims)]
lines(plotX,plotY, col = "black", lwd = 2)



################ polygenic


library(scales)
plot(c(1,30),c(0,15),type="n",ylab=expression(italic(""*L*"")),xlab="Generation",cex.lab=1.3)
for(i in 1:nrow(keepPolyLims)){
  lines(1:80,keepPolyLims[i,],col=alpha("orange",alpha=0.01))
}
lines(1:80,colMeans(keepPolyLims,na.rm=TRUE),col="orange",lwd=3)


plot(allPolyLims,allPolyNs,col=alpha("orange",alpha=0.05),xlab=expression(italic(""*L*"")[0]),ylab="Persistence",cex.lab=1.3)
#logistic model
logMod <- glm(formula = allPolyNs ~ allPolyLims, family = binomial(logit))
summary(logMod)
yPreds <- predict(logMod,type="response")
plotX <- allPolyLims[order(allPolyLims)]
plotY <- yPreds[order(allPolyLims)]
lines(plotX,plotY, col = "black", lwd = 2)





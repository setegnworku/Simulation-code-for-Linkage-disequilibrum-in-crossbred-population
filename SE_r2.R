
###############################################################################
# Function to estimate the SE(r2) from n_rep replicates
# reads a scenario in terms of nrep,Nab,Ncd,pA,pB,Pc,PD,D,rab,rcd
# makes input matrix M with nrep replications of the input parameters
# Uses the function simCross_r2 which gets r2 for a single replicate
###############################################################################

SE_r2 <- function(input) {
  nrep <- as.numeric(input[1]);
  Nab <-as.numeric(input[2]); Ncd <- as.numeric(input[3]);
  pA  <- as.numeric(input[4]); pB  <- as.numeric(input[5]);
  pC  <- as.numeric(input[6]); pD  <- as.numeric(input[7]);
 rab <- as.numeric(input[8]); rcd <- as.numeric(input[9]) #r<- input[10]
  M<-cbind(rep(nrep,nrep),rep(Nab,nrep),rep(Ncd,nrep),
           rep(pA,nrep),rep(pB,nrep),
           rep(pC,nrep),rep(pD,nrep),
           rep(rab,nrep),rep(rcd,nrep))
  
  results<-apply(M,1,simCross_r2)
  res<- as.data.frame(results)
  #note: this results in a matrix with 6 rows instead of columns.
  mu<-c(mean(results[1,]),mean(results[2,]),mean(results[3,]))
  sigma<-c(sd(results[1,]),sd(results[2,]),sd(results[3,]))
  return(c(mu,sigma))

}
## Example for one r2 vale with sample size of 900.
r=0.2
N_alt=100
pAvals<- rep( seq(0.05,0.95,0.1),each=10)
pBvals<- rep( seq(0.05,0.95,0.1),each=10)

rab<-rcd<- rep(r,100)
pCvals<-rep( seq(0.05,0.95,0.10),10)
pDvals<- rep( seq(0.05,0.95,0.10),10)

datpar<- data.frame(pAvals,pBvals,pCvals,pDvals,rab,rcd)


datfilter=datpar[datpar$pAvals<=datpar$pCvals, ]
nrep<- rep(2,nrow(datfilter))
NAB<- NCD<- rep( N_alt,nrow(datfilter))
datwithrep<- cbind(nrep,NAB,NCD,datfilter)
library(dplyr)
dat=datwithrep
M<- dat[,1:9]## The first 9 columns are the input parameters
M<- as.matrix(M)

results<-apply(M,1,SE_r2)### Apply r
results<-t(results)# transpose the result
colnames(results)<- c('mr_2_True','mr_2_hap','mr_2_geno','sdr_2_True','sdr_2_hap','sdr_2_geno')###
results<-cbind(dat,results)## scenarios with parameter estimate with standard error.

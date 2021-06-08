

simCross_r2 <- function(input) {
  Nrep<- input[1]
  Nab <- as.numeric(input[2]); Ncd <-as.numeric(input[3]); # Nab size of  line1 and Ncd size of line 2
  pA  <- input[4]; pB  <- input[5];
  pC  <- input[6]; pD  <- input[7];
  rab <- input[8]; rcd <- input[9]

  
  ##breed 1
  pa <- 1-pA
  pb <- 1-pB
  ##breed 2
  pc <- 1-pC
  pd <- 1-pD
  Dab=sqrt(rab*pA*pa*pB*pb)
  Dcd=sqrt(rcd*pC*pc*pD*pd)
  
  r12 <- Dab^2/(pA*pa*pB*pb)
  ##breed 2 r2
  r22 <- Dcd^2/(pC*pD*pc*pd)
  
  rcom<- ((Dab+Dcd)^2)/((pa*pA+pC*pc)*(pB*pb+pD*pd)) ### The combined r2for cross breed (the true r squared)
  ##breed 1
  ## Once we have the LD level and the allele freqeuncy  we computed the haplotype frequency using the equations fab=pa*pb+Dab, we 
  ## sampled the two haplotype from the multinominal distribution for the required population size. This applies for the
  ## two lines
  
  ### Haplotype frequency for line AB
  fAB = pA*pB+Dab
  fAb = pA*pb-Dab
  faB = pa*pB-Dab
  fab = pa*pb+Dab
  
  ##Haplotype freqeuncy for line CD
  fCD = pC*pD+Dcd
  fCd = pC*pd-Dcd
  fcD = pc*pD-Dcd
  fcd = pc*pd+Dcd
  
  
  ## Sample the haplotype for line AB 
  yab <- rmultinom(Nab,2,c(fAB,fAb,faB,fab))  #sample two haplotypes for each individual from line AB
  xab<-t(yab) #dim n*4 is the haplotypes
  #------------------------

  x <- xab
  N <- Nab
  
  #------------------------
  h<-matrix(nrow=4,ncol=N)  #rows and colums this way because R fills by column
  ## The homozygous at haplotype level 
  h[,x[,1]==2]<-c(1,1,1,1)        #make halpotype dataset; first haplotype second halpotype, count A and B as one.
  h[,x[,2]==2]<-c(1,0,1,0)  #Ab Ab
  h[,x[,3]==2]<-c(0,1,0,1)
  h[,x[,4]==2]<-c(0,0,0,0)
  
  ##  Hetrozygous at haplotype levels
  #--- the four columns in x are simulated based on the order of fAB,fAb,faB,fab, so are h[,x[,col]==2]
  h[,x[,1]==1&x[,2]==1]<-c(1,1,1,0) #AB Ab
  h[,x[,1]==1&x[,3]==1]<-c(1,1,0,1)
  h[,x[,1]==1&x[,4]==1]<-c(1,1,0,0)
  h[,x[,2]==1&x[,3]==1]<-c(1,0,0,1)
  h[,x[,2]==1&x[,4]==1]<-c(1,0,0,0)
  h[,x[,3]==1&x[,4]==1]<-c(0,1,0,0) #aB ab
  h<-t(h)
  hab <- h[1:Nab,]## 
  
  ycd <- rmultinom(N,2,c(fCD,fCd,fcD,fcd))  #sample two haplotypes for each individual for  line CD
  xcd<-t(ycd)#dim 100*4 is the haplotypes
  #------------------------

  x <- xcd
  N <- as.numeric(Ncd)## NOte Nab Ncd are the same so it doesnot matter 
  #------------------------
  h<-matrix(nrow=4,ncol=N)  #rows and columns this way because R fills by column
  h[,x[,1]==2]<-c(1,1,1,1)        #make halpotype dataset; 
  h[,x[,2]==2]<-c(1,0,1,0)  #Cd Cd
  h[,x[,3]==2]<-c(0,1,0,1)
  h[,x[,4]==2]<-c(0,0,0,0)
  
  #--- the four columns in x are simulated based on the order of fCD,fCd,fcD,fcd so are h[,x[,col]==2]
  h[,x[,1]==1&x[,2]==1]<-c(1,1,1,0) #CD Cd
  h[,x[,1]==1&x[,3]==1]<-c(1,1,0,1)
  h[,x[,1]==1&x[,4]==1]<-c(1,1,0,0)
  h[,x[,2]==1&x[,3]==1]<-c(1,0,0,1)
  h[,x[,2]==1&x[,4]==1]<-c(1,0,0,0)
  h[,x[,3]==1&x[,4]==1]<-c(0,1,0,0) #cD cd
  h<-t(h)
  ## Nab eqauls Ncd 
  hcd <- h[1:Nab,]
  
  
  ## Here sample the haplotype for the cross breed ( sample from line AB and line CD for each of the haplotype)
  h <- matrix(NA,ncol=4,nrow=Nab) 
  ## Note I aribitrary use line  AB or ab , and for the second breed too line CD or cd 
  habcdc<-habp<-hcdp<- matrix(NA,ncol=4,nrow=Nab)## habcdc just the cross breed from line ab and cd, habp  store parental haplotype from breed ab, hcdp  store parental haplotype from breed cd
  habp[,1]<-paste0(hab[,1],hab[,2])# First hapolotype for the first line AB
  habp[,2]<-paste0(hab[,3],hab[,4])## second haplotype for the first line  AB
  
  hcdp[,1]<-paste0(hcd[,1],hcd[,2])# First haplotype for the second line CD
  hcdp[,2]<-paste0(hcd[,3],hcd[,4])# Second haplotype for the second line  CD
  
  habcdc[,1] <- apply(habp[,c(1,2)],1,function(x)sample(x,size=1)) # x- user-defined function; 1 means rowwise; sample one  haplotype from breed AB
  habcdc[,2] <- apply(hcdp[,c(1,2)],1,function(x)sample(x,size=1))# sample the other haplotype from breed CD
  m1<- apply(habcdc[,c(1,2)],1,function(x)substring(x,1,1)) # break the haplotype to form the allele
  m2<- apply(habcdc[,c(1,2)],1,function(x)substring(x,2,2))
  
  h[1:Nab,1]<- as.numeric(m1[1,])#  first locus first allele
  h[1:Nab,2]<- as.numeric(m1[2,])### second locus first allele
  h[1:Nab,3]<- as.numeric(m2[1,]) ## first locus second allele
  h[1:Nab,4]<- as.numeric(m2[2,])## secondlocus second allele
  
  g<-matrix(nrow=N,ncol=2)  #genotypes
  g[,1]=h[,1]+h[,2] #g is genotype column1 is the first allele at locus 1 and coulum3 is the second allele at locus 1
  g[,2]=h[,3]+h[,4] # column2 is the first allele at locus 2 and coulum4 is the second allele at locus 2
  cov1<-cov(h[,1],h[,3])
  
  cov2<-cov(h[,2],h[,4])
  
  cov_avg<-(cov1+cov2)/2
  
  var1<-var(h[,1])
  
  var2<-var(h[,2])
  
  var3<-var(h[,3])
  
  var4<-var(h[,4])
  
  varL1<-(var1+var2)/2  ## Variance of line 1
  
  varL2<-(var3+var4)/2 ## variance of line 2
  
  r2_haplo<-cov_avg^2/(varL1*varL2)
  
  
  #-----------------------
  r_geno<-cor(g[,1],g[,2])
  r2_geno<-r_geno^2
  #------------------------
  r2_geno
  #------------------------
  R2<-c( rcom,r2_haplo,r2_geno)## True r2, r2 using haplotype data and r2 using genotype data
}
#------------------------
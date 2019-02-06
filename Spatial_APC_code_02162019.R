###########################################################################
#           by: Pavel Chernyavskiy; last modified 01312019                # 
#                                                                         #
#                                                                         #
###########################################################################

od<-read.csv("Opioid_5yr_15_64_2002_2016.csv") #left censored = -1
table(od$Age_grp, od$Year_grp) #quick check

#must sort the data by State then year then age
dat.s<-od[order(od$State,od$year,od$age),]

#*** APC Design matrix function ***#
#For more info see: Rosenberg et al (2014), PMID: 25146089
designmatrix <- function(PVP){
  N <- nrow(PVP$D$DATA)
  J <- matrix(1, nrow = N)
  a <- matrix(PVP$D$DATA[,1])
  avals <- PVP$D$a
  aM <- length(avals)
  p <- matrix(PVP$D$DATA[,2])
  pvals <- PVP$D$p
  pN <- length(pvals)
  c <- matrix(PVP$D$DATA[,3])
  cvals <- PVP$D$c
  cK <- length(cvals)
  
  # Age 
  Ad <- matrix(NaN,nrow=N,ncol=aM)
  for (i in 1:aM)
    Ad[,i] <- a==avals[i]
  # Per
  Pd <- matrix(NaN,nrow=N,ncol=pN)
  for (i in 1:pN)
    Pd[,i] <- p==pvals[i]
  # Coh
  Cd <- matrix(NaN,nrow=N,ncol=cK)
  for (i in 1:cK)
    Cd[,i] <- c==cvals[i]
  abar <- mean(PVP$D$a)
  pbar <- mean(PVP$D$p)
  cbar <- mean(PVP$D$c)
  a0 <- a - abar
  p0 <- p - pbar
  c0 <- c - cbar
  Xa <- cbind(J, a0)
  Ra <- solve(t(Xa)%*%Xa,t(Xa))
  XAD <- diag(N) - Xa%*%Ra
  Xp <- cbind(J, p0)
  Rp <- solve(t(Xp)%*%Xp,t(Xp))
  XPD <- diag(N) - Xp%*%Rp
  Xc <- cbind(J, c0)
  Rc <- solve(t(Xc)%*%Xc,t(Xc))
  XCD <- diag(N) - Xc%*%Rc
  Ad0 <- XAD%*%Ad
  Pd0 <- XPD%*%Pd
  Cd0 <- XCD%*%Cd
  X <- cbind(J, a0, c0, Ad0[,2:(aM-1)], Pd0[,2:(pN-1)], Cd0[,2:(cK-1)])
  D <- list(X = X)
}
#build APC design matrix by state
X.by.geog<-function(use){
  XMAT<-cbind(use$age,use$year,use$cohort,use$deaths,use$Pop)
  colnames(XMAT)<-c("Age","Period","Cohort","Events","Offset")
  rownames(XMAT)<-NULL
  PVP<-list(D=list(DATA=XMAT, a=unique(XMAT[,1]),p=unique(XMAT[,2]),c=unique(XMAT[,3])))
  mat.list<-designmatrix(PVP)
  rm(use,XMAT,PVP)
  return(mat.list$X)
}
library(foreach)
#could be done in parallel (using %dopar%) if number of regions is large
ugeog<-unique(dat.s$State)
#foreach stacks the APC design matrix for each state
XM<-foreach(i=1:length(ugeog),.combine='rbind') %do% X.by.geog(subset(dat.s,dat.s$State==ugeog[i])) 

###################################################
#partition X matrix for easy post-processing
#Would need to be adjusted for a different dataset
#Column 1 contains the Intercept
A<-XM[,2] #centered Age
C<-XM[,3] #centered Cohort
A_dev<-XM[,4:11] #10-2 = 8 age curvature parms (10 total age grps)
colnames(A_dev)<-paste("A_dev",seq(1,ncol(A_dev)),sep="")
P_dev<-XM[,12] #3-2 = 1 period curvature parm (3 period grps)
C_dev<-XM[,13:ncol(XM)] #12-2 = 10 cohort curv parms (12 cohort grps)
colnames(C_dev)<-paste("C_dev",seq(1,ncol(C_dev)),sep="")


##################################################################################################
#*** Fit standard model with independent random effects (Chernyavskiy et al, 2017) using brms ***#
library(brms)
library(rstan)
rstan_options(auto_write=TRUE)
#define weakly-informative priors
my_priors<-c(
  set_prior("normal(0,2.5)", class = "b"),    
  set_prior("normal(0,1)",class="sd"), 
  set_prior("lkj(2)", class = "cor"))
#make a dataset to use in brms
brms.dat<-cbind(dat.s,A,C,A_dev,P_dev,C_dev)

###########################################################################
#specify formula (standard model with left-censoring)
f1<-bf(deaths|cens(censor) ~ A + C + (A + C|State) + offset(log(Pop)) +  
         A_dev1 + A_dev2 + A_dev3 + A_dev4 + A_dev5 + A_dev6 + A_dev7 + A_dev8 + 
         P_dev + 
         C_dev1 + C_dev2 + C_dev3 + C_dev4 + C_dev5 + C_dev6 + C_dev7 + C_dev8 + C_dev9 + C_dev10)

#run standard model
bm1<-brm(formula=f1, data=brms.dat, family=negbinomial(),  
         prior = my_priors,
         control = list(adapt_delta = 0.8, max_treedepth=15), seed=10111985,
         cores=3, iter=5000,warmup=1000, chains=3, thin=10, init_r=0.3) 
bm1 #<-view output


##################################################
#*** produce adjacency matrix for spatial APC ***#
#library(devtools)
#devtools::install_github("hrbrmstr/albersusa")
library(albersusa)
library(spdep)
library(sp)
us_poly<-usa_composite("laea")
#MUST SORT us_poly SAME WAY AS THE DATA ARE SORTED#
us_poly<-us_poly[order(us_poly$name),] 
nb.usa<-poly2nb(us_poly,queen=TRUE,row.names = us_poly$name)
Adj.usa<-nb2mat(nb.usa,style='B',zero.policy = TRUE)
colnames(Adj.usa)<-rownames(Adj.usa)
Dn<-diag(rowSums(Adj.usa)) #matrix with number of neighbors diagonal#

#extract formatted standata from brms model (bm1)
standat<-standata(bm1)
#append other matrices into standat
standat$W<-as.matrix(Adj.usa)
standat$Q<-as.matrix(Dn - Adj.usa)
standat$In<-as.matrix(diag(nrow(Dn)))
str(standat) #<-check 

###############################
#fit spatial APC model in Stan#
#change model_code argument for place first (spatial_APC_place) and people first (spatial_APC_ppl)#
sm1_5yr<-stan(model_code = spatial_APC_place, 
              model_name = "Opioid NegBin Spatial APC",
              data = standat, 
              control = list(adapt_delta = 0.8, max_treedepth=15), seed=10111985, 
              pars = c("b_Intercept", "b","SD", "rho", "eta",       
                       "psi","log_lik","mu_fit","shape",
                       "region_RR","region_LAT","region_drift"),          
              cores=3, iter=5000, warmup = 1000, chains=3, thin = 10, init_r=0.3) 
#init_r of 0.5, 0.75 also works
sm1_5yr #<-view output




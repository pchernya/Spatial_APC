################################################################################################################################
#*** Spatialy-Independent random effects***#
library(mvtnorm)
library(matrixcalc)
library(lme4)
NSIMS<-50 #set number of simulated datasets

#fit a base glmer model for its matrices
m0<-glmer.nb(deaths ~  A + C + A_dev + P_dev + C_dev + (1 + A + C|State) + offset(Pop),
          data=brms.dat,family="poisson",nAGQ=1,
          control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))
summary(m0)

#beta vector similar to final spatial model
beta<-c(-9.36,0.08,0.07,0.98,1.24,1.30,1.31,1.24,1.09,0.80,0.44,-0.01,-0.33,-0.75,
        -1.15,-1.35,-1.47,-1.68,-2.13,-2.76,-3.46,-3.66)

#Age and Coh random slopes + correlated
#no other strong correlations
Rmat<-matrix(c(1,-0.1,0.1,
               -0.1,1,0.7,
                0.1,0.7,1),byrow=TRUE,nrow=3,ncol=3)
Sigma<-diag(c(0.2,0.03,0.02))%*%Rmat%*%diag(c(0.2,0.03,0.02))

#draw a set of random effects
set.seed(10111985)
u<-rmvnorm(n=51, sigma=Sigma, method="chol")
#Get the matrices & offset
X0<-as.matrix(getME(m0,"X"))
Z0<-t(as.matrix(getME(m0,"Ztlist")[[1]]))
Za<-t(as.matrix(getME(m0,"Ztlist")[[2]]))
Zc<-t(as.matrix(getME(m0,"Ztlist")[[3]]))
O<-getME(m0,"offset") #already on log scale
#make mu
mu_ind<-exp(O + X0%*%beta + Z0%*%u[,1] + Za%*%u[,2] + Zc%*%u[,3])

#loop to generate counts using mu
Y_ind<-matrix(NA,nrow=length(O),ncol=NSIMS)
set.seed(101119851) 
for(i in 1:NSIMS){
  Y_ind[,i]<-MASS::rnegbin(mu_ind,theta=40) 
}
# pick off random effects if needed
#u0<-u[,1]
#ua<-u[,2]
#uc<-u[,3]

#write Y_ind
#saveRDS(Y_ind, "Simulation//Y_ind.rds")


###############################################################################################################################
#*** Spatially-correlated random effects***#
#use same beta vector
#use between state centroid distance
st_coords<-read.csv("Simulation//state_coords.csv")
st_Dist<-as.matrix(dist(st_coords[,-1]))

#realistically strong spatial correlation: eff range of 3/0.06 = 50, approx 90th pctl of distances
Sigma_0<-(0.2^2)*exp(-0.06*st_Dist)
Sigma_a<-(0.03^2)*exp(-0.06*st_Dist)
Sigma_c<-(0.02^2)*exp(-0.06*st_Dist)

#draw a set of random effects
#Age and Coh random slopes + correlated
#no other strong correlations
set.seed(101119854) 
u0<-rmvnorm(n=1, sigma=Sigma_0, method="chol")
ua<-rmvnorm(n=1, mean=-0.04*u0, sigma=Sigma_a, method="chol")
uc<-rmvnorm(n=1, mean= 0.04*u0 + 0.7*ua, sigma=Sigma_c, method="chol")
cor(cbind(t(u0),t(ua),t(uc)))

#make mu_sp
mu_sp<-exp(O + X0%*%beta + Z0%*%t(u0) + Za%*%t(ua) + Zc%*%t(uc))

#loop to generate counts using mu_sp
Y_sp<-matrix(NA,nrow=length(O),ncol=NSIMS)
set.seed(101119851)
for(i in 1:NSIMS){
  Y_sp[,i]<-MASS::rnegbin(mu_sp,theta=40) #real shape estimate was 37.32 
}

#write Y_sp
#saveRDS(Y_sp, "Simulation//Y_sp.rds")

#############################################################################################################################
#############################################################################################################################
#*** SIM CODE ***# 
library(albersusa)
library(spdep)
library(sp)
library(rstan)
library(brms)
rstan_options(auto_write = TRUE)
#load in standata template#
#extract standata from base brms model#
bm1<-readRDS("Application models//5-yr models//Op_NB_APC_5yr.rds")
standat<-standata(bm1) #switch out Y during loop

us_poly<-usa_composite("laea")
#SORT SAME WAY AS THE DATA ARE SORTED#
us_poly<-us_poly[order(us_poly$name),] 
nb.usa<-poly2nb(us_poly,queen=TRUE,row.names = us_poly$name)
Adj.usa<-nb2mat(nb.usa,style='B',zero.policy = TRUE)
colnames(Adj.usa)<-rownames(Adj.usa)
Dn<-diag(rowSums(Adj.usa)) #diagonal matrix with number of neighbors on diagonal#

#append spatial weights matrix - Dn
standat$W<-as.matrix(Adj.usa)
standat$Q<-as.matrix(Dn - Adj.usa)
standat$In<-as.matrix(diag(nrow(Dn)))

#Spatial APC loop#
sim_func_sp<-function(Y.sim,control.parms,iter,warmup,thin){
  
  #apply censoring
  Y.use<-Y.sim
  cens.use<-c(rep(0,standat$N))
  cens.use[which(Y.use < 10)] <- -1
  Y.use[which(Y.use < 10)]<- 10
  #replace simulated data and cens vectors in standata
  standat$cens <- cens.use
  standat$Y <- Y.use

  #fit the spatial APC model
  sp.model<-stan(model_code = spat_APC_txt_9_ppl, model_name = "NegBin Spatial APC",
                 data = standat, control = control.parms, seed=10111985, 
                 pars = c("b_Intercept", "b","SD", "rho", "eta", "psi","log_lik","shape",
                          "region_RR","region_LAT","region_drift"),
                 cores=3, chains=3, iter=iter, warmup = warmup, thin = thin, init_r=0.2)
  
  #pre-process estimates (parms should have ~1200 samples)
  temp<-rstan::extract(sp.model)
  #hmc_diag<-rstan::check_hmc_diagnostics(sp.model)
  #str(temp)
  
  #specify parameters to record
  waic.mod<-loo::waic(extract_log_lik(sp.model))$estimates[3,1] #compute the WAIC
  shape<-mean(temp$shape)
  Int<-mean(temp$b_Intercept)
  LAT<-mean(temp$b[,1])
  NetDrift<-mean(temp$b[,2])
  SD.0<-mean(temp$SD[,1])
  SD.a<-mean(temp$SD[,2])
  SD.c<-mean(temp$SD[,3])
  rho.0<-mean(temp$rho[,1])
  rho.a<-mean(temp$rho[,2])
  rho.c<-mean(temp$rho[,3])
  ###
  eta.a0<-mean(temp$eta[,1]) #place-first specific
  eta.c0<-mean(temp$eta[,2]) #place-first specific
  eta.ac<-mean(temp$eta[,3]) #place-first specific
  psi.a0<-mean(temp$psi[,1]) #place-first specific
  psi.c0<-mean(temp$psi[,2]) #place-first specific
  psi.ac<-mean(temp$psi[,3]) #place-first specific
  ###
#  eta.a0<-mean(temp$eta[,2]) #people-first specific
#  eta.c0<-mean(temp$eta[,3]) #people-first specific
#  eta.ac<-mean(temp$eta[,1]) #people-first specific
#  psi.a0<-mean(temp$psi[,2]) #people-first specific
#  psi.c0<-mean(temp$psi[,3]) #people-first specific
#  psi.ac<-mean(temp$psi[,1]) #people-first specific
  ###
  region.RR<-colMeans(temp$region_RR)
  region.LAT<-colMeans(temp$region_LAT)
  region.drift<-colMeans(temp$region_drift)
  
  #make a list to output
  out.list<-list(waic=waic.mod,
                 pop.est=c(Int,LAT,NetDrift,SD.0,SD.a,SD.c,rho.0,rho.a,rho.c,
                           eta.a0,eta.c0,eta.ac,psi.a0,psi.c0,psi.ac,shape),
                 reg.est=matrix(cbind(region.RR,region.LAT,region.drift),nrow=51,ncol=3)
                 )
  return(out.list)
}

#*run the loop*#
library(foreach)
library(doParallel)
cl <- makeCluster(5) #<-can try up to 15 = 3 chains x 5 workers#
registerDoParallel(cl)
#stopCluster(cl) #can use to remove parallel cores when done#

#strong spatial corr - place first#
res.1<-foreach(i=1:NSIMS, .errorhandling=c('stop'), .export=c("standat"),
               .packages=c('rstan','brms','loo')
               ) %dopar% sim_func_sp(Y.sim=Y_sp[,i], 
                                     control.parms=list(adapt_delta = 0.8, max_treedepth=15),
                                     iter=800,warmup=400,thin=3)
saveRDS(res.1,"Simulation//strong_cor_sp_mod_place.rds")
#res.1<-readRDS("Simulation//strong_cor_sp_mod_place.rds")

#independent data - place first#
res.2<-foreach(i=1:NSIMS, .errorhandling=c('stop'), .export=c("standat"),
               .packages=c('rstan','brms','loo')
) %dopar% sim_func_sp(Y.sim=Y_ind[,i], 
                      control.parms=list(adapt_delta = 0.8, max_treedepth=15),
                      iter=800,warmup=400,thin=3)
saveRDS(res.2,"Simulation//indep_sp_mod_place.rds")
#res.2<-readRDS("Simulation//indep_sp_mod_place.rds")

#strong spatial corr - people first#
res.3<-foreach(i=1:NSIMS, .errorhandling=c('stop'), .export=c("standat"),
               .packages=c('rstan','brms','loo')
) %dopar% sim_func_sp(Y.sim=Y_sp[,i], 
                      control.parms=list(adapt_delta = 0.8, max_treedepth=15),
                      iter=800,warmup=400,thin=3)
saveRDS(res.3,"Simulation//strong_cor_sp_mod_people.rds")
#res.3<-readRDS("Simulation//strong_cor_sp_mod_people.rds")

#indepedent data - people first
res.4<-foreach(i=1:NSIMS, .errorhandling=c('stop'), .export=c("standat"),
               .packages=c('rstan','brms','loo')
) %dopar% sim_func_sp(Y.sim=Y_ind[,i], 
                      control.parms=list(adapt_delta = 0.8, max_treedepth=15),
                      iter=800,warmup=400,thin=3)
saveRDS(res.4,"Simulation//indep_sp_mod_people.rds")
#res.4<-readRDS("Simulation//indep_sp_mod_people.rds")



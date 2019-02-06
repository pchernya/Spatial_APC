###########################################################################
#              *** Stan code for spatial APC models ***#                  #
###########################################################################

spatial_APC_place<-"
// *** Standard model Stan code created using the brms package (Buerkner, 2017) *** //
// *** Spatial APC code by Pavel Chernyavskiy for the manuscript 'Spatially-varying age-period-cohort analysis with 
//     application to US mortality, 2002-2016', latest version 01/30/2019  *** //

data { 
int<lower=1> N;                 // total number of observations 
int Y[N];                       // response variable 
int<lower=1> K;                 // number of population-level effects 
matrix[N, K] X;                 // population-level design matrix 
vector[N] offset;               // has already been logged ie log(Pop) 
int<lower=-1,upper=2> cens[N];  // indicates left/right/no censoring 

// set up random effects dimensions & random effects vecs
int<lower=1> J_1[N];
int<lower=1> N_1; // total unique pop strata (ie geographies)
int<lower=1> M_1; // 3 ran eff vecs here (Int, LAT, drift)
// random effects design vectors (N x 1)
vector[N] Z_1_1;
vector[N] Z_1_2;
vector[N] Z_1_3;

// W is (0/1) spatial adjacency matrix
matrix[N_1,N_1] W;
// Dn is diag matrix with number of 1st order neighbors
// define Q = W-Dn and pass into Stan from R
matrix[N_1, N_1] Q;
// Identity matrix to be passed in from R
matrix[N_1, N_1] In;
}

//centers the X matrix: X -> Xc
//removes the Int: K -> Kc
transformed data {
int Kc = K - 1; 
matrix[N, K - 1] Xc;  
vector[K - 1] means_X;   
for (i in 2:K) { 
means_X[i - 1] = mean(X[, i]); 
Xc[, i - 1] = X[, i] - means_X[i - 1]; 
}} 

parameters { 
// pop-level fixed effects and temp Int
vector[Kc] b; 
real temp_Intercept;   

// ranef SDs
// SD[1]=int, SD[2]=LAT, SD[3]=drift
vector<lower=0>[M_1] SD;  

// spatial smoothing parameters
// rho[1]=int; rho[2]=LAT; rho[3]=drift
vector<lower=0, upper=1>[M_1] rho;

// spatial cross-corr parameters
vector[M_1] eta;
vector[M_1] psi;

// random effects vectors
vector[N_1] r_1_1; //random int
vector[N_1] r_1_2; //random LAT
vector[N_1] r_1_3; //random drift

// neg-binomial shape parameter kappa
real<lower=0> shape;  
}

transformed parameters{
vector<lower=0>[M_1] tau;
matrix[N_1,N_1] Ci_int;
matrix[N_1,N_1] Ci_LAT;
matrix[N_1,N_1] Ci_drift;
vector[N_1] mu_int;
vector[N_1] mu_LAT;
vector[N_1] mu_drift;

// compute inverse variances
for(i in 1:M_1){ tau[i] = 1/(SD[i]^2); }

// compute inverse-variance covariance matrices (precision matrices)
Ci_int =   tau[1]*(rho[1]*Q + (1-rho[1])*In);
Ci_LAT =   tau[2]*(rho[2]*Q + (1-rho[2])*In);
Ci_drift = tau[3]*(rho[3]*Q + (1-rho[3])*In);

//***********************************************************************//
// order 1: (int|age|coh) PLACE-FIRST
//          p(int)*p(age|int)*p(coh|age,int)
// eta[1] = eta_a0; eta[2] = eta_c0; eta[3] = eta_ac
// psi[1] = psi_a0; psi[2] = psi_c0; psi[3] = psi_ac

// random int MVN mean
mu_int = rep_vector(0, N_1); 
// random LAT MVN mean
mu_LAT = (eta[1]*r_1_1 + psi[1]*(W*r_1_1)); 
// random drift MVN mean                    
mu_drift = (eta[2]*r_1_1 + psi[2]*(W*r_1_1)) +
           (eta[3]*r_1_2 + psi[3]*(W*r_1_2));
}

model {
vector[N] mu;

target += multi_normal_prec_lpdf(r_1_1| mu_int, Ci_int);
target += multi_normal_prec_lpdf(r_1_2| mu_LAT, Ci_LAT);
target += multi_normal_prec_lpdf(r_1_3| mu_drift, Ci_drift);

//fixed effects contribution to mean
//offset already logged in standat
mu = temp_Intercept + Xc*b + offset; 
//modifies fixed effect mean by ran effects using +=
for (n in 1:N) { 
mu[n] += r_1_1[J_1[n]]*Z_1_1[n] + r_1_2[J_1[n]]*Z_1_2[n] + r_1_3[J_1[n]]*Z_1_3[n];
mu[n] = exp(mu[n]); 
}
// priors including all constants (weakly informative) 
target += normal_lpdf(sum(r_1_1) | 0, 0.001*N_1); // soft sum-to-zero constraints 
target += normal_lpdf(sum(r_1_2) | 0, 0.001*N_1); 
target += normal_lpdf(sum(r_1_3) | 0, 0.001*N_1); 

target += normal_lpdf(temp_Intercept | 0, 10);   
target += normal_lpdf(b | 0, 2); 
target += normal_lpdf(SD | 0, 1); 
target += beta_lpdf(rho   | 0.5, 0.5); 
target += normal_lpdf(eta | 0, 2); 
target += normal_lpdf(psi | 0, 2); 
target += gamma_lpdf(shape | 0.01, 0.01); 

// data likelihood with censored data 
for (n in 1:N) {
if (cens[n] == 0) {
target += neg_binomial_2_lpmf(Y[n] | mu[n], shape);
} else if (cens[n] == -1) {
target += neg_binomial_2_lcdf(Y[n] | mu[n], shape);
}}
}

generated quantities { 
vector[N_1] region_RR;
vector[N_1] region_LAT;
vector[N_1] region_drift;
vector[N] mu_fit;
vector[N] log_lik;
// actual population-level intercept 
real b_Intercept = temp_Intercept - dot_product(means_X, b); 
// region-specific RR
region_RR = exp(r_1_1);
// region-specific LAT %/yr
region_LAT = 100*(exp(b[1] + r_1_2)-1);
// region-specific net drift %/yr
region_drift = 100*(exp(b[2] + r_1_3)-1);

// fitted mean, log-likelihood
mu_fit = b_Intercept + Xc*b + offset; 
for(n in 1:N){
mu_fit[n] += r_1_1[J_1[n]]*Z_1_1[n] + r_1_2[J_1[n]]*Z_1_2[n] + r_1_3[J_1[n]]*Z_1_3[n];
mu_fit[n] = exp(mu_fit[n]); 
// log likelihood with censored data
if (cens[n] == 0) {
log_lik[n] = neg_binomial_2_lpmf(Y[n] | mu_fit[n], shape);
} else if (cens[n] == -1) {
log_lik[n] = neg_binomial_2_lcdf(Y[n] | mu_fit[n], shape);
}}
}


"
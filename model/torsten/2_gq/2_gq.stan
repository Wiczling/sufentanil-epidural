data {
  int<lower=1> nt;
  int<lower=1> nObs;
  int<lower=1> nSubjects;
  array[nObs] int<lower=1> iObs;
  array[nSubjects] int<lower=1> start;
  array[nSubjects] int<lower=1> end;
  array[nt] int<lower=1> cmt;
  array[nt] int evid;
  array[nt] int addl;
  array[nt] int ss;
  array[nt] real amt;
  array[nt] real time;
  array[nt] real rate;
  array[nt] real ii;
  vector<lower=0>[nObs] cObs;
  int<lower=0, upper=1> runestimation; //   a switch to evaluate the likelihood
}

transformed data {
  vector[nObs] logCObs = log(cObs);
  int nTheta = 5;
  int nCmt = 3;
  int nIIV = 5;
  array[nCmt] real biovar = rep_array(1.0, nCmt);
  array[nCmt] real tlag= rep_array(0.0, nCmt);
}

parameters {

  real<lower=0, upper=500> CLHat;
  real<lower=0, upper=500> QHat;
  real<lower=0, upper=3500> V1Hat;
  real<lower=0, upper=3500> V2Hat;
  real<lower=0, upper=10> KAHat;
  real<lower=0> sigma;
  real<lower=3> nu; // normality constant
  // Inter-Individual variability
  vector<lower=0.01, upper=2>[nIIV] omega;
  vector<lower=0, upper=10>[4] fr;
  matrix[nIIV, nSubjects] etaStd;
  cholesky_factor_corr[nIIV] L;
}

transformed parameters {
  vector<lower=0>[nIIV] thetaHat;
  matrix<lower=0>[nSubjects, nIIV] etaM; // variable required for Matt's trick
  array[nTheta] real<lower=0> theta;
  matrix<lower=0>[nCmt, nt] x;
  row_vector<lower=0>[nt] cHat;
  row_vector<lower=0>[nObs] cHatObs;
  
  thetaHat[1] = CLHat;
  thetaHat[2] = QHat;
  thetaHat[3] = V1Hat;
  thetaHat[4] = V2Hat;
  thetaHat[5] = KAHat;
  
 // Matt's trick to use unit scale 
  etaM =  exp(diag_pre_multiply(omega, L * etaStd))'; 
  
  for(j in 1:nSubjects)
  {
    theta[1] = thetaHat[1] * fr[1] * etaM[j, 1] ; // CL
    theta[2] = thetaHat[2] * fr[2] * etaM[j, 2] ; // Q
    theta[3] = thetaHat[3] * fr[3] * etaM[j, 3] ; // V1
    theta[4] = thetaHat[4] * fr[4] * etaM[j, 4] ; // V2
    theta[5] = thetaHat[5] * etaM[j, 5] ; // KA
    
    x[,start[j]:end[j]] = pmx_solve_twocpt(time[start[j]:end[j]], 
                                       amt[start[j]:end[j]],
                                       rate[start[j]:end[j]],
                                       ii[start[j]:end[j]],
                                       evid[start[j]:end[j]],
                                       cmt[start[j]:end[j]],
                                       addl[start[j]:end[j]],
                                       ss[start[j]:end[j]],
                                       theta, biovar, tlag);

    cHat[start[j]:end[j]] = x[2,start[j]:end[j]] ./ theta[3] .*1000; // divide by V1
  }

  cHatObs  = cHat[iObs];
}

model{
  //Informative Priors
  CLHat ~ lognormal(log(45.3),0.25);
  QHat  ~ lognormal(log(38.3),0.25);
  V1Hat ~ lognormal(log(7.90),0.25);
  V2Hat ~ lognormal(log(481),0.25);
  KAHat ~ lognormal(log(1),0.5);
  L~lkj_corr_cholesky(10);
  nu ~ gamma(2,0.1);
  fr ~ lognormal(log(1), 0.25);
 // Inter-individual variability (see transformed parameters block
 // for translation to PK parameters)
  to_vector(etaStd) ~ normal(0, 1);
  omega ~ lognormal(log(0.4),0.25);
  sigma ~ lognormal(log(0.20), 0.25);

  if(runestimation==1){
    logCObs ~ student_t(nu,log(cHatObs), sigma);
  }
}

generated quantities{
  vector[nObs] log_lik;
  for(i in 1:nObs){
   log_lik[i] = student_t_lpdf(cObs[i] | nu, log(fmax(machine_precision(),cHatObs[i])), sigma);
 }
}

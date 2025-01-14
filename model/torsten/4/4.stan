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
  int nTheta = 7;
  int nCmt = 4;
  int nIIV = 7;
  array[nCmt] real biovar = rep_array(1.0, nCmt);
  array[nCmt] real tlag= rep_array(0.0, nCmt);
}

parameters {

  real<lower=0, upper=500> CLHat;
  real<lower=0, upper=500> Q1Hat;
  real<lower=0, upper=500> Q2Hat;
  real<lower=0, upper=3500> V1Hat;
  real<lower=0, upper=3500> V2Hat;
  real<lower=0, upper=3500> V3Hat;
  real<lower=0, upper=20> KAHat;
  real<lower=0> sigma;
  real<lower=3> nu; // normality constant
  // Inter-Individual variability
  vector<lower=0.01, upper=2>[nIIV] omega;
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
  matrix[nCmt, nCmt] K;
  real ka;
  real k10;
  real k12;
  real k21;
  real k13;
  real k31;
  thetaHat[1] = CLHat;
  thetaHat[2] = Q1Hat;
  thetaHat[3] = Q2Hat;
  thetaHat[4] = V1Hat;
  thetaHat[5] = V2Hat;
  thetaHat[6] = V3Hat;
  thetaHat[7] = KAHat;
  
 // Matt's trick to use unit scale 
  etaM =  exp(diag_pre_multiply(omega, L * etaStd))'; 
  
  for(j in 1:nSubjects)
  {
    theta[1] = thetaHat[1] * etaM[j, 1] ; // CL
    theta[2] = thetaHat[2] * etaM[j, 2] ; // Q1
    theta[3] = thetaHat[3] * etaM[j, 3] ; // Q2
    theta[4] = thetaHat[4] * etaM[j, 4] ; // V1
    theta[5] = thetaHat[5] * etaM[j, 5] ; // V2
    theta[6] = thetaHat[6] * etaM[j, 6] ; // V3
    theta[7] = thetaHat[7] * etaM[j, 7] ; // KA
    
    k20 = theta[1]/theta[4];
    k23 = theta[2]/theta[4];
    k32 = theta[2]/theta[5];
    k24 = theta[3]/theta[4];
    k42 = theta[3]/theta[6];
    ka    = theta[7]; 
    
    K = rep_matrix(0, nCmt, nCmt);

    K[1, 1] = -ka;
    K[2, 1] =  ka;
    K[2, 2] = -(k20 + k23 + k24);
    K[2, 3] = k32;
    K[2, 4] = k42;
    K[3, 2] =  k23;
    K[3, 3] = -k32;
    K[4, 2] =  k24;
    K[4, 4] = -k42;
    
    x[,start[j]:end[j]] = pmx_solve_linode(time[start[j]:end[j]], 
                                       amt[start[j]:end[j]],
                                       rate[start[j]:end[j]],
                                       ii[start[j]:end[j]],
                                       evid[start[j]:end[j]],
                                       cmt[start[j]:end[j]],
                                       addl[start[j]:end[j]],
                                       ss[start[j]:end[j]],
                                       K, biovar, tlag);
 
    cHat[start[j]:end[j]] = x[2,start[j]:end[j]] ./ theta[3] .*1000; // divide by V1
  }

  cHatObs  = cHat[iObs];
}

model{
  //Informative Priors
  // https://pubmed.ncbi.nlm.nih.gov/8533912/
  // VC = 14.6, VP1 = 66, VP2 = 608, Cl = 52.8 L/h, and Q = 102 L/h i 40.8 L/h

  CLHat ~ lognormal(log(52.8),0.25);
  Q1Hat ~ lognormal(log(102),0.25);
  Q2Hat ~ lognormal(log(40.8),0.25);
  V1Hat ~ lognormal(log(14.6),0.25);
  V2Hat ~ lognormal(log(66),0.25);
  V3Hat ~ lognormal(log(608),0.25);
  KAHat ~ lognormal(log(1),1);
  L~lkj_corr_cholesky(10);
  nu ~ gamma(2,0.1);
  to_vector(etaStd) ~ normal(0, 1);
  omega ~ lognormal(log(0.4),0.5);
  sigma ~ lognormal(log(0.20), 0.5);

  if(runestimation==1){
    logCObs ~ student_t(nu,log(cHatObs), sigma);
  }
}

generated quantities{
  array[nt] real cCond;
  array[nt] real cPred;
  row_vector[nt] cHatPred;
  matrix[nCmt, nt] xPred;
  matrix[nIIV, nSubjects] etaStdPred;
  matrix<lower=0>[nSubjects, nIIV] etaPredM;
  corr_matrix[nIIV] rho;
  array[nTheta] real<lower=0> thetaPred;
  matrix[nCmt, nCmt] KPred;
  real kaPred;
  real k20Pred;
  real k23Pred;
  real k32Pred;
  real k24Pred;
  real k42Pred;
  vector[nObs] log_lik;
  vector[nObs] errors;
  
  
  rho = L * L';

    for(i in 1:nSubjects){
      for(j in 1:nIIV){ 
        etaStdPred[j, i] = normal_rng(0, 1);
      }
    }

    etaPredM = exp(diag_pre_multiply(omega, L * etaStdPred))';

    for(j in 1:nSubjects){
     
    thetaPred[1] = thetaHat[1] * etaPredM[j, 1] ; // CL
    thetaPred[2] = thetaHat[2] * etaPredM[j, 2] ; // Q1
    thetaPred[3] = thetaHat[3] * etaPredM[j, 3] ; // Q2
    thetaPred[4] = thetaHat[4] * etaPredM[j, 4] ; // V1
    thetaPred[5] = thetaHat[5] * etaPredM[j, 5] ; // V2
    thetaPred[6] = thetaHat[6] * etaPredM[j, 6] ; // V3
    thetaPred[7] = thetaHat[7] * etaPredM[j, 7] ; // KA
      
    k20Pred = thetaPred[1]/thetaPred[4];
    k23Pred = thetaPred[2]/thetaPred[4];
    k32Pred = thetaPred[2]/thetaPred[5];
    k24Pred = thetaPred[3]/thetaPred[4];
    k42Pred = thetaPred[3]/thetaPred[6];
    kaPred    = thetaPred[7]; 
    
    KPred = rep_matrix(0, nCmt, nCmt);

    K[1, 1] = -kaPred;
    K[2, 1] =  kaPred;
    K[2, 2] = -(k20Pred + k23Pred + k24Pred);
    K[2, 3] = k32Pred;
    K[2, 4] = k42Pred;
    K[3, 2] =  k23Pred;
    K[3, 3] = -k32Pred;
    K[4, 2] =  k24Pred;
    K[4, 4] = -k42Pred;
    
    xPred[,start[j]:end[j]] = pmx_solve_linode(time[start[j]:end[j]], 
                                       amt[start[j]:end[j]],
                                       rate[start[j]:end[j]],
                                       ii[start[j]:end[j]],
                                       evid[start[j]:end[j]],
                                       cmt[start[j]:end[j]],
                                       addl[start[j]:end[j]],
                                       ss[start[j]:end[j]],
                                       KPred, biovar, tlag);
    
                                       
    cHatPred[start[j]:end[j]] = xPred[2,start[j]:end[j]] ./ thetaPred[3]*1000; // divide by V1
  }

  for(i in 1:nt){
      cCond[i] = exp(student_t_rng(nu,log(fmax(machine_precision(),cHat[i])), sigma));     // individual predictions
      cPred[i] = exp(student_t_rng(nu,log(fmax(machine_precision(),cHatPred[i])), sigma)); // population predictions
 }
 
  for(i in 1:nObs){
   errors[i] = logCObs[i]-log(fmax(machine_precision(),cHatObs[i]));
   log_lik[i] = student_t_lpdf(logCObs[i] | nu, log(fmax(machine_precision(),cHatObs[i])), sigma);
 }
 
}

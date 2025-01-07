functions {
  real partial_sum(array[] int ind, int start, int end, 
                   vector logDV, 
                   array[] int startx, 
                   array[] int endx,
                   array[] real time, 
                   array[] real amt,
                   array[] real rate,
                   array[] real ii,
                   array[] int evid,
                   array[] int cmt,
                   array[] int addl,
                   array[] int ss,
                   array[] matrix K, 
                   array[] vector itheta,
                   array[] real biovar, 
                   array[] real tlag,
                   real sigma,
                   real nu) {
                     
    real lp = 0;

    for (z in start : end) {
      
    int nti=endx[z]-startx[z]+1;
    array[nti] int include = evid[startx[z]:endx[z]];
    array[nti] int includex = rep_array(0,nti);
    
    for (n in 1:nti) {
      if (include[n]>0)
      includex[n]=1;
    }
    
    array[sum(includex)] int ns1;
    
    int n1 = 1;
    for (n in 1:nti) {
      if (includex[n]==1) {
        ns1[n1] = n;
        n1 += 1;
        }}
    
    matrix[4,nti] x = pmx_solve_linode(time[startx[z]:endx[z]], 
                                       amt[startx[z]:endx[z]],
                                       rate[startx[z]:endx[z]],
                                       ii[startx[z]:endx[z]],
                                       evid[startx[z]:endx[z]],
                                       cmt[startx[z]:endx[z]],
                                       addl[startx[z]:endx[z]],
                                       ss[startx[z]:endx[z]],
                                       K[z], biovar, tlag);
 
    row_vector[nti] cHat = x[2,:] ./ itheta[z,3] .*1000; // divide by V1

    lp = lp + student_t_lpdf(logDV[ns1] | nu, log(cHat[ns1]),  sigma);
      
    }
    return lp;
  }
}
data {
  int<lower=1> nt;
  int<lower=1> nObs;
  int<lower=1> nSubjects;
  array[nObs] int<lower=1> iObs;
  array[nSubjects] int<lower=1> startx;
  array[nSubjects] int<lower=1> endx;
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
  
  int grainsize = 1;
  array[nObs] int ind = rep_array(1, nObs);
  vector[nt] logDV = rep_vector(0,nt);
  logDV[iObs]=logCObs;
}

parameters {

  real<lower=0, upper=500> CLHat;
  real<lower=0, upper=500> QHat;
  real<lower=0, upper=3500> V1Hat;
  real<lower=0, upper=3500> V2Hat;
  real<lower=0, upper=10> KAHat;
    real<lower=0, upper=10> KA14Hat;
      real<lower=0, upper=10> KA41Hat;
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
  array[nSubjects] matrix[nCmt, nCmt] K;
  array[nSubjects] vector[nIIV] itheta;
  
  real ka;
  real ka14;
  real ka41;
  real k10;
  real k12;
  real k21;
  
  thetaHat[1] = CLHat;
  thetaHat[2] = QHat;
  thetaHat[3] = V1Hat;
  thetaHat[4] = V2Hat;
  thetaHat[5] = KAHat;
  thetaHat[6] = KA14Hat;
  thetaHat[7] = KA41Hat;
  
 // Matt's trick to use unit scale 
  etaM =  exp(diag_pre_multiply(omega, L * etaStd))'; 
  
  for(j in 1:nSubjects){
    {
      
    vector[nTheta] theta;
    matrix[nCmt, nCmt] Ki;
        
    theta[1] = thetaHat[1] * etaM[j, 1] ; // CL
    theta[2] = thetaHat[2] * etaM[j, 2] ; // Q
    theta[3] = thetaHat[3] * etaM[j, 3] ; // V1
    theta[4] = thetaHat[4] * etaM[j, 4] ; // V2
    theta[5] = thetaHat[5] * etaM[j, 5] ; // KA
    theta[6] = thetaHat[6] * etaM[j, 6] ; // KA14
    theta[7] = thetaHat[7] * etaM[j, 7] ; // KA41
    
    itheta[j] = theta;
    
    k10 = theta[1]/theta[3];
    k12 = theta[2]/theta[3];
    k21 = theta[2]/theta[4];
    ka    = theta[5]; 
    ka14  = theta[6];
    ka41  = theta[7];
    
    Ki = rep_matrix(0, nCmt, nCmt);
    Ki[1, 1] = -(ka+ka14);
    Ki[1, 4] = ka41;
    Ki[2, 1] = ka;
    Ki[2, 2] = -(k10 + k12);
    Ki[2, 3] = k21;
    Ki[3, 2] = k12;
    Ki[3, 3] = -k21;
    Ki[4, 1] = ka14;
    Ki[4, 4] = -ka41;
    
    K[j] = Ki;
    }
  }
}   

model{
  //Informative Priors
  CLHat ~ lognormal(log(45.3),0.25);
  QHat  ~ lognormal(log(38.3),0.25);
  V1Hat ~ lognormal(log(7.90),0.25);
  V2Hat ~ lognormal(log(481),0.25);
  KAHat ~ lognormal(log(1),0.5);
  KA14Hat ~ lognormal(log(4.85),0.5);
  KA41Hat ~ lognormal(log(0.080),0.5);
  L~lkj_corr_cholesky(10);
  nu ~ gamma(2,0.1);
  to_vector(etaStd) ~ normal(0, 1);
  omega ~ lognormal(log(0.4),0.25);
  sigma ~ lognormal(log(0.20), 0.25);

  if(runestimation==1){
    // logCObs ~ student_t(nu,log(cHatObs), sigma);
    target += reduce_sum(partial_sum, ind, grainsize, 
                         logDV, startx, endx,
                         time, 
                         amt,
                         rate,
                         ii,
                         evid,
                         cmt,
                         addl,
                         ss,
                         K, itheta, biovar, tlag,
                         sigma, nu);
  }
}

generated quantities{
  array[nt] real cCond;
  array[nt] real cPred;
  row_vector[nt] cHatPred;
  row_vector[nt] cHatCond;
  matrix[nCmt, nt] xPred;
  matrix[nCmt, nt] xCond;
  matrix[nIIV, nSubjects] etaStdPred;
  matrix<lower=0>[nSubjects, nIIV] etaPredM;
  corr_matrix[nIIV] rho;
  array[nTheta] real<lower=0> thetaPred;
  matrix[nCmt, nCmt] KPred;
  real kaPred;
  real ka14Pred;
  real ka41Pred;
  real k10Pred;
  real k12Pred;
  real k21Pred;
  
  rho = L * L';

    for(i in 1:nSubjects){
      for(j in 1:nIIV){ 
        etaStdPred[j, i] = normal_rng(0, 1);
      }
    }

    etaPredM = exp(diag_pre_multiply(omega, L * etaStdPred))';

    for(j in 1:nSubjects){
     
    thetaPred[1] = thetaHat[1] * etaPredM[j, 1] ; // CL
    thetaPred[2] = thetaHat[2] * etaPredM[j, 2] ; // Q
    thetaPred[3] = thetaHat[3] * etaPredM[j, 3] ; // V1
    thetaPred[4] = thetaHat[4] * etaPredM[j, 4] ; // V2
    thetaPred[5] = thetaHat[5] * etaPredM[j, 5] ; // KA
    thetaPred[6] = thetaHat[6] * etaPredM[j, 6] ; // KA14
    thetaPred[7] = thetaHat[7] * etaPredM[j, 7] ; // KA41 
      
    k10Pred = thetaPred[1]/thetaPred[3];
    k12Pred = thetaPred[2]/thetaPred[3];
    k21Pred = thetaPred[2]/thetaPred[4];
    kaPred    = thetaPred[5]; 
    ka14Pred  = thetaPred[6];
    ka41Pred  = thetaPred[7];
    
    KPred = rep_matrix(0, nCmt, nCmt);

    KPred[1, 1] = -(kaPred+ka14Pred);
    KPred[1, 4] = ka41Pred;
    KPred[2, 1] = kaPred;
    KPred[2, 2] = -(k10Pred + k12Pred);
    KPred[2, 3] = k21Pred;
    KPred[3, 2] = k12Pred;
    KPred[3, 3] = -k21Pred;
    KPred[4, 1] = ka14Pred;
    KPred[4, 4] = -ka41Pred;
    
    xCond[,startx[j]:endx[j]] = pmx_solve_linode(time[startx[j]:endx[j]], 
                                       amt[startx[j]:endx[j]],
                                       rate[startx[j]:endx[j]],
                                       ii[startx[j]:endx[j]],
                                       evid[startx[j]:endx[j]],
                                       cmt[startx[j]:endx[j]],
                                       addl[startx[j]:endx[j]],
                                       ss[startx[j]:endx[j]],
                                       K[j], biovar, tlag);
                                       
    xPred[,startx[j]:endx[j]] = pmx_solve_linode(time[startx[j]:endx[j]], 
                                       amt[startx[j]:endx[j]],
                                       rate[startx[j]:endx[j]],
                                       ii[startx[j]:endx[j]],
                                       evid[startx[j]:endx[j]],
                                       cmt[startx[j]:endx[j]],
                                       addl[startx[j]:endx[j]],
                                       ss[startx[j]:endx[j]],
                                       KPred, biovar, tlag);
    
                                       
    cHatCond[startx[j]:endx[j]] = xCond[2,startx[j]:endx[j]] ./ itheta[j,3]*1000; // divide by V1
    cHatPred[startx[j]:endx[j]] = xPred[2,startx[j]:endx[j]] ./ thetaPred[3]*1000; // divide by V1
  }

  for(i in 1:nt){
      cCond[i] = exp(student_t_rng(nu,log(fmax(machine_precision(),cHatCond[i])), sigma));     // individual predictions
      cPred[i] = exp(student_t_rng(nu,log(fmax(machine_precision(),cHatPred[i])), sigma)); // population predictions

 }
}

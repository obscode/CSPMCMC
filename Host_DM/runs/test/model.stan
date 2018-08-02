
functions {
   real modelCeph(real M0, real DM, real betaP, real betaVI, real betaOH, real logP, 
                  real VI, real OH) {
      return ( M0 + DM + betaP*logP + betaVI*VI + betaOH*(OH-9.5));
   }
    
   real toflux(real mag, real zp) {
      return pow(10, -0.4*(mag - zp));
   }

   real fluxerr (real mag, real e_mag, real zp) {
      real f;
      f = toflux(mag, zp);
      return (f*e_mag/1.087);
   }
}

data {
   int<lower=4> N;    // Total number of data points in all sets
   int<lower=1> S;    // Number of data sets
   
   real P[N];         // the 1st predictor
   real VI[N];        // the 2nd predictor
   real OH[N];        // the 3rd predictor
   real mag[N];         // the outcome
   real e_mag[N];
   int ID[N];         // Which set does the data point belong to?
   int in_mag;          // ==1 if fit in magnitudes, ==0 if fit in flux 

   int<lower=1> NMW;        // Number of MW Cepheids
   real pi[NMW];             // parallax in arc-sec
   real e_pi[NMW];       // error in parallax
   real magMW[NMW];
   real e_magMW[NMW];           //error in V magnitude for MW cepheids
   real VI_MW[NMW];
   real P_MW[NMW];
   real OH_MW;
   int LK_prior;            // Apply Lutz-Kelker prior?
   
   int <lower=1> NLMC;      // number of data in LMC
   real P_LMC[NLMC];
   real VI_LMC[NLMC];
   real OH_LMC;
   real magLMC[NLMC];
   real e_magLMC[NLMC];

   int<lower=1> N4258;          //Number of data in NGC4258
   real P_4258[N4258];   
   real VI_4258[N4258];
   real OH_4258[N4258];
   real mag4258[N4258];
   real e_mag4258[N4258];
}

parameters {
   real <lower=10.0, upper=40.0> DM[S];           // intercept
   real betaP;           // Period-slope
   real M;               // Absolute mag of Cepheids
   real betaOH;           // slope of metallicity dependence
   real<lower=0> eps;    // dispersion
   real<lower=0.001> lowP;
   real<lower=0.001> highP;
   real<lower=0> betaVI;
real <lower=0> var1;   //additive 

   real<lower=0.001> eps_MW;
   real dzp_MW;
   real<lower=0.000125, upper=1.0> pi_true[NMW]; // true parallax, no farther 
                                                 // than 8kpc

   real DM_LMC;
   real<lower=0> eps_LMC;
   real dzp_LMC;

   real DM_4258;
   real<lower=0> eps_4258;
}

model {
   real modl;
   real dist;
   real edist;
   real zp;
   real e;


   zp = 25.0;
   for (i in 1:N) {
      modl = modelCeph(M, DM[ID[i]], betaP, betaVI, betaOH, P[i], VI[i], OH[i]);
      e = eps + exp(-(P[i]-0.4)/lowP) + exp(-(2.1-P[i])/highP);
      if (in_mag==0){
         modl = toflux(modl, zp);
         mag[i] ~ normal(modl, sqrt(var1 + e*e*mag[i]*mag[i] + e_mag[i]*e_mag[i]));
      }
      else {
         mag[i] ~ normal(modl, sqrt(e*e +e_mag[i]*e_mag[i]));
      }
   }

   zp = 3.0;
   dzp_MW ~ normal(0, 0.03);
   for (i in 1:NMW) 
      pi[i] ~ normal(pi_true[i], e_pi[i]);
   if (LK_prior == 1) {
      for (i in 1:NMW)
         increment_log_prob(-17.2812 - 3*log(pi_true[i]));       
   }
   for (i in 1:NMW) {
      dist = -5*log10(pi_true[i]) - 5;
      modl = modelCeph(M, dist, betaP, betaVI, betaOH, P_MW[i], VI_MW[i], OH_MW) +
              dzp_MW;
      e = eps_MW + exp(-(P_MW[i] - 0.4)/lowP) + exp(-(2.1-P_MW[i])/highP);
      if (in_mag == 0){
         modl = toflux(modl, zp);
         magMW[i] ~ normal(modl, sqrt(e*e*magMW[i]*magMW[i] + 
                                       e_magMW[i]*e_magMW[i]));
      }
      else {
         magMW[i] ~ normal(modl, sqrt(e*e + e_magMW[i]*e_magMW[i]));
      }
   }

   zp = 12.0;
   DM_LMC ~ normal(18.49,0.05);
   dzp_LMC ~ normal(0.0, 0.03);
   for (i in 1:NLMC) { 
      modl = modelCeph(M, DM_LMC, betaP, betaVI, betaOH, P_LMC[i], 
              VI_LMC[i], OH_LMC) + dzp_LMC;
      e = eps_LMC + exp(-(P_LMC[i]-0.4)/lowP) + exp(-(2.1-P_LMC[i])/highP);
      if (in_mag ==0){
         modl = toflux(modl, zp);
         magLMC[i] ~ normal(modl ,
              sqrt(e*e*magLMC[i]*magLMC[i]+ e_magLMC[i]*e_magLMC[i])); 
      } else {
         magLMC[i] ~ normal(modl ,sqrt(e*e+ e_magLMC[i]*e_magLMC[i]));   
      }
   }

   zp = 23.0;
   DM_4258 ~ normal(29.387000,0.056800);
   for (i in 1:N4258) {
      modl = modelCeph(M, DM_4258, betaP, betaVI, betaOH, P_4258[i], VI_4258[i],OH_4258[i]);
      e = eps_4258 + exp(-(P_4258[i]-0.4)/lowP) + exp(-(2.1-P_4258[i])/highP);
      if (in_mag == 0){
         modl <- toflux(modl, zp);
         mag4258[i] ~ normal(modl, sqrt(e*e*mag4258[i]*mag4258[i]+
              e_mag4258[i]*e_mag4258[i]));
      }
      else {
         mag4258[i] ~ normal(modl,sqrt(e*e+ e_mag4258[i]*e_mag4258[i]));
      }
   }
}
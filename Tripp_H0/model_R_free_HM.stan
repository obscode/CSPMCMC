data {
   int<lower=1> C;             // Number of Cepheid hosts
   cov_matrix[C] Ccov;         // Covariance matrix for Cepheids
   vector[C] DMCeph;           // mean cepheid distances
   int<lower=1> NSNe;          // Number of SNe
   int <lower=-1> host[NSNe];  // Host galaxy index
   int<lower=1> NObs;          // Number of observations
   int<lower=1> NFs;          // Number of filters
   real ms[NObs];               // magnitudes
   real vms[NObs];              // error in magnitudes
   real cs[NSNe];              // colors
   real vcs[NSNe];             // variance in color
   real st[NSNe];              // s_BV
   real est[NSNe];             // error in s_BV
   real zcmb[NSNe];            // CMB redshift
   real zhel[NSNe];            // heliocentric redshift
   real K[NSNe];            //  log10(Host mass)
   int sindex[NObs];             // SN index
   int findex[NObs];           // filter index
   int findex1;           // filter index of first color filter
   int findex2;           // filter index of second color filter
   int Nphotsys;          // Number of photometric systems
   int photsys[NSNe];     // index of photometric system
   vector[NFs] zperr[Nphotsys];   //zero-point errors
   real M0;                      // host mass zero-point
   real sigmaK[NSNe];            // error in Host Mass estimate
}


transformed data {
   vector[2] mc[NObs];
   for (i in 1:NObs) {
      mc[i,1] = ms[i];
      mc[i,2] = cs[sindex[i]];
   }
}


parameters {
   real <lower=-5.0,upper=5.0>a[NFs];
   real <lower=-10.0,upper=10.0>b[NFs];
   real <lower=-10.0,upper=10.0>c[NFs];
   real <lower=0.0,upper=10.0>Rl[NFs];
   real <lower=0.0,upper=1.0>evar[NFs];
   real <lower=0.0,upper=10.0>vpec;
   real <lower=40, upper=100> H0;
   real <lower=-10, upper=10> ct[NSNe];
   real <lower=-1, upper=1> alpha;
   vector[C] DMhost;
   vector<lower=-5, upper=5>[NFs] zpoff[Nphotsys];   // zero-point offsets
}

transformed parameters {
   real DM[NSNe];
   cov_matrix[2] covar[NObs];
   for (i in 1:NSNe) {
      if (host[i] < 0) {
         DM[i] = 5.0*log10(((1+zhel[i])/(1+zcmb[i]))*(300000.0/H0)*
                 (zcmb[i] + zcmb[i]*zcmb[i]*0.79)) + 25.0;
      } else {
         DM[i] = DMhost[host[i]];
      }
   }
   for (i in 1:NObs) {
      covar[i][1,1] = vms[i] + evar[findex[i]] +
                               square(sigmaK[sindex[i]]*alpha);
      if (host[sindex[i]] < 0) {
         covar[i][1,1] = covar[i][1,1] + pow(0.000723*vpec/zcmb[sindex[i]],2);
      }
      covar[i][2,2] = vcs[sindex[i]];
      covar[i][2,1] = 0;
      if (findex[i] == findex1) covar[i][2,1] = covar[i][2,1] + vms[i];
      if (findex[i] == findex2) covar[i][2,1] = covar[i][2,1] - vms[i];
      covar[i][1,2] = covar[i][2,1];
   }
}

model {
   vector[2] mct;
   real mass;
   
   DMhost ~ multi_normal(DMCeph, Ccov);

   for (i in 1:NObs) {
      mass = -0.4*(K[sindex[i]] - DM[sindex[i]]) + 1.04
      mct[1] = (a[findex[i]] - 19) + b[findex[i]]*(st[sindex[i]] - 1) + 
               c[findex[i]]*(st[sindex[i]]-1)*(st[sindex[i]]-1) + 
               Rl[findex[i]]*ct[sindex[i]] + 
               alpha*(mass - M0) + DM[sindex[i]];
      if ( photsys[sindex[i]] > 0) {
         mct[1] = mct[1] + zpoff[photsys[sindex[i]],findex[i]];
      }
      mct[2] = ct[sindex[i]];

      mc[i] ~ multi_normal(mct, covar[i]);
   }
   for (i in 1:Nphotsys) {
      for (j in 1:NFs) {
         zpoff[i,j] ~ normal(0, zperr[i,j]);
      }
   }
}

generated quantities {
   real g_modl[NObs]; 
   real mass;

   for (i in 1:NObs) {
      mass = -0.4*(K[sindex[i]] - DM[sindex[i]]) + 1.04;
      g_modl[i] = (a[findex[i]] - 19) + b[findex[i]]*(st[sindex[i]] - 1) +
                   c[findex[i]]*(st[sindex[i]]-1)*(st[sindex[i]]-1) +
                   Rl[findex[i]]*ct[sindex[i]] + 
                   alpha*(mass - M0) + DM[sindex[i]];
   }
}

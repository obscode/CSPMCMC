functions {
   real A_lamb(vector Al_coef, real Rv, real EBV, int order) {
      int ii;   // 1D index
      real Al;  // the answer
      ii = 1;
      Al = 0;
      for (j in 0:order) {
         for (i in 0:order) {
            if (j+i <= order) {
               Al = Al + Al_coef[ii]*pow(Rv,i)*pow(EBV,j+1);
               ii = ii+1;
            }
         }
      }
      return Al;
   }
}

data {
   int<lower=1> C;             // Number of Cepheid hosts
   cov_matrix[C] Ccov;         // Covariance matrix for Cepheids
   vector[C] DMCeph;           // mean cepheid distances
   int<lower=1> NSNe;          // Number of SNe
   int <lower=-1> host[NSNe];  // Host galaxy index
   int<lower=1> NObs;          // Number of observations
   int<lower=1> NFs;           // Number of filters
   int<lower=1> NFs2;          // Number of filters in Alamb matrix
   real ms[NObs];              // magnitudes
   real vms[NObs];             // error in magnitudes
   real cs[NSNe];              // colors
   real vcs[NSNe];             // variance in color
   real st[NSNe];              // s_BV
   real est[NSNe];             // error in s_BV
   real zcmb[NSNe];            // CMB redshift
   real zhel[NSNe];            // heliocentric redshift
   int sindex[NObs];           // SN index
   int findex[NObs];           // filter index
   int findex1;                // filter index of first color filter
   int findex2;                // filter index of second color filter
   int Al_order[NFs2];         // actual order of poly
   int N_coef;                 // size of array
   vector[N_coef] Amat[NFs2];  // the coefficents
   int Nphotsys;          // Number of photometric systems
   int photsys[NSNe];     // index of photometric system
   vector[NFs] zperr[Nphotsys];   //zero-point errors
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
   real <lower=0.0,upper=10.0>R;
   real <lower=0.0,upper=1.0>evar[NFs];
   real <lower=0.0,upper=10.0>vpec;
   real <lower=40, upper=100> H0;
   real <lower=-10, upper=10> ct[NSNe];
   vector<lower=-5, upper=5>[NFs] zpoff[Nphotsys];   // zero-point offsets
}

transformed parameters {
   cov_matrix[2] covar[NObs];
   real Rl[NFs];
   real DM[NSNe];
   for (i in 1:NSNe) {
      if (host[i] < 0) {
         DM[i] = 5.0*log10(((1+zhel[i])/(1+zcmb[i]))*(300000.0/H0)*
                 (zcmb[i] + zcmb[i]*zcmb[i]*0.79)) + 25.0;
      } else {
         DM[i] = DMCeph[host[i]];
      }
   }
   for (i in 1:NObs) {
      covar[i][1,1] = vms[i] + evar[findex[i]];
      if (host[sindex[i]] < 0) {
         covar[i][1,1] = covar[i][1,1] + pow(0.000723*vpec/zcmb[sindex[i]],2);
      }
      covar[i][2,2] = vcs[sindex[i]];
      covar[i][2,1] = 0;
      if (findex[i] == findex1) covar[i][2,1] = covar[i][2,1] + vms[i];
      if (findex[i] == findex2) covar[i][2,1] = covar[i][2,1] - vms[i];
      covar[i][1,2] = covar[i][2,1];
   }
   for (i in 1:NFs) {
      Rl[i] = A_lamb(Amat[i], R, 0.1, Al_order[i])/
           (A_lamb(Amat[findex1], R, 0.1, Al_order[findex1]) - 
            A_lamb(Amat[findex2], R, 0.1, Al_order[findex2]));
   }
}

model {
   vector[2] mct;
   
   for (i in 1:NObs) {
      mct[1] = (a[findex[i]] - 19) + b[findex[i]]*(st[sindex[i]] - 1) + 
               c[findex[i]]*(st[sindex[i]]-1)*(st[sindex[i]]-1) + 
               Rl[findex[i]]*ct[sindex[i]] + DM[sindex[i]];
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


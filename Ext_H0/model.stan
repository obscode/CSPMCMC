
functions {
   real A_lamb(vector Al_coef, real EBV, real Rv, int Al_order){
      int ii;  // 1-D index
      real Al;  // the answer
      ii =1;
      Al =0;
      for (j in 0:Al_order){
         for (i in 0:Al_order){
            if(j+i <= Al_order){
               Al = Al + Al_coef[ii]*pow(Rv,i)*pow(EBV,j+1);
               ii = ii+1;
            }       
         }
      }
      return Al;
   }

   real vA_lamb(vector Al_coef, real EBV, real Rv, real eEBV, real eRv,
                real cov, int Al_order){
      int ii;  // 1-D index
      real dAldE; // derivative of Al w.r.t E(B-V)
      real dAldR; // derivative of Al w.r.t Rv
      ii =1;
      dAldE =0;
      dAldR =0;
      for (j in 0:Al_order){
         for (i in 0:Al_order){
            if(j+i <= Al_order){
               dAldE = dAldE + Al_coef[ii]*(j+1)*pow(Rv,i)*pow(EBV,j);
               dAldR = dAldR + Al_coef[ii]*i*pow(Rv,i-1)*pow(EBV,j+1);
               ii = ii+1;
            }       
         }
      }
      return dAldE*dAldE*eEBV*eEBV + dAldE*dAldR*cov + dAldR*dAldR*eRv*eRv;
   }

   real modelSN(real DM, vector a, vector B, real Al, real alpha, real K,
       real M0) {
       real mass;
       mass = 0.4*(DM - K) + 1.04;
       return ( dot_product(B, a) - 19 + DM + Al + alpha*(mass-M0));
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
   int<lower=1> S;    // Number of Cepheid hosts
   cov_matrix[S] Cov;   // Covariance matrix
   vector[S] DMobs;      // mean

   // SUPERNOVAE
   int in_mag;
   int <lower=1> Nsn;   // number of SN
   real st[Nsn];
   real e_st[Nsn];
   
   real EBV[Nsn];      //E(B-V)
   real Rv[Nsn];       //Rv reddening
   matrix[2,2] ERprec[Nsn]; // E(B-V)/Rv precision matrices
   int Nphotsys;            // Number of photometric systems
   int photsys[Nsn];        // index of photometric system

   int <lower=-1> host[Nsn];      //Indexes the cepheid host
   real zcmb[Nsn];          // redshift CMB
   real zhelio[Nsn];        // redshift Heliocentric 

   int <lower=1> Nbasis;    // number of basis functions
   vector [Nbasis] Bs[Nsn]; // stretch basis functions
   real K[Nsn];             // Host K-band mag"
   real sigmaK[Nsn];        //  error therein"
   real M0;                 //  Host mass zero-point"

   real alpha;

   int <lower=1> NSNobs;  // number of observations
   int <lower=1> Nfilt; // number of filters
   int filt[NSNobs];      // which filter
   int oid[NSNobs];       // which object
   real m_sn[NSNobs];      // Bmax, Vmax, Imax , ect..
   real em_sn[NSNobs];             
   vector[Nfilt] zperr[Nphotsys]; // error in zero-points

   int Al_order[Nfilt];          // order of the polynomial
   int N_coef;
   vector[N_coef] Al_coef[Nfilt];  //coefficients of poly
}

parameters {
   real H0;              // Hubble constant
   real<lower=0, upper=100> pec_vel;  // peculiar velocity in units of 100km/s

   vector<lower=0, upper=100>[S] DM;    // The true DMs

    real<lower=0, upper=10> eps_sn[Nfilt];    // intrinsic dispersion
    vector<lower=-10, upper=10>[2] EBV_Rv[Nsn];          // E(B-V) and RV
    vector<lower=-5, upper=5>[Nfilt] zpoff[Nphotsys];  // zero-point errors
    vector<lower=-10, upper=10>[Nbasis] a[Nfilt];        // basis coefficients
}

transformed parameters {
   real dist[Nsn];    // DM for hosts, distant and cross-validation
   for (i in 1:Nsn) {
      if ( host[i] > 0) {
         dist[i] = DM[host[i]];
      } else {
         dist[i] = 5.0*log10((1+zhelio[i])/(1+zcmb[i])*300000.0/H0*
                   (zcmb[i]+zcmb[i]*zcmb[i]*0.79))+ 25.0;
      }
   }
}


model {
   real modl;
   real edist;
   real Al;
   vector[2] mu;
   matrix[2,2] prec;
   real zp;
   real varm;

   // The distance modulii
   DM ~ multi_normal(DMobs, Cov);
  //SUPERNOVAE
  zp = 16.0;
  for (i in 1:Nsn) {
     mu[1] = EBV[i];
     mu[2] = Rv[i];
     prec = ERprec[i];
     EBV_Rv[i] ~ multi_normal_prec(mu, prec);
  }

  for (i in 1:NSNobs){
     Al = A_lamb(Al_coef[filt[i]], EBV_Rv[oid[i]][1], EBV_Rv[oid[i]][2], 
                  Al_order[filt[i]]);
     modl = modelSN(dist[oid[i]], a[filt[i]], Bs[oid[i]], Al, alpha,
                    K[oid[i]], M0);
     if ( photsys[oid[i]] > 0) {
        modl = modl + zpoff[photsys[oid[i]],filt[i]];
     }
     if (host[oid[i]] > 0) {
        if ( in_mag == 0) {
           modl = toflux(modl, zp);
           m_sn[i] ~normal(modl, 
               sqrt(square(eps_sn[filt[i]]*m_sn[i]) +
                    square(em_sn[i]) + square(alpha*sigmaK[oid[i]]*m_sn[i]))) ;
        } else {
           m_sn[i] ~normal(modl, sqrt(square(eps_sn[filt[i]]) + 
                                   square(em_sn[i]) + 
                                   square(alpha*sigmaK[oid[i]]))) ;
        }
     }
     else{ 
        if ( in_mag == 0) {
           modl = toflux(modl, zp);
           edist = 0.000723*pec_vel/zcmb[oid[i]]*modl/1.087;
           m_sn[i] ~normal(modl, 
                     sqrt(square(eps_sn[filt[i]]*m_sn[i]) + 
                                square(edist) + 
                                square(em_sn[oid[i]]) + 
                                square(alpha*sigmaK[oid[i]]*m_sn[i])));
        } 
        else {
           edist = 0.000723*pec_vel/zcmb[oid[i]];
           m_sn[i] ~normal(modl, sqrt(square(eps_sn[filt[i]]) + 
                           square(edist) + square(em_sn[i]) +
                           square(alpha*sigmaK[oid[i]])));
        }
         
     }
   }
   for (i in 1:Nphotsys) {
      for (j in 1:Nfilt) {
         zpoff[i,j] ~ normal(0, zperr[i,j]);
      }
   }
}
generated quantities {
      real g_modl[NSNobs];
      real g_Al[NSNobs];

      for (i in 1:NSNobs) {
         g_Al[i] = A_lamb(Al_coef[filt[i]], EBV_Rv[oid[i]][1], 
            EBV_Rv[oid[i]][2], Al_order[filt[i]]);
         g_modl[i] = modelSN(dist[oid[i]], a[filt[i]], Bs[oid[i]], g_Al[i],
                             alpha, K[oid[i]], M0);
      }
}

# This part of the code is responsible for putting together the STAN code.
# I'm outsourcing this to its own module because some of the options were
# impossible to do in the STAN language with if/else blocks. 

from hashlib import md5
import pystan
import pickle
import os
from numpy import random, array
from numpy import zeros,linalg

# First, the static stuff.  These are blocks of code that will be assembled
# by the generate_STAN function as needed.
STAN_functions = '''
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
'''

#########################
## Data blocks
#########################
data_ceph = '''
data {
   int<lower=1> S;    // Number of Cepheid hosts
   cov_matrix[S] Cov;   // Covariance matrix
   vector[S] DMobs;      // mean
'''

trans_data_ceph = '''
transformed data {
   vector[S] DM;
   for (i in 1:S) {
      DM[i] = DMobs[i];
   }
}'''

data_sn = '''
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
'''

data_sn_1filt = '''
   real m_sn[Nsn];          // Bmax, Vmax, Imax , ect..
   real em_sn[Nsn];             

   int Al_order;            // order of the polynomial
   int N_coef;
   vector[N_coef] Al_coef;  //coefficients of poly
   real zperr[Nphotsys];    // error in zero-points of photsys
'''

data_sn_Nfilt = '''
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
'''

data_cv = '''
   int <lower=1> NCV;   // number of cross-validations
   int <lower=0> CVids[Nsn];  // ID into DMCV
'''

#########################
## Parameter blocks
#########################

param_common_fixed = '''
parameters {
   real H0;              // Hubble constant
   real<lower=0, upper=100> pec_vel;  // peculiar velocity in units of 100km/s
'''
param_common = param_common_fixed + \
   "\n   vector<lower=0, upper=100>[S] DM;    // The true DMs\n"


param_sn_1filt = '''
   real<lower=0, upper=10> eps_sn;   // intrinsic dispersion
   vector<lower=-10, upper=10>[2] EBV_Rv[Nsn];          // E(B-V) and RV
   real<lower=-5, upper=5> zpoff[Nphotsys];   // zero-point offset
   vector<lower=-10, upper=10>[Nbasis] a;       // basis coefficients
'''

param_sn_Nfilt = '''
    real<lower=0, upper=10> eps_sn[Nfilt];    // intrinsic dispersion
    vector<lower=-10, upper=10>[2] EBV_Rv[Nsn];          // E(B-V) and RV
    vector<lower=-5, upper=5>[Nfilt] zpoff[Nphotsys];  // zero-point errors
    vector<lower=-10, upper=10>[Nbasis] a[Nfilt];        // basis coefficients
'''


transparam = '''
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
'''

transparam_cv = '''
transformed parameters {
   real dist[Nsn];    // DM for hosts, distant and cross-validation
   for (i in 1:Nsn) {
      if ( host[i] > 0) {
         dist[i] = DM[host[i]];
      } else {
         if (CVids[i] > 0) {
            dist[i] = DMCV[CVids[i]];
         } else {
            dist[i] = 5.0*log10((1+zhelio[i])/(1+zcmb[i])*300000.0/H0*
                   (zcmb[i]+zcmb[i]*zcmb[i]*0.79))+ 25.0;
         }
      }
   }
}
'''         
         
#########################
## Model blocks
#########################

model_common = '''

model {
   real modl;
   real edist;
   real Al;
   vector[2] mu;
   matrix[2,2] prec;
   real zp;
   real varm;
'''
#model_common = model_common_fixed +\
#      "\n   // The distance modulii\n   DM ~ multi_normal(DMobs, Cov);"


model_sn_1filt = '''
   //SUPERNOVAE
   zp = 16.0;
   for (i in 1:Nsn){
      mu[1] = EBV[i];
      mu[2] = Rv[i];
      prec = ERprec[i];
      EBV_Rv[i] ~ multi_normal_prec(mu, prec);
      Al = A_lamb(Al_coef, EBV_Rv[i][1], EBV_Rv[i][2], Al_order);
      modl = modelSN(dist[i], a, Bs[i], Al, alpha, K[i], M0);
      if (photsys[i] > 0) {
         modl = modl + zpoff[photsys[i]];
      }
      if (host[i] > 0) {
         if ( in_mag == 0) {
            modl = toflux(modl, zp);
            m_sn[i] ~normal(modl, sqrt(square(eps_sn*m_sn[i]) + 
                            + square(m_sn[i]*alpha*sigmaK[i]) +
                            square(em_sn[i]))) ;
         } else {
            m_sn[i] ~normal(modl, sqrt(square(eps_sn) + 
                            + square(alpha*sigmaK[i]) +
                            square(em_sn[i]))) ;
         }
      }
      else { 
         if ( in_mag == 0) {
            modl = toflux(modl, zp);
            edist = 0.000723*pec_vel/zcmb[i]*modl/1.087;
            m_sn[i] ~normal(modl, sqrt(square(eps_sn*m_sn[i]) +
                            square(edist) + square(em_sn[i]) +
                            square(m_sn[i]*alpha*sigmaK[i])));
         } else {
            edist = 0.000723*pec_vel/zcmb[i];
            m_sn[i] ~normal(modl, sqrt(square(eps_sn) + square(edist) + 
                            square(em_sn[i]) + square(alpha*sigmaK[i])));
         }
            
      }
   }
   for ( i in 1:Nphotsys) {
      zpoff[i] ~ normal(0, zperr[i]);
   }
}

generated quantities {
      real g_modl[Nsn];
      real g_Al[Nsn];

      for (i in 1:Nsn) {
         g_Al[i] = A_lamb(Al_coef, EBV_Rv[i][1], EBV_Rv[i][2], Al_order);
         g_modl[i] = modelSN(dist[i], a, Bs[i], g_Al[i], alpha, K[i], M0);
      }
}
'''

model_sn_Nfilt = '''
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
'''

# This next bit is a little trickery that allows us to make pickles of
# STAN models based on an MD5 hash of the model. That way, we can re-use
# models and save the compilation time.
def cache_loc():
   ''' find a good place to store the STAN models.'''
   pd = os.environ.get('PYTHONDATA', None)
   if pd is not None:
      if not os.path.isdir(os.path.join(pd, 'STAN_models')):
         os.mkdir(os.path.join(pd, 'STAN_models'))
      return os.path.join(pd, 'STAN_models')
   else:
      return '.'

def generate_STAN(cfg, cacheloc=None, outfile=None, codefile=None):
   if codefile is None:
      code = generate_STAN_code(cfg)
   else:
      f = open(codefile)
      code = f.read()
      f.close()

   # Look for a pre-existing pickle. This next line maps a STAN model to a 
   #  unique hash code
   code_hash = md5(code.encode('ascii')).hexdigest()
   
   if outfile is not None:
      with open(outfile, 'w') as f:
         f.write(code)

   if cacheloc is None:
      cacheloc = cache_loc()
   cache_fn = 'script_H0-model-{}.pkl'.format(code_hash)
   cache_fn = os.path.join(cacheloc, cache_fn)
   try:
      # This will fail if we have inconsistent STAN versions, e.g.
      sm = pickle.load(open(cache_fn, 'rb'))
   except:
      sm = pystan.StanModel(model_code=code, verbose=True)
      with open(cache_fn, 'wb') as f:
         pickle.dump(sm, f)
   return sm

def generate_init_dict(cfg, data):
   '''Generate a dictionary of initial values to improve convergence
   rate (or at least avoid nan's).'''
   init = {}

   # First off, the fixed parameters
   init['H0'] = random.uniform(50,100)
   if not cfg.model.fixed_DM:
      init['DM'] = random.uniform(20,40, size=data['S'])

   init['a'] = random.uniform(-5,5, size=(data['Nfilt'],data['Nbasis']))
   init['eps_sn'] = random.uniform(0,1, size=len(cfg.data.sn_filt))
   init['pec_vel'] = random.uniform(0,3)
   init['zpoff'] = random.normal(0, data['zperr'])

   # Conditionals
   init['EBV_Rv'] = zeros((data['Nsn'], 2))
   for i in range(data['Nsn']):
      #C = array([[data['e_EBV'][i]**2, data['cov_EBV_Rv'][i]],
      #           [data['cov_EBV_Rv'][i],data['e_Rv'][i]**2]])
      mu = array([data['EBV'][i], data['Rv'][i]])
      init['EBV_Rv'][i,:] = random.multivariate_normal(mu, 
                            linalg.inv(data['ERprec'][i]))
      #init['EBV_Rv'][i,0] = random.normal(data['EBV'][i],data['e_EBV'][i])
      #init['EBV_Rv'][i,1] = random.normal(data['Rv'][i],data['e_Rv'][i])

   if cfg.model.cv is not None:
      init['DMCV'] = random.uniform(20,40, size=data['NCV'])

   if cfg.model.HostMass:
      init['alpha'] = random.uniform(-1,1)

   return init

def generate_STAN_code(cfg):
   # Functions and data
   model = STAN_functions + data_ceph
   model += data_sn
   if not cfg.model.HostMass:
      model += "\n   real alpha;\n"

   model += data_sn_Nfilt
   if cfg.model.cv is not None:
      model += data_cv
   model += "}\n"
   if cfg.model.fixed_DM:
      model += trans_data_ceph

   # Parameters
   if cfg.model.fixed_DM:
      model += param_common_fixed
   else:
      model += param_common
   model += param_sn_Nfilt
   if cfg.model.cv is not None:
      model += "   real DMCV[NCV];\n"
   if cfg.model.HostMass:
      model += "   real <lower=-2, upper=2> alpha;\n"
   model += "}\n"
   if cfg.model.cv:
      model += transparam_cv
   else:
      model += transparam

   # Model part
   model += model_common
   if not cfg.model.fixed_DM:
      model += "\n   // The distance modulii\n   DM ~ multi_normal(DMobs, Cov);"
   model += model_sn_Nfilt

   return model


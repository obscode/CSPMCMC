# This part of the code is responsible for putting together the STAN code.
# I'm outsourcing this to its own module because some of the options were
# impossible to do in the STAN language with if/else blocks. 

from hashlib import md5
import pystan
import pickle
import os
from numpy import random, array
from numpy import zeros

# First, the static stuff.  These are blocks of code that will be assembled
# by the generate_STAN function as needed.
STAN_functions = '''
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
'''

#########################
## Data blocks
#########################
data_ceph = '''
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
'''
data_4258 = '''
   int<lower=1> N4258;          //Number of data in NGC4258
   real P_4258[N4258];   
   real VI_4258[N4258];
   real OH_4258[N4258];
   real mag4258[N4258];
   real e_mag4258[N4258];
'''
data_LMC = '''   
   int <lower=1> NLMC;      // number of data in LMC
   real P_LMC[NLMC];
   real VI_LMC[NLMC];
   real OH_LMC;
   real magLMC[NLMC];
   real e_magLMC[NLMC];
'''
data_MW = '''
   int<lower=1> NMW;        // Number of MW Cepheids
   real pi[NMW];             // parallax in arc-sec
   real e_pi[NMW];       // error in parallax
   real magMW[NMW];
   real e_magMW[NMW];           //error in V magnitude for MW cepheids
   real VI_MW[NMW];
   real P_MW[NMW];
   real OH_MW;
   int LK_prior;            // Apply Lutz-Kelker prior?
'''


#########################
## Parameter blocks
#########################

param_common = '''
parameters {
   real <lower=10.0, upper=40.0> DM[S];           // intercept
   real betaP;           // Period-slope
   real M;               // Absolute mag of Cepheids
   real betaOH;           // slope of metallicity dependence
   real<lower=0> eps;    // dispersion
   real<lower=0.001> lowP;
   real<lower=0.001> highP;
'''

param_MW = '''
   real<lower=0.001> eps_MW;
   real dzp_MW;
   real<lower=0.000125, upper=1.0> pi_true[NMW]; // true parallax, no farther 
                                                 // than 8kpc
'''
param_LMC = '''
   real DM_LMC;
   real<lower=0> eps_LMC;
   real dzp_LMC;
'''
param_4258 = '''
   real DM_4258;
   real<lower=0> eps_4258;
'''   

param_flux = '''
'''


#########################
## Model blocks
#########################


model_common = '''
model {
   real modl;
   real dist;
   real edist;
   real zp;
   real e;

'''

model_hosts_normal = '''
   zp = 25.0;
   for (i in 1:N) {
      modl = modelCeph(M, DM[ID[i]], betaP, betaVI%s, betaOH, P[i], VI[i], OH[i]);
      e = eps + exp(-(P[i]-0.4)/lowP) + exp(-(2.1-P[i])/highP);
      if (in_mag==0){
         modl = toflux(modl, zp);
         mag[i] ~ normal(modl, sqrt(var1 + e*e*mag[i]*mag[i] + e_mag[i]*e_mag[i]));
      }
      else {
         mag[i] ~ normal(modl, sqrt(e*e +e_mag[i]*e_mag[i]));
      }
   }
'''   
model_hosts_lognormal = '''
   zp = 25.0;
   for (i in 1:N) {
      modl = modelCeph(M, DM[ID[i]], betaP, betaVI%s, betaOH, P[i], VI[i], OH[i]);
      e = eps + exp(-(P[i]-0.4)/lowP) + exp(-(2.1-P[i])/highP);
      if (in_mag==0){
         modl = 0.921034*(z - modl);
         mag[i] ~ lognormal(modl, 0.921034*sqrt(e*e + e_mag[i]*e_mag[i]));
      }
      else {
         mag[i] ~ normal(modl, sqrt(e*e +e_mag[i]*e_mag[i]));
      }
   }
'''   

model_MW = '''
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
      modl = modelCeph(M, dist, betaP, %s, betaOH, P_MW[i], VI_MW[i], OH_MW) +
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
'''

model_LMC = '''
   zp = 12.0;
   DM_LMC ~ normal(18.49,0.05);
   dzp_LMC ~ normal(0.0, 0.03);
   for (i in 1:NLMC) { 
      modl = modelCeph(M, DM_LMC, betaP, %s, betaOH, P_LMC[i], 
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
'''

model_4258_normal = '''
   zp = 23.0;
   DM_4258 ~ normal(%f,%f);
   for (i in 1:N4258) {
      modl = modelCeph(M, DM_4258, betaP, %s, betaOH, P_4258[i], VI_4258[i],OH_4258[i]);
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
'''

model_4258_lognormal = '''
   zp = 23.0;
   DM_4258 ~ normal(%f,%f);
   for (i in 1:N4258) {
      modl = modelCeph(M, DM_4258, betaP, %s, betaOH, P_4258[i], VI_4258[i],OH_4258[i]);
      e = eps_4258 + exp(-(P_4258[i]-0.4)/lowP) + exp(-(2.1-P_4258[i])/highP);
      if (in_mag == 0){
         modl = 0.921034*(zp - modl);
         mag4258[i] ~ lognormal(modl, 0.921034*sqrt(e*e + e_mag4258[i]*e_mag4258[i]));
      }
      else {
         mag4258[i] ~ normal(modl,sqrt(e*e+ e_mag4258[i]*e_mag4258[i]));
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

DMs = {'LMC':18.459, '4258':29.344, 'N1365':31.265, 'M101':29.137, 
      'N4639':31.458, 'N5917':32.22, 'N1015':32.436, 'N3370':32.022, 
      'N3447':31.855, 'N1448':31.252, 'U9391':32.852, 'N3982':31.691,
      'N5584':31.723, 'N3972':31.578, 'N4536':30.862, 'N2442':31.540,
      'N4424':31.06, 'N4038':31.24, 'N7250':31.454, 'N3021':32.45, 
      'N1309':32.461} 

def generate_init_dict(cfg, data, cephlist):
   '''Generate a dictionary of initial values to improve convergence
   rate (or at least avoid nan's).'''
   init = {}

   # First off, the fixed parameters
   init['M'] = -2.96 + random.normal(0,0.5)
   init['betaP'] = -3.07 + random.normal(0.1)
   if cfg.model.betaVI_tied:
      init['betaVI'] = 0.3 + random.normal(0,0.01)
   else:
      init['m_betaVI'] = 0.3 + random.normal(0,0.01)
      init['s_betaVI'] = random.normal(0,1)
      init['betaVI'] = 0.3 + random.normal(0, 0.01)
   init['betaOH'] = -0.2 + random.normal(0, 0.1)
   init['eps'] = random.uniform(0, 0.5)

   init['DM'] = array([DMs[name] for name in cephlist]) + \
         random.normal(0, 0.5, size=len(cephlist))

   init['lowP'] = random.uniform(0.001, 1.0)
   init['highP'] = random.uniform(0.001, 1.0)

   # Conditionals
   if not cfg.model.in_mag:
      init['var1'] = random.uniform(0, 0.05)
   if cfg.model.use_MW:
      init['pi_true'] = random.normal(data['pi'], data['e_pi'])
      if not cfg.model.betaVI_tied:
         init['betaVI_MW'] = random.uniform(0,10)
      init['eps_MW'] = random.uniform(0,0.5)
      init['dzp_MW'] = random.uniform(0,0.03)
   if cfg.model.use_LMC:
      init['DM_LMC'] = random.normal(18.49, 0.05)
      init['dzp_LMC'] = random.normal(0, 0.03)
      if not cfg.model.betaVI_tied:
         init['betaVI_LMC'] = 0.3 + random.normal(0,0.1)
      init['eps_LMC'] = random.uniform(0,0.03)
   if cfg.model.use_4258:
      init['DM_4258'] = random.normal(29.40, 0.23)
      if not cfg.model.betaVI_tied:
         init['betaVI_4258'] = 0.3 + random.normal(0,0.1)
      init['eps_4258'] = random.uniform(0,0.5)

   return init

def generate_STAN_code(cfg):
   # Functions and data
   model = STAN_functions + data_ceph
   if cfg.model.use_MW:
      model += data_MW
   if cfg.model.use_LMC:
      model += data_LMC
   if cfg.model.use_4258:
      model += data_4258
   model += "}\n"

   # Parameters
   model += param_common
   if cfg.model.betaVI_tied:
      model += "   real<lower=0> betaVI;\n"
   else:
      model += "   real<lower=0> m_betaVI;\n"
      model += "   real<lower=0> s_betaVI;\n"
      model += "   real<lower=0> betaVI[S];\n"
   if not cfg.model.in_mag:
      model += "real <lower=0> var1;   //additive \n"
   if cfg.model.use_MW:
      model += param_MW
      if not cfg.model.betaVI_tied:
         model += "   real<lower=0> betaVI_MW;"
   if cfg.model.use_LMC:
      model += param_LMC
      if not cfg.model.betaVI_tied:
         model += "   real<lower=0> betaVI_LMC;"
   if cfg.model.use_4258:
      model += param_4258
      if not cfg.model.betaVI_tied:
         model += "   real<lower=0> betaVI_4258;"
   model += "}\n"

   # Model part
   model += model_common
   if cfg.model.in_mag and not cfg.model.lognormal:
      model += "   real var1;\n\n   var1=0.0;\n"
   if cfg.model.betaVI_tied:
      if cfg.model.lognormal:
         model += model_hosts_lognormal % ""
      else:
         model += model_hosts_normal % ""
   else:
      model += "   betaVI ~ normal(m_betaVI, s_betaVI);\n"
      if cfg.model.lognormal:
         model += model_hosts_lognormal % "[ID[i]]"
      else:
         model += model_hosts % "[ID[i]]"
   if cfg.model.use_MW:
      if cfg.model.betaVI_tied:
         model += model_MW % "betaVI"
      else:
         model += "   betaVI_MW ~ normal(m_betaVI, s_betaVI);\n"
         model += model_MW % "betaVI_MW"

   if cfg.model.use_LMC:
      if cfg.model.betaVI_tied:
         model += model_LMC % "betaVI"
      else:
         model += "   betaVI_LMC ~ normal(m_betaVI, s_betaVI);\n"
         model += model_LMC % "betaVI_LMC"

   if cfg.model.use_4258:
      if cfg.model.betaVI_tied:
         if cfg.model.lognormal:
            model += model_4258_lognormal % (cfg.model.mu_4258,cfg.model.sig_4258,"betaVI")
         else:
            model += model_4258_normal % (cfg.model.mu_4258,cfg.model.sig_4258,"betaVI")
      else:
         model += "   betaVI_4258 ~ normal(m_betaVI, s_betaVI);\n"
         if cfg.model.lognormal:
            model += model_4258_lognormal % \
               (cfg.model.mu_4258,cfg.model.sig_4258,"betaVI_4258")
         else:
            model += model_4258_normal % \
               (cfg.model.mu_4258,cfg.model.sig_4258,"betaVI_4258")

   model += "}"
   return model


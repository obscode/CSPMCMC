'''A module that assembles info from the config file and generates
a STAN model string.'''
import variable
import priors
import pystan, pickle
from hashlib import md5
import os
import numpy as np

def cache_loc():
   #find a good place to keep STAN model pickles
   pd = os.environ.get('PYTHONDATA', None)
   if pd is not None:
      if not os.path.isdir(os.path.join(pd,'STAN_models')):
         os.mkdir(os.path.join(pd, 'STAN_models'))
      return os.path.join(pd, 'STAN_models')
   else:
      return '.'

stypes = {
     np.float16:'real',
     np.float32:'real',
     np.float64:'real',
     np.float128:'real',
     np.int8:'int',
     np.int16:'int',
     np.int32:'int',
     np.int64:'int'}

functions = '''
functions {
   real A_lamb(vector Al_coef, real Rv, real EBV, int order) {
         int ii;   // 1D index
         real Al;  // the answer
         ii <- 1;
         Al <- 0;
         for (j in 0:order) {
            for (i in 0:order) {
               if (j+i <= order) {
                  Al <- Al + Al_coef[ii]*pow(Rv,i)*pow(EBV,j+1);
                  ii <- ii+1;
               }
            }
         }
         return Al;
   }

   real expon_log(real x, real x0, real tau) {
      return(-log(tau) + (x0 - x)/tau);
   }

    real NGauss_log(real[] y, vector theta, real[] mu, real[] sigma) {
       int N;
       int K;
       real P;
       real P0;
       N <- size(y);
       K <- dims(theta)[1];
       P <- 0;
       for (n in 1:N) {
          P0 <- 0;
          for (k in 1:K) {
             P0 <- P0 + theta[k]*exp(normal_log(y[n], mu[k], sigma[k]));
          }
          P <- P + log(P0);
       }
       return(P);
    }

    real BGauss_log(real[] y, int[] bins, real[] mu, real[] sigma) {
       int N;
       real P;
       N <- size(y);
       P <- 0;
       for (n in 1:N) {
          P <- P + normal_log(y[n], mu[bins[n]], sigma[bins[n]]);
       }
       return (P);
    }
}

'''
data0 = '''
data {
   int<lower=1> Nf;         // Number of filters
   int<lower=1> NSNe;       // Number of SNe
   int<lower=1> Nobs;       // Number of observations
   int<lower=1> Ncoef;      // Number of coefficients
   int<lower=1> Nknots;     // Number of B-spline knots
   int<lower=1> f0;         // Normalizing filter
   real m[Nobs];            // magnitudes
   real vm[Nobs];           // error in magnitudes
   real st[NSNe];           // s_BV
   real vst[NSNe];          // variance in s_BV
   int findex[Nobs];        // filter index
   int findex0[NSNe];       // index into m[] for normalizing filter
   int sindex[Nobs];        // SN index
   int bindex[NSNe];        // R_V bin
   vector[Ncoef] Amat[Nf];  //  A_lamb matrix
   int Al_order[Nf];        // F99 A_lamb matrix
   vector[Nknots] Bs[NSNe]; // B-spline design matrix
   vector[Nknots] dBs[NSNe]; // B-spline derivatives
'''
data0_labs = ['Nf','NSNe','Nobs','Ncoef','Nknots','f0','m','vm','st','vst','findex',
      'findex0','sindex','bindex','Amat','Al_order','Bs','dBs']

# The following model does not change based on cfg file.
the_model = '''
   for (i in 1:Nobs) {
      if (findex[i] != f0) {
         A0 <- A_lamb(Amat[f0], R_V[sindex[i]], EBV[sindex[i]], Al_order[f0]);
         Al <- A_lamb(Amat[findex[i]], R_V[sindex[i]], EBV[sindex[i]], 
                      Al_order[findex[i]]);
         mod <- m[findex0[sindex[i]]] - 
                dot_product(Bs[sindex[i]],a[findex[i]-1])
                + Al - A0;
         m[i] ~ normal(mod, sqrt(vm[i] + vm[findex0[sindex[i]]] + 
                  square(dot_product(dBs[sindex[i]],a[findex[i]-1]))*vst[sindex[i]] +
                  evar[findex[i]-1]));
      }
   }
'''

def prior_to_string(prior):
   if isinstance(prior, priors.Uniform):
      # dealt with in bounds
      pstr = None
   elif isinstance(prior, priors.Normal):
      pstr = "normal(%f,%f)" % (prior.mu,prior.sigma)
   elif isinstance(prior, priors.Exponential):
      pstr = "expon(%f,%f)" % (prior.mu,prior.tau);
   elif isinstance(prior, priors.InverseGamma):
      pstr = "inv_gamma(%f,%f)" % (prior.alpha,prior.beta)
   elif isinstance(prior, priors.Dirichlet):
      pstr = "dirichlet(["
      for i in range(len(prior.alpha)-1):
         pstr += "%f," % prior.alpha[i]
      pstr += "%f])" % prior.alpha[-1]
   else:
      pstr = None
   return pstr

def PriorString(var):
   '''Take a prior instance and map it to a STAN string.'''
   if isinstance(var.prior, priors.Uniform) \
         or isinstance(var.prior,priors.Uninformative) \
         or isinstance(var.prior,priors.Dirichlet):
      return ""
   if isinstance(var.prior, priors.STAN):
      return "   "+var.prior.string+"\n"
   if var.shape == 0:
      vector = False
   elif var.name in ['EBV']:
      vector = True
      upper = 'NSNe'
   else:
      vector = True
      upper = 'Nf'

   if not vector:
      return "   %s ~ %s" % (var.name, prior_to_string(var.prior))

   # Dealing with a vector type
   if type(var.prior) is type([]) or var.cindex is not None:
      if type(var.prior) is type([]):
         pstr = [prior_to_string(p) for p in var.prior]
      else:
         pstr = [prior_to_string(var.prior)]*var.shape[0]
      if var.cindex is not None:
         ii = 0
         for id,val in zip(var.cindex,var.cvalue):
            pstr[id] = "normal(%f, %f)" % (val, max(0.001, abs(val)*0.001))
      retstr = ""
      for i,pstr in enumerate(pstr):
         if pstr is not None:
            retstr += "   %s[%d] ~ %s;\n" % (var.name, i+1, pstr)
   else:
      pstr = prior_to_string(var.prior)
      if pstr is not None:
         retstr = "   for (i in 1:%s) %s[i] ~ %s;\n" % (upper, var.name, pstr)
      else:
         retstr = ""
   return retstr


def generate_stan_code(vinfo, idata):

   # Given the vinfo, configuration and data, build the STAN
   # model

   # First the parameters and Data sections
   par = 'parameters {\n   vector[Nknots] a[Nf-1];\n'
   data = data0*1
   varies = []
   for var in vinfo:
      if var.name == 'a': continue
      if var.vary:
         varies.append(var)
         # Add to the parameters
         done = False
         if isinstance(var.prior, priors.Uniform):
            line = "   real <lower=%f,upper=%f>%s" % \
                   (var.prior.lower,var.prior.upper,var.name)
         elif var.lower is not None and var.upper is not None:
            line = "   real <lower=%f,upper=%f>%s" % \
                   (var.lower,var.upper,var.name)
         elif var.lower is not None:
            line = "   real <lower=%f>%s" % \
                   (var.lower,var.name)
         elif var.upper is not None:
            line = "   real <upper=%f>%s" % \
                   (var.upper,var.name)
         elif isinstance(var.prior, priors.Exponential):
            line = "   real <lower=0>%s" % (var.name)
         elif var.positive:
            line = "   real <lower=0>%s" % (var.name)
         elif isinstance(var.prior, priors.Dirichlet):
            line = "   simplex[%d] %s" % (var.shape[0], var.name)
            done = True
         else:
            line = "   real %s" % (var.name)

         if var.shape == 0 or done:
            line += ";\n"
         elif var.name in ['EBV']:
            line += "[NSNe];\n"
         elif var.name in ['R_V']:
            if var.shape == 0:
               line += ";\n"
            else:
               line += "[NSNe];\n"
         elif var.name in ['R0','eR0']:
            line += "[%d];\n" % (var.shape[0])
         else:
            line += "[Nf-1];\n"
         par += line
      else:
         # Held constant, so put it in the data
         if var.shape == 0:
            data += "   real %s;\n" % (var.name)
         elif var.name in ['EBV','Rv']:
            data += "   real %s[NSNe];\n" % (var.name)
         elif var.name in ['theta']:
            data += "   simplex[%d] %s;\n" % (var.shape, var.name)
         else:
            data += "   real %s[Nf-1];\n" % (var.name)
         idata[var.name] = var.value
            
   # close off data and parameters
   par += "}\n"
   # Add any extra data in
   varnames = [var.name for var in vinfo]
   for key in idata:
      if key not in data0_labs and key not in varnames:
         var = idata[key]
         if type(var) is np.ndarray:
            tp = stypes[var.dtype.type]
            shp = '['+','.join(map(str, var.shape))+'];\n'
            data += "   " + tp  +" " + key + shp
         else:
            if type(var) is type(0):
               data += "   int %s;\n" % (key)
            else:
               data += "   real %s;\n" % (key)

   data += "}\n"
   
   # now we work on the model section
   model = "model {\n   real mod;\n   real Al;\n   real A0;\n"
   # Start with priors
   for var in varies:
      model += PriorString(var)

   model += the_model + "\n}\n"
   model += '''
\n\ngenerated quantities {
   real mcorr[Nobs];

   for (i in 1:Nobs) {
      mcorr[i] <- m[i] - A_lamb(Amat[findex[i]], R_V[sindex[i]], EBV[sindex[i]],
                                Al_order[findex[i]]);
   }
}\n'''

   return (functions + data + par + model),idata

def generate_stan(vinfo, idata, cacheloc=None, outfile=None, stanfile=None, **args):

   # get the code and dict
   code,idata = generate_stan_code(vinfo, idata)
   # Now look for pickle
   if stanfile is not None:
      fin = open(stanfile, 'r')
      code = fin.read()
      fin.close()
   code_hash = md5(code.encode('ascii')).hexdigest()

   if outfile is not None:
      with open(outfile, 'w') as f:
         f.write(code)

   if cacheloc is None:
      cacheloc = cache_loc()
   cache_fn = 'color_edge-model-{}.pkl'.format(code_hash)
   cache_fn = os.path.join(cacheloc,cache_fn)

   try:
      sm = pickle.load(open(cache_fn, 'rb'))
   except:
      sm = pystan.StanModel(model_code=code)
      with open(cache_fn, 'wb') as f:
         pickle.dump(sm, f)
   return sm.sampling(data=idata,**args)

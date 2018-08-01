#!/usr/bin/env pyhton

'''Script that runs EBV in parallel using STAN.'''

import sys,os,string
from numpy import *
import pickle
import time
import priors
import pystan
import config
#import LMCMCMC
import EBV
import make_stan_model
import bspline2
from scipy.interpolate import splrep,splev

if len(sys.argv) < 2:
   print "Usage: EBV_stan config-file [stan-file]"
   sys.exit(1)

stanfile = None
if len(sys.argv) > 2:
   stanfile = sys.argv[2]
   if not os.path.isfile(stanfile):  stanfile = None


if not os.path.isfile(sys.argv[1]):
   print "Error, can't find config-file ",sys.argv[1]
   sys.exit(1)

#nchains=int(sys.argv[2])

# poly coefficients from CCM+O'Donnel 
f = open('Ia_A_poly.pickle')
d = pickle.load(f)
f.close()

t0 = time.time()
cf = EBV.get_config(sys.argv[1])
data = EBV.get_data(cf)
vinfo = EBV.setup_varinfo(cf, data)
nchains = cf.Sampler.chains

coefs = d[cf.Model.redlaw]

Amat = array([coefs[f] for f in data.filters])
Aorder = array([coefs['order'][f] for f in data.filters])
Ncoef = Amat.shape[1]
Nf = len(data.filters)

data['knots'] = array(cf.Model.knots)
#data['knots'] = [data.st.min()]*4 + [median(data.st)] + [data.st.max()*1.001]*4
#bs = bspline.Bspline(knot_vector=data.knots, order=3)
bs = bspline2.bspline_basis(data.knots, data.st, 3, gradient=cf.Model.gradient)
#data['Bs'] = array([bs(st) for st in data.st])
data['Bs'] = bs
print data['Bs'].shape
# first derivative of splines (for computing errors)
sts = linspace(data.st.min(),data.st.max(),100)
bss = bspline2.bspline_basis(data['knots'], sts, 3, gradient=cf.Model.gradient)
tcks = [splrep(sts, bss[:,i], s=0, k=3) for i in range(bss.shape[1])]
data['dBs'] = array([splev(data.st, tck, der=1) for tck in tcks]).T
print data['dBs'].shape
#data['dBs'] = array([bs.d(st) for st in data.st])

data_in = dict(
      Nf = Nf,
      NSNe = data.Nobj,
      Nobs = len(data.ms),
      Ncoef = Ncoef,
      Nknots = data.Bs.shape[1],
      Bs = data.Bs,
      dBs = data.dBs,
      m = data.ms,
      vm = data.vms,
      st = data.st,
      vst = data.vst,# + 0.06**2,   # RMS of dm15-s_BV relation
      findex = data.fids+1,
      sindex = data.oids+1,
      bindex = data.bids+1,
      f0 = data.f0 + 1,
      findex0 = data.findex0 + 1,
      Amat = Amat,
      Al_order = Aorder)

if 'extra_data' in cf.sections.keys():
   for key in cf.extra_data.options:
      val = getattr(cf.extra_data, key)
      if shape(val) > 0:
         val = array(val)
      data_in[key] = val


# Initial guesses
init = []
for i in range(nchains):
   init.append({})
   for var in vinfo.varnames:
      if vinfo[var].vary:
         #if var == 'a':
         #   val = random.uniform(-2,2, size=Nf-1)
         #elif var == 'b':
         #   val = random.uniform(-1,1, size=Nf-1)
         #elif var == 'c':
         #   val = random.uniform(-1,1, size=Nf-1)
         if var == 'a':
            val = random.uniform(-1,1, size=(Nf-1,data.Bs.shape[1]))
         elif var == 'R_V':
            if cf.Priors.Rv_global:
               val = random.uniform(1.0,4.0)
            else:
               val = random.uniform(1.0,4.0, size=data.Nobj)
         elif var == 'evar':
            val = random.uniform(0.0001, 0.25, size=Nf-1)
         elif var == 'EBV':
            val = random.exponential(0.2, size=data.Nobj)
         elif var == 'tau':
            val = random.uniform(0.01, 1)
         elif var == 'R0':
            if cf.Priors.Rv_binned:
               val = random.uniform(1.5,5.0, size=len(cf.Priors.Rv_bins))
            else:
               val = random.uniform(1.5,5.0, size=cf.Priors.NGauss)
         elif var == 'eR0':
            if cf.Priors.Rv_binned:
               val = random.uniform(0, 2.0, size=len(cf.Priors.Rv_bins))
            else:
               val = random.uniform(0, 2.0, size=cf.Priors.NGauss)
         elif var == 'theta':
            val = random.dirichlet(ones(cf.Priors.NGauss,))
            print val
         elif var == 'muR':
            val = random.uniform(1.0,4.0)
         elif var == 'sigR':
            val = random.uniform(0.0,2.0)
         else:
            raise ValueError, "Unknown var %s" % var
         init[-1][var] = val

fit = make_stan_model.generate_stan(vinfo, data_in, outfile="test.stan",
      iter=cf.Sampler.N_final, warmup=cf.Sampler.burn_final, chains=nchains,
      init=init, stanfile=stanfile)

print "Finished sampling"

samples = fit.extract()
#d = dict(samples=samples,
#         data=data)
samples['data'] = data
fout = open(cf.Sampler.outfile, 'w')
pickle.dump(samples, fout)
fout.close()


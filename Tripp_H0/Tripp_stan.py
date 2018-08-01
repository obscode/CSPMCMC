#!/usr/bin/env pyhton

'''Script that runs Tripp model in parallel using STAN.'''

import sys,os,string
from numpy import *
import pickle
import time
from hashlib import md5
import config
import get_data
import pystan

def get_stan_model(stanfile, cacheloc=None, **args):
   f = open(stanfile)
   code = f.read()
   f.close()
   code_hash = md5(code.encode('ascii')).hexdigest()

   if cacheloc is None:
      pd = os.environ.get('PYTHONDATA', None)
      if pd is not None:
         if not os.path.isdir(os.path.join(pd, 'STAN_models')):
            os.mkdir(os.path.join(pd, 'STAN_models'))
         cacheloc = os.path.join(pd, 'STAN_models')
      else:
         cacheloc = '.'
   cache_fn = 'Trip_H0-model-{}.pkl'.format(code_hash)
   cache_fn = os.path.join(cacheloc,cache_fn)

   try:
      sm = pickle.load(open(cache_fn, 'rb'))
   except:
      sm = pystan.StanModel(model_code=code)
      with open(cache_fn, 'wb') as f:
         pickle.dump(sm, f)
   return sm


t0 = time.time()
cf = config.config(sys.argv[1])
data = get_data.get_data(cf)
nchains = cf.Sampler.chains

needed_data = ['C','Ccov','DMCeph','NSNe','host','NObs','NFs','ms','vms',
      'cs','vcs','st','est','zcmb','zhel','sindex','findex','findex1',
      'findex2','Nphotsys','photsys','zperr']
if cf.Model.Rv_global:
   needed_data += ['Amat','N_coef','Al_order','NFs2']
if cf.Model.HostMass:
   needed_data.append('K')
   needed_data.append('M0')
   needed_data.append('sigmaK')

data_in = {}
for key in needed_data:
   data_in[key] = data[key]
data_extra = {}
for key in data.keys():
   if key not in data_in:
      data_extra[key] = data[key]

# Initial guesses
init = []
for i in range(nchains):
   init.append({
      'a':random.uniform(-2,2, size=data.NFs),
      'b':random.uniform(-1,1, size=data.NFs),
      'c':random.uniform(-1,1, size=data.NFs),
      'evar':random.uniform(0.0001, 0.25, size=data.NFs),
      'vpec':random.uniform(0,1),
      'H0':random.uniform(50,100),
      'ct':random.normal(data_in['cs'], sqrt(data_in['vcs'])),
      'zpoff':random.normal(0, 0.1, size=(data_in['Nphotsys'],data_in['NFs'])),
      'DMhost':random.normal(data_in['DMCeph'], sqrt(diag(data_in['Ccov'])))
   })
   if cf.Model.Rv_global:
      init[-1]['R'] = random.uniform(1.0,4.0)
   else:
      init[-1]['Rl'] = random.uniform(0.0,7.0, size=data.NFs)
   if cf.Model.HostMass:
      init[-1]['alpha'] = random.uniform(-1,1)

if cf.Model.Rv_global:
   if cf.Model.fixed_DM:
      model_file = 'model_R_tied_fixed.stan'
   else:
      model_file = 'model_R_tied.stan'
else:
   if cf.Model.fixed_DM:
      model_file = 'model_R_free_fixed.stan'
   else:
      model_file = 'model_R_free.stan'
if cf.Model.HostMass:
   model_file = model_file.replace('.stan','_HM.stan')

sm = get_stan_model(model_file)
fit = sm.sampling(data=data_in, iter=cf.Sampler.N_final, 
   warmup=cf.Sampler.burn_final, chains=nchains, init=init)

if cf.Sampler.out is not None:
   f = open(cf.Sampler.out, 'w')
else:
   f = open(['sampler.out','sampler_fixed.out'][cf.Model.fixed_DM], 'w')
f.write(str(fit))
f.close()

print "Finished sampling"

samples = fit.extract(permuted=False)
d = dict(data=data_in, extra=data_extra, samples=samples, 
      flatnames=fit.flatnames)

fout = open(cf.Sampler.outfile, 'w')
pickle.dump(d, fout)
fout.close()


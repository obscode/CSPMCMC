import numpy as np
from numpy import*
from astropy.io import ascii
from astropy.table import Table , Column
import pystan
import pickle
import config
import sys , os , string
import generate_STAN
import get_data

# go in /data/Cepheids/runs
# creates a pickle file and a table 

cfg = config.config(sys.argv[1])
if len(sys.argv) > 2:
   codefile = sys.argv[2]
else:
   codefile = None

# PYSTAN MODEL
model = generate_STAN.generate_STAN(cfg, outfile='model.stan', codefile=codefile)

# Data
dat,extras = get_data.get_data(cfg)
if not cfg.model.HostMass:
   dat['alpha'] = 0

# Initial guess for parameters
samplefile = cfg.sampler.sample0
if samplefile is not None:
   import STANstats
   c = STANstats.STANchains(samplefile)
   d = []
   for i in range(cfg.sampler.chains):
      d0 = generate_STAN.generate_init_dict(cfg, dat)
      d.append(d0)
      for key in d0:
         if d0 in c.params:
            d[-1][key] = random.normal(c.median(key), c.std(key))
else:
   d = [generate_STAN.generate_init_dict(cfg,dat) \
         for i in range(cfg.sampler.chains)]

if __name__ == "__main__":
   
   #___________________________________________________________________________
   #  FIT 
   fit2 = model.sampling(data=dat, iter=cfg.sampler.iter,
         warmup=cfg.sampler.burn, chains=cfg.sampler.chains,
         init=d)
   
   
   if cfg.sampler.summary is not None:
      of = open(cfg.sampler.summary, 'w')
   else:
      of = open('summary.txt', 'w')
   try:
      of.write(str(fit2))
   except:
      of.write('printing summary failed, likely non-converged chains\n')
   of.close()
   #________________________________________
   #Make pickle file
   
   
   samples = fit2.extract(permuted=cfg.sampler.permuted)
   if cfg.sampler.output is not None:
      filename = cfg.sampler.output
      # Now we can add extra data to dat before saving to pickle file.
      for key in extras:
         dat[key] = extras[key]
      if not cfg.sampler.permuted:
         d = dict(data=dat, samples=samples, flatnames=fit2.flatnames)
      else:
         d = samples
         d['data'] = dat
      fout = open(filename, 'w')
      pickle.dump(d, fout)
      fout.close()

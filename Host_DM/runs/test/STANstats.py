'''Pymc statistics module.
Author:  Chris Burns
         Carnegie Observatories, cburns@obs.carnegiescience.edu

Contents:  This module has a class, chains, that reads in chains produced by
           pystan, assigns the chains to easily-accessible member variables,
           allows for doing statistics and plotting of histograms, confidence
           regions, etc.
Requires:  - pickle file with samples from STAN
           - NumPy (developed under version 1.6.0)
           - scipy (developed under verison 0.9.0)
           - matplotlib (developed under version 1.0.1)
           - pymc (if you want gp interpolation)
'''

from numpy import *
from glob import glob
import os,string,copy,pickle, re
from MCMCstats import MCMCchains

vec_pat = re.compile(r'(.*)\[([0-9,]+)\]')

def parse_flatnames(flatnames):
   params = []
   slices = {}
   shapes = {}
   for i,var in enumerate(flatnames):
      var = str(var)
      res = vec_pat.search(var)
      if res:
         # The parameter has an index, possibly more than one
         vname,ind = res.groups()
         inds = map(int, ind.split(','))
         
         if vname not in params:
            params.append(vname)
            slices[vname] = slice(i,i+1)
            shapes[vname] = [ind+1 for ind in inds[::-1]]
         else:
            slices[vname] = slice(slices[vname].start, i+1)
            for i,ind in enumerate(inds[::-1]):
               if ind+1 > shapes[vname][i]:
                  shapes[vname][i] = ind+1
      else:
         params.append(var)
         slices[var] = slice(i,i+1)
   for vname in shapes:
      shapes[vname] = tuple(shapes[vname])

   return params,slices,shapes

def parse_vinfo(vinfo):

   params = [var.name for var in vinfo.vars if var.vary]
   slices = {}
   shapes = {}
   for param in params:
      if vinfo[param].shape == 0:
         slices[param] = slice(vinfo[param].slice,
                                    vinfo[param].slice+1)
      else:
         slices[param] = vinfo[param].slice
         shapes[param] = vinfo[param].shape
   return params,slices,shapes

class STANchains(MCMCchains):

   def __init__(self, filename=None, d=None, chains=None, flatnames=None, 
         vinfo=None, nchains=1,verbose=False):

      self.d = d
      self.nchains = nchains
      MCMCchains.__init__(self,filename, chains, verbose, flatnames=flatnames,
            vinfo=vinfo)

   def load_chains(self, filename=None, chains=None, flatnames=None, vinfo=None):
      if filename is None and chains is None and self.d is None:
         raise ValueError, "One of filename or chains  or d must be defined"
      if chains is not None:
         Ndim = len(shape(chains))
         if not 0 < Ndim < 4:
            raise ValueError, "chains should be 1, 2, or 3D"
         elif Ndim == 1:
            # assume one variable with iterations
            self.chains = self.chains.reshape((1,self.chains.shape[0],1))
         elif Ndim == 2:
            # assume only one chain
            self.chains = self.chains.reshape((1,self.chains.shape[0],
                                                self.chains.shape[1]))
         self.chains = transpose(chains, (1,0,2))
         if flatnames is not None:
            self.params,self.slices,self.shapes = parse_flatnames(flatnames)
         elif vinfo is not None:
            self.params,self.slices = parse_vinfo(vinfo)
         else:
            self.params = ['p%d' % i for i in range(self.chains.shape[2])]
            for i,par in enumarate(self.params):
               self.slices[par] = slice(i,i+1)
            
      else:
         if filename is not None:
            if not os.path.isfile(filename):
               raise IOError, "Error, cannot find file %s" % filename
            f = open(filename)
            d = pickle.load(f)
            f.close()
         else:
            d = self.d
         if 'lp__' in d:
            # Saved with permuted=True, so just one chain
            self.params = [var for var in d if var not in ['lp__','data']]
            self.chains = []
            self.slices = {}
            self.shapes = {}
            i = 0
            for var in self.params:
               if len(d[var].shape) == 1:
                  self.chains.append(d[var].reshape((-1,1)))
                  self.slices[var] = slice(i, i+1)
                  self.shapes[var] = (1,)
                  i += 1
               elif len(d[var].shape) == 2:
                  self.chains.append(d[var])
                  size = d[var].shape[1]
                  self.slices[var] = slice(i,i+size)
                  self.shapes[var] = (size,)
                  i += size
               else:
                  self.chains.append(
                      reshape(d[var],(-1,prod(d[var].shape[1:]))))
                  size = prod(d[var].shape[1:])
                  self.slices[var] = slice(i, i+size)
                  self.shapes[var] = d[var].shape[1:]
                  i += prod(d[var].shape[1:])
            self.chains = hstack(self.chains)
            self.chains.shape = (self.nchains,self.chains.shape[0]/self.nchains,
                                 self.chains.shape[1])
         elif 'samples' in d:
            # saved with permuted = False
            self.chains = transpose(d['samples'], axes=(1,0,2))
            if 'flatnames' in d:
               self.params,self.slices,self.shapes = parse_flatnames(d['flatnames'])
            elif 'vinfo' in d:
               self.params,self.slices,self.shapes = parse_vinfo(d['vinfo'])
            else:
               self.params = ['p%d' % i for i in range(self.chains.shape[2])]
               for i,par in enumerate(self.params):
                  self.slices[par] = slice(i,i+1)
         if 'data' in d:
            self.data = d['data']
         if 'extra' in d:
            self.extra = d['extra']

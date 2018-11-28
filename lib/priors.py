'''A module that contains a prior class that can be used with emcee
(it has to be pickleable)'''

import numpy as np
from scipy.special import gamma as Gamma
import os

class Uninformative:
   '''Uninformative prior.  Trivial.'''
   def __call__(self, value):
      return 0.0
class STAN:
   '''A prior encoded in STAN language. Only works with STAN, obviously.'''
   def __init__(self, string):
      self.string = string

   def __call__(self, value):
      raise RuntimeError, "This is a STAN prior and should only be used with pystan"

class Uniform:
   '''Uniform (uninformative over fixed range) prior.'''
   def __init__(self, lower, upper):
      self.lower = lower
      self.upper = upper

   def __call__(self, value):
      if len(np.shape(value)) == 0:
         value = np.atleast_1d(value)

      if np.sometrue(np.less(value,self.lower)) or \
            np.sometrue(np.greater(value, self.upper)):
         return -np.inf

      return 0.0

class Normal:
   '''Normal distribution with mean \mu and std deviation \sigma'''
   def __init__(self, mu, sigma):
      self.mu = mu
      self.sigma = sigma
      self.var = np.power(sigma,2)

   def __call__(self, value):
      return np.sum(-0.5*(np.power(value-self.mu,2)/self.var+np.log(self.var)+\
            np.log(2*np.pi)))

class Exponential:
   '''Exponentially declining prior with position parameter \mu
   and scale \tau.'''
   def __init__(self, mu, tau):
      self.mu = mu
      self.tau = tau
      self.beta = np.power(self.tau,-1)
      self.lbeta = np.log(self.beta)

   def __call__(self, value):
      if np.sometrue(np.less(value, 0)):
         return -np.inf
      return np.sum(self.lbeta - self.beta*(value-self.mu))

class InverseGamma:
   '''Inverse Gamma distribution.  \Gamma(\alpha,\beta).'''
   def __init__(self, alpha, beta):
      self.alpha = alpha
      self.beta = beta
      self.const = np.log(np.power(beta, alpha)/Gamma(alpha))

   def __call__(self, value):
      if np.sometrue(np.less(value, 0)):
         return -np.inf
      return np.sum(self.const - (self.alpha+1)*np.log(value) - self.beta/value)

class Cauchy:
   '''Cauchy distribution.'''
   def __init__(self, mu, tau):
      self.mu = mu
      self.tau = tau

   def __call__(self, value):
      return np.sum(np.log(np.pi*self.tau) - \
            np.log(1 + np.power((value-self.mu)/self.tau,2)))

class HalfCauchy:
   '''Half-Cauchy distribution.'''
   def __init__(self, mu, tau):
      self.mu = mu
      self.tau = tau

   def __call__(self, value):
      if sometrue(np.less(value, mu)):
         return -np.inf
      return np.sum(np.log(np.pi*self.tau) - \
            np.log(1 + np.power((value-self.mu)/self.tau,2)))

class Dirichlet:
   
   def __init__(self, alpha):
      self.alpha = alpha
      self.dim = len(alpha)
      self.const = np.log(Gamma(np.sum(alpha))) - np.sum(np.log(Gamma(alpha)))

   def __call__(self, value):
      if sometrue(np.less_equal(value, 0)):
         return -np.inf
      if sometrue(np.greater_equal(value, 1)):
         return -np.inf
      return self.const + np.sum((self.alpha-1)*np.log(value))






class MWZP:
   '''Prior based on the Milky-Way zero-points from Fouke et al. (2007) and
   from Monson et al. (2012)'''
   ZP_MW = {'B':(-3.225, 0.207), 'V':(-3.953,0.173), 'R':(-4.405,0.180),
        'I':(-4.706, 0.168), 'J':(-5.258,0.155), 'H':(-5.543,0.146),
        'K':(-5.647, 0.144),'[3.6]':(-5.8,0.10), '[4.5]':(-5.77,0.10)}

   def __init__(self, filters):
      self.filters = filters
      self.shape = len(filters)
      self.mus = np.array([self.ZP_MW[filt][0] for filt in self.filters])
      self.sigs = np.array([self.ZP_MW[filt][1] for filt in self.filters])
      self.vars = np.power(self.sigs,2)

   def __call__(self, value):
      if len(np.shape(value)) != 1:
         raise ValueError, "value must be array-like"
      if value.shape[0] != self.shape:
         raise ValueError, "value must have shape (%d,)" % self.shape

      return np.sum(-0.5*(np.power(value-self.mus,2)/self.vars+np.log(self.vars)+\
            np.log(2*np.pi)))

'''Here follows some utility functions.'''
str2prior = {'U':Uniform,
          'N':Normal,
          'E':Exponential,
          'IG':InverseGamma,
          'Dir':Dirichlet}


def get_value(seq, data):
   '''parse a string and extract values.'''
   try:
      # try integers first
      args = map(int, seq)
   except:
      try:
         # next try floats
         args = map(float, seq)
      except:
         # okay, they're strings
         args = []
         for item in seq:
            if item in data.keys():
               args.append(data[item])
            elif os.path.isfile(item):
               try:
                  vals = loadtxt(item, unpack=True)
                  if len(vals.shape) == 1:
                     args.append(vals)
                  else:
                     args += list(vals)
               except:
                  raise ValueError, "failed to load values from %s" % item
            else:
               args.append(item)
   return np.asarray(args)

def get_prior(cf, res, data):
   '''Convert a string from configuration file into prior class.'''
   p = res[0]
   if p == 'C':
      vary = False
      prior = None
      value = get_value(res[1:], data)
   elif p.lower() == 'stan':
      vary = True
      prior = STAN(",".join(res[1:]))
      value = None
   elif p.lower() == 'dir':
      vary = True
      value = get_value(res[1:], data)
      prior = Dirichlet(value)
   else:
      if p not in str2prior:
         raise ValueError, "Unknown prior string %s" % p
      vary = True
      prior = str2prior[p](*get_value(res[1:], data))
      value = None

   return prior,vary,value



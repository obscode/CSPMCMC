#!/usr/bin/env pyhton

import sys,os,string
from scipy.optimize import minimize
from numpy import *
import config
from get_data import get_data
import priors
import variable
import re
import pickle
#import Variable

# a pattern that matches item assignment
assign_pat = r'\[([^\}]+)\]'

def get_config(file):
   '''Get a configuration file and create configuration object.'''
   cf = config.config(file)
   return cf

def setup_varinfo(cf,data):
   '''Given a configuration file, setup information about which variables
   are free or fixed and the slices into the parameter list.'''
   Nf = len(cf.data.filters)
   No = data.Nobj
   varinfo = variable.VarInfo()
   #i = 0
   vars = ['a','R_V','EBV','evar','tau','s0']
   for var in ['R0','eR0','theta','muR','sigR']:
      if getattr(cf.Priors,var): vars.append(var)

   for var in vars:
      res = getattr(cf.Priors, var,None)
      if res is None:
         # Assume a completely uninformative prior
         base_prior = priors.Uninformative()
         vary = True
         ids = []
      else:
         ids = []
         vals = []
         # Scan the options for direct assignment:
         for opt in cf.Priors.options:
            s = re.search(var+assign_pat, opt)
            if s is not None:
               id = s.group(1)
               if re.match('[\+\-]?\d+', id):
                  ids.append(int(id))
               elif id in cf.data.filters:
                  ids.append(cf.data.filters.index(id))
               else:
                  raise ValueError, "config: Invalid index on variable %s" % var
               vals.append(cf.Priors[opt])

         base_prior,vary,value = priors.get_prior(cf, res, data)

      if var in ['tau','muR','sigR','s0']:
         # variables that are scalars
         if vary:
            value = 0
         else:
            value = value[0]
      elif var in ['R_V']:
         if cf.Priors.Rv_global:
            if vary:
               value = 0
            else:
               value = value[0]
         else:
            if vary:
               value = zeros((No,))
            else:
               if value.shape[0] == 1:
                  value = value*ones((No,))
               elif value.shape[0] != No:
                  raise ValueError, "Parameter %s needs to have shape (%d)" % \
                        (var,No)

      elif var in ['EBV']:
         # variables that are indexed by object
         if vary:
            value = zeros((No,))
         else:
            if value.shape[0] == 1:
               value = value*ones((No,))
            elif value.shape[0] != No:
               raise ValueError, "Parameter %s needs to have shape (%d)" % \
                     (var,No)
      elif var in ['R0','eR0','theta']:
         # Variables that have dynamic size based on Ngauss, Rv_bins, etc
         if vary:
            if cf.Priors.Rv_binned:
               value = zeros((len(cf.Priors.Rv_bins),))
            else:
               value = zeros((cf.Priors.NGauss,))+1.0/cf.Priors.NGauss
         else:
            if cf.Priors.Rv_binned:
               if value.shape[0] != len(cf.Priors.Rv_bins):
                  raise ValueError, "Parameter %s needs to have shape (%d)" %\
                        (var,len(cf.Priors.Rv_bins))
            else:
               if value.shape[0] != cf.Priors.NGauss:
                  raise ValueError, "Parameter %s needs to have shape (%d)" %\
                        (var,cf.Priors.NGauss)
      elif var in ['a']:
         # An array of vectors:  shape should be (No,Nf)
         if vary:
            value = zeros((No,Nf-1))
         else:
            if value.shape[0] == 1:
               value = value*ones((No,Nf-1))
            else:
               raise ValueError, "I don't know what I'm doing"

      else:
         # variables that are indexed by filter
         if vary:
            value = zeros((Nf,))
         else:
            if value.shape[0] == 1:
               value = value*ones((Nf,))
            elif value.shape[0] != Nf:
               raise ValueError, "Parameter %s needs to have shape (%d)" % \
                     (var,Nf)
      v = variable.Variable(var, value, vary, base_prior)
      if (var == 'EBV' and cf.Priors.EBVpos) \
         or (var in ['R_V','R0','muR'] and cf.Priors.R_Vpos) \
         or var in ['evar','eR0','sigR']:
         v.positive=True
      else:
         v.positive=False

      if ids and len(shape(value)) > 0:
         prior = []
         cindex = []
         cvalue = []
         for i in range(len(ravel(value))):
            if i in ids:
               p,vary,value = priors.get_prior(cf, vals[ids.index(i)],data)
               if vary:
                  prior.append(p)
               else:
                  cindex.append(i)
                  cvalue.append(value[0])
                  prior.append(p)
            else:
               prior.append(base_prior)
         if cindex:
            v.cindex = array(cindex)
            v.cvalue = array(cvalue)
         v.prior = prior

      varinfo.add_variable(v)

   # Do some checking
   if cf.Priors.EBVpos:
      varinfo.EBV.positive=True
   varinfo.R_V.lower = cf.Priors.Rv_lower
   varinfo.R_V.upper = cf.Priors.Rv_upper

   return varinfo

def guess(d, varinfo):
   '''provides a guess of the initial point to start the MCMC, given the 
   data object.'''
   if varinfo.M_l.vary:
      varinfo.M_l.value = array([-19 for f in d.filters])
   if varinfo.b.vary:
      varinfo.b.value = array([-0.5 for f in d.filters])
   if varinfo.R_V.vary: varinfo.R_V.value = 2.0
   if varinfo.EBV.vary: varinfo.EBV.value = array([0.0]*d.Nobj)
   if varinfo.s0.ary: varinfo.s0.value = 1.0
   p0 = varinfo.getParVector()
   # Now find a-priori maximum
   nll = lambda *args: -lnprob(*args)
   result = minimize(nll, p0, args=(d,varinfo), method='Nelder-Mead')
   if varinfo.EBV0.vary: 
      result['x'][varinfo.EBV0.slice] = \
            max(result['x'][varinfo.EBV0.slice], 0.01)
   if varinfo.EBV.vary: 
      result['x'][varinfo.EBV.slice] = \
            where(result['x'][varinfo.EBV.slice] < 0.01, 0.01,
                  result['x'][varinfo.EBV.slice])
   return result['x']


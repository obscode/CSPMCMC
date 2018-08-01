#!/usr/bin/env pyhton

import sys,os,string
from numpy import *
import config
import priors
import variable
import re
import get_data as gd

# a pattern that matches item assignment
assign_pat = r'\[([^\}]+)\]'

def get_config(file):
   '''Get a configuration file and create configuration object.'''
   cf = config.config(file)
   return cf

def get_data(cf):
   return gd.get_data(cf)

def setup_varinfo(cf,data):
   '''Given a configuration file, setup information about which variables
   are free or fixed and the slices into the parameter list.'''
   Nf = len(cf.Data.fcombs)       # Number of filters combs
   No = data.Nobj                 # Number of SNe
   Nm = data.ms.shape[0]          # dimension of mag arrays
   Nc = data.m1ids.shape[0]       # number of contsraints
   varinfo = variable.VarInfo()
   #i = 0
   vars = ['a','b','c','R','evar','H0','vpec','mt']

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
               elif id in cf.Data.filters:
                  ids.append(cf.Data.filters.index(id))
               else:
                  raise ValueError, "config: Invalid index on variable %s" % var
               vals.append(cf.Priors[opt])

         base_prior,vary,value = priors.get_prior(cf, res, data)

      if var in ['vpec','H0']:
         # variables that are scalars
         if vary:
            value = 0
         else:
            value = value[0]
      elif var in ['R']:
         # special case, pepending on now R is interpreted
         if cf.Model.Rv_global:
            if vary:
               value = 0
            else:
               value = value[0]
         else:
            if vary:
               value = zeros((Nf,))
            else:
               if value.shape[0] == 1:
                  value = value*ones((Nf,))
               elif value.shape[0] != Nf:
                  raise ValueError, "Parameter %s needs to have shape (%d)" % \
                        (var,Nf)

      elif var in ['mt']:
         # variables that are indexed by observation
         if vary:
            value = zeros((Nm,))
         else:
            if value.shape[0] == 1:
               value = value*ones((Nm,))
            elif value.shape[0] != Nm:
               raise ValueError, "Parameter %s needs to have shape (%d)" % \
                     (var,Nm)
      else:
         # variables that are indexed by filter combination
         if vary:
            value = zeros((Nf,))
         else:
            if value.shape[0] == 1:
               value = value*ones((Nf,))
            elif value.shape[0] != Nf:
               raise ValueError, "Parameter %s needs to have shape (%d)" % \
                     (var,Nf)
      v = variable.Variable(var, value, vary, base_prior)

      if ids and len(shape(value)) > 0:
         prior = []
         cindex = []
         cvalue = []
         for i in range(len(value)):
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

   return varinfo


'''Simple module for loading the data for Cepheid Metallicity fitter.'''
from numpy import *
from collect_data_s import collect_data
import pickle
from astropy.io import ascii

zpt_errors = {'u':0.091,
              'g':0.024,
              'r':0.020,
              'i':0.022,
              'B':0.031,
              'V':0.020,
              'Y':0.040,
              'J':0.040,
              'H':0.022}


class dataContainer(dict):
   '''A class derived from dictionary type that has a 
   __getattr__ function, allowing access by dot-notation.'''
   def __init__(self, *args, **kwargs):
      dict.__init__(self, *args, **kwargs)

   def __getattr__(self, key):
      try:
         return self.__getitem__(key)
      except:
         raise AttributeError, "No such attribute %s" % (key)

def get_data(cf):
   '''Import data and put it in a DataContainer objects for easy retrieval later.
   single argument is a configuration objects (config.py) that defines the data
   file and any filters that need to be applied.'''

   d = dataContainer()

   # Cepheid host data and covariance matrix
   f = open(cf.Data.hostdata)
   lines = f.readlines()
   cephlist = lines[0].strip().split()
   d['DMCeph'] = array(map(float, lines[1].strip().split()))
   d['C'] = len(cephlist)
   lines = lines[2:]
   if len(lines) != d.C:
      raise RuntimeError, "Problem with covariance matrix"
   d['Ccov'] = array([map(float, line.strip().split()) for line in lines])
   f.close()

   # Assemble the SN data
   d['filters'] = cf.Data.fcombs
   if type(d.filters) is not type([]):
      d['filters'] = [d.filters]
   f1,f2 = cf.Data.color

   d['all_filters'] = [f for f in d.filters]
   if f1 not in d.filters:
      d['all_filters'].append(f1)
   if f2 not in d.filters:
      d['all_filters'].append(f2)

   zhel,zcmb,hosts,M,eM,st,est,EBVgal,eEBVgal,t0s,names,sources = \
         collect_data(d.all_filters, cf.Data.filename)   
      
   if cf.Filters.omit is None:
      d['omit'] = []
   else:
      d['omit'] = cf.Filters.omit

   gids = array([name not in d.omit for name in names])
   if cf.Filters.cmax:
      gids = gids*less(M['B'] - M['V'], cf.Filters.cmax)
   gids = gids*greater_equal(st,cf.Filters.stmin)
   gids = gids*less_equal(st,cf.Filters.stmax)
   gids = gids*greater_equal(zcmb,cf.Filters.zmin)
   # make sure we have the Cepheid hosts
   for i in range(len(gids)):
      if names[i] in hosts:
         if hosts[names[i]] in cephlist:
            gids[i] = True
   # ensure we have valid colors
   gids = gids*less(M[f1], 90)
   gids = gids*less(M[f2], 90)

   # Gather Host masses
   if cf.Model.HostMass:
      d['M0'] = 11.
      K = []
      sigmaK = []
      tab = ascii.read(cf.Data.hostprops)
      hnames = list(tab['SN'])
      for i,name in enumerate(names):
         if name in hosts and hosts[name] in cephlist:
            idx = cephlist.index(hosts[name])
            mu = d['DMCeph'][idx]
         else:
            mu = 5*log10(3e5*zcmb[i]/73.) + 25
         if name in hnames:
            idx = hnames.index(name)
            if tab['K'][idx] > 0:
               K.append(tab['K'][idx])
               sigmaK.append(cf.Model.HostMassSigma)
            elif tab['M1'][idx] > 0:
               K.append(-2.5*(tab['M1'][idx]-1.04) + mu)
               sigmaK.append(cf.Model.HostMassSigma)
            elif tab['M2'][idx] > 0:
               K.append(-2.5*(tab['M2'][idx]-1.04) + mu)
               sigmaK.append(cf.Model.HostMassSigma)
            else:
               if cf.Model.ImputeMissingMass:
                  K.append(-2.5*(cf.Model.MissingMass - 1.04) + mu)
                  sigmaK.append(cf.Model.MissingMassSigma)
               else:
                  K.append(-1)
                  sigmaK.append(-1)
         else:
            K.append(-1)
            sigmaK.append(-1)
      gids = gids*greater(K,0)

   d['zhel']  = zhel[gids]
   d['zcmb'] = zcmb[gids]
   d['st'] = st[gids]
   d['est'] = est[gids]
   d['cs'] = M[f1][gids] - M[f2][gids]
   d['vcs'] = power(eM[f1][gids],2) + power(eM[f2][gids],2)
   d['zhel'] = zhel[gids]
   d['zcmb'] = zcmb[gids]
   d['names'] = [names[i] for i in range(len(gids)) if gids[i]]
   d['NFs'] = len(d.filters)
   d['NSNe'] = len(d.names)
   if cf.Model.HostMass:
      d['K'] = array(K)[gids]
      d['sigmaK'] = array(sigmaK)[gids]

   non_csp = []
   photsysids = []
   for name in d['names']:
      if name in sources:
         if sources[name] not in non_csp:
            non_csp.append(sources[name])
         photsysids.append(non_csp.index(sources[name])+1)
      else:
         photsysids.append(0)
   d['Nphotsys'] = len(non_csp)
   d['photsys'] = array(photsysids)

   d['host'] = []
   for i in range(len(d.names)):
      if d.names[i] in hosts and hosts[d.names[i]] in cephlist:
         d['host'].append(cephlist.index(hosts[d.names[i]])+1)
      else:
         d['host'].append(-1)

   if f1 in d.filters:
      d['findex1'] = d.filters.index(f1) + 1
   else:
      d['findex1'] = d.all_filters.index(f1) + 1
   if f2 in d.filters:
      d['findex2'] = d.filters.index(f2) + 1
   else:
      d['findex2'] = d.all_filters.index(f2) + 1

   d['ms'] = []
   d['vms'] = []
   #d['t0'] = []
   d['findex'] = []   # index array of which filter
   d['sindex'] = []
   d['zperr'] = zeros((d['Nphotsys'],d['NFs']))
   for i,f in enumerate(d.filters):
      ggids = greater(M[f][gids], 0)
      d.ms.append(M[f][gids][ggids])
      d.vms.append(power(eM[f][gids][ggids],2))
      #d.t0.append(t0s[f][gids][ggids])
      d.findex.append(array([i+1]*len(d.ms[-1])))
      d.sindex.append(arange(M[f][gids].shape[0])[ggids]+1)
      for j in range(d['Nphotsys']):
         d['zperr'][j,i] = zpt_errors[f]
   d['ms'] = concatenate(d.ms)
   d['NObs'] = d.ms.shape[0]
   d['vms'] = concatenate(d.vms)
   #d['t0'] = concatenate(d.t0)
   d['findex'] = concatenate(d.findex)
   d['sindex'] = concatenate(d.sindex)

   if cf.Model.Rv_global:
      f = open('Ia_A_poly.pickle')
      dd = pickle.load(f)
      f.close()
   
      rl = cf.Model.redlaw
      d['Amat'] = array([dd[rl][filt] for filt in d.all_filters])
      d['Al_order'] = array([dd[rl]['order'][filt] for filt in d.all_filters])
      d['N_coef'] = len(d.Amat[0])
      d['NFs2'] = len(d.all_filters)


   return d

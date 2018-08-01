'''Simple module for loading the data for Cepheid Metallicity fitter.'''
from numpy import *
from collect_data_s import collect_data
from A_lamb import A_lamb

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

   file = cf.data.filename
   filters = cf.data.filters
   if type(filters) is not type([]):
      if len(filters) > 1:
         filters = [f for f in filters]
      else:
         raise ValueError, "Need more than 1 filter to fit colors!"

   d['filters'] = filters
   d['all_filters'] = [f for f in ['u','B','V','g','r','i','Y','J','H'] \
         if f in d.filters + ['B','V']]

   zhel,zcmb,hosts,M,eM,st,est,EBVgal,eEBVgal,t0s,names,sources = \
         collect_data(d.all_filters, cf.data.filename)   
      
   if cf.Filters.omit is None:
      d['omit'] = []
   else:
      d['omit'] = cf.Filters.omit
      if type(d['omit']) is not type([]):
         d['omit'] = [d['omit']]

   BmV = M['B'] - M['V']

   gids = array([name not in d.omit for name in names])
   if cf.Filters.cmax:
      gids = gids*less(M['B'] - M['V'], cf.Filters.cmax)
   if cf.Filters.smin:
      gids = gids*greater_equal(st,cf.Filters.smin)
   if cf.Filters.smax:
      gids = gids*less_equal(st,cf.Filters.smax)

   d['st'] = st[gids]
   d['est'] = est[gids]
   d['vst'] = power(est[gids],2)
   d['EBVgal'] = EBVgal[gids]
   d['eEBVgal'] = eEBVgal[gids]
   d['names'] = [names[i] for i in range(len(gids)) if gids[i]]
   d['NFs'] = len(d.filters)
   d['Nobj'] = len(d.names)
   BmV = BmV[gids]

   d['ms'] = []
   d['ems'] = []
   d['vms'] = []
   d['fids'] = []   # index array of which filter
   dids = []
   for i,f in enumerate(d.filters):
      ggids = greater(M[f][gids], 0)
      #d.ms.append(M[f][gids][ggids] - A_lamb(f, EBVgal[gids][ggids], 3.1,
      #   redlaw=cf.Model.mw_redlaw))
      d.ms.append(M[f][gids][ggids])
      m1 = A_lamb(f, EBVgal[gids][ggids]-eEBVgal[gids][ggids], 3.1,
            redlaw=cf.Model.mw_redlaw)
      m2 = A_lamb(f, EBVgal[gids][ggids]+eEBVgal[gids][ggids], 3.1,
            redlaw=cf.Model.mw_redlaw)
      err = absolute(m2-m1)/2
      d.vms.append(power(eM[f][gids][ggids],2)+power(err,2))
      d.ems.append(sqrt(d.vms[-1]))
      d.fids.append(array([i]*len(d.ms[-1])))
      dids.append(arange(M[f][gids].shape[0])[ggids])
   d['ms'] = concatenate(d.ms)
   d['ems'] = concatenate(d.ems)
   d['vms'] = concatenate(d.vms)
   d['fids'] = concatenate(d.fids)
   d['oids'] = concatenate(dids)
   d['f0'] = 0
   d['findex0'] = nonzero(equal(d['fids'], 0))[0]

   # Binning based on inital guess of EBVhost
   EBV0 = BmV - (-0.082 + 0.048*(d['st']-1) + 1.39*power(d['st']-1,2))
   d['bids'] = searchsorted(cf.Priors.Rv_bins, EBV0)

   return d

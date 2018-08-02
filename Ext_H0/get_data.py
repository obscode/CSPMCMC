'''This module is responsible for getting the data needed for the STAN
runs. It must have a function called get_datget_data, which takes the configuration
instance as an argument and returns a dictionary with the following key/data:

S:           Number of Cepheid hosts (int)
DMobs:       Array of Cepheid distance modulii (float array[C])
Cov:         DM Covariance matrix for Cepheid hosts (float array[C,C])
Nsn:         Number of SNeIa (int)
zcmb:        CMB frame redshifts of each SN (float array[Nsn])
zhelio:      heliocentric redshift of each SN (float array[Nsn])
host:        Host galaxy index for each SN. -1 indicates no Cepheid host,
             otherwise, the index should correspond to value in DMobs
             indexed from 1 -> C (int array[Nsn])
K:           K-band magnitude of host galaxy (float array[NSNe])
M0:          host mass zero-point (float)
sigmaK:      error in host mass estimate (dex) (float array[NSNe])
Nphotsys:    Number of photometric systems being considered
photsys:     index of which photometric system, indexed from 1. Value of
             0 or less indicates natural system (no syst. error)
             (int array[Nsn])
st:          stretch/dm15 (float array[Nsn])
e_st:        error in st (float array[Nsn])
EBV:         E(B-V) for each SN (float array[Nsn])
Rv:          Rv or each SN (float array[Nsn])
ERprec       E(B-V)/Rv Precision matrix for each SN (float array[Nns,2,2])
             ERprec[i,j,k] = inv(Cov[i,j,k]), where Cov is covariance matrix
             of [E(B-V),Rv] for SN i.
NSNobs:      Number of SN observations (int)
Nfilt:       Number of filters used in observations (int)
m_sn:        magnitudes at maximum (float array[NSNobs])
em_sn:       error of m_sn (float array[NObs])
oid:         index array indicating which SN in m_sn/em_sn, indexed from 1
             (int array[NSNobs])
filt:        index array indicating which filter in ms/vms, indexed from 1
             (int array[NSNobs])
zperr:       zero-point error for filter,photsys combination
             (float array[Nphotsys,Nfilt])
in_mag:      Work in magnitudes? (int, boolean)
Nbasis:      Number of basis functions (int)
Bs:          Basis functions evaluated at st (float array[Nsn,Nbasis])
Al_order:    Order of polynomial fit to extinction vs (E(B-V),Rv) (int)
N_coef:      Dimension of coefficient matrix (int)
Al_coef:     Extinction vs (E(B-V),Rv) coefficients (float array[Nfilt,N_coef])
NCV:         Number of objects to remove for cross-validation. (int)
CVids:       Indeces of SNe to remove for cross-valication (int array[NCV])

'''

from numpy import*
from numpy import linalg
from astropy.io import ascii
import pickle
import config
import sys , os , string
from glob import glob
import bspline2
from collect_data_s import collect_data

zpt_errors = {'u':0.091,
              'g':0.024,
              'r':0.020,
              'i':0.022,
              'B':0.031,
              'V':0.020,
              'Y':0.040,
              'J':0.040,
              'H':0.022}

def toflux(mag, err, zp=22):
    f = power(10, -0.4*(mag - zp))
    ef = f*err/1.087
    return f,ef

def parse_ext(filename, names):
   '''Parse the filename for 6-tuple floats separated by whitespace,
   which are taken as name,E(B-V),+/-,Rv,+/-,cov(E(B-V),Rv). Return
   the data assocated with given names.'''

   with open(filename) as fin:
      lines = fin.readlines()
   lines = [line.strip().split() for line in lines]
   lines = [line for line in lines if len(line) == 6]
   data = {}
   for line in lines:
      try:
         data[line[0]] = map(float, line[1:])
      except:
          print "Warning:  ",line[0],"cannot be parsed"
   gids = array([name in data for name in names])
   data = array([data.get(name, [-1,-1,-1,-1,-1]) for name in names])
   return gids,data[:,0],data[:,1],data[:,2],data[:,3],data[:,4]


def get_data(cfg):
   ''' Load all the data based on the contents of the config file cfg.'''
   f = open(cfg.data.hostdata)
   lines = f.readlines()
   cephlist = lines[0].strip().split()
   DMs = array(map(float, lines[1].strip().split()))
   S = len(cephlist)
   lines = lines[2:]
   if len(lines) != S:
      raise RuntimeError, "Problem with covariance matrix"
   Cov = array([map(float, line.strip().split()) for line in lines])
   f.close()

   in_mag=cfg.model.in_mag
   # The data to be returned
   data = {'S':S, 'DMobs':DMs, 'Cov':Cov, 'in_mag':int(in_mag)}

   zhel,zcmb,hosts,M,eM,st,est,EBVgal,eEBVgal,t0s,sn,sources =\
          collect_data(cfg.data.sn_filt, cfg.data.sndata) 
   print sources
   # eids tells us which SNe have extinction estimates.
   eids,EBV,e_EBV,Rv,e_Rv,covER = parse_ext(cfg.data.extinctionData,sn)

   # Go through conditions and build up a yes/no array of which SNe
   # to include
   ceph_ids = array([hosts[name] in cephlist for name in sn])

   # Now, deal with filters:
   gids = greater(st, cfg.filter.st_min)
   gids = gids*less(st, cfg.filter.st_max)
   gids = gids*greater(zhel, cfg.filter.z_min)
   gids = gids*less(zhel, cfg.filter.z_max)
   gids = gids*less(EBV, cfg.filter.EBV_max)
   if cfg.filter.omit is not None:
      if type(cfg.filter.omit) is not type([]):
         cfg.filter.omit = [cfg.filter.omit]
      gids = gids*array([sname not in cfg.filter.omit for sname in sn])
   if cfg.filter.sources:
      if not (type(cfg.filter.sources) is type([])):
         cfg.filter.sources = [cfg.filter.sources]
      gids = gids*array([src in cfg.filter.sources for src in sources])
   gids = gids*eids

   data['M0'] = 11.
   K = []
   sigmaK = []
   tab = ascii.read(cfg.data.hostprops)
   hnames = list(tab['SN'])
   for i,name in enumerate(sn):
      if hosts[name] in cephlist:
         idx = cephlist.index(hosts[name])
         mu = DMs[idx]
      else:
         mu = 5*log10(3e5*zcmb[i]/73.) + 25
      if name in hnames:
         idx = hnames.index(name)
         if tab['K'][idx] > 0:
            K.append(tab['K'][idx])
            sigmaK.append(cfg.model.HostMassSigma)
         elif tab['M1'][idx] > 0:
            K.append(-2.5*(tab['M1'][idx]-1.04) + mu)
            sigmaK.append(cfg.model.HostMassSigma)
         elif tab['M2'][idx] > 0:
            K.append(-2.5*(tab['M2'][idx]-1.04) + mu)
            sigmaK.append(cfg.model.HostMassSigma)
         else:
            if cfg.model.ImputeMissingMass:
               K.append(-2.5*(cfg.model.MissingMass - 1.04) + mu)
               sigmaK.append(cfg.model.MissingMassSigma)
            else:
               K.append(-1)
               sigmaK.append(-1)
      else:
         K.append(-1)
         sigmaK.append(-1)
   if cfg.model.HostMass:
      gids = gids*greater(K,0)

   # Move this to main script, as it is procedural.
   #if not cfg.model.HostMass:
   #   data['alpha'] = 0
   # Now add ceph_ids so that they don't get thrown out with filtering
   #gids = gids + ceph_ids
   Nfilt = len(cfg.data.sn_filt)

   # extra error for objects whose photometric systems is not CSP
   all_sources = ['CSP']+list(set(sources.values())-set(['CSP']))
   data['Nphotsys'] = len(all_sources)-1
   print all_sources

   f=open('Ia_A_poly.pickle')
   d=pickle.load(f)
   f.close()

   # Ensure a list, even if a single filter
   if type(cfg.data.sn_filt) is not type([]):
      cfg.data.sn_filt = [cfg.data.sn_filt]

   Al_coef=array([d['fm'][filt] for filt in cfg.data.sn_filt])
   Al_order=array([d['fm']['order'][filt] for filt in cfg.data.sn_filt])
   N_coef=len(Al_coef[0])

   # Now check for good data. Note that even ceph_ids are thrown out if
   #   there is no data! So first, we figure out if any objects are
   #   completely missing data in requested filters
   nfilts = sum([less(M[f],90) for f in cfg.data.sn_filt], axis=0)
   gids = gids*greater(nfilts, 0)
   Nsn = sum(gids)
   oids0 = arange(1,Nsn+1)    # STAN indexes from 1
   
   m_sn = []
   em_sn = []
   oid = []
   filtid = []
   zperr = []
   for i,filt in enumerate(cfg.data.sn_filt):
      ggids = less(M[filt][gids], 90)
      m_sn.append(M[filt][gids][ggids])
      em_sn.append(eM[filt][gids][ggids])
      # Remember, STAN indexes from 1
      oid.append(oids0[ggids])
      filtid.append([i+1]*sum(ggids))
      zperr.append([zpt_errors[filt]]*data['Nphotsys'])
   m_sn = concatenate(m_sn)
   em_sn = concatenate(em_sn)
   oid = concatenate(oid)
   filtid = concatenate(filtid)
   zperr = array(zperr).T
   NSNobs = m_sn.shape[0]
 
   data['oid'] = oid
   data['filt'] = filtid
   data['NSNobs'] = NSNobs
   data['Nfilt'] = Nfilt
   data['zperr'] = zperr

   if cfg.model.in_mag==0:
       data['m_sn'],data['em_sn'] = toflux(m_sn,em_sn, 16.0)
   else:
      data['m_sn'] = m_sn
      data['em_sn'] = em_sn
   data['Nsn'] = Nsn
   data['Al_coef'] = Al_coef
   data['N_coef'] = N_coef
   data['Al_order'] = Al_order

   data['st'] = st[gids]
   data['e_st'] = est[gids]
   data['EBV'] = EBV[gids]
   data['Rv'] = Rv[gids]
   data['zcmb'] = zcmb[gids]
   data['zhelio'] = zhel[gids]
   data['K'] = array(K)[gids]
   data['sigmaK'] = array(sigmaK)[gids]
   snnames=[sn[i] for i in range(len(gids)) if gids[i]]

   # precision matrix = inv(covariance matrix)
   ERprec = []
   for i in range(data['EBV'].shape[0]):
      mat = array([[e_EBV[gids][i]**2, covER[gids][i]],
                   [covER[gids][i], e_Rv[gids][i]**2]])
      ERprec.append(linalg.inv(mat))
   data['ERprec'] = array(ERprec)
   
   data['host'] = array([cephlist.index(hosts[sn[i]])+1 \
                    if hosts[sn[i]] in cephlist else 0 \
                    for i in range(len(sn))], dtype=int)[gids]

   non_csp = []
   photsysids = []
   for name in snnames:
      if name in sources:
         if sources[name] not in non_csp:
            non_csp.append(sources[name])
         photsysids.append(non_csp.index(sources[name])+1)
      else:
         photsysids.append(0)
   data['Nphotsys'] = len(non_csp)
   data['photsys'] = array(photsysids)

   # BAsis functions
   if cfg.model.basis == 'poly':
      data['Nbasis'] = cfg.model.order+1
      bs = []
      for j in range(data['Nbasis']):
         bs.append(power(data['st']-1,j))
      bs = array(bs).T
   else:
      bs = bspline2.bspline_basis(cfg.model.knots, data['st'], cfg.model.order,
            gradient=cfg.model.gradient)
      data['Nbasis'] = bs.shape[1]
   data['Bs'] = bs

   # some extra stuff that STAN doesn't like but would be nice to have
   # (for plotting)
   extras = dict(cephlist=cephlist, snnames=snnames,
                 source=[sources.get(name,'CSP') for name in snnames])
                 
   if cfg.model.cv is not None:
      k = 1
      if type(cfg.model.cv) is type([]):
         cvnames = cfg.model.cv
      else:
         cvnames = [cfg.model.cv]
      CVids = []
      for name in extras['snnames']:
         if name in cvnames:
            CVids.append(k)
            k = k + 1
         else:
            CVids.append(0)
      if not CVids:
         raise "Error, no cv names found in sample"
      data['CVids'] = array(CVids)
      data['NCV'] = len(nonzero(data['CVids'])[0])

   return data,extras

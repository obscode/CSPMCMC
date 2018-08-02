# This module is responsible for loading the data from tables and generating
# the data dictionary that STAN will use. 

from numpy import*
from astropy.io import ascii
import pickle
import config
import sys , os , string
from glob import glob

# H-band data from Riess et al. (2016)
R16tab = ascii.read('Riess+2016tab4.mrt.dat')
SNhosts = {}
SNhosts['H'] = [host for host in list(set(R16tab['Field'])) 
                if host not in ['N4258','M31']]

def toflux(mag, err, zp=22):
    f = power(10, -0.4*(mag - zp))
    ef = f*err/1.087
    return f,ef
   

def get_data(cfg):
   ''' Load all the data based on the contents of the config file cfg.'''
   if cfg.filter.exclude_hosts is not None:
      cephlist = [host for host in SNhosts['H'] \
                  if host not in cfg.filter.exclude_hosts]
   else:
      cephlist = SNhosts['H']
   in_mag=cfg.model.in_mag
   
   # offset to the V-I color
   VI0 = cfg.model.VI0
   if VI0 is None:
      VI0 = 0
   
   Nset=len(cephlist)
   zp = 22
   
   VI=[]
   logP=[]
   mag=[]
   e_mag=[]
   Nobs=[]
   OH=[]
   N=[]
   ID=[]
   for i,ceph in enumerate(cephlist):
      #data=ascii.read(ceph+'_table.dat')
      data = R16tab[R16tab['Field'] == ceph]
      LogP = log10(data['Per'])

      gids = greater_equal(LogP,cfg.filter.logP_min)
      gids = gids*less_equal(LogP,cfg.filter.logP_max)
      if cfg.model.with_rej and 'Flag' in data.colnames:
         gids = gids*-(data['Flag'] == 'rej')
      mag.append(data['F160W'][gids])
      e_mag.append(data['{sigma}tot'][gids])
      OH.append(data['[O/H]'][gids])
      VI.append(data['F555W-F814W'][gids])
      if cfg.model.in_mag==0:
         f,ef = toflux(mag[-1],e_mag[-1], 25.0)
         mag[-1] = f
         if not cfg.model.lognormal:
            e_mag[-1] = ef
            
      logP.append(LogP[gids])
      N = sum(gids)
      ID.append((ones(N, dtype='int8'))*(i+1))
       
   # The data to be returned
   ID = concatenate(ID)
   data = {'N':len(ID),'S':Nset,'P':concatenate(logP),'VI':concatenate(VI),
         'mag':concatenate(mag),'ID':ID,
           'e_mag':concatenate(e_mag), 'in_mag':in_mag}
   data['OH'] = concatenate(OH)
   
   #____NGC 4258
   data4258=R16tab[R16tab['Field'] == 'N4258']
   LogP = log10(data4258['Per'])

   gids = greater_equal(LogP,cfg.filter.logP_min)
   gids = gids*less_equal(LogP,cfg.filter.logP_max)
   mag4258=data4258['F160W'][gids]
   e_mag4258=data4258['{sigma}tot'][gids]
   if cfg.model.in_mag==0:
      f,ef = toflux(mag4258,e_mag4258, 23.0)
      mag4258 = f
      if not cfg.model.lognormal:
         e_mag4258 = ef
   VI_4258=data4258['F555W-F814W'][gids] - VI0
   OH_4258=data4258['[O/H]'][gids]
   logP_4258=LogP[gids]

   data['N4258'] = len(mag4258)
   data['P_4258'] = logP_4258
   data['VI_4258'] = VI_4258 
   data['mag4258'] = mag4258
   data['e_mag4258'] = e_mag4258
   data['OH_4258'] = OH_4258

   #____M31
   dataM31=R16tab[R16tab['Field'] == 'M31']
   LogP = log10(dataM31['Per'])

   gids = greater_equal(LogP,cfg.filter.logP_min)
   gids = gids*less_equal(LogP,cfg.filter.logP_max)
   magM31=dataM31['F160W'][gids]
   e_magM31=dataM31['{sigma}tot'][gids]
   if cfg.model.in_mag==0:
      magM31,e_magM31 = toflux(magM31,e_magM31, 23.0)
   VI_M31=dataM31['F555W-F814W'][gids] - VI0
   logP_M31=LogP[gids]

   data['NM31'] = len(magM31)
   data['P_M31'] = logP_M31
   data['VI_M31'] = VI_M31 
   data['magM31'] = magM31
   data['e_magM31'] = e_magM31
   data['OH_M31'] = dataM31['[O/H]'][gids]
   
   #_____MW parallax/cepheids
   dataMW=ascii.read('MW_parallax_table.dat')
   gids = greater_equal(dataMW['logP'],cfg.filter.logP_min_MW)
   gids = gids*less_equal(dataMW['logP'],cfg.filter.logP_max_MW)
   magMW=[]
   e_magMW=[]
   magMW=dataMW['H'][gids]
   e_magMW=linspace(0.005,0.005,len(dataMW['H'][gids]))
   if cfg.model.in_mag==0:
       magMW,e_magMW = toflux(magMW,e_magMW, 3.0)
   VI_MW=dataMW['V-i'][gids] - VI0
   #VI_MW=dataMW['V-i'][gids]
   logP_MW=dataMW['logP'][gids]
   NobsMW=len(dataMW['logP'][gids])
   pi=dataMW['pi'][gids]
   e_pi=dataMW['e_pi'][gids]

   data['P_MW'] = logP_MW
   data['pi'] = pi 
   data['magMW'] = magMW 
   data['VI_MW'] = VI_MW 
   data['OH_MW'] = cfg.model.OH_MW
   data['NMW'] = NobsMW
   data['e_pi'] = e_pi
   data['e_magMW'] = e_magMW
   data['LK_prior'] = cfg.model.LK_prior
   
   #_____LMC
   dataLMC=ascii.read('LMC_opt+NIR_withoutVmissing.dat')
   gids = greater_equal(dataLMC['LogP'],cfg.filter.logP_min)
   gids = gids*less_equal(dataLMC['LogP'],cfg.filter.logP_max)
   magLMC=dataLMC['H'][gids]
   e_magLMC=linspace(0.005,0.005,len(dataLMC['H'][gids]))
   if cfg.model.in_mag==0:
      magLMC,e_magLMC = toflux(magLMC,e_magLMC, 12.0)
   VI_LMC=dataLMC['V-i'][gids] - VI0
   #VI_LMC=dataLMC['V-i'][gids]
   logP_LMC=dataLMC['LogP'][gids]
   NobsLMC=len(dataLMC['LogP'][gids])

   data['NLMC'] = NobsLMC
   data['P_LMC'] = logP_LMC
   data['magLMC'] = magLMC
   data['e_magLMC'] = e_magLMC
   data['VI_LMC'] = VI_LMC 
   data['OH_LMC'] = cfg.model.OH_LMC

   # some extra stuff that STAN doesn't like but would be nice to have
   # (for plotting)
   extras = dict(cephlist=cephlist)

   return data,extras

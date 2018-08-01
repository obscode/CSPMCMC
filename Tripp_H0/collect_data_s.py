from numpy import *
from glob import glob
import string
import types

def parse_file(files):
   data = {}
   all_names = []
   hosts = {}

   if type(files) is not types.ListType:
      files = [files]

   for file in files:
      f = open(file)
      lines = f.readlines()
      lines = map(string.strip, lines)
      lines = map(string.split, lines)
      
      for line in lines:
         if line[0][0] == "#":  continue
         if len(line) == 11:
            name,z,filter,dm15,edm15,Mmax,eMmax,c_m_dm15,DM,eDM,method = line
            EBVgal = 0
            eEBVgal = 0
         elif len(line) == 14:
            name,z,zcmb,filter,dm15,edm15,Mmax,eMmax,t0,EBVgal,eEBVgal,c_m_dm15,host,source = line

         elif len(line) == 15:
            name,z,zcmb,filter,dm15,edm15,Mmax,eMmax,t0,EBVgal,eEBVgal,c_m_dm15,DM,eDM,method = line
         else:
            raise ValueError, "Could not parse line from input file"
         if name not in all_names:  
            all_names.append(name)
            hosts[name] = host
      
         if filter not in data:
            data[filter] = {}
         data[filter][name] = map(float,[z,zcmb,dm15,edm15,Mmax,eMmax,EBVgal,eEBVgal,c_m_dm15,t0])
         data[filter][name].append(source)
   
   all_names.sort()
   return data,all_names,hosts

###############################  Data Collection ###############################
def collect_data(fs, filename):
   Ms = {}
   eMs = {}
   t0s = {}
   sources = {}
   
   names = []
   zhel = []
   zcmb = []
   hosts = []
   dm15s = []
   edm15s = []
   EBVgals = []
   eEBVgals = []

   data,all_names,hosts = parse_file(filename)
   
   for sn in all_names:
      for filt in data:
         if sn in data[filt]: break
      z,zc,dm15,edm15,Bmax,eBmax,EBVgal,eEBVgal,c_B_dm15,t0,source = \
         data[filt][sn]
      names.append(sn)
      zhel.append(z)
      zcmb.append(zc)
      dm15s.append(dm15)
      edm15s.append(edm15)
      EBVgals.append(EBVgal)
      eEBVgals.append(eEBVgal)
      if source != 'CSP':
         sources[sn] = source
      for f in data:
         if f not in Ms:
            Ms[f] = []
            eMs[f] = []
            t0s[f] = []
         if sn in data[f]:
            Ms[f].append(data[f][sn][4])
            eMs[f].append(data[f][sn][5])
            t0s[f].append(data[f][sn][9])
         else:
            Ms[f].append(-1)
            eMs[f].append(-1)
            t0s[f].append(-1)
       
   zhel = array(zhel)
   zcmb = array(zcmb)
   dm15s = array(dm15s)
   edm15s = array(edm15s)
   EBVgals = array(EBVgals)
   eEBVgals = array(eEBVgals)
   for f in Ms:
      Ms[f] = array(Ms[f])
      eMs[f] = array(eMs[f])
      t0s[f] = array(t0s[f])
   return(zhel,zcmb,hosts,Ms,eMs,dm15s,edm15s,EBVgals,eEBVgals,t0s,names,
          sources)

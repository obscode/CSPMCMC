from numpy import *
from glob import glob
import string
import types

def parse_file(files):
   data = {}
   all_names = []

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
         elif len(line) == 13:
            name,z,filter,dm15,edm15,Mmax,eMmax,EBVgal,eEBVgal,c_m_dm15,DM,eDM,method = line
         elif len(line) == 15:
            name,z,zcmb,filter,dm15,edm15,Mmax,eMmax,t0,EBVgal,eEBVgal,c_m_dm15,DM,eDM,method = line
         else:
            raise ValueError, "Could not parse line from input file"
         if name not in all_names:  all_names.append(name)
      
         if filter not in data:
            data[filter] = {}
         data[filter][name] = map(float,[z,dm15,edm15,Mmax,eMmax,EBVgal,eEBVgal,c_m_dm15,DM,eDM])
         data[filter][name].append(method)
   
   all_names.sort()
   return data,all_names

###############################  Data Collection ###############################
def collect_data(fs, filename='/home/cburns/data/MCMC/Tripp/st_template_fits_20121210.dat'):
   Ms = {}
   eMs = {}
   #for f in fs:
   #   Ms[f] = []
   #   eMs[f] = []
   
   names = []
   dm15s = []
   edm15s = []
   EBVgals = []
   eEBVgals = []

   data,all_names = parse_file(filename)
   
   for sn in all_names:
      #has_all_filters = True
      #for f in fs:
      #   if sn not in data[f]:
      #      has_all_filters = False
      #      break
      #if not has_all_filters:  continue
      if sn not in data['B'] or sn not in data['V']:
         continue

      z,dm15,edm15,Bmax,eBmax,EBVgal,eEBVgal,c_B_dm15,DM,eDM,method = data['B'][sn]
      names.append(sn)
      dm15s.append(dm15)
      edm15s.append(edm15)
      EBVgals.append(EBVgal)
      eEBVgals.append(eEBVgal)
      for f in data:
      #for f in fs:
         if f not in Ms:
            Ms[f] = []
            eMs[f] = []
         if sn in data[f]:
            Ms[f].append(data[f][sn][3])
            eMs[f].append(data[f][sn][4])
         else:
            Ms[f].append(-1)
            eMs[f].append(-1)
       
   dm15s = array(dm15s)
   edm15s = array(edm15s)
   EBVgals = array(EBVgals)
   eEBVgals = array(eEBVgals)
   for f in Ms:
      Ms[f] = array(Ms[f])
      eMs[f] = array(eMs[f])
   return(Ms,eMs,dm15s,edm15s,EBVgals,eEBVgals,names)

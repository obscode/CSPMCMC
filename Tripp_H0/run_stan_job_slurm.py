#!/usr/bin/env python

import sys,string,os
sys.path.append('../lib/')
from subprocess import Popen,PIPE,STDOUT
import time
import ConfigParser
import config
import glob

depends = ['Tripp_stan.py','get_data.py','collect_data_s.py',
      'A_lamb.py','Ia_A_poly.pickle','plot_Tripp.py']
depends += glob.glob('../lib/*.py')

if len(sys.argv) < 2 or '-h' in sys.argv:
   print "Usage:   runjob.py config-file"
   sys.exit(1)

cfgfile = sys.argv[1]
if not os.path.isfile(cfgfile):
   print "Error:  %s not found" % cfgfile
   sys.exit(1)

#nchains = int(sys.argv[2])

cf = config.config(cfgfile)
datafiles = cf.Data.filename
if type(datafiles) is not type([]):
   datafiles = [datafiles]
depends += datafiles
depends.append(cf.Data.hostdata)
if cf.Model.HostMass:
   depends.append(cf.Data.hostprops)
if cf.Model.Rv_global:
   if cf.Model.fixed_DM:
      depends.append('model_R_tied_fixed.stan')
   else:
      depends.append('model_R_tied.stan')
else:
   if cf.Model.fixed_DM:
      depends.append('model_R_free_fixed.stan')
   else:
      depends.append('model_R_free.stan')
if cf.Model.HostMass:
   depends[-1] = depends[-1].replace('.stan','_HM.stan')

nchains = cf.Sampler.chains
if nchains > 12:
   print "Error: nchains must be at most 12"
   sys.exit(1)

cfgfile = os.path.realpath(cfgfile)
dest = os.path.dirname(cfgfile)
bcfgfile = os.path.basename(cfgfile)

submit = '''#!/bin/bash
#SBATCH -p OBS
#SBATCH -N 1
#SBATCH -n {0}
#SBATCH -e {1}/stderr
#SBATCH -o {1}/stdout
export PATH=/share/apps/obs/cburns/anaconda2/bin:${{PATH}}
source activate pystan

DEST={1}
SDIR={2}
LHOST=$SLURM_JOB_NODELIST
echo "WORKING at $DEST on $LHOST"
echo "OUTPUT going to $DEST"

mkdir -p $DEST
mkdir -p $DEST/object_plots
if [ ! -d $DEST ]
then
   echo $DEST not created
   exit
fi
'''.format(nchains,dest,os.getcwd())
for file in depends:
   submit += "\ncp ${SDIR}/%s $DEST" % (file)

submit += "\ncd $DEST\npython Tripp_stan.py %s" % (bcfgfile)
if not cf.Model.fixed_DM:
   submit += "\npython plot_Tripp.py %s" % (bcfgfile)

fout = open(os.path.join(dest, 'submit.sh'), 'w')
fout.write(submit)
fout.close()

p = Popen(['sbatch','-J','Tripp_stan'], stdout=PIPE, stdin=PIPE, stderr=STDOUT)
out = p.communicate(input=submit)
print out

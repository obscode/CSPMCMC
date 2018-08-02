#!/usr/bin/env python

import sys,string,os
sys.path.append('../lib/')
from subprocess import Popen,PIPE,STDOUT
import time
import config
import glob

depends = ['run_stan.py',
           'get_data.py','generate_STAN.py','plot_results.py']
depends += ['LMC_opt+NIR_withoutVmissing.dat']
depends += ['MW_parallax_table.dat']
depends += ['Riess+2016tab4.mrt.dat','Riess+2016tab5.dat']
depends += glob.glob('../lib/*.py')

submit_str = '''#!/bin/bash
#SBATCH -p OBS
#SBATCH -N 1
#SBATCH -n {0}
#SBATCH --mem=32768
#SBATCH -e {1}/stderr     # stderr
#SBATCH -o {1}/stdout     # stdout

export PATH=/share/apps/obs/cburns/anaconda2/bin:${{PATH}}
#export PYTHONPATH=/home/cburns/CarPy_local/lib_local
source activate pystan

DEST={1}
SDIR={2}
LHOST=$SLURM_JOB_NODELIST
echo "WORKING at $DEST on $LHOST"
echo "OUTPUT going to $DEST"

mkdir -p $DEST
if [ ! -d $DEST ]
then
   echo $DEST not created
   exit
fi
'''

if len(sys.argv) < 2 or '-h' in sys.argv:
   print "Usage:   run_stan_job_slurm.py config-files"
   sys.exit(1)

cfgfiles = sys.argv[1:]
for cfgfile in cfgfiles:
   if not os.path.isfile(cfgfile):
      print "Error:  %s not found" % cfgfile
      sys.exit(1)
   
   cfg = config.config(cfgfile)
   nchains = cfg.sampler.chains
   if nchains > 16:
      print "Error: nchains must be at most 16"
      sys.exit(1)
   
   cfgfile = os.path.realpath(cfgfile)
   dest = os.path.dirname(cfgfile)
   bcfgfile = os.path.basename(cfgfile)

   submit = submit_str.format(nchains,dest,os.getcwd())
   
   for file in depends:
      submit += "\ncp ${SDIR}/%s $DEST" % (file)
   
   submit += "\ncd $DEST\npython run_stan.py %s" % \
         (bcfgfile)
   submit += "\npython plot_results.py %s" % (bcfgfile)
   
   f = open(os.path.join(dest,'submit.sh'), 'w')
   f.write(submit)
   f.close()
   p = Popen(['sbatch','-J','DM_Hosts'], stdout=PIPE, stdin=PIPE, 
       stderr=STDOUT)
   out = p.communicate(input=submit)

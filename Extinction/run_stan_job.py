#!/usr/bin/env python

import sys,string,os
sys.path.append('../lib/')
from subprocess import Popen,PIPE,STDOUT
import time
import ConfigParser
import config
import glob

depends = ['EBV.py','EBV_stan.py','make_stan_model.py',
      'get_data.py','collect_data_s.py','plot_EBV.py','priors.py',
      'variable.py', 'A_lamb.py','Ia_A_poly.pickle','../lib/bspline2.py']
depends += glob.glob('../lib/*.py')

if len(sys.argv) < 2 or '-h' in sys.argv:
   print "Usage:   runjob.py config-file chains"
   sys.exit(1)

cfgfile = sys.argv[1]
if not os.path.isfile(cfgfile):
   print "Error:  %s not found" % cfgfile
   sys.exit(1)

#nchains = int(sys.argv[2])

cf = config.config(cfgfile)
datafiles = cf.data.filename
if type(datafiles) is not type([]):
   datafiles = [datafiles]
depends += datafiles

nchains = cf.Sampler.chains
if nchains > 12:
   print "Error: nchains must be at most 12"
   sys.exit(1)

cfgfile = os.path.realpath(cfgfile)
dest = os.path.dirname(cfgfile)
bcfgfile = os.path.basename(cfgfile)

submit = '''
#!/bin/bash
#PBS -S /bin/bash
#PBS -l nodes=1:ppn={0}
#PBS -l walltime=10:00:00
#PBS -e {1}/stderr
#PBS -o {1}/stdout
export PATH=/home/cburns/Enthought/Canopy_64bit/User/bin:${{PATH}}
export PYTHONPATH=${{PYTHONPATH}}:/home/cburns/CarPy/dist/lib_local

DEST={1}
SDIR={2}
LHOST=`cat $PBS_NODEFILE`
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

submit += "\ncd $DEST\npython EBV_stan.py %s" % (bcfgfile)
submit += "\npython plot_EBV.py %s" % (bcfgfile)

p = Popen(['qsub','-N','color_edge'], stdout=PIPE, stdin=PIPE, stderr=STDOUT)
out = p.communicate(input=submit)

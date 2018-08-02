#!/bin/bash
#PBS -S /bin/bash           # use bash dammit
#PBS -l nodes=1:ppn=4     # Memory
#PBS -l walltime=10:00:00   # maximum 10-hour runtime
#PBS -e /data002/cburns/MCMC/H0/runs/lognormal/all/B/stderr.${PBS_JOBID}          # stderr
#PBS -o /data002/cburns/MCMC/H0/runs/lognormal/all/B/stdout.${PBS_JOBID}          # stdout

export PATH=/data002/cburns/anaconda2/bin:${PATH}
source activate pystan
export PYTHONPATH=/home/cburns/CarPy_local/lib:${PYTHONPATH}

DEST=/data002/cburns/MCMC/H0/runs/lognormal/all/B
SDIR=/data002/cburns/MCMC/H0
LHOST=`cat $PBS_NODEFILE`
echo "WORKING at $DEST on $LHOST"
echo "OUTPUT going to $DEST"

mkdir -p $DEST
if [ ! -d $DEST ]
then
   echo $DEST not created
   exit
fi

cp ${SDIR}/script.H0_pystan_extinction.py $DEST
cp ${SDIR}/bspline2.py $DEST
cp ${SDIR}/get_data.py $DEST
cp ${SDIR}/generate_STAN.py $DEST
cp ${SDIR}/plot_results.py $DEST
cp ${SDIR}/Ia_A_poly.pickle $DEST
cp ${SDIR}/supernovae_DR3_extinction.dat $DEST
cp ${SDIR}/CSPI_masses.dat $DEST
cd $DEST
rm -f traces_fixed.pickle
python script.H0_pystan_extinction.py base_fixed.cfg
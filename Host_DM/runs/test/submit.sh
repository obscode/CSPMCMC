#!/bin/bash
#SBATCH -p OBS
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=32768
#SBATCH -e /home/cburns/CSPMCMC/Host_DM/runs/test/stderr     # stderr
#SBATCH -o /home/cburns/CSPMCMC/Host_DM/runs/test/stdout     # stdout

export PATH=/share/apps/obs/cburns/anaconda2/bin:${PATH}
#export PYTHONPATH=/home/cburns/CarPy_local/lib_local
source activate pystan

DEST=/home/cburns/CSPMCMC/Host_DM/runs/test
SDIR=/home/cburns/CSPMCMC/Host_DM
LHOST=$SLURM_JOB_NODELIST
echo "WORKING at $DEST on $LHOST"
echo "OUTPUT going to $DEST"

mkdir -p $DEST
if [ ! -d $DEST ]
then
   echo $DEST not created
   exit
fi

cp ${SDIR}/run_stan.py $DEST
cp ${SDIR}/get_data.py $DEST
cp ${SDIR}/generate_STAN.py $DEST
cp ${SDIR}/plot_results.py $DEST
cp ${SDIR}/LMC_opt+NIR_withoutVmissing.dat $DEST
cp ${SDIR}/MW_parallax_table.dat $DEST
cp ${SDIR}/Riess+2016tab4.mrt.dat $DEST
cp ${SDIR}/Riess+2016tab5.dat $DEST
cp ${SDIR}/../lib/gelman_rubin.py $DEST
cp ${SDIR}/../lib/myplotlib.py $DEST
cp ${SDIR}/../lib/bspline2.py $DEST
cp ${SDIR}/../lib/MCMCstats.py $DEST
cp ${SDIR}/../lib/fit_poly.py $DEST
cp ${SDIR}/../lib/STANstats.py $DEST
cp ${SDIR}/../lib/sigfig.py $DEST
cp ${SDIR}/../lib/config.py $DEST
cp ${SDIR}/../lib/columns.py $DEST
cd $DEST
python run_stan.py base.cfg
python plot_results.py base.cfg
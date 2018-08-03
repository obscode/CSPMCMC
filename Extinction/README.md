# Extinction Estimates based on SNIa Colors

The code in this folder implements the methods outlined in  
[Burns et al. 2014)(http://adsabs.harvard.edu/abs/2014ApJ...789...32B).
It uses data obtained by fitting SNIa light curves, resulting in 
magnitudes in each filter at maximum, and a measure of the light-curve (LC)
width (stretch, dm15, x1, etc) for each SN. The code then uses these
magnitudes and LC widths to determine intrinsic colors and extinction
properties for the sample.

## Data Format
The input data is the same as the other two codes (Tripp_H0 and Ext_H0) to
keep things simple. Each line in the data file represents a single SN and filter
combination and each field should be white-space delimited. The fields are:

1- SN name (string, no white space)
2- Heliocentric redshift (float)
3- CMB-frame redshift (float)
4
 hosted SNeIa). The result is a file, DM_covar.dat, which contains the distances and their covariance matrix.
This file is used as input for the other codes (Tripp_H0 and Ext_H0) to determined the Hubble constant using
more distanct SNe Ia.

The code can be run simply by specifying a configuration file. An example (`example.cfg`) is provided here.
The configuration controls key parameters of the fitting procedure and model. The code is then run as
follows:
  
  `python run_stan.py example.cfg`

If you want to modify the code to work on another data-set, you can easily do so by editing get_data.py. The
comments therein describe how the data need to be prepared to send to the STAN code.

There are two scripts, `run_stan_job.py` and `run_stan_job_slurm.py` that can be used to submit jobs to
PBS and SLURM queue managers, respectively. Make a working folder somewhere and put a configuration file there.
Run either script from this folder, specifying the path to the configuration file. The scripts will take
care of generating submit scripts and submitting them to the job scheduler. The submit scripts will also
copy over all the needed files to the working folder.

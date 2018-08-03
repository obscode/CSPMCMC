# Determining the Hubble Constant using the Extinction Model

The code in this folder implements the methods outlined in  Burns et al. (2018)
(submitted for publication).  It uses data obtained by fitting SNIa light
curves, resulting in magnitudes in each filter at maximum, and a measure of the
light-curve (LC) width (stretch, dm15, x1, etc) for each SN. It also requires
extinction parameters (E(B-V), Rv, and covariance) for each SN. These are
provided by the `Extinction` code and that should be run first. The code then
uses these magnitudes, LC widths, extinctions as well as extra meta-data for
each SN ( redshifts, host galaxies, etc) to compute distances. Given a certain
number of SNe that occurred in hosts with distances determined with the
`Host_DM` code, the Hubble diagram can be constructed and calibrated to obtain
the Hubble constant. Along with the SN data and extinction data you will need
the `DM_cov.dat` file produced by the `Host_DM` code.

## Data Format
The input SN data is the same as the other two codes (`Extinction` and Tripp_H0)
to keep things simple. Each line in the data file represents a single SN and
filter combination and each field should be white-space delimited. The fields
are:

1. SN name (string, no white space)
2. Heliocentric redshift (float)
3. CMB-frame redshift (float)
4. filter (string)
5. LC width (float)
6. error in LC width (float)
7. magnitude at maximum, K-corrected to rest frame and corrected for MW
   extintion (float)
8. Error in magnitude at maximum
9. phase of first observation (not currently used) (float)
10. E(B-V) from Milky-Way (not currently used) (float)
11. error in MW E(B-V) (not currently used) (float)
12. covariance between magnitude at maximum and LC width (float)
13. Host galaxy name (string)
14. Source of photometry (string)

The data can be split up over multiple files. Note that the host galaxy name
is only important for matching a SN to a host that has a Cepheid distance.
For those objects, the galaxy host name must match one of the names in
`DM_cov.dat`.

The code can be run simply by specifying a configuration file. An example 
(`example.cfg`) is provided here.  The configuration controls key parameters 
of the fitting procedure and model (see comments). The code is then run as
follows:
  
  `python Ext_H0.py example.cfg`

If you want to modify the code to work on another data-set, you can easily do 
so by editing get_data.py. The comments therein describe how the data need to
be prepared to send to the STAN code.

There is a post-run script, `plot_Tripp.py` that should be run to examine the
results and produce a table of results, including the Hubble constant and
its error.

  `python plot_results.py extample.cfg`

There are two more scripts, `run_stan_job.py` and `run_stan_job_slurm.py` that
can be used to submit jobs to PBS and SLURM queue managers, respectively. Make
a working folder somewhere and put a configuration file there.  Run either
script from this folder, specifying the path to the configuration file. The
scripts will take care of generating submit scripts and submitting them to the
job scheduler. The submit scripts will also copy over all the needed files to
the working folder.

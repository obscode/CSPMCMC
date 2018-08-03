# CSPMCMC
MCMC codes used by the Carnegie Supernova Project

These are the python codes used for the analysis of type Ia supernovae.
Currently, there are four sets of code:

* `Host_DM`: Determine the distances to hosts with Cepheid variables. 
* `Extinction`:  Compute the extinction properties of a sample of SNe Ia.
* `Tripp_H0`:  Determine the Hubble constant using the Tripp method
* `Ext_H0`:  Determine the Hubble constant using the extinctions determined
       from `Extinction`.

Each set of codes runs independently, though `Tripp_H0` and `Ext_H0` require
the distance modulus output from `Host_DM` and `Ext_H0` requires the
extinction properties output from `Extinction`.

The following python modules/packages are required to run the code, all of which
are available through `anaconda`. The versions are  those used when 
developing the code.

* Python 2.7
* Numpy 1.11.3
* Scipy 0.19.0
* matplotlib 2.0.2 (for graphical output)
* pystan 2.14.0.0
* astropy 1.3.3

There are modules used in common between all the codes. These are located in
the `lib` folder and should be made available in python's search path
(e.g., by setting `$PYTHONPATH`). If you use the provided PBS/SLURM submission
scripts in each folder, they will take care of copying the needed python
modules/scripts to the working folder.

Each folder has a readme with instructions on usage.

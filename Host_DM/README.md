# Cepheid Host Galaxy Distances

The code in this folder uses Cepheid variable data to compute the distance to their host galaxies (which also
hosted SNeIa). The result is a file, DM_covar.dat, which contains the distances and their covariance matrix.
This file is used as input for the other codes (Tripp_H0 and Ext_H0) to determined the Hubble constant using
more distanct SNe Ia.

The code can be run simply by specifying a configuration file. An example (`example.cfg`) is provided here.
The configuration controls key parameters of the fitting procedure and model. The code is then run as
follows:
  
  python run_stan.py example.cfg

If you want to modify the code to work on another data-set, you can easily do so by editing get_data.py. The
comments therein describe how the data need to be prepared to send to the STAN code.

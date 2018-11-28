'''MCMC statistics module.
Author:  Chris Burns
         Carnegie Observatories, cburns@obs.carnegiescience.edu

Contents:  This module has a class, chains, that reads in chains produced by
           an MCMC algorithm. It is assumed there will be C parallel chains
           of length N iterations for P parameters. The chains array will
           therefore have shape (C,N,P). The member functions and data
           allows for doing statistics and plotting of histograms, confidence
           regions, etc.
Requires:  - MCMC module in question (emcee, PyMC, pySTAN, etc)
           - NumPy (developed under version 1.6.0)
           - scipy (developed under verison 0.9.0)
           - matplotlib (developed under version 1.0.1)
           - gelman_rubin.py (module <<URL>>)
           - fit_poly.py (module <<URL>>)
           - pymc (if you want gp interpolation)
'''

from numpy import *
from scipy.optimize import brentq,brent
from scipy import histogram2d
from scipy.special import chdtr
from scipy.interpolate import bisplrep,bisplev
from scipy import integrate
from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib import mlab
from matplotlib.colors import Normalize
from matplotlib.patches import Ellipse
from matplotlib.transforms import IdentityTransform as it
from glob import glob
import os,string,copy,pickle, re
try:
   from pymc import gp
except:
   gp = None
import gelman_rubin
import fit_poly
if 'OMP_NUM_THREADS' not in os.environ:
   os.environ['OMP_NUM_THREADS'] = '4'
try:
   import triangle
except:
   try:
      import corner as triangle
   except:
      triangle = None

norm = Normalize(0,1)

'''The data is going to be stored using the following scheme.  This gets
complicated because we have 1) multiple chains that were run in parallel,
2) possibly multiple chains in the db, 3) parameters with more than one
dimension.  The data is stored as:

   data[i,j,k]

   i = chain number (for parallel chains)
   j = in-chain index (iterations along the MCMC chain)
   k = index into the parameter

   There will be a list of parameters:  self.pars.  For each of these, there
   will be a self.slice into k for parameter vectors.
'''

def wrap_M(x, M, lower, upper, log=False):
   if len(shape(x)) == 0:
      scalar = True
      x = array([x])
   else:
      scalar = False
   res = M(x)
   res = where(less(x,lower), 0, res)
   res = where(greater(x,upper), 0, res)
   if log:  res = exp(res)
   if scalar:
      return res[0]
   else:
      return res

def wrap_poly(x, x0, pars, lower, upper, log=False):
   if len(shape(x)) == 0:
      scalar = True
      x = array([x])
   else:
      scalar = False
   res = fit_poly.poly(x, x0, pars)
   res = where(less(x, lower), 0, res)
   res = where(greater(x, upper), 0, res)
   if log:  res = exp(res)
   if scalar:
      return res[0]
   else:
      return res


def plus_minus(arr, bins=30, conf=0.68, xrange=None, func='poly',
      fit_log=True, order=7, debug=False, zero_pad=False,
      end_tol=[None,None]):
   hist0,bins = histogram(arr, bins=bins, range=xrange)
   xb = (bins[1:] + bins[:-1])/2
   if fit_log:
      gids = greater(hist0, 0)
      xb = xb[gids]
      var = 1./hist0[gids]
      hist = log(hist0[gids])
   else:
      var = hist0*1
      hist = hist0*1
   if xrange is None:
      xrange = (bins[0],bins[-1])
   xplot = linspace(xrange[0]*0.9,xrange[1]*1.1,101)
   if debug:
      fig = plt.figure()
      if fit_log:
         y1 = hist
         y2 = exp(hist)
      else:
         y1 = log(hist)
         y2 = hist
      ax1 = fig.add_subplot(211)
      ax1.plot(xb, y1, 'o')
      ax2 = fig.add_subplot(212)
      ax2.plot(xb, y2, 'o')

   if func == 'gp' or func == 'poly':
      if func == 'gp':
         if not gp:
            raise RuntimeError, "To use GP interpolation, you need to install pymc"
         scale = xb.max() - xb.min()
         M = gp.Mean(lambda x:  zeros(x.shape[0], 
            dtype=float32)+ median(hist))
         C = gp.Covariance(gp.matern.euclidean, diff_degree=3, 
            scale=scale*0.5, amp=std(hist))
 
         # Pad with zeros
         if zero_pad and not fit_log:
            obs_mesh = concatenate([xb.min()+(xb-xb.max())[:-1],xb,
               xb.max()+(xb-xb.min())[1:]])
            obs = concatenate([hist[1:]*0,hist,hist[1:]*0])
            var = concatenate([hist[1:]*0,var,hist[1:]*0])
         else:
            obs_mesh = xb
            obs = hist
         gp.observe(M,C,obs_mesh=obs_mesh,obs_vals=obs, obs_V=var)
 
         func = lambda x:  wrap_M(x, M, xb[0], xb[-1], log=fit_log)
 
      else:
         x0 = xb[argmax(hist)]
         pars,epars = fit_poly.fitpoly(xb,hist,w=1./var,x0=x0,k=order)
         func = lambda x:  wrap_poly(x, x0, pars, xb[0], xb[-1], 
                                     log=fit_log)
 
      if debug:
         ax1.plot(xplot, log(func(xplot)), '-')
         ax2.plot(xplot, func(xplot), '-')
      oneside=False
      if argmax(hist) == 0:
         mod = xb[0]
         oneside=True
      elif argmax(hist) == len(xb)-1:
         mod = xb[-1]
         oneside=True
      else:
         mod0 = xb[argmax(hist)]
         try:
            mod = brent(lambda x:  -func(x),brack=(xb.min(),mod0,xb.max()))
         except: 
            # monotonic.  Take extremum
            oneside=True
            if func(xb[0]) > func(xb[-1]):
               mod = xb[0]
            else:
               mod = xb[-1]

      fac = integrate.quad(func, xb[0], xb[-1])[0]
      prob = lambda x:  func(x)/fac

      #end tolerance  if requested
      lower_limit = False
      upper_limit = False
      if end_tol[0] is not None and float(hist0[0])/hist0.max() > end_tol[0]:
         lower_limit = True
      if end_tol[1] is not None and float(hist0[-1])/hist0.max() > end_tol[1]:
         upper_limit = True
      if lower_limit and upper_limit:
            #  too flat, return mode, but no limits
            return mod,nan,nan
      elif lower_limit and not upper_limit:
         # one-sided
         tail = (1-conf)
         upper = brentq(\
               lambda x: integrate.quad(prob, x, xplot[-1])[0]-tail, 
                  mod, xplot[-1])
         return mod,nan,upper
      elif upper_limit and not lower_limit:
         tail = (1-conf)
         lower = brentq(\
               lambda x: integrate.quad(prob, xplot[0], x)[0]-tail, 
               xplot[0], xplot[-1])
         return mod,lower,nan

      if debug:
         ax1.axvline(mod, color='red')
         ax2.axvline(mod, color='red')
      
      if oneside:
         tail = (1-conf)
      else:
         tail = (1-conf)/2
      if integrate.quad(prob, xplot[0], mod)[0] < tail:
         # No lower bound
         minus = nan
      else:
         lower = brentq(\
               lambda x: integrate.quad(prob, xplot[0], x)[0]-tail, 
               xplot[0], mod)
         minus = mod - lower
         if debug:
            ax1.axvline(lower, color='orange')
            ax2.axvline(lower, color='orange')
      #test for upper bound
      if integrate.quad(prob, mod, xplot[-1])[0] < tail:
         # No upper bound
         plus = nan
      else:
         upper = brentq(\
               lambda x: integrate.quad(prob, x, xplot[-1])[0]-tail, 
                  mod, xplot[-1])
         plus = upper-mod
         if debug:
            ax1.axvline(upper, color='orange')
            ax2.axvline(upper, color='orange')

   else:
      hist = hist*1.0/sum(hist)
      mid = argmax(hist)
      mod = xb[mid]
      if debug:
         ax1.axvline(mod, color='red')
         ax2.axvline(mod, color='red')
      i0 = 0
      i1 = len(hist)-1
      prob = 0
      while (prob < (1-conf)/2):
         if i0 < mid:
            i0 += 1
         else:
            break
         prob = sum(hist[0:i0])
      if i0 == 0:
         lower = None
      else:
         lower = xb[i0]
         if debug:
            ax1.axvline(lower, color='orange')
            ax2.axvline(lower, color='orange')
      while(prob < 1-conf):
         if i1 > mid:
            i1 -= 1
         else:
            break
         prob = sum(hist[0:i0]) + sum(hist[i1:])
      if i1 == len(xb)-1:
         upper = None
      else:
         upper = xb[i1]
         if debug:
            ax1.axvline(upper, color='orange')
            ax2.axvline(upper, color='orange')
      if upper is not None:
         plus = upper-mod
      else:
         plus = nan
      if lower is not None:
         minus = mod-lower
      else:
         minus = nan
   return mod, minus, plus



class MCMCchains:

   def __init__(self, filename=None, chains=None, verbose=False, **kwargs):

      self.verbose = verbose
      self.kwargs = kwargs.copy()

      self.chains = None    # The raw data as an ND-array
      self.params = []      # the parameter names
      self.slices = {}      # The slices into self.chains for each parameter
      self.shapes = {}      # The shape of the parameter
      self.parlabels = []   # parameter labels for plots
      self.indexlabels = [] # index labels for plots
      self.load_chains(filename, chains, **kwargs)

      self.parlabels = []
      self.indexlabels = []
      for par in self.params:
         npar = self.get_trace(par).shape[2]
         self.parlabels += [par]*npar
         self.indexlabels += range(npar)

   def load_chains(self, filename, **kwargs):
      '''This function should be defined by the subclass. It is responsible
      for loading the chains into the instance. The following must be
      defined:
         self.params:  list of strings that represent the P parameters
         self.chains[C,N,P]:  the chains C: chain #, N: iteration, P: param
         self.slices[P]:  slices into self.chains[:,:,:] for P
         self.parlabels[i]:  parameter label for parameter i
         self.indexlabels[i]:  parameter index for parameter i
      '''

   def __contains__(self, item):
      return item in self.params

   def __getitem__(self, key):
      if type(key) is type(""):
         try:
            return self.get_trace(key)
         except:
            raise KeyError, "Key not found:  "+key
      else:
         raise KeyError, "Key not found: "+str(key)

   def __getattr__(self, name):
      if 'params' not in self.__dict__:
         raise AttributeError
      if name in self.__dict__['params']:
         return self.get_trace(name)
      elif name == 'flatchains':
         return self.chains[:,:,:].reshape((-1,self.chains.shape[-1]))
      else:
         raise AttributeError, "No such attribute " + str(name)

   def get_trace(self,key, index=None, merge=False, do_reshape=True):
      '''Given a key name, we look for the proper trace.  If there
      are multiple chains and [merge]=True, concatenate the
      chains into one longer chain. '''
      if index is None:
         # See if an explicit index is asked for:
         if key.find(':') > 0:
            key,index = key.split(':')
            try:
               index = int(index)
            except:
               raise ValueError, "only integers can follow ':'"
            index = slice(index, index+1)
         else:
            index = slice(None)
      else:
         index = slice(index, index+1)

      if key not in self.params and key.split(':')[0] not in self.params:
         raise ValueError, "parameter not found: %s" % (key)
      if merge:
         # Now indexed by [iter,param]
         data = self.flatchains[:,self.slices[key]]
         if key in self.shapes and len(self.shapes[key]) > 1 and do_reshape:
            newshape = (data.shape[0],) + self.shapes[key]
            data = reshape(data, newshape)
            return data
         else:
            return data[:,index]
      else:
         # Now indexed by [chain,iter,param]
         data = self.chains[:,:,self.slices[key]]
         if key in self.shapes and len(self.shapes[key]) > 1 and do_reshape:
            newshape = (data.shape[0:2]) + self.shapes[key]
            data = reshape(data, newshape)
            return data
         else:
            return data[:,:,index]

   def get_trace0(self, key, index=None):
      '''Like get_trace, but imposes merge=True and squezes all 
      dimentions.'''
      tr = self.get_trace(key, index, merge=True)
      return squeeze(tr)

   def gelman_rubin(self, var, burn=0, imax=None, refined=True):
      '''compute the Gelman-Rubin statistic on variable [var].  Optionally,
      specify a burn-in [burn] and maximum index [imax].  You can choose the
      unrefined GR statistic by setting refined=False.'''
      arr = self.get_trace(var)
      if arr.shape[0] == 1:
         raise ValueError, "Error:  must have more than one chain"
      if arr.shape[2] == 1:
         return gelman_rubin.R(arr, caxis=0, iaxis=1)[0]
      return gelman_rubin.R(arr, caxis=0, iaxis=1)
      

   def mean(self, param, merge=True):
      '''Compute the weighted mean of the parameter.  If merge=True, then
      concatenate the chains first, otherwise return mean for each
      chain independently.'''
      ar = self.get_trace(param)
      if merge:
         ar = concatenate([ar[i] for i in range(ar.shape[0])])
         return squeeze(mean(ar, axis=0))
      else:
         return squeeze(mean(ar, axis=1))

   def variance(self, param, merge=True):
      '''compute the variance of the parameter.  If merge=True, then
      concatenate the chains first, otherwise return variance for each
      chain independently.'''
      ar = self.__getitem__(param)
      if merge:
         ar = concatenate([ar[i] for i in range(ar.shape[0])])
         return squeeze(var(ar, axis=0))
      else:
         return squeeze(var(ar, axis=1))

   def std(self, param, merge=True):
      return sqrt(self.variance(param, merge=merge))

   def mad(self, param, mode='median', merge=True):
      '''Compute the median absolute deviation of the parameter. If
      mrege=True, then concatenate chains first, otherwise return 
      mad for each chain independently.'''
      if mode == 'median':
         m = self.median(param, merge=merge)
      else:
         m = self.mean(param, merge=merge)
      m = atleast_1d(m)
      ar = self.get_trace(param, merge=merge)
      if merge:
         return median(absolute(ar - m[newaxis,:]), axis=0)
      else:
         return median(absolute(ar - m[:,newaxis,:]), axis=1)

   def median(self, param, merge=True):
      '''Compute the (weighted) median.  If merge is True, then concatenate
      the chains first, otherwise return median of each chain independently.'''
      ar = self.__getitem__(param)
      if merge:
         ar = concatenate([ar[i] for i in range(ar.shape[0])])
         return squeeze(median(ar, axis=0))
      else:
         return squeeze(median(ar, axis=1))

   def mode_plus_minus(self, param, index=None, conf=0.68, bins=30, 
         xrange=None, merge=True, func='poly', fit_log=True, 
         end_tol=[None,None], order=5,
         verbose=False):
      '''Compute the mode (maximum likelihood) of the posterior for
      [param] and also compute the interval (lower,upper) that contains
      the fraction [conf] of the probability.  [bins] controls how
      finely you want the PDF sampled.   Use [xrange] to restrict the
      range of the PDF.  If the lower or upper bounds
      cannot be found, they are set to None.'''
      ars = self.__getitem__(param)
      if ars.shape[2] == 1:
         scalar = True
      else:
         scalar = False

      mods = []
      pluss = []
      minuss = []

      if index is not None:
         ran = [index]
      else:
         ran = range(ars.shape[2])
      # Loop over variables
      for i in ran:
         mods.append([])
         pluss.append([])
         minuss.append([])
         subar = ars[:,:,i]
         if merge:
            subar = concatenate([subar[j] for j in range(subar.shape[0])])
            subar = subar.reshape((1,subar.shape[0]))
         for j in range(subar.shape[0]):
            mod,minus,plus = plus_minus(subar[j], bins, conf, xrange, 
                  func=func, fit_log=fit_log, end_tol=end_tol, debug=verbose,
                  order=order)
            mods[-1].append(mod)
            minuss[-1].append(minus)
            pluss[-1].append(plus)

      mods = array(mods)
      pluss = array(pluss)
      minuss = array(minuss)
      if scalar and merge:
         return (mods[0,0],pluss[0,0],minuss[0,0])
      elif scalar:
         return (mods[0,:],pluss[0,:],minuss[0,:])
      elif merge:
         return (mods[:,0],pluss[:,0],minuss[:,0])
      else:
         return (mods,pluss,minuss)

   def covar(self, params, center='mean', merge=True, chain=0):
      ''' Compute the covariance between [param] (list of names).  Specify 
      whether you want to use the mean or median for the central values by 
      specifying parameter [center].'''
      arrs = []
      for param in params:
         if not merge:
            arrs.append(self.get_trace(param)[chain])
         else:
            arrs.append(self.get_trace(param, merge=True))
              
      arrs = concatenate(arrs, axis=1)  # Now shape iter,param
      npar = arrs.shape[1]

      if center == 'mean':
         arm = mean(arrs, axis=0)
      else:
         arm = median(arrs, axis=0)

      darr = arrs - arm[newaxis,:]
      covar = zeros((npar,npar), float64)
      for i in range(npar):
         for j in range(i,npar):
            covar[i,j] = mean(darr[:,i]*darr[:,j])
            covar[j,i] = covar[i,j]
      return covar

   def correlation(self, params, center='mean', merge=True, chain=0):
      ''' Compute the correlation between [param] (list of names).  Specify 
      whether you want to use the mean or median for the central values by 
      specifying parameter [center].'''
      arrs = []
      for param in params:
         if not merge:
            arrs.append(self.get_trace(param)[chain])
         else:
            arrs.append(self.get_trace(param, merge=True))
              
      arrs = concatenate(arrs, axis=1)  # Now shape iter,param
      npar = arrs.shape[1]

      if center == 'mean':
         arm = mean(arrs, axis=0)
      else:
         arm = median(arrs, axis=0)

      N = arrs.shape[0]
      stds = sqrt(sum(power(arrs - arm[newaxis,:],2), axis=0)/(N-1))

      darr = arrs - arm[newaxis,:]
      corr = zeros((npar,npar), float64)
      for i in range(npar):
         for j in range(i,npar):
            corr[i,j] = mean(darr[:,i]*darr[:,j])/(stds[i]*stds[j])
            corr[j,i] = corr[i,j]
      return corr

   def corClick(self, event):
      fig = plt.figure(101)
      fig.clear()
      ax = fig.add_subplot(111)
      i = int(event.xdata)
      j = int(event.ydata)
      par1 = self.parlabels[i] + ":%d" % self.indexlabels[i]
      par2 = self.parlabels[j] + ":%d" % self.indexlabels[j]
      self.plot2dscatter(par1, par2, levels=[0.68], falpha=0.1)
      plt.draw()

   def plot_cor_matrix(self, params, center='mean', merge=True, chain=0,
         cut=None):

      boundaries = [0]
      midpts = []
      for par in params:
         npar = self.get_trace(par).shape[2]
         midpts.append(boundaries[-1]+npar/2)
         boundaries.append(boundaries[-1]+npar)
      corr = self.correlation(params, center=center, merge=merge, chain=chain)
      if cut:
         corr = where(greater(absolute(corr), cut), corr, 0)
      f = plt.figure(figsize=(8,8))
      ax = f.add_subplot(111)
      im = ax.imshow(corr, origin='lower left', interpolation='nearest')
      f.colorbar(im, ax=ax)
      for i in range(len(boundaries)):
         ax.axvline(boundaries[i]-0.5, color='k')
         ax.axhline(boundaries[i]-0.5, color='k')
      ax.set_xticks(midpts)
      ax.set_xticklabels(params, rotation='vertical')
      ax.set_yticks(midpts)
      ax.set_yticklabels(params)

      ax.format_coord = lambda x,y: "cov(%s[%d],%s[%d]) = %g" % \
            (self.parlabels[int(x)], self.indexlabels[int(x)],
             self.parlabels[int(y)], self.indexlabels[int(y)],
             corr[int(x),int(y)])
      plt.draw() 
      self.cid = f.canvas.mpl_connect('button_press_event', self.corClick)



   def plot_param(self, param, index=None, burn=0, thin=1, outfile=None, 
         close=False, bins=30, merge=True, chain=None):
      '''Plot the trace and histogram of the [param].  Optionally
      burn and thin by amounts given.  If outfile is given, save
      the figure as [outfile].  If [close], close the figure once
      it is saved.'''

      fig = plt.figure(figsize=(10,6))
      ax_trace = fig.add_axes([0.1,0.5,0.4,0.4])
      ax_acorr = fig.add_axes([0.1,0.1,0.4,0.3])
      ax_hist = fig.add_axes([0.6,0.1,0.3,0.8])
      ax_trace.set_xlabel('iteration', size=8)
      ax_trace.set_ylabel("%s" % param)
      ax_acorr.set_xlabel('lag')
      ax_acorr.set_ylabel('auto-correlation')

      arr = self.get_trace(param, index=index)
      if merge:
         arr = concatenate([arr[i] for i in range(arr.shape[0])])
         [ax_hist.hist(arr[burn::thin,i], bins=bins, histtype='step') \
               for i in range(arr.shape[1])]
         [ax_acorr.acorr(arr[burn::thin,i], maxlags=100, 
            detrend=mlab.detrend_mean) for i in range(arr.shape[1])]
         ax_trace.plot(arr[burn::thin])
      elif chain is not None:
         arr = arr[chain]
         [ax_hist.hist(arr[burn::thin,i], bins=bins, histtype='step') \
               for i in range(arr.shape[1])]
         [ax_acorr.acorr(arr[burn::thin,i], maxlags=100) \
               for i in range(arr.shape[1])]
         ax_trace.plot(arr[burn::thin])
      else:
         if arr.shape[2] != 1:
            raise ValueError, "Cannot plot multiple traces of a multi-valued variable"
         for i in range(arr.shape[0]):
            ax_trace.plot(arr[i,burn::thin,0])
            ax_hist.hist(arr[i,burn::thin,0], bins=bins, histtype='step')
            ax_acorr.acorr(arr[i, burn::thin,0], maxlags=100)

      #ax_hist.set_ylim(ax_trace.get_ylim())
      #ax_hist.set_xticks([])
      ax_hist.set_yticks([])
      
      if outfile is not None:
         fig.savefig(outfile)
      plt.draw()
      if close:
         plt.close(fig)


   def plot2dsurf(self, param1, param2, ax=None, xrange=None, yrange=None,
         bins=30, smooth=False, bfac=2, sfac=1., dd=3, cmap=cm.gray_r, 
         levels=[], ccolor='red', fill=False, ccmap=None, falpha=1.0,
         outfile=None, zorder=None):
      '''Plot up a 2D binned paramter plot for [param1] and [param2].
      if [ax] is supplied, use it to plot, otherwise, open up a new figure
      and axes.  You can specify [xrange] and [yrange].  [bins] will be
      passed to histogram2d.  If [smooth], the binned surface is smoothed
      using either a bivariate spline or a Gaussian Process (if pymc.gp is
      available).  If [cmap] is None, no image is drawn.  If [levels] is
      specified as fractions (0.68, 0.95, etc), draw the contours that
      enclose this fraction of the data.'''
      if ax is None:
         fig = plt.figure()
         ax = fig.add_subplot(111)
         own_ax = True
      else:
         own_ax = False

      #if ccmap is not None and ccolor is not None:
      #   # Cmap takes precedence
      #   ccolor = None

      tr1 = self.get_trace0(param1)
      tr2 = self.get_trace0(param2)
      if len(tr1.shape) != 1 or len(tr2.shape) != 1:
         raise RuntimeError, "Error, variables must be scalars, try using ':' notation"
      #tr1 = tr1[:,0]
      #tr2 = tr2[:,0]
      range = [[tr1.min(), tr1.max()],
               [tr2.min(), tr2.max()]]
      if xrange is not None:
         range[0] = list(xrange)
      if yrange is not None:
         range[1] = list(yrange)

      # first, bin up the data (all of it)
      grid,xs,ys = histogram2d(tr1,tr1, bins=bins,
            range=range)
      grid = grid.T*1.0
      xplot = linspace(xs[0], xs[-1], 101)
      yplot = linspace(ys[0], ys[-1], 101)
      extent = [xs[0], xs[-1], ys[0], ys[-1]]

      xs = (xs[1:] + xs[:-1])/2
      ys = (ys[1:] + ys[:-1])/2
      
      x,y = meshgrid(xs,ys)
      tx = xs[::bfac]
      ty = ys[::bfac]
      if smooth and not gp:
         tck = bisplrep(ravel(x), ravel(y), ravel(grid), task=-1, tx=tx, ty=ty)
         x = linspace(xs[0], xs[-1], 501)
         y = linspace(ys[0], ys[-1], 501)
         grid = bisplev(x, y, tck).T
      elif smooth and gp:
         M = gp.Mean(lambda x:  zeros(x.shape[:-1], dtype=float)+median(grid))
         scalerat = (tr2.max() - tr2.min())/(tr1.max()-tr1.min())
         C = gp.Covariance(gp.matern.aniso_euclidean, diff_degree=dd,
               scale=(tr1.max()-tr1.min())*sfac, amp=std(grid),
               scalerat=scalerat)
         x,y = meshgrid(xs,ys)
         mesh = vstack((ravel(x), ravel(y))).T
         gp.observe(M, C, obs_mesh=mesh, obs_vals=ravel(grid), obs_V=ravel(grid))
         dplot = dstack(meshgrid(xplot,yplot))
         grid,Vsurf = gp.point_eval(M,C,dplot)

      grid = where(grid < 0, 0, grid)

      if cmap:
         ax.imshow(grid, extent=extent,
               origin='lower', aspect='auto', interpolation='nearest',
               cmap=cmap)
      if levels:
         prob = ravel(grid)/sum(grid)
         sprob = sort(prob)
         cprob = 1.0 - cumsum(sprob)
         clevels = []
         for l in levels:
            id = nonzero(greater(cprob-l,0))[0][-1]
            clevels.append(sprob[id])
         prob.shape = grid.shape

         clevels.sort()
         norm = Normalize(clevels[0]*0.5,clevels[-1]*1.3)
         if fill:
            ax.contourf(prob, levels=clevels+[1],
               extent=extent, origin='lower',
               alpha=falpha, cmap=ccmap, norm=norm, zorder=zorder)
         ax.contour(prob, levels=clevels, colors=ccolor,
               extent=extent, origin='lower', linewidths=2, zorder=zorder)

      if own_ax:
         ax.set_xlabel("$%s$" % param1)
         ax.set_ylabel("$%s$" % param2)
         if xrange is not None:
            ax.set_xlim(xrange[0],xrange[1])
         if yrange is not None:
            ax.set_ylim(yrange[0],yrange[1])
         plt.draw()
         if outfile is not None:
            fig.savefig(outfile)
         return fig

   def plot2dscatter(self, param1, param2, param3=None, ax=None, xrange=None, 
         yrange=None, zrange=None, cmap=None, levels=[], ccolor='red', 
         falpha=1.0, outfile=None, zorder=None):
      '''Plot up a 2D scatter plot for [param1] and [param2], optinally colored
      by [param3] with colormap [cmap].  if [ax] is supplied, use it to plot,
      otherwise, open up a new figure and axes.  You can specify [xrange] and
      [yrange].  If [levels] is specified as fractions (0.68, 0.95, etc), draw
      the contours that enclose this fraction of the data.''' 
      
      if ax is None:
         fig = plt.figure() 
         ax = fig.add_subplot(111) 
         own_ax = True 
      else: 
         own_ax = False

      tr1 = squeeze(self.get_trace(param1, merge=True))
      tr2 = squeeze(self.get_trace(param2, merge=True))
      if param3 is not None:
         tr3 = self.get_trace(param3, merge=True)
      if len(tr1.shape) != 1 or len(tr2.shape) != 1:
         raise RuntimeError, "Error, variables must be scalars, try using ':' notation"
      range = [[tr1.min(), tr1.max()],
               [tr2.min(), tr2.max()]]
      if xrange is not None:
         range[0] = list(xrange)
      if yrange is not None:
         range[1] = list(yrange)

      if param3 is not None:
         if zrange is None:
            zrange = [tr3.min(), tr3.max()]
         ax.scatter(tr1, tr2, c=tr3, marker='.', linewidths=0, vmin=zrange[0],
               vmax=zrange[1], alpha=falpha)
         plt.colorbar(ax.collections[0])
      else:
         ax.plot(tr1, tr2, ',', color='k', alpha=falpha)

      if levels:
         for level in levels:
            i0,j0,major,minor,pa = error_ellipse(tr1, tr2, center='median',
                  level=level)
            ell = Ellipse(xy=(i0,j0), width=2*major,
               height=2*minor, angle=pa, facecolor='none', edgecolor=ccolor,
               zorder=100)
            ax.add_artist(ell)

      if own_ax:
         ax.set_xlabel("$%s$" % param1)
         ax.set_ylabel("$%s$" % param2)
         if xrange is not None:
            ax.set_xlim(xrange[0],xrange[1])
         if yrange is not None:
            ax.set_ylim(yrange[0],yrange[1])
         plt.draw()
         if outfile is not None:
            fig.savefig(outfile)
         return fig
   
   def triangle_plot(self, pars=None, **args):
      '''Plot a triangle.corner plot for list of parameters pars (or all of them
      if pars is None.'''
      if triangle is None:
         print "Error:  you need triangle_plot module to use this feature"
         return None
      vlabels = []
      arr = []
      if pars is None:
         pars = self.params
      for par in pars:
         x = self.get_trace(par, merge=True)
         arr.append(x)
         if x.shape[1] == 1:
            vlabels.append(par)
         else:
            vlabels = labels + [par+str(i) for i in range(x.shape[1])]
      arr = hstack(arr)
      if args.get('labels',None) is None:
         args['labels'] = vlabels
      plt = triangle.corner(arr, truths=median(arr, axis=0), **args)
      return plt



def error_ellipse(x, y, center='mean', level=0.68):
   '''compute the position, major axis, minor axis, and position angle
   of an ellipse that represents and [level] confidence interval.'''

   if center == 'mean':
      p1 = mean(x)
      p2 = mean(y)
   else:
      p1 = median(x)
      p2 = median(y)

   covar = covarc([x,y], center)


   eigval,eigvec = linalg.eig(covar)
   eigval = map(abs, eigval)

   dchi2 = brentq(lambda x:  chdtr(2,x) - level, 0, 100)

   if eigval[0] > eigval[1]:
      major = sqrt(eigval[0]*dchi2)
      minor = sqrt(eigval[1]*dchi2)
      pa = arctan2(eigvec[1,0], eigvec[0,0])*180.0/pi
   else:
      major = sqrt(eigval[1]*dchi2)
      minor = sqrt(eigval[0]*dchi2)
      pa = arctan2(eigvec[1,1], eigvec[0,1])*180.0/pi

   return (p1,p2,major,minor,pa)

def covarc(arrs, center):
      arrs = asarray(arrs).T  # Now shape iter,param
      npar = arrs.shape[1]

      if center == 'mean':
         arm = mean(arrs, axis=0)
      else:
         arm = median(arrs, axis=0)

      darr = arrs - arm[newaxis,:]
      covar = zeros((npar,npar), float64)
      for i in range(npar):
         for j in range(i,npar):
            covar[i,j] = mean(darr[:,i]*darr[:,j])
            covar[j,i] = covar[i,j]
      return covar

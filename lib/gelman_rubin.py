#!/usr/bin/env python
# This module computes the Gelman-Rubin diagnostic for a group of
# chains.
from numpy import *

def R(chains, caxis=0, iaxis=1, burn=0, imax=None, refined=True):
   '''Given [chains], a numpy array with at least 2 dimensions, compute the R
   statistic.  Chain index should be set with [caxis].  Chain iteration
   index should be set with [iaxis]'''
   chains = asarray(chains)
   if len(chains.shape) < 2:
      raise ValueError, "chains must be at least a 2D array"

   if burn > 0 or imax is not None:
      if burn > chains.shape[iaxis]:
         raise ValueError, "burn can't be more than %d" % chains.shape[iaxis]
      sl = [slice(None)]*len(shape(chains))
      sl[iaxis] = slice(burn,imax,None)
      sl = tuple(sl)
      chains = chains[sl]
   M = chains.shape[caxis]
   N = chains.shape[iaxis]


   s2_m = var(chains, axis=iaxis)  # Now 1-less Dim with chain as 
   theta_bar_m = mean(chains, axis=iaxis)
   if caxis > iaxis:
      # reduction in shape, the caxis is now one dimension over
      caxis -= 1
   theta_bar = mean(theta_bar_m, axis=caxis)
   # need a more complicated slice tuple with newaxis at caxis
   sl = [slice(None)]*len(shape(s2_m))
   sl[caxis] = newaxis
   sl = tuple(sl)
   W = mean(s2_m, axis=caxis)
   d_theta_bar = theta_bar_m-theta_bar[sl]
   B = N/(M-1)*sum(power(d_theta_bar,2), axis=caxis)
   Vhat = (N-1.0)/N*W + (M+1.0)/(N*M)*B
   
   if not refined:
      return Vhat/W

   theta_bar2 = mean(power(theta_bar_m,2),axis=caxis)
   cov1 = 1.0/(M - 1.0)*sum((s2_m - W[sl])*\
                            (power(theta_bar_m,2)-theta_bar2[sl]),axis=caxis)
   cov2 = 1.0/(M - 1.0)*sum((s2_m - W[sl])*d_theta_bar, axis=caxis)
   var_Vhat = ((N-1.0)/N)**2/M*var(s2_m,axis=caxis) + \
         2*((M+1.)/(N*M))**2/(M-1.0)*B**2
   var_Vhat += 2*(M+1.0)*(N-1.0)/N**2/M**2*(cov1 - 2*theta_bar*cov2)
   dhat = 2*Vhat**2/var_Vhat

   return sqrt((dhat+3)*Vhat/(dhat+1)/W)


def gelman_rubin(traces, vars=[], chain=-1, burn=0, imax=None, refined=True):
   '''Given a list of trace objects, compute the Gelman-Rubin diagnostic
   for each variable.  This is returned as a dictionary indexed by 
   the variable name.  Specify specific variables using the list [vars].'''

   if len(vars) == 0:
      vars = [v for v in traces[0].trace_names[chain] \
            if v.find('adaptive_scale_factor') < 0]
   else:
      for v in vars:
         if v not in traces[0].trace_names[chain]:
            raise ValueError, "variable %s not found in chain" % (v)
   Rs = {}
   for v in vars:
      chs = [tr.trace(v).gettrace(chain=chain) for tr in traces]
      try:
         chs = array(chs)
      except:
         raise ValueError, "Chains for %s don't have the same size" % v
      Rs[v] = R(chs, caxis=0, iaxis=1, burn=burn, imax=imax, refined=refined)
   return Rs

def gelman_rubin_trace(traces, N=100, vars=[], chain=-1, refined=True):
   '''Given a list of trace objects, compute the Gelman-Rubin diagnostic
   for each variable as a function of burn-in.  Specify the number of points 
   as [N].  [N] burn-in values spaced between 0 and iter/2 will be tried.  
   R values are retruned returned as a dictionary indexed by 
   the variable name.  Specify specific variables using the list [vars].'''

   if len(vars) == 0:
      vars = [v for v in traces[0].trace_names[chain] \
            if v.find('adaptive_scale_factor') < 0]
   else:
      for v in vars:
         if v not in traces[0].trace_names[chain]:
            raise ValueError, "variable %s not found in chain" % (v)
   Niter = traces[0].trace(vars[0]).gettrace(chain=chain).shape[0]
   imaxs = linspace(Niter/10,Niter,N)
   Rs = {}
   for v in vars:
      chs = [tr.trace(v).gettrace(burn=0, chain=chain) for tr in traces]
      try:
         chs = array(chs)
      except:
         raise ValueError, "Chains for %s don't have the same size" % v
      Rs[v] = array([R(chs, caxis=0, iaxis=1,imax=imax, refined=refined)\
            for imax in imaxs])
      #if len(shape(Rs[v])) == 2:
      #   Rs[v] = Rs[v].T
   return Rs,imaxs

if __name__ == "__main__":
   from pymc import database
   import sys
   from rdarg import rdarg
   try:
      from matplotlib import pyplot as plt
   except:
      plt = None

   if len(sys.argv) ==1 or '-h' in sys.argv:
      print '''
 Usage:  gelman_rubin [-chain chain] [-vars vars] 
              [-mode {report|burn|plot}] [-over over]
              [-Rmax Rmax] database-files
   compute the Gelman-Rubin statistic for parallel chains.
   where:  -chain:  specify which chain in multi-chain files
           -vars:   comma-separated list of variables to analyze
                    (default all)
           -mode:   report:  print out R-values for each var
                    burn:    compute how much burn needed to reach
                             R < Rmax
                    plot:    plot out R as a function of iteration
           -over:   print out only those vars with R > Rmax
           -Rmax:   Maximum value of R permissible (default 1.3)
      '''
      sys.exit(1)
           

   argv = sys.argv[1:]
   argv,chain = rdarg(argv, '-chain', int, -1)
   argv,vars = rdarg(argv, '-vars', None, None)
   argv,mode = rdarg(argv, '-mode', None, 'report')
   argv,over = rdarg(argv, '-over', int, 0, single=1)
   argv,Rmax = rdarg(argv, '-Rmax', float, 1.3)
   if vars is not None:
      vars = vars.split(',')

   trs = [database.hdf5.load(f) for f in argv]

   varlist = trs[0].trace_names[chain]
   for v in vars:
      if v not in varlist:
         id = vars.index(v)
         del vars[id]
         i = 0
         while (v+str(i) in varlist):
            vars.insert(id, v+str(i))
            i += 1

   if mode == 'report':
      Rs = gelman_rubin(trs, vars=vars, chain=chain)
      for v in Rs:
         if len(shape(Rs[v])) == 0:
            if over and Rs[v] > Rmax:
               print v + ":\t" + "%.4f" % Rs[v]
            else:
               print v + ":\t" + "%.4f" % Rs[v]
         else:
            if over:
               gids = greater(Rs[v], Rmax)
            else:
               gids = greater(Rs[v], 0)
            if sometrue(gids):
               print v + ":\t",
               for i in range(Rs[v].shape[0]):
                  if gids[i]:  print " %.4f" % Rs[v][i],
               print ''
   elif mode == 'burn':
      # find the minimum burn for each var
      minburn = 0
      Rs = gelman_rubin(trs, vars=vars, chain=chain)
      for v in Rs:
         R = atleast_1d(Rs[v])
         if alltrue(less(R, Rmax)):
            continue
         else:
            R,burns = gelman_rubin_trace(trs, vars=[v], chain=chain)[v]
            if shape(R) > 0:
               for i in range(R.shape[1]):
                  ids = nonzero(less(R[:,i], Rmax))[0]
                  if len(ids) == 0:
                     mR = R[:,i].min()
                     print "Warning!  %s[%d] has min R = %.4f" % (v,i,mR)
                  else:
                     minburn = max(minburn,burns[ids[0]])
                     print "%s[%d] requires burn > %d" % (v,i,ids[0])
            else:
               ids = nonzero(less(R,Rmax))[0]
               if len(ids) == 0:
                  print "Warning!  %s has min R = %.4f" % (v,R.min())
               else:
                  minburn = max(minburn,burns[ids[0]])
                  print "%s requires burn > %d" % (v,ids[0])
      print "Minimum burn-in needed:",minburn
   elif mode == 'plot':
      
      if plt is None:
         raise RuntimeError, "need matplotlib for plot mode"
      fig = plt.figure()
      ax = fig.add_subplot(111)
      R,burns = gelman_rubin_trace(trs, vars=vars, chain=chain, refined=False)
      for var in R: 
         ax.plot(burns,R[var])
      ax.axhline(1.03, color='red')
      ax.set_xlabel('burn')
      ax.set_ylabel('R')
      ax.set_title('Gelman-Rubin Statistic')

      plt.show()
         
      







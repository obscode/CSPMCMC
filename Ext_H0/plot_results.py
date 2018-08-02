#!/usr/bin/env python

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from numpy import *
import pickle
import sys,os,string
from myplotlib import PanelPlot,rcParams
import config
import STANstats
import get_data
import sigfig
import bspline2
import corner

def tomag(flux,eflux,zp):
   m = -2.5*log10(flux) + zp
   dm = eflux/flux*1.087
   return m,dm

rcParams['font.size'] = 16

cfg = config.config(sys.argv[1])
c = STANstats.STANchains(cfg.sampler.output)
data = c.data

# MCMC parameters
for var in c.params:
   locals()[var] = c.median(var)
   locals()['e_'+var] = c.std(var)
# Data from pickle file
for var in data:
   locals()[var] = data[var]

if not cfg.model.in_mag:
   m_sn,em_sn = tomag(m_sn, em_sn, 16.0)

if not cfg.model.HostMass:
   e_alpha = 0

# compute H0-derived distance modulus, but replace host-value distance modulus
#  where host = -1
DMz = where(less(host-1,0),
      5.0*log10(((1+zhelio)/(1+zcmb))*(3e5/H0)*(zcmb + zcmb**2*0.79))+25,
      DM[host-1])
if len(cfg.data.sn_filt) == 1:
   filt = ones(DMz.shape[0], dtype=int32)
   oid = arange(1,DMz.shape[0]+1).astype(int16)
   #M_sn = array([M_sn])
   #beta_st = array([beta_st])
   #gamma_st = array([gamma_st])
   a = array([a]).T


sts = linspace(st.min(), st.max(), 100)
if cfg.model.basis == 'poly':
   bs = array([power(sts-1,i) for i in range(cfg.model.order+1)]).T
else:
   bs = bspline2.bspline_basis(cfg.model.knots, sts, cfg.model.order, 
         gradient=cfg.model.gradient)

for i in range(len(cfg.data.sn_filt)):
   fids = equal(filt, i+1)
   oids = oid[fids] - 1

   model = g_modl[oids]
   #model = dot(Bs[oids,:], a[:,i]) - 19 + DMz[oids] + g_Al[fids]
   #model = M_sn[i] + DMz[oids] + beta_st[i]*(st[oids] - 1) + \
   #      gamma_st[i]*(st[oids] - 1)**2 + g_Al[fids]
   resids = m_sn[fids] - model
   of = open('resids_%s.dat' % (cfg.data.sn_filt[i]), 'w')
   print >> of, "#       SN   zcmb    st    mod   m_sn    res     Al"
   for j in range(len(model)):
      print >> of, "%10s %.4f %5.3f %6.3f %6.3f %6.3f %6.3f" % \
            (snnames[oids[j]],zcmb[oids[j]],st[oids[j]],model[j],
             m_sn[fids][j],resids[j], g_Al[fids][j])
   of.close()

   cids = greater(host[oids], 0)
   csp_ids = array([source[j].find('CSP')==0 for j in oids
      ])
   plt.plot(st[oids], resids, '.', color='k')
   plt.plot(st[oids][csp_ids], resids[csp_ids], 'o', mfc='white', mec='k')
   plt.plot(st[oids][cids], resids[cids], 'o', mfc='red', mec='none',
         zorder=50)
   plt.xlabel('stretch')
   plt.ylabel('residuals')
   plt.ylim(-1.5,1.5)
   plt.tight_layout()
   plt.savefig('resids_%s.eps' % (cfg.data.sn_filt[i]))

   plt.figure()
   zs = linspace(0.0001, zcmb.max(),100)
   if len(cfg.data.sn_filt) == 1:
      zerr = sqrt(eps_sn**2 + power(2.17*0.001/zs,2))
   else:
      zerr = sqrt(eps_sn[i]**2 + power(2.17*0.001/zs,2))
   plt.plot(zcmb[oids], resids, '.', color='k')
   plt.plot(zcmb[oids][csp_ids], resids[csp_ids], 'o', mfc='white', mec='k')
   plt.plot(zcmb[oids][cids], resids[cids], 'o', mfc='red', mec='none', 
         zorder=50)
   plt.plot(zs, zerr, '-', color='red')
   plt.plot(zs, -zerr, '-', color='red')
   plt.ylim(-3,3)
   plt.xlabel('z')
   plt.ylabel('residuals')
   plt.ylim(-1.5,1.5)
   plt.tight_layout()
   plt.savefig('resids_z_%s.eps' % (cfg.data.sn_filt[i]))
   plt.close()

   plt.figure()
   masses = 0.4*(dist[oids] - K[oids]) + 1.04
   resids = resids + alpha*(masses - M0)

   xx = linspace(masses.min(), masses.max(), 100)
   plt.errorbar(masses, resids, fmt='o', xerr=cfg.model.HostMassSigma,
         yerr=sqrt(e_g_modl[oids]**2 + em_sn[fids]**2), capsize=0)
   plt.scatter(masses, resids, marker='o', c=zcmb[oids],zorder=100)
   plt.plot(masses[cids], resids[cids], 'o', mec='k', mfc='red', zorder=200)
   plt.plot(xx, alpha*(xx-M0), '-', color='red', zorder=200)
   plt.plot(xx, (alpha+e_alpha)*(xx-M0), '--', color='red', zorder=200)
   plt.plot(xx, (alpha-e_alpha)*(xx-M0), '--', color='red', zorder=200)
   plt.colorbar()
   plt.xlabel('$\log_{10}(Mass/M_\odot)$')
   plt.ylabel('residuals')
   plt.tight_layout()
   plt.ylim(-1.5, 1.5)
   plt.savefig('resids_mass_%s.eps' % (cfg.data.sn_filt[i]))


   plt.figure()
   abs_m = m_sn[fids] - DMz[oids]
   m_corr = abs_m - g_Al[fids]
   if len(cfg.data.sn_filt) > 1:
      models = dot(c.get_trace('a', merge=True)[:,:,i], bs.T) - 19
   else:
      models = dot(c.get_trace('a', merge=True), bs.T) - 19
   model = median(models, axis=0)
   emodel = std(models, axis=0)
   plt.plot(st[oids], m_corr, '.', color='k')
   plt.plot(st[oids][csp_ids], m_corr[csp_ids], 'o', color='blue')
   plt.plot(st[oids][cids], m_corr[cids], 'o', mec='k', mfc='red')
   #plt.plot(st[oids], abs_m, '.', color='red')
   #plt.plot(st[oids][csp_ids], abs_m[csp_ids], 'o', color='red')
   plt.plot(sts, model, '-', color='red')
   plt.plot(sts, model+emodel, '--', color='red')
   plt.plot(sts, model-emodel, '--', color='red')
   plt.gca().invert_yaxis()
   plt.xlabel('$s_{BV}$')
   plt.ylabel('$M_{corr}$')
   plt.tight_layout()
   plt.savefig('Phillips_%s.eps' % (cfg.data.sn_filt[i]))
   plt.close()


a_s = c.get_trace('a', merge=True)
eps_sn_s = c.get_trace('eps_sn', merge=True)
H0_s = c.get_trace('H0', merge=True)
# SN parameters
if len(cfg.data.sn_filt) > 1:
   for i in range(len(cfg.data.sn_filt)):
      arr = [a_s[:,j,i] for j in range(a_s.shape[1])]
      arr += [eps_sn_s[:,i], H0_s[:,0]]
      arr = array(arr).T

      labs = ['$a_%d$' % j for j in range(a_s.shape[1])]
      labs = labs + ['$\epsilon_{sn}$','$H_0$']
      #tp1 = c.triangle_plot(['M_sn:%d' % i, 'beta_st:%d' % i, 'gamma_st:%d' % i,
      #                       'eps_sn:%d' % i, 'H0'])
      tp1 = corner.corner(arr, labels=labs, truths=median(arr, axis=0))
      tp1.savefig('SN_triangle_%s.eps' % cfg.data.sn_filt[i])
      plt.close(tp1)
else:
   arr = [a_s[:,j] for j in range(a_s.shape[1])]
   arr = arr + [eps_sn_s[:,0], H0_s[:,0]]
   arr = array(arr).T
   labs = ['$a_%d$' % j for j in range(a_s.shape[1])]
   labs = labs + ['$\epsilon_{sn}$','$H_0$']
   tp1 = corner.corner(arr, labels=labs, truths=median(arr,axis=0))
   tp1.savefig('SN_triangle.eps')
   plt.close(tp1)

# Now we output some tables.
fout = open('results_table.txt','w')
fout.write("SNe\n")
fout.write("---\n")
fout.write('v_pec:  %s +/- %s\n' % sigfig.round_sig_error(pec_vel,e_pec_vel,2))
headers = ['filter']
headers = headers + ["a[%d]" % i for i in range(Nbasis)]
headers = headers + ["eps"]
if len(cfg.data.sn_filt) == 1:
   cols = [[a[i]] for i in range(Nbasis)] + [[eps_sn]]
   errors = [[e_a[i]] for i in range(Nbasis)] + [[ e_eps_sn ]]
   labels = [cfg.data.sn_filt]
else:
   cols = [a[i,:] for i in range(Nbasis)] + [eps_sn]
   errors = [e_a[i,:] for i in range(Nbasis)] + [e_eps_sn]
   labels = cfg.data.sn_filt
lines = sigfig.format_table(cols=cols, errors=errors,
                       n=2, labels=labels, headers=headers)
[fout.write(line+"\n") for line in lines]

# The answer!
fout.write("H_0:  %s  +/-  %s\n" % sigfig.round_sig_error(H0,e_H0,2))

fout.close()


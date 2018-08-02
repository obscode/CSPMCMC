#!/usr/bin/env python

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from numpy import *
import pickle
import sys,os,string
from astropy.io import ascii
from myplotlib import PanelPlot
import config
import STANstats
import get_data
import sigfig
try:
   import corner
except:
   corner = None

def MAD(a):
   '''return the median absolute deviation *1.48'''
   return 1.48*median(absolute(a-median(a)))

def RMS(a):
   '''Return the root-mean-square'''
   return sqrt(mean(power(a - median(a),2)))

def tomag(flux,eflux,zp):
   m = -2.5*log10(flux) + zp
   dm = eflux/flux*1.087
   return m,dm

def toflux(mag,emag,zp):
   flux = power(10, -0.4*(mag-zp))
   eflux = emag*flux/1.087
   return flux,eflux


cfg = config.config(sys.argv[1])
with open(cfg.sampler.output) as f:
   d = pickle.load(f)
c = STANstats.STANchains(chains=d['samples'], flatnames=d['flatnames'])
if not cfg.model.NGauss:
   cfg.model.NGauss = 1

# MCMC parameters
for var in c.params:
   locals()[var] = c.median(var)
   locals()['e_'+var] = c.std(var)
# Data from pickle file
for var in d['data']:
   locals()[var] = d['data'][var]

if not cfg.model.in_mag:
   flux,e_flux = mag,e_mag
   mag,e_mag = tomag(mag, e_mag, 25.0)
   flux4258,e_flux4258 = mag4258, e_mag4258
   mag4258,e_mag4258 = tomag(mag4258, e_mag4258, 23.0)
   fluxLMC,e_fluxLMC = magLMC,e_magLMC
   magLMC,e_magLMC = tomag(magLMC, e_magLMC, 12.0)
   fluxMW,e_fluxMW = magMW, e_magMW
   magMW,e_magMW = tomag(magMW, e_magMW, 3.0)

w_P = []
w_VI = []
w_OH = []
w_res = []
s_res = []
w_labs = []

Pmin,Pmax = inf,-inf
VImin,VImax = inf,-inf
OHmin,OHmax = inf,-inf
fig1 = PanelPlot(1,2, pwidths=[1], pheights=[0.2, 0.8], figsize=(10,6))
fig2 = PanelPlot(1,2, pwidths=[1], pheights=[0.2, 0.8], figsize=(10,6))
fig3 = PanelPlot(1,2, pwidths=[1], pheights=[0.2, 0.8], figsize=(10,6))

aresids = []
alabels = []

if cfg.model.use_MW:
   if 'betaVI_MW' not in c.params:
      betaVI_MW = betaVI
      e_betaVI_MW = e_betaVI
   dist = -5*log10(pi_true) - 5
   model = dist + betaVI*VI_MW + betaP*P_MW + betaOH*(OH_MW-9.5) + M
   resids = magMW - model
   if cfg.model.in_mag:
      aresids.append(resids)
   else:
      r = fluxMW - toflux(model, 0, 3.0)[0]
      aresids.append(r/RMS(r))
   alabels.append('MW')
   w_res.append(median(resids))
   s_res.append(std(resids))
   OH_MW = array([OH_MW]*len(VI_MW))
   w_P.append(median(P_MW))
   w_VI.append(median(VI_MW))
   w_OH.append(median(OH_MW))
   w_labs.append('MW')
   fig1.axes[0].plot(P_MW, resids, 'o', color='blue')
   fig2.axes[0].plot(VI_MW,resids, 'o', color='blue')
   fig3.axes[0].plot([OH_MW]*len(VI_MW),resids, 'o', color='blue')
   fig1.axes[1].plot(P_MW, magMW - dist - betaVI_MW*VI_MW - betaOH*(OH_MW-9.5), 
                     'o', color='blue', label='MW')
   fig2.axes[1].plot(VI_MW, magMW - dist - betaP*P_MW - betaOH*(OH_MW-9.5), 
                     'o', color='blue', label='MW')
   fig3.axes[1].plot(OH_MW, magMW - dist - betaP*P_MW - betaVI_MW*VI_MW, 
                     'o', color='blue', label='MW')
   fig1.axes[0].axhline(eps_MW, linestyle='--', color='blue')
   fig1.axes[0].axhline(-eps_MW, linestyle='--', color='blue')
   Pmin = min(Pmin, P_MW.min())
   Pmax = max(Pmax, P_MW.max())
   VImin = min(VImin, VI_MW.min())
   VImax = max(VImax, VI_MW.max())
   OHmin = min(OHmin, OH_MW.min())
   OHmax = max(OHmax, OH_MW.max())

if cfg.model.use_LMC:
   if 'betaVI_LMC' not in c.params:
      betaVI_LMC = betaVI
      e_betaVI_LMC = e_betaVI
   model = DM_LMC + betaVI_LMC*VI_LMC + betaP*P_LMC + \
            betaOH*(OH_LMC-9.5) + M
   resids = magLMC - model
   w_res.append(median(resids))
   s_res.append(std((resids)))
   if cfg.model.in_mag:
      aresids.append(resids)
   else:
      r = fluxLMC - toflux(model,0.0, 12.0)[0]
      aresids.append(r/RMS(r))
   alabels.append('LMC')
   OH_LMC = array([OH_LMC]*len(P_LMC))
   w_P.append(median(P_LMC))
   w_VI.append(median(VI_LMC))
   w_OH.append(median(OH_LMC))
   w_labs.append('LMC')
   fig1.axes[0].plot(P_LMC, resids, 's', color='k')
   fig2.axes[0].plot(VI_LMC, resids, 's', color='k')
   fig3.axes[0].plot(OH_LMC, resids, 's', color='k')
   fig1.axes[1].plot(P_LMC, resids + betaP*P_LMC + M, 's', color='k', 
                     label='LMC')
   fig2.axes[1].plot(VI_LMC, resids + betaVI_LMC*VI_LMC + M, 's', color='k',
         label='LMC')
   fig3.axes[1].plot(OH_LMC, resids + betaOH*(OH_LMC-9.5) + M, 's', color='k',
         label='LMC')
   fig1.axes[0].axhline(eps_LMC, linestyle='--', color='k')
   fig1.axes[0].axhline(-eps_LMC, linestyle='--', color='k')
   Pmin = min(Pmin, P_LMC.min())
   Pmax = max(Pmax, P_LMC.max())
   VImin = min(VImin, VI_LMC.min())
   VImax = max(VImax, VI_LMC.max())
   OHmin = min(OHmin, OH_LMC.min())
   OHmax = max(OHmax, OH_LMC.max())

if cfg.model.use_4258:
   if 'betaVI_4258' not in c.params:
      betaVI_4258 = betaVI
      e_betaVI_4258 = e_betaVI
   model = DM_4258 + betaVI_4258*VI_4258 + betaP*P_4258 + \
            betaOH*(OH_4258-9.5) + M
   resids = mag4258 - model
   if cfg.model.in_mag:
      aresids.append(resids)
   else:
      r = flux4258 - toflux(model, 0.0, 23.0)[0]
      aresids.append(r/RMS(r))
   alabels.append('N4258')
   w_res.append(median(resids))
   s_res.append(std((resids)))
   w_P.append(median(P_4258))
   w_VI.append(median(VI_4258))
   w_OH.append(median(OH_4258))
   w_labs.append('4258')
   fig1.axes[0].plot(P_4258, resids, '^', color='red')
   fig2.axes[0].plot(VI_4258, resids, '^', color='red')
   fig3.axes[0].plot(OH_4258, resids, '^', color='red')
   fig1.axes[1].plot(P_4258, resids + betaP*P_4258 + M, '^', color='red',
         label='4258')
   fig2.axes[1].plot(VI_4258, resids + betaVI_4258*VI_4258+M, '^', color='red',
         label='4258')
   fig3.axes[1].plot(OH_4258, resids + betaOH*(OH_4258-9.5)+M, '^', color='red',
         label='4258')
   if cfg.model.NGauss > 1:
      #mid = argmax(c.median('theta_4258'))
      mid = argmax(c.median('theta'))
      #fig1.axes[0].axhline(eps_4258[mid], linestyle='--', color='red')
      #fig1.axes[0].axhline(-eps_4258[mid], linestyle='--', color='red')
      fig1.axes[0].axhline(eps[mid], linestyle='--', color='red')
      fig1.axes[0].axhline(-eps[mid], linestyle='--', color='red')
   else:
      #fig1.axes[0].axhline(eps_4258, linestyle='--', color='red')
      #fig1.axes[0].axhline(-eps_4258, linestyle='--', color='red')
      fig1.axes[0].axhline(eps, linestyle='--', color='red')
      fig1.axes[0].axhline(-eps, linestyle='--', color='red')
   Pmin = min(Pmin, P.min())
   Pmax = max(Pmax, P.max())
   VImin = min(VImin, VI_4258.min())
   VImax = max(VImax, VI_4258.max())
   OHmin = min(OHmin, OH_4258.min())
   OHmax = max(OHmax, OH_4258.max())

#xx1 = array([Pmin,Pmax])
xx1 = linspace(Pmin,Pmax, 100)
xx2 = array([VImin,VImax])
xx3 = array([OHmin,OHmax])
fig1.axes[1].plot(xx1, M + betaP*xx1, '-', color='k')
fig3.axes[1].plot(xx3, M + betaOH*(xx3-9.5), '-', color='k')
if cfg.model.use_MW:
   fig2.axes[1].plot(xx2, M + betaVI_MW*xx2, '-', color='blue')
if cfg.model.use_LMC:
   fig2.axes[1].plot(xx2, M + betaVI_LMC*xx2, '-', color='k')
if cfg.model.use_4258:
   fig2.axes[1].plot(xx2, M + betaVI_4258*xx2, '-', color='red')
fig1.axes[0].axhline(0, linestyle='-', color='k')
fig2.axes[0].axhline(0, linestyle='-', color='k')
fig3.axes[0].axhline(0, linestyle='-', color='k')
fig1.axes[0].set_xlabel(r'$\log_{10}\left(P\right)$')
fig2.axes[0].set_xlabel('$V-I$')
fig3.axes[0].set_xlabel('$[O/H]$')
fig1.axes[0].set_ylabel('resids')
fig1.axes[1].set_ylabel('corrected mag')
fig2.axes[0].set_ylabel('resids')
fig2.axes[1].set_ylabel('corrected mag')
fig3.axes[0].set_ylabel('resids')
fig3.axes[1].set_ylabel('corrected mag')
fig1.axes[1].legend()
fig2.axes[1].legend()
fig3.axes[1].legend()
plt.draw()
fig1.set_limits()
fig1.draw()
fig2.set_limits()
fig2.draw()
fig3.set_limits()
fig3.draw()
fig1.fig.savefig('anchors_P.pdf')
fig2.fig.savefig('anchors_VI.pdf')
fig3.fig.savefig('anchors_OH.pdf')
plt.close(fig1.fig)
plt.close(fig2.fig)
plt.close(fig3.fig)


symbs = ['o','s','^','d','p','v']*5
cols = ['k']*6+['red']*6+['blue']*6+['green']*6+['orange']*6
Pmin,Pmax = inf,-inf
VImin,VImax = inf,-inf
OHmin,OHmax = inf,-inf
fig1 = PanelPlot(1,2, pwidths=[1], pheights=[0.2, 0.8], figsize=(10,6))
fig2 = PanelPlot(1,2, pwidths=[1], pheights=[0.2, 0.8], figsize=(10,6))
fig3 = PanelPlot(1,2, pwidths=[1], pheights=[0.2, 0.8], figsize=(10,6))

if len(shape(betaVI)) == 0:
   betaVI = ones((S,))*betaVI
   e_betaVI = ones((S,))*e_betaVI

sresids = []

for i in range(S):
   figi = PanelPlot(1,2, pwidths=[1], pheights=[0.2,0.8], figsize=(10,6))
   gids = equal(ID, i+1)
   model = DM[i] + betaVI[i]*VI[gids] + betaP*P[gids] + \
            betaOH*(OH[gids]-9.5) + M
   resids = mag[gids] - model
   if cfg.model.in_mag:
      sresids.append(resids)
   else:
      r = flux[gids] - toflux(model, 0.0, 25.0)[0]
      sresids.append(r/RMS(r))
   w_res.append(median(resids))
   s_res.append(std((resids)))
   w_P.append(median(P[gids]))
   w_VI.append(median(VI[gids]))
   w_OH.append(median(OH[gids]))
   w_labs.append(cephlist[i])
   fig1.axes[0].plot(P[gids],resids, symbs[i], color=cols[i])
   figi.axes[0].plot(P[gids],resids, 'o', color='k')
   fig2.axes[0].plot(VI[gids],resids, symbs[i], color=cols[i])
   fig3.axes[0].plot(OH[gids],resids, symbs[i], color=cols[i])
   fig1.axes[1].plot(P[gids], resids + betaP*P[gids] + M, symbs[i], 
         color=cols[i], label=cephlist[i])
   figi.axes[1].plot(P[gids], resids + betaP*P[gids] + M, 'o', 
         color='k', label=cephlist[i])
   reals = c.get_trace('M', merge=True)[newaxis,:,0] + \
          (c.get_trace('DM', merge=True)[newaxis,:,i]-DM[i]) + \
          c.get_trace('betaP', merge=True)[newaxis,:,0]*xx1[:,newaxis]
   mreals = median(reals, axis=1)
   sreals = std(reals, axis=1)
   fig2.axes[1].plot(VI[gids], resids + betaVI[i]*VI[gids] + M, symbs[i], 
         color=cols[i], label=cephlist[i])
   fig3.axes[1].plot(OH[gids], resids + betaOH*(OH[gids]-9.5) + M, symbs[i], 
         color=cols[i], label=cephlist[i])
   figi.axes[0].axhline(0, linestyle='-', color='red')
   figi.axes[0].plot(xx1, sreals, '--', color='red')
   figi.axes[0].plot(xx1, -sreals, '--', color='red')
   #figi.axes[0].axhline(eps, linestyle='--', color='red')
   #figi.axes[1].plot(xx1, M + betaP*xx1, '-', color='red')
   figi.axes[1].plot(xx1, mreals, '-', color='red')
   figi.axes[1].plot(xx1, mreals+sreals, '--', color='red')
   figi.axes[1].plot(xx1, mreals-sreals, '--', color='red')
   figi.axes[0].set_xlabel(r'$\log_{10}\left(P\right)$')
   figi.axes[0].set_ylabel('resids')
   figi.axes[1].set_ylabel('corrected mag')
   figi.axes[1].legend(fontsize=8)
   figi.set_limits()
   figi.draw()
   figi.fig.savefig('SN_hosts_P_%s.pdf' % cephlist[i])
   plt.close(figi.fig)
   Pmin = min(Pmin, P.min())
   Pmax = max(Pmax, P.max())
   VImin = min(VImin, VI.min())
   VImax = max(VImax, VI.max())
   OHmin = min(OHmin, OH.min())
   OHmax = max(OHmax, OH.max())
xx1 = array([Pmin,Pmax])
xx2 = array([VImin,VImax])
xx3 = array([OHmin,OHmax])

if 'm_betaVI' not in d:
   m_betaVI = mean(betaVI)

fig1.axes[1].plot(xx1, M + betaP*xx1, '-', color='k')
fig2.axes[1].plot(xx2, M + m_betaVI*xx2, '-', color='k')
fig3.axes[1].plot(xx3, M + betaOH*(xx3-9.5), '-', color='k')
fig1.axes[0].axhline(0, linestyle='-', color='k')
fig2.axes[0].axhline(0, linestyle='-', color='k')
fig2.axes[0].axhline(0, linestyle='-', color='k')
if cfg.model.NGauss > 1:
   aeps = eps[argmax(c.median('theta'))]
else:
   aeps = eps
fig1.axes[0].axhline(aeps, linestyle='--', color='k')
fig1.axes[0].axhline(-aeps, linestyle='--', color='k')
fig2.axes[0].axhline(aeps, linestyle='--', color='k')
fig2.axes[0].axhline(-aeps, linestyle='--', color='k')
fig3.axes[0].axhline(aeps, linestyle='--', color='k')
fig3.axes[0].axhline(-aeps, linestyle='--', color='k')
fig1.axes[0].set_xlabel(r'$\log_{10}\left(P\right)$')
fig2.axes[0].set_xlabel('$V-I$')
fig3.axes[0].set_xlabel('$[O/H]$')
fig1.axes[0].set_ylabel('resids')
fig1.axes[1].set_ylabel('corrected mag')
fig2.axes[0].set_ylabel('resids')
fig2.axes[1].set_ylabel('corrected mag')
fig3.axes[0].set_ylabel('resids')
fig3.axes[1].set_ylabel('corrected mag')
fig1.axes[1].legend(fontsize=8)
fig2.axes[1].legend(fontsize=8)
fig3.axes[1].legend(fontsize=8)
plt.draw()
fig1.set_limits()
fig1.draw()
fig2.set_limits()
fig2.draw()
fig3.set_limits()
fig3.draw()
fig1.fig.savefig('SN_hosts_P.pdf')
fig2.fig.savefig('SN_hosts_VI.pdf')
fig3.fig.savefig('SN_hosts_OH.pdf')
plt.close(fig1.fig)
plt.close(fig2.fig)
plt.close(fig3.fig)

# Redidual histogram
sresids = concatenate(sresids)
aresids.append(sresids)
alabels.append('SN Hosts')

fig = plt.figure()
ax = fig.add_subplot(111)
ax.hist(aresids, label=alabels, histtype='step', stacked=True, normed=True,
        bins=100, linewidth=2)
ax.set_xlabel('Model residuals', fontsize=16)
ax.set_xlim(-5,5)
ax.legend()
plt.tight_layout()
fig.savefig('resids_hist.pdf')

fig3 = plt.figure()         # to told weighted average of residuals
fig4 = plt.figure()         # to told weighted average of residuals
fig5 = plt.figure()         # to told weighted average of residuals
ax3 = fig3.add_subplot(111)
ax4 = fig4.add_subplot(111)
ax5 = fig5.add_subplot(111)
for i in range(len(w_res)):
   ax3.errorbar([w_P[i]], [w_res[i]], fmt=symbs[i], color=cols[i], 
         label=w_labs[i], yerr=s_res[i], capsize=0, ms=10)
   ax4.errorbar([w_VI[i]], [w_res[i]], fmt=symbs[i], color=cols[i], 
         label=w_labs[i], yerr=s_res[i], capsize=0, ms=10)
   ax5.errorbar([w_OH[i]], [w_res[i]], fmt=symbs[i], color=cols[i], 
         label=w_labs[i], yerr=s_res[i], capsize=0, ms=10)
lgd1 = ax3.legend(fontsize=8, loc=3, ncol=4, bbox_to_anchor=(0.,1.02,1.,0.102),
      mode='expand')
lgd2 = ax4.legend(fontsize=8, loc=3, ncol=4, bbox_to_anchor=(0.,1.02,1.,0.102), 
      mode='expand')
lgd3 = ax5.legend(fontsize=8, loc=3, ncol=4, bbox_to_anchor=(0.,1.02,1.,0.102), 
      mode='expand')
ax3.set_xlabel(r'$\log_{10}\left(P\right)$')
ax4.set_xlabel(r'$V-I$')
ax5.set_xlabel(r'$[O/H]$')
ax3.set_ylabel(r'median residuals')
ax4.set_ylabel(r'median residuals')
ax5.set_ylabel(r'median residuals')
fig3.savefig('ceph_res_comb_P.pdf', bbox_extra_artists=(lgd1,), bbox_inches='tight')
fig4.savefig('ceph_res_comb_VI.pdf', bbox_extra_artists=(lgd2,), bbox_inches='tight')
fig5.savefig('ceph_res_comb_OH.pdf', bbox_extra_artists=(lgd3,), bbox_inches='tight')
plt.close(fig3)
plt.close(fig4)
plt.close(fig5)

# Now, let's do some triangle plots and output the parameters of interest.
# Cepheid parameters:
if corner is not None:
   tp1 = c.triangle_plot(['M','betaP','betaVI','betaOH'])
   tp1.savefig('Ceph_triangle.pdf')
else:
   print "Warning:  corner is not installed, so no triangle plots. To install:"
   print "pip install corner"


# Now we output some tables.
fout = open('results_table.txt','w')
fout.write("Cepheids\n")
fout.write("--------\n")
fout.write('M:      %s +/- %s\n' % sigfig.round_sig_error(M,e_M,2))
fout.write('betaP:  %s +/- %s\n' % sigfig.round_sig_error(betaP,e_betaP,2))
fout.write('betaOH:  %s +/- %s\n' % sigfig.round_sig_error(betaOH,e_betaOH,2))

hosts = []
headers = ['Host','DM','betaVI']
headers += ['eps%d' % i for i in range(cfg.model.NGauss)]
cols = [[],[]]
ecols = [[],[]]
for i in range(cfg.model.NGauss):
   cols.append([])
   ecols.append([])

if cfg.model.use_MW:
   hosts.append('MW')
   cols[0].append(-1); ecols[0].append(-1)
   cols[1].append(betaVI_MW); ecols[1].append(e_betaVI_MW)
   cols[2].append(eps_MW); ecols[2].append(e_eps_MW)
   for i in range(1,cfg.model.NGauss):
      cols[i+2].append(-1)
      ecols[i+2].append(-1)
if cfg.model.use_LMC:
   hosts.append('LMC')
   cols[0].append(DM_LMC);     ecols[0].append(e_DM_LMC)
   cols[1].append(betaVI_LMC); ecols[1].append(e_betaVI_LMC)
   cols[2].append(eps_LMC);    ecols[2].append(e_eps_LMC)
   for i in range(1,cfg.model.NGauss):
      cols[i+2].append(-1)
      ecols[i+2].append(-1)
if cfg.model.use_4258:
   hosts.append('4258')
   cols[0].append(DM_4258);     ecols[0].append(e_DM_4258)
   cols[1].append(betaVI_4258); ecols[1].append(e_betaVI_4258)
   if cfg.model.NGauss > 1:
      for i in range(cfg.model.NGauss):
         #cols[2+i].append(eps_4258[i]);
         #ecols[2+i].append(e_eps_4258[i])
         cols[2+i].append(eps[i]);
         ecols[2+i].append(e_eps[i])
   else:
      #cols[2].append(eps_4258);
      #ecols[2].append(e_eps_4258)
      cols[2].append(eps);
      ecols[2].append(e_eps)

hosts += cephlist
cols[0] = concatenate([cols[0], DM]); ecols[0] = concatenate([ecols[0], e_DM]) 
cols[1] = concatenate([cols[1], betaVI]); ecols[1] = concatenate([ecols[1], 
                                                     e_betaVI]) 
if cfg.model.NGauss > 1:
   for i in range(cfg.model.NGauss):
      cols[2+i] = concatenate([cols[2+i], [eps[i]]*S])
      ecols[2+i] = concatenate([ecols[2+i], [e_eps[i]]*S]) 
else:
   cols[2] = concatenate([cols[2], [eps]*S])
   ecols[2] = concatenate([ecols[2], [e_eps]*S])

lines = sigfig.format_table(cols=cols, errors=ecols, n=2, headers=headers, 
        labels=hosts)
[fout.write(line+"\n") for line in lines]

fout.close()

# Final covariance matrix. We need to be a bit more robust here
DMs = c.get_trace('DM', merge=True)
devs = absolute(DMs - c.median('DM', merge=True)[newaxis,:])
# Do a 5-sigma clip to get rid fo really diviant points for NGC 4424
gids = less(devs, 5*1.4826*c.mad('DM', merge=True)[newaxis,:])
gids = product(gids, axis=1).astype(bool)
C = cov(DMs[gids,:].T)

DMs = c.median('DM')
eDMs = c.std('DM')
f = open('DM_cov.dat','w')
[f.write(ceph+" ") for ceph in cephlist]
f.write('\n')
[f.write("%f " % DM) for DM in DMs]
f.write('\n')
for i in range(C.shape[0]):
   for j in range(C.shape[1]):
      f.write("%f " % C[i,j])
   f.write('\n')

# Now make a nice figure
fig = plt.figure()
ax = fig.add_subplot(111)
sids = argsort(diag(C))
names = [cephlist[i] for i in sids]
CC = zeros(C.shape)
for i in range(C.shape[0]):
   for j in range(C.shape[1]):
      CC[i,j] = C[sids[i],sids[j]]
img = ax.imshow(CC, interpolation='nearest', origin='lower', vmin=0, vmax=0.005)
plt.colorbar(img)
plt.xticks(arange(CC.shape[0]), names, rotation='vertical', fontsize=14)
plt.yticks(arange(CC.shape[1]), names, rotation='horizontal', fontsize=14)
plt.tight_layout()
fig.savefig('Hosts_covar.pdf')
plt.close(fig)

Rdata = ascii.read('Riess+2016tab5.dat')
sids = array([list(Rdata['Host']).index(name) for name in cephlist])
fig = plt.figure()
ax = fig.add_subplot(111)
ax.errorbar(Rdata['mu_ceph'][sids], DMs-Rdata['mu_ceph'][sids],
      fmt='o', capsize=0, xerr=Rdata['sigma2'][sids], yerr=eDMs)
ax.set_xlabel('DM(Riess)')
ax.set_ylabel('DM(MCMC)-DM(Riess)')
ax.axhline(0)
ax.set_ylim(-0.5,0.5)
fig.savefig('Delta-DMs.pdf')
plt.close(fig)

Adata = ascii.read('Riess+2016tab4.mrt.dat')
Adata_g = Adata.group_by('Field')
metals = Adata_g['[O/H]'].groups.aggregate(mean)
sids2 = array([list(Adata_g.groups.keys['Field']).index(name) for name in cephlist])
fig = plt.figure()
ax = fig.add_subplot(111)
ax.errorbar(metals[sids2], DMs-Rdata['mu_ceph'][sids],
      fmt='o', capsize=0, yerr=eDMs)
ax.set_xlabel('12+log(O/H)')
ax.set_ylabel('DM(MCMC)-DM(Riess)')
ax.axhline(0)
fig.savefig('Delta-DMs_metal.pdf')
plt.close(fig)

VIs = Adata_g['F555W-F814W'].groups.aggregate(mean)
fig = plt.figure()
ax = fig.add_subplot(111)
ax.errorbar(VIs[sids2], DMs-Rdata['mu_ceph'][sids],
      fmt='o', capsize=0, yerr=eDMs)
ax.set_xlabel('V-I')
ax.set_ylabel('DM(MCMC)-DM(Riess)')
ax.axhline(0)
fig.savefig('Delta-DMs_VI.pdf')
plt.close(fig)

Ps = Adata_g['Per'].groups.aggregate(mean)
fig = plt.figure()
ax = fig.add_subplot(111)
ax.errorbar(Ps[sids2], DMs-Rdata['mu_ceph'][sids],
      fmt='o', capsize=0, yerr=eDMs)
ax.set_xlabel('Period (days)')
ax.set_ylabel('DM(MCMC)-DM(Riess)')
ax.axhline(0)
fig.savefig('Delta-DMs_P.pdf')
plt.close(fig)

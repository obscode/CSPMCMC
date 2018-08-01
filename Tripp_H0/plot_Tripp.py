#!/usr/bin/env python
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt, rcParams
from matplotlib.ticker import MaxNLocator
import myplotlib
import get_data
import config
from numpy import *
import pickle
import sys,os,string
from scipy.optimize import leastsq
import sigfig
from astropy.io.ascii import read
from A_lamb import A_lamb
import STANstats


plt.style.use('serif')
rcParams['font.size'] = 14

cfgfile = sys.argv[1]
base = os.path.dirname(cfgfile)

cf = config.config(cfgfile)
data = get_data.get_data(cf)
pfile = cf.Sampler.outfile

ch = STANstats.STANchains(pfile)

names = data.names
Np = len(names)
d0 = 1.0

a = ch.get_trace('a', merge=True) - 19
b = ch.get_trace('b', merge=True)
c = ch.get_trace('c', merge=True)

evar = atleast_1d(ch.get_trace('evar', merge=True))
Rl = atleast_1d(ch.get_trace('Rl', merge=True))
dz = ch.median('vpec')*100/3e5

ma = atleast_1d(ch.median('a')) - 19
mb = atleast_1d(ch.median('b'))
mc = atleast_1d(ch.median('c'))
sa = atleast_1d(ch.std('a'))
sb = atleast_1d(ch.std('b'))
sc = atleast_1d(ch.std('c'))
mevar = atleast_1d(ch.median('evar'))
mRl = atleast_1d(ch.median('Rl'))
eRl = atleast_1d(ch.std('Rl'))
Rlambs = ch.get_trace('Rl', merge=True)
if cf.Model.Rv_global:
   R = ch.get_trace('R', merge=True)
H0 = ch.median('H0')
eH0 = ch.std('H0')
if 'alpha' in ch.params:
   malpha = ch.median('alpha')
   salpha = ch.std('alpha')
else:
   malpha = 0
   salpha = 0


#DMcorr = 5*log10(1.0*H0/72)

# OUTPUT a table of the results
f = open("fit_results.txt", 'w')
outcols = []
outcols2 = []
for i in range(len(data.filters)):
   the_a = a[:,i]
   the_b = b[:,i]
   the_c = c[:,i]
   the_R = Rlambs[:,i]
   the_sig = evar[:,i]
   covmat = cov([the_a,the_b,the_c,the_R])

   med_a = median(the_a)
   med_b = median(the_b)
   med_c = median(the_c)
   med_R = median(the_R)
   e_a = sqrt(covmat[0,0])
   e_b = sqrt(covmat[1,1])
   e_c = sqrt(covmat[2,2])
   e_R = sqrt(covmat[3,3])
   sig = sqrt(median(the_sig))

   l = ['%s(%s)' % (data.filters[i],"-".join(cf.Data.color))]
   l += list(sigfig.round_sig_error(med_a, e_a, 2))
   l += list(sigfig.round_sig_error(med_b, e_b, 2))
   l += list(sigfig.round_sig_error(med_c, e_c, 2))
   l += list(sigfig.round_sig_error(med_R, e_R, 2))
   l += [sigfig.round_sig(sig, 2)]
   outcols.append(l)

outcols2 = [[sigfig.round_sig(covmat[i,j],2) for i in range(3)] \
      for j in range(3)]

format = ""                  
for j in range(len(outcols[0])):
   maxlen = max([len(outcols[i][j]) for i in range(len(outcols))])
   format += " %%%ds" % maxlen

print >>f, "H0 = %.3f +\- %.3f" % (H0, eH0)
print >>f, format % ('#f','a','+/-','b','+/-','c','+/-','R','+/-','sig')
for outcol in outcols:
   print >>f, format % tuple(outcol)

print >>f, "#Covariance Matrix"
maxlen = 0
for i in range(3):
   for j in range(3):
      maxlen = max(maxlen, len(outcols2[i][j]))
format = ("%%%ds " % maxlen)*3
for outcol in outcols2:
   print >>f, format % tuple(outcol)

if cf.Model.Rv_global:
   print >>f, "# Reddening law"
   mR,eR = sigfig.round_sig_error(median(R), std(R), 2)
   print >>f, "Rv = %s +/- %s" % (mR,eR)

f.close()

fcombs = [(data.filters[i],cf.Data.color[0],cf.Data.color[1]) \
      for i in range(len(data.filters))]
print fcombs
if ('B','B','V') in fcombs and ('H','B','V') in fcombs:
   fig2 = myplotlib.PanelPlot(2,2, rect=[0,0,0.93,1])
else:
   fig2 = None

resids = []
if 'SN2006mr' in names:
   id06mr = names.index('SN2006mr')
else:
   id06mr = None

for i in range(len(data.filters)):
   f1 = data.filters[i]
   f2,f3 = cf.Data.color
   if (f1,f2,f3) == ('B','B','V') and fig2 is not None:
      ax21 = fig2.axes[2]
      ax22 = fig2.axes[3]
   elif (f1,f2,f3) == ('H','B','V') and fig2 is not None:
      ax21 = fig2.axes[0]
      ax22 = fig2.axes[1]
      ax22.xaxis.set_major_locator(MaxNLocator(5))
   else:
      ax21 = None
      ax22 = None
   cids = equal(data.findex, i+1)  # findex is from 1
   oids = data.sindex[cids]-1      # sindex is from 1
   if id06mr is not None:
      o6id = nonzero(equal(oids,id06mr))[0][0]

   m1 = data.ms[cids]
   vm1 = data.vms[cids]
   cs = data.cs[oids]
   vcs = data.vcs[oids]
   Rs = mRl[i]
   eRs = eRl[i]
   zcmb = data.zcmb[oids]
   zhel = data.zhel[oids]
   DMs = ch.median('DM')[oids]# + DMcorr
   if id06mr is not None:
      DMs[o6id] = 31.25
   # calibrating hosts
   cids = greater(array(data.host)[oids], 0)

   if cf.Model.HostMass:
      hmass = 0.4*(DMs - data.K[oids]) + 1.04
      ehmass = data.sigmaK[oids]

   zids = ~isnan(zcmb)
   dm15s = data.st[oids]
   edm15s = data.est[oids]
   #d0 = self.d0

   #fig = plt.figure()
   fig = myplotlib.PanelPlot(2,1,pwidths=(0.4,0.4))
   fig3 = plt.figure()
   fig.axes[1].xaxis.set_major_locator(MaxNLocator(5))
   #fig.right_pad = 0.1
   ax = fig.axes[0]
   ax.invert_yaxis()
   if ax21: ax21.invert_yaxis()
   fig.axes[1].set_xlabel('$m_%s - m_%s$' % (f2,f3))
   if ax22:  ax22.set_xlabel('$m_%s - m_%s$' % (f2,f3))
   fig.axes[0].set_xlabel('$s_{BV}$')
   if ax21: ax21.set_xlabel('$s_{BV}$')

   ax.set_ylabel('$m_%s - R_%s(m_%s - m_%s)$' % (f1,f1,f2, f3))
   if ax21:
      ax21.set_ylabel('$m_%s - R_%s(m_%s - m_%s)$' % (f1,f1,f2, f3))
   x1 = dm15s
   y1 = m1 - Rs*cs - DMs
   ey = sqrt(vm1 + Rs**2*vcs + power(2.17*zids*dz/zcmb,2))
   for axis in [ax, ax21]:
      if not axis: continue
      axis.errorbar(x1, y1, xerr=edm15s, yerr=ey, fmt='o', capsize=0,
            color='k', mfc='white', ecolor='k')
      if id06mr is not None:
         axis.plot([x1[o6id]],[y1[o6id]], '*', mec='k', mfc='yellow', 
               ms=15, zorder=1000)
      axis.plot(x1[cids], y1[cids], 'o', mfc='red', mec='none', zorder=1000)

   rs = y1 - ma[i] - mb[i]*(x1-d0) - mc[i]*(x1-d0)**2
   sids1 = argsort(absolute(rs)/ey)
   resids.append(rs)

   sig = sqrt(mevar[i])

   ax3 = fig3.add_subplot(111)
   ax3.plot(zcmb, rs, 'o', mec='k', mfc='white',color='k')
   ax3.plot(zcmb[cids], rs[cids], 'o', mec='none', mfc='red',zorder=50)
   xx = linspace(0, zcmb.max(), 100)
   yy = sqrt(sig**2 + power(2.17*dz/xx,2))
   ax3.axhline(0, color='0.5')
   ax3.plot(xx, yy, '-', color='red')
   ax3.plot(xx, -yy, '-', color='red')
   ax3.set_xlabel('$z_{CMB}$')
   ax3.set_ylabel('residuals (mag)')
   ax3.invert_yaxis()
   ax3.set_ylim(-1.5,1.5)
   plt.tight_layout()
   fig3.savefig('residuals_%s_z.eps' % (f1+f2+f3))

   if cf.Model.HostMass:
      fig4 = plt.figure()
      ax4 = fig4.add_subplot(111)
      ax4.errorbar(hmass, rs, fmt='o', xerr=ehmass, yerr=sig, capsize=0)
      res = ax4.scatter(hmass, rs, marker='o', c=zcmb, zorder=20)
      ax4.plot(hmass[cids], rs[cids], 'o', mec='k', mfc='red', zorder=50)
      xx = linspace(hmass.min(), hmass.max(), 100)
      ax4.plot(xx, (xx-data['M0'])*malpha, '-', color='red')
      ax4.plot(xx, (xx-data['M0'])*(malpha+salpha), '--', color='red')
      ax4.plot(xx, (xx-data['M0'])*(malpha-salpha), '--', color='red')
      plt.colorbar(res)
      ax4.set_xlabel('$\log_{10}(Mass/M_\odot)$')
      ax4.set_ylabel('residuals (mag)')
      ax4.set_ylim(-1.5,1.5)
      plt.tight_layout()
      fig4.savefig('residuals_%s_mass.eps' % (f1+f2+f3))


   xx = arange(101)/100.*(data.st.max() - data.st.min()) + data.st.min()
   yy = ma[i] + mb[i]*(xx-d0) + mc[i]*(xx-d0)**2
   print "ma=%f, mb=%f, mc=%f" % (ma[i], mb[i], mc[i])
   ax.plot(xx,yy,color='red', zorder=50)
   if ax21:  ax21.plot(xx, yy, color='red', zorder=50)
   x0 = 0.95
   ha = 'right'
   ax.text(x0, 0.05, 
      "$P^0 = %.2f \\pm %.2f$\n$P^1 = %.2f \\pm %.2f$\n$P^2 = %.2f \\pm %.2f$\n$\\sigma=%.2f$" %\
      (ma[i],sa[i],mb[i],sb[i],mc[i],sc[i],sig), transform=ax.transAxes, 
      va='bottom', ha=ha, fontdict={'size':12})
   if ax21:
      ax21.text(x0, 0.05, 
         "$P^0 = %.2f \\pm %.2f$\n$P^1 = %.2f \\pm %.2f$\n$P^2 = %.2f \\pm %.2f$\n$\\sigma=%.2f$" %\
         (ma[i],sa[i],mb[i],sb[i],mc[i],sc[i],sig), transform=ax21.transAxes, 
         va='bottom', ha=ha, fontdict={'size':12})
   
   # Now let's make some realizations of the function
   N = len(a[:,i])
   n = N/100
   if n < 1:  n = 1
   funcs = a[:,i,newaxis] + b[:,i,newaxis]*(xx[newaxis,:]-d0) +\
           c[:,i,newaxis]*(xx[newaxis,:]-d0)**2
   std_f = std(funcs, axis=0)
   ax.plot(xx,yy+std_f,'--', color='red', zorder=50)
   ax.plot(xx,yy-std_f,'--', color='red', zorder=50)
   ax.fill_between(xx, yy-sig,yy+sig, facecolor="0.7", edgecolor="0.7",
         zorder=0)
   if ax21:
      ax21.plot(xx,yy+std_f,'--', color='red', zorder=50)
      ax21.plot(xx,yy-std_f,'--', color='red', zorder=50)
      ax21.fill_between(xx, yy-sig,yy+sig, facecolor="0.7", edgecolor="0.7",
            zorder=0)
   plt.draw()

   # Second panel
   ax = fig.axes[1]
   ax.text(1.1, 0.5, '$m_%s - P(s-1)$' % (f1),
         ha='left', va='center', rotation=90, transform=ax.transAxes)
   if ax22:
      ax22.text(1.1, 0.5, '$m_%s - P(s-1)$' % (f1),
            ha='left', va='center', rotation=90, transform=ax22.transAxes)
   x2 = cs
   ex = sqrt(vcs)
   y2 = m1 - mb[i]*(dm15s-d0) -mc[i]*(dm15s-d0)**2 -DMs
   ey = sqrt(vm1 + mb[i]**2*edm15s**2 + 2*mc[i]*dm15s*edm15s**2 +
            power(2.17*zids*dz/zcmb,2))
   sids2 = argsort(absolute(y2 - ma[i] - Rs*cs)/ey)
   for axis in [ax,ax22]:
      if not axis: continue
      axis.errorbar(x2, y2, xerr=ex, yerr=ey, fmt='o', capsize=0,
            color='k', mfc='white', ecolor='k')
      axis.plot(x2[cids],y2[cids], 'o', mec='k', mfc='red', zorder=1000)
   sig = sqrt(mevar[i])
   xx = linspace(cs.min(), cs.max(), 100)
   yy = ma[i] + Rs*xx # + mc[i]*(x-d0)**2
   ax.plot(xx,yy,color='red', zorder=50)
   if ax22:  ax22.plot(xx,yy,color='red', zorder=50)
   x0 = 0.05
   ha = 'left'
   
   # Now let's make some realizations of the function
   N = len(a[:,i])
   n = N/100
   if n < 1:  n = 1
   funcs = a[:,i,newaxis] + Rlambs[:,i,newaxis]*xx[newaxis,:]
   std_f = std(funcs, axis=0)
   ax.plot(xx,yy+std_f,'--', color='red', zorder=50)
   ax.plot(xx,yy-std_f,'--', color='red', zorder=50)
   ax.fill_between(xx, yy-sig,yy+sig, facecolor="0.7", edgecolor="0.7",
         zorder=0)
   if ax22:
      ax22.plot(xx,yy+std_f,'--', color='red', zorder=50)
      ax22.plot(xx,yy-std_f,'--', color='red', zorder=50)
      ax22.fill_between(xx, yy-sig,yy+sig, facecolor="0.7", edgecolor="0.7",
            zorder=0)
   #plt.draw()
   (ymin,ymax) = ax.get_ylim()
   #ax.set_ylim(ymax,ymin)
   ax.text(0.05, 0.05, 
      "$R_%s = %.2f \\pm %.2f$" %\
      (f1,Rs,eRs), transform=ax.transAxes, 
      va='bottom', ha='left', fontdict={'size':12})
   if ax22:
      ax22.text(0.05, 0.05, 
         "$R_%s = %.2f \\pm %.2f$" %\
         (f1,Rs,eRs), transform=ax22.transAxes, 
         va='bottom', ha='left', fontdict={'size':12})


   fig.set_limits()
   fig.draw()
   fig.fig.savefig("%s_fit.eps" % (f1+f2+f3))

   names = [data.names[k] for k in oids]
   # Now label 10 the most outlying objects
   for id in sids1[-10:]:
      fig.axes[0].text(x1[id], y1[id], names[id]+" >", va='top', ha='right',
            fontdict={'size':8})
   for id in sids2[-10:]:
      fig.axes[1].text(x2[id], y2[id], names[id]+" >", va='top', ha='right',
            fontdict={'size':8})
   fig.draw()
   fig.fig.savefig("%s_fit_ann.eps" % (f1+f2+f3))
   #plt.draw()
   #fig.savefig(ofile+"_%s_fit_ann.eps" % (f1+f2+f3))

   fig.close()
if fig2 is not None:
   fig2.set_limits(pad=0.04)
   ymax = max([ax.get_ylim()[0] for ax in fig2.axes])
   ymin = min([ax.get_ylim()[1] for ax in fig2.axes])
   [ax.set_ylim(ymax,ymin) for ax in fig2.axes]
   fig2.axes[2].set_xlabel('')
   fig2.axes[3].set_xlabel('')
   fig2.draw()
   fig2.fig.savefig(os.path.join(base,'BBV_HBV_compare.eps'))
   fig2.close()

f = open("residuals.dat", 'w')
mat = zeros((len(data.names),len(resids)), dtype=float32) + 99
for i in range(len(resids)):
   ids = data.sindex[data.findex==i+1]-1   # (f,s)index indexed from 1
   mat[ids,i] = resids[i]
fmt1 = "%-10s " + "%7s "*len(resids)
fmt2 = "%-10s " + "%7.3f "*len(resids)
strs = [data.filters[i]+cf.Data.color[0]+cf.Data.color[1] \
      for i in range(len(resids))]
print >> f, fmt1 % tuple(["#Name"] + strs)
for i in range(len(data.names)):
   vals = [data.names[i]]
   vals += mat[i,:].tolist()
   print >> f, fmt2 % tuple(vals)
f.close()

if cf.Model.Rv_global:
   f = plt.figure()
   ax = f.add_subplot(111)
   ax.set_xlabel('$R_V$')
   ax.set_ylabel('$N$')
   ax.hist(R, bins=30)
   mR = median(R)
   eR = std(R)
   ax.text(0.95, 0.95, "$R_V = %.2f \pm %.2f$" % (mR, eR),
      transform=ax.transAxes, ha='right', va='top',
      fontdict={'size':10})
   plt.draw()
   f.savefig("Rv_hist.eps")
   plt.close(f)

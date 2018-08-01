#!/usr/bin/env python
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib import rcParams
from numpy import *
import pickle
import sys,os,string
import sigfig
import EBV
import STANstats
from A_lamb import A_lamb
#import bspline
import bspline2

d0 = 1

rcParams['font.size'] = 16
rcParams['font.family'] = 'serif'

def pick_handler(event):
   me = event.mouseevent
   art = event.artist
   xdata = art.get_xdata()
   ydata = art.get_ydata()
   ind = event.ind
   print "SN %s at %.3f,%.3f" % (names[ind], xdata[ind], ydata[ind])

def get_color(s, f1, f2):
   # Damn this is getting complicated!
   f1ids = equal(s.fids, s.filters.index(f1))
   f2ids = equal(s.fids, s.filters.index(f2))
   d1ids = s.oids[f1ids]
   d2ids = s.oids[f2ids]
   dids = list(set(d1ids.tolist())&set(d2ids.tolist()))
   g1ids = array([id in dids for id in d1ids])
   g2ids = array([id in dids for id in d2ids])
   return s.ms[f1ids][g1ids] - s.ms[f2ids][g2ids],\
          s.vms[f1ids][g1ids] + s.vms[f2ids][g2ids],\
          d1ids[g1ids]

argv = sys.argv[1:]

cfg = EBV.get_config(argv[0])
data = EBV.get_data(cfg)
vinfo = EBV.setup_varinfo(cfg, data)

base = os.path.dirname(argv[0])
ofile = os.path.join(base, "results")

R_global = cfg.Priors.Rv_global
blue = cfg.Priors.blue
prior = cfg.Priors.red_prior

trs = STANstats.STANchains(filename=cfg.Sampler.outfile)
with open(cfg.Sampler.outfile) as fin:
   d = pickle.load(fin)
   data = d['data']

names = data.names
Np = len(names)

# basis coefficients
a = trs.get_trace('a', merge=True)
ma = trs.median('a')
sa = trs.std('a')

taus = trs.get_trace('tau', merge=True)
reds = trs.get_trace('EBV', merge=True)
Rv = trs.get_trace('R_V', merge=True)

#Indexed by [fid,iter,oid]
redlaw = cfg.Model.redlaw
A_lambdas = array([A_lamb(f,reds,Rv,redlaw) for f in data['filters']])

evar = trs.median('evar')
if len(shape(evar)) == 0:
   evar = ones((len(colors),))*evar

# OUTPUT a table of the results
f = open(ofile + ".txt", 'w')
wanted_colors = [('u','B'),('B','V'),('g','r'),('V','r'),('V','i'),('r','i'),
                 ('V','Y'),('Y','J'),('J','H'),('V','J'),('V','H')]
f.write('# Bspline knot points:\n')
f.write('[%.3f' % data.knots[0])
[f.write(', %.3f' % data.knots[i]) for i in range(1,len(data.knots)-1)]
f.write(', %.3f]\n' % data.knots[-1])
outcols = []
for cols in wanted_colors:
   if cols[0] not in data.filters or cols[1] not in data.filters:
      continue
   f1 = data.filters.index(cols[0])-1
   f2 = data.filters.index(cols[1])-1
   if f1 == -1:
      # B - something
      the_a = a[:,f2,:]
   elif f2 == -1:
      # something - B
      the_a = -a[:,f1,:]
   else:
      the_a = a[:,f2,:] - a[:,f1,:]

   med_a = median(the_a, axis=0)
   covmat = cov(the_a.T)

   e_a = sqrt(diag(covmat))

   varc = (evar[f1] + evar[f2])/2
   this_col = ['%s-%s' % (cols[0],cols[1])]
   for i in range(med_a.shape[0]):
      this_col += list(sigfig.round_sig_error(med_a[i], e_a[i], 2))
   this_col += [sigfig.round_sig(sqrt(varc),2)]

   outcols.append(this_col)

format = ""                  
for j in range(len(outcols[0])):
   maxlen = max([len(outcols[i][j]) for i in range(len(outcols))])
   format += " %%%ds" % maxlen

labs = ['#clr']
for i in range(med_a.shape[0]):
   labs += ['a[%d]' % i,'+/-']
labs += ['sig']
labs = tuple(labs)

print >>f, format % labs
for outcol in outcols:
   print >>f, format % tuple(outcol)

# Now for the scalars
med_tau = median(taus)
std_tau = std(taus)

print >> f, "tau  = %s +/- %s" % sigfig.round_sig_error(med_tau,std_tau,2)

for i in range(len(data.names)):
   C = cov(reds[:,i],Rv[:,i])
   print >> f, "%10s  %5.3f %5.3f %4.2f %4.2f %f" % \
           (data.names[i],median(reds[:,i]),sqrt(C[0,0]),median(Rv[:,i]),
            sqrt(C[1,1]),C[0,1]) 
f.close()

#bs = bspline.Bspline(knot_vector=data.knots, order=3)
x = arange(101)/100.*(data.st.max() - data.st.min()) + data.st.min()
#bx = array([bs(xi) for xi in x])
bx = bspline2.bspline_basis(data.knots, x, 3, gradient=cfg.Model.gradient)

colors = [(data.filters[0],filt) for filt in data.filters[1:]]
for i,(f1,f2) in enumerate(colors):
   fig = plt.figure()
   ax = fig.add_subplot(111)
   ax.set_xlabel('$s_{BV}$')
   ax.set_ylabel('$%s - %s$' % (f1, f2))
   cs,vcs,dids = get_color(data, f1,f2)
   ax.errorbar(data.st[dids], cs, xerr=data.est[dids], 
         yerr=sqrt(vcs), fmt='o', capsize=0, color='k')
   rids = greater(data.EBVgal[dids], 0.1)
   ax.plot(data.st[dids][rids], cs[rids], 'o', color='red')
   sig = sqrt(evar[i])
   #y = ma[i] + mb[i]*(x-d0) + mc[i]*power(x-d0,2)
   y = dot(bx, ma[i])
   ax.plot(x,y,color='red')
   
   # Now let's make some realizations of the function
   funcs = array([dot(bx,a[j,i]) for j in range(a.shape[0])])
   std_f = std(funcs, axis=0)
   ax.plot(x,y+std_f,'--', color='red')
   ax.plot(x,y-std_f,'--', color='red')
   ax.fill_between(x, y-sig,y+sig, facecolor="0.7", edgecolor="0.7",
         zorder=0)
   plt.draw()
   fig.savefig(ofile+"_%s%s_fit.eps" % (f1,f2))

   for j in range(len(dids)):
      name = data.names[dids[j]]
      xst = data.st[dids[j]]
      yst = cs[j]
      ax.text(xst,yst, name+"_", ha='right',va='top', fontdict={'size':8})

   fig.savefig(ofile+"_%s%s_label.eps" % (f1,f2))
   plt.close(fig)

   # now do the correction
   fig = plt.figure()
   ax = fig.add_subplot(111)
   ax.set_xlabel('$s_{BV}$')
   ax.set_ylabel('$%s - %s$ (corrected)' % (f1, f2))
   gids = greater(data.fids, 0)
   cids = data.fids - 1
   id1 = data.filters.index(f1)
   id2 = data.filters.index(f2)
   cmod = median(A_lambdas[id1,:,dids]-A_lambdas[id2,:,dids], axis=1)
   ax.errorbar(data.st[dids], cs - cmod, xerr=data.est[dids], 
         yerr=sqrt(vcs), fmt='o', capsize=0, color='k', mfc='white')
   ax.plot(data.st[dids][rids], cs[rids]-cmod[rids], 'o', color='red')
   sig = sqrt(evar[i])
   ax.plot(x,y,color='red')
   
   # Output the residuals
   fout = open(ofile+"_%s%s_resids.dat" % (f1,f2), 'w')
   print >> fout, "#%9s %5s %5s %6s %5s %5s" % \
         ("SN","st","est","res","err","#sig")
   resids = cs - cmod - dot(data.Bs[dids,:],ma[i])
   nsigmas = resids/sqrt(vcs + evar[i])
   sids = argsort(absolute(nsigmas))
   for k in sids[::-1]:
      print >> fout, "%-10s %5.3f %5.3f %6.3f %5.3f %5.1f" % \
            (data.names[dids[k]],data.st[dids][k], data.est[dids][k],resids[k],
             sqrt(vcs[k]),nsigmas[k])
   fout.close()


   # Now let's make some realizations of the function
   ax.plot(x,y+std_f,'--', color='red')
   ax.plot(x,y-std_f,'--', color='red')
   ax.fill_between(x, y-sig,y+sig, facecolor="0.7", edgecolor="0.7",
         zorder=0)
   plt.draw()
   fig.savefig(ofile+"_%s%s_corr.eps" % (f1,f2))

   for j in range(len(dids)):
      name = data.names[dids[j]]
      xst = data.st[dids[j]]
      yst = cs[j] - cmod[j]
      ax.text(xst,yst, name+"_", ha='right',va='top', fontdict={'size':8})
   fig.savefig(ofile+"_%s%s_corr_label.eps" % (f1,f2))
   plt.close(fig)

   fig = plt.figure()
   ax = fig.add_subplot(111)
   ax.set_xlabel('$E(B-V)_{MW}$')
   ax.set_ylabel('residuals')
   cmod += median(dot(data.Bs[dids,:],a[:,i,:].T), axis=1) 
   ax.errorbar(data.EBVgal[dids], cs - cmod, xerr=data.eEBVgal[dids], 
         yerr=sqrt(vcs), fmt='o', capsize=0, color='k', mfc='white')
   ax.axhline(0, linestyle='-', color='red')
   ax.axhline(sig, linestyle='--', color='red')
   ax.axhline(-sig, linestyle='--', color='red')
   fig.savefig(ofile+"_%s%s_resids_EBVgal.eps" % (f1,f2))
   plt.close(fig) 

   fig = plt.figure()
   ax = fig.add_subplot(111)
   ax.set_xlabel('$s_{BV}$')
   ax.set_ylabel('residuals')
   ax.errorbar(data.st[dids], cs - cmod, xerr=data.est[dids], 
         yerr=sqrt(vcs), fmt='o', capsize=0, color='k', mfc='white')
   ax.axhline(0, linestyle='-', color='red')
   ax.axhline(sig, linestyle='--', color='red')
   ax.axhline(-sig, linestyle='--', color='red')
   ax.set_ylim(-1,1)
   fig.savefig(ofile+"_%s%s_resids.eps" % (f1,f2))
   plt.close(fig) 

   fig = plt.figure()
   ax = fig.add_subplot(111)
   ax.set_xlabel('$E(%s - %s)$' % (f1, f2))
   ax.set_ylabel('Number')
   EBVs = cs - dot(data.Bs[dids,:],ma[i])
   freq,bins,patches = ax.hist(EBVs, bins=20)
   if f1 == 'B' and f2 == 'V':
      xx = arange(101)/100.*bins.max()+0.01
      dx = bins[1] - bins[0]
      if prior == 'Cauchy':
         bet = median(taus)
         ebet = std(taus)
         #fac = sum(freq*dx)/bet/arctan(10./bet) 
         fac = (sum(freq)*dx)/bet/arctan(10./bet) 
         yy = fac*power(1+power(xx/bet,2),-1)
         lab = "$\\beta = %.3f \pm %.3f$" % (bet,ebet)
      elif prior == 'Exp':
         tau = median(taus)
         etau = std(taus)
         fac = sum(freq)*dx/tau
         yy = fac*exp(-xx/tau)
         lab = "$\\tau = %.3f \pm %.3f$" % (tau,etau)
      elif prior == 'disk':
         tau = median(taus)
         etau = std(taus)
         tmin = median(tmins)
         fact = sum(freq)*dx*(pi/2-tmin)*arccos(sin(tmin))
         u = (xx/tau + 1)
         yy = fact*power(tau*(pi/2-tmin)*u*sqrt(u**2 - 1),-1)
         lab = "$\\tau = %.3f \pm %.3f$" % (tau,etau)
      else:
         lab = ""
         yy = None
  
      if yy is not None:
         ax.plot(xx,yy,'-',color='red')
         ax.text(0.95,0.95,lab, ha='right', va='top', transform=ax.transAxes)
   ax.set_ylim(0,freq.max()*1.05)
   plt.draw()
   fig.savefig(ofile+"_%s%s_hist.eps" % (f1,f2))
   plt.close(fig)

if 'B' in data.filters and 'V' in data.filters:
   cs,vcs,dids = get_color(data, 'B','V')
   f = plt.figure()
   ax = f.add_subplot(111)
   ax.set_xlabel('$E(B-V)$ (parameter)')
   ax.set_ylabel('$E(B-V)$ (observed)')
   x = median(reds, axis=0)
   dx = std(reds, axis=0)
   f2 = data.filters.index('V')
   y = cs - dot(data.Bs[dids,:], ma[f2-1])
   dy = sqrt(vcs)
   ax.errorbar(x[dids], y, xerr=dx[dids], yerr=dy, fmt='o', color='k',
         mfc='white')
   xx = array([x.min(), x.max()])
   ax.plot(xx, xx, color='red')
   plt.draw()
   f.savefig(ofile+"_EBV_EBV.eps")
   plt.close(f)

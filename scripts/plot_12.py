#!/usr/bin/env python

import numpy as np
from sklearn.metrics import mean_squared_error
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter)


###
fig, axs = plt.subplots(6,2,figsize=(6.4,12))
plt.rcParams['axes.titlepad'] = 0 
params = {'mathtext.default': 'regular' }          
plt.rcParams.update(params)
min = 3
max = 11 
for i, ax in enumerate(axs.flat):
   #if i > 10: break
   ax.tick_params(size=6,width=1)
   ax.tick_params(which='minor',width=1, size=4)
   ax.tick_params(labelsize=9)
   ax.set_xlim(min,max)
   ax.set_ylim(min,max)
   ax.xaxis.set_major_locator(MultipleLocator(2))
   ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
   ax.xaxis.set_minor_locator(MultipleLocator(1))
   ax.yaxis.set_major_locator(MultipleLocator(2))
   ax.yaxis.set_major_formatter(FormatStrFormatter('%d'))
   ax.yaxis.set_minor_locator(MultipleLocator(1))

   lims = np.array([min,max])
   ax.plot(lims, lims, 'k--', lw=1,alpha=0.75, zorder=10)
   ax.fill_between(lims, lims+0.5,lims-0.5, facecolor='grey',alpha=0.5, zorder=10)
   ax.fill_between(lims, lims+1,lims+0.5, facecolor='grey',alpha=0.25, zorder=10)
   ax.fill_between(lims, lims-1,lims-0.5, facecolor='grey',alpha=0.25, zorder=10)
###


path = {
       'cdk5/r_predt.smi':'CDK5',  
       'jak1/r_predt.smi':'JAK1',
       'mapk14/r_predt.smi':'MAPK14',
       'pik3ca/r_predt.smi':'PIK3CA',
       'tyk2/r_predt.smi':'TYK2',
       'f7/r_predt.smi':'F7',
       'prss2/r_predt.smi':'PRSS2',
       'mmp2/r_predt.smi':'MMP2',
       'mgll/r_predt.smi':'MGLL',
       'il4/r_predt.smi':'IL4',
       'gnrhr/r_predt.smi':'GNRHR',
       'drd2/r_predt.smi':'DRD2',
   }

for i,para in enumerate(path):
    exp=[]
    prd=[]
    std=[]
    with open(para) as f:
         for line in f:
             cells = line.split()
             exp.append(float(cells[1]))
             prd.append(float(cells[2]))
             std.append(float(cells[3]))

    slope, intercept, r, p, std_err = stats.linregress(exp, prd)
    print(path[para],'RMSE %.2f'%mean_squared_error(exp,prd),'R2 %.2f'%(r*r))

    axs.flat[i].scatter(exp, prd,s=20 , alpha=0.80, zorder=20, label=path[para]) #+' RMSE: %.2f'%rmse)
    axs.flat[i].errorbar(exp, prd, std, lw=3,linestyle='', alpha=0.80, zorder=20, capsize=1,elinewidth=1)
    axs.flat[i].legend(frameon=False,loc='upper left',fontsize=9)


for i in range(2):
    axs.flat[10+i].set_xlabel(r'Expt. pK$_{\rm i}$', fontsize=9)
for i in range(6):
    axs.flat[i*2].set_ylabel(r'Pred. pK$_{\rm i}$', fontsize=9)

#plt.delaxes(axs[5,1])

plt.savefig('si.jpg',bbox_inches='tight',dpi=300) #,transparent=True)
plt.show()


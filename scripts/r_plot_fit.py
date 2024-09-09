#!/usr/bin/env python

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)


fig, ax = plt.subplots()

exp = []
prd = []

with open('r_train.smi') as f:
     for line in f:
         cells = line.strip().split()
         exp.append(float(cells[1]))
         ps = []
         for p in cells[2:]:
             ps.append(float(p))
         prd.append(np.mean(ps))

x = np.array(exp)
y = np.array(prd)
xy = np.vstack([x,y])
z = gaussian_kde(xy)(xy)
idx = z.argsort()
x, y, z = x[idx], y[idx], z[idx]
ax.scatter(x, y, c=z, s=50,alpha=0.5,cmap='jet', edgecolor='none')

ax.xaxis.set_major_locator(MultipleLocator(1))
ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax.xaxis.set_minor_locator(MultipleLocator(.5))
ax.yaxis.set_major_locator(MultipleLocator(1))
ax.yaxis.set_major_formatter(FormatStrFormatter('%d'))
ax.yaxis.set_minor_locator(MultipleLocator(.5))
plt.tick_params(size=7,width=1)
plt.tick_params(which='minor',width=1, size=4)
plt.tick_params(labelsize=15)
plt.xlabel('Exp. pIC$_{50}$', fontsize=18)
plt.ylabel('Pred. pIC$_{50}$', fontsize=18)

ax.set_xlim(4,12)
ax.set_ylim(4,12)

lims = [
    np.min([plt.xlim(), plt.ylim()]),  # min of both axes
    np.max([plt.xlim(), plt.ylim()]),  # max of both axes
]

# now plot both limits against eachother
plt.plot(lims, lims, 'k--', alpha=0.75, zorder=0)

plt.savefig('fit_mean.png',bbox_inches='tight',dpi=300,transparent=True)
plt.show()


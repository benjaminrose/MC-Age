"""plot_PCASignificance.py -- plot the PCA significance

* Benjamin Rose
* University of Notre Dame
* benjamin.rose@me.com
* Python 3.6
* 2018-05-30
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.mlab as mlab

RHO = np.genfromtxt('../resources/PCASignificance.csv', delimiter=',')
MEASURED_S = 0.64
X_LIMIT = 0.71

sns.set(context='talk', style='ticks', font='serif', color_codes=True)

#set axes ticks and gridlines
ax = plt.gca()
ax.tick_params(axis='both', top='on', right='on', direction='in')
ax.tick_params(axis='x', pad=10)
ax.grid(which='major', axis='both', color='0.90', linestyle='-')
ax.set_axisbelow(True)
plt.xlim(0, X_LIMIT)

#plot distribution and measurement line
sns.distplot(RHO, kde=False, hist_kws={"linewidth": 0})
# sns.kdeplot(RHO, cumulative=True)
plt.axvline(x=MEASURED_S, color='r', lw=1.5)

#plot a guassian
x = np.linspace(0, X_LIMIT, 100)
plt.plot(x, 2200*mlab.normpdf(x, 0, 0.1), linestyle='--', c='k', alpha=0.65)

plt.xlabel("Spearman's Correlation Coefficient")
plt.ylabel('Occurrences')

fig = plt.gcf()
fig.set_tight_layout({'pad': 1.5})
plt.savefig('bootstrap.pdf', bbox_inches='tight')
plt.show()
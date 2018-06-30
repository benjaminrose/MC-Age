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

RHO = np.genfromtxt('../resources/PCASignificance.csv', delimiter=',')

sns.set(context='talk', style='ticks', font='serif', color_codes=True)

#set axes ticks and gridlines
ax = plt.gca()
ax.tick_params(axis='both', top='on', right='on', direction='in')
ax.grid(which='major', axis='both', color='0.90', linestyle='-')
ax.set_axisbelow(True)
plt.xlim(0, 0.45)

sns.distplot(RHO, kde=False)
# sns.kdeplot(RHO, cumulative=True)

plt.xlabel("Spearman's Correlation Coefficient")
plt.ylabel('Occurrences')

plt.savefig('bootstrap.pdf', bbox_inches='tight')
plt.show()
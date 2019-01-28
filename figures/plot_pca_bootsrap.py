"""plot_pca_bootstrap.py
"""
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# sns.set(context='talk', style='ticks', palette='colorblind', font='serif', color_codes=True)
sns.set(context='talk', style='ticks', font='serif', color_codes=True)
bootresult = np.genfromtxt('../resources/pca_bootstrap.csv', delimiter=',')

plt.figure('PCA bootstrap', figsize=(7.5, 5))
plt.xticks(np.arange(-1, 1.25, step=0.25))
plt.xlim([-0.9, 0.9])
sns.distplot(bootresult[:, 3], kde=False, hist_kws={'linewidth': 0}, label=r'age')
sns.distplot(bootresult[:, 2], kde=False, hist_kws={'linewidth': 0}, label=r'mass')
sns.distplot(bootresult[:, 1], kde=False, hist_kws={'linewidth': 0}, label=r'$c$')
sns.distplot(bootresult[:, 0], kde=False, hist_kws={'linewidth': 0}, label=r'$x_1$')
plt.legend(loc=9, frameon=False, ncol=2)
sns.despine(left=True)  # removes box, but not labels & ticks
plt.gca().get_yaxis().set_visible(False)
plt.savefig('PCA_bootstrap.pdf')
plt.show()




# f, axes = plt.subplots(2, 2, sharex=True)
# sns.distplot(bootresult[:, 0], kde=False, ax=axes[0, 0], label=r'$x_1$')
# sns.distplot(bootresult[:, 1], kde=False, ax=axes[0, 1], label=r'$c$')
# sns.distplot(bootresult[:, 2], kde=False, ax=axes[1, 0], label=r'mass')
# sns.distplot(bootresult[:, 3], kde=False, ax=axes[1, 1], label=r'age')
# axes[0, 0].set_title(r'$x_1$')
# axes[0, 1].set_title(r'$c$')
# axes[1, 0].set_title(r'mass')
# axes[1, 1].set_title(r'age')
# sns.despine(left=True)  # removes box, but not labels & ticks
# axes[0, 0].get_yaxis().set_visible(False)
# axes[0, 1].get_yaxis().set_visible(False)
# axes[1, 0].get_yaxis().set_visible(False)
# axes[1, 1].get_yaxis().set_visible(False)
# plt.show()
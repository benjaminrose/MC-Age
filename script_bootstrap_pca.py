""" script_bootstrap_pca.py - bootstrap uncertainties on PCA coefficients

Benjamin Rose, University of Notre Dame & STScI, benjamin.rose@me.com, 2019-01-16
"""

import numpy as np
from astropy.table import Table
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from astropy.stats import bootstrap

import matplotlib.pyplot as plt
import seaborn as sns

# TODO move to util.py
################
# copied from PCASignificance.py
################
Z_CUT = 0.2
HR_CUT = 0.7

# import data for PCA analysis
## HR
HR = pd.read_csv('data/Campbell_local.tsv', sep='\\t',
                 usecols=['SNID', 'redshift', 'hr', 'err_mu'], index_col='SNID')
HR.rename(columns={'err_mu': 'hr uncert'}, inplace=True)
HR = HR[HR['redshift'] < Z_CUT]
HR = HR[HR['hr'] < HR_CUT]
# print('Hubble Residual:')
# print(HR.describe())

## SALT2 parameters (x_1 & c)
t = Table.read('data/SDSS_Photometric_SNe_Ia.fits')
salt = t['CID', 'Z', 'X1', 'X1_ERR', 'COLOR', 'COLOR_ERR'].to_pandas()
salt.columns = salt.columns.str.lower()
salt.rename(columns={'cid': 'SNID', 'z': 'redshift'}, inplace=True)
salt.set_index('SNID', inplace=True)
# print('\nSALT2 parameters:')
# print(salt.describe())

## stellar mass
galaxy = pd.read_csv('resources/kcorrect_stellarmass.csv',
                     usecols=['GAL', 'redshift', 'stellarmass'], index_col='GAL')
galaxy.rename(columns={'redshift': 'gal redshift', 'stellarmass': 'stellar mass'}, inplace=True)
# print('\nGalaxy Stellar Mass:')
# print(galaxy.describe())

## age
age = pd.read_csv('resources/ages_campbell.tsv', sep='\\t', skiprows=[1],
                  usecols=['# sn id', 'age'], dtype={'age': np.float64, '# sn id': np.int})
age.rename(columns={'# sn id': 'SNID'}, inplace=True)
age.set_index('SNID', inplace=True)
# print('\nLocal Envir. Ages:')
# print(age.describe())


# combine into on array
data = pd.concat([HR, salt, galaxy, age], axis=1)
data.dropna(inplace=True)
data['stellar mass'] = np.log10(data['stellar mass'])
# print('\nAnalysis Data:')
# print(data.describe())

################
#end copied from PCASignificance.py
################

# PCA preprocessing
pca = PCA(n_components=4)
FEATURES = ['x1', 'color', 'stellar mass', 'age']
#x = data.loc[:, FEATURES].values
scaled_data = StandardScaler().fit_transform(data.loc[:, FEATURES].values)

# run PCA
#principal_components = pca.fit_transform(x)
pca.fit(scaled_data)
original_components = pca.components_
print(original_components)

# Perform bootstrap
def bootfunc(x):
    """Stats to be performed on each bootstrap re-sample.

    This function performs PCA, gets PC1, then converts to same
    handedness as on the original data set.
    """
    pc1 = pca.fit(x).components_[0]
    if original_components[0].dot(pc1) < 0:
        pc1 = -pc1
    return pc1
# bootresult = bootstrap(scaled_data, 3)
# print('test:', bootresult.shape, scaled_data.shape,
#     bootresult[0].shape, bootfunc(scaled_data).shape)
bootresult = bootstrap(scaled_data, 100000, bootfunc=bootfunc)
print(bootresult[:5])
# convert so all bootstrap results are using the same handed coordinate system
# if pca.components_[0].dot()

np.savetxt('resources/pca_bootstrap.csv', bootresult, delimiter=',')

print('bootstrap results')
print(bootresult[:5])
print('mean: ', np.mean(bootresult, axis=0))
# print('std: ', np.std(bootresult[:100], axis=0))
# print('std: ', np.std(bootresult[:500], axis=0))
print('std: ', np.std(bootresult, axis=0))
print('min: ', np.min(bootresult, axis=0))
print('max: ', np.max(bootresult, axis=0))

f, axes = plt.subplots(2, 2, sharex=True)
sns.distplot(bootresult[:, 0], kde=False, ax=axes[0, 0], label=r'$x_1$')
sns.distplot(bootresult[:, 1], kde=False, ax=axes[0, 1], label=r'$c$')
sns.distplot(bootresult[:, 2], kde=False, ax=axes[1, 0], label=r'mass')
sns.distplot(bootresult[:, 3], kde=False, ax=axes[1, 1], label=r'age')
axes[0, 0].set_title(r'$x_1$')
axes[0, 1].set_title(r'$c$')
axes[1, 0].set_title(r'mass')
axes[1, 1].set_title(r'age')
sns.despine(left=True)  # removes box, but not labels & ticks
axes[0, 0].get_yaxis().set_visible(False)
axes[0, 1].get_yaxis().set_visible(False)
axes[1, 0].get_yaxis().set_visible(False)
axes[1, 1].get_yaxis().set_visible(False)
plt.show()

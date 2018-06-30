"""PCASignificance.py -- testing the significance of our PCA results.

* Benjamin Rose
* University of Notre Dame
* benjamin.rose@me.com
* Python 3.6
* 2018-05-21
"""
import numpy as np
from astropy.table import Table
import pandas as pd
from scipy.stats import spearmanr
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA


SPEARMAN_CUT = 0.626    # cutoff for a "good" Spearman value. Is the value calculated for HRvPC1
Z_CUT = 0.2
HR_CUT = 0.7
ITERATIONS = 32000
# with our 103x4 array, we can only iterate a max of 103*102*101*100 = 106,110,600
# assuming we don't repeat
# Is this the number of data sets we can make, or the number of sn?
# if it is the number of sn then the number of iterations needs to max out at number of SN/103.


# import data for PCA analysis
## HR
HR = pd.read_csv('data/CampbellHoltzman_mb.tsv', sep='\\t',
                 usecols=['SNID', 'redshift', 'hr', 'err_mu'], index_col='SNID')
HR.rename(columns={'err_mu': 'hr uncert'}, inplace=True)
HR = HR[HR['redshift'] < Z_CUT]
HR = HR[HR['hr'] < HR_CUT]
print('Hubble Residual:')
print(HR.describe())

## SALT2 parameters (x_1 & c)
t = Table.read('data/SDSS_Photometric_SNe_Ia.fits')
salt = t['CID', 'Z', 'X1', 'X1_ERR', 'COLOR', 'COLOR_ERR'].to_pandas()
salt.columns = salt.columns.str.lower()
salt.rename(columns={'cid': 'SNID', 'z': 'redshift'}, inplace=True)
salt.set_index('SNID', inplace=True)
print('\nSALT2 parameters:')
print(salt.describe())

## stellar mass
galaxy = pd.read_csv('resources/kcorrect_stellarmass.csv',
                     usecols=['GAL', 'redshift', 'stellarmass'], index_col='GAL')
galaxy.rename(columns={'redshift': 'gal redshift', 'stellarmass': 'stellar mass'}, inplace=True)
print('\nGalaxy Stellar Mass:')
print(galaxy.describe())

## age
age = pd.read_csv('resources/ages_campbell.tsv', sep='\\t', skiprows=[1],
                  usecols=['# sn id', 'age'], dtype={'age': np.float64, '# sn id': np.int})
age.rename(columns={'# sn id': 'SNID'}, inplace=True)
age.set_index('SNID', inplace=True)
print('\nLocal Envir. Ages:')
print(age.describe())


# combine into on array
data = pd.concat([HR, salt, galaxy, age], axis=1)
data.dropna(inplace=True)
data['stellar mass'] = np.log10(data['stellar mass'])
print('\nAnalysis Data:')
print(data.describe())


# PCA prep: extraction and standardization
pca = PCA(n_components=4)
FEATURES = ['x1', 'color', 'stellar mass', 'age']
x = data.loc[:, FEATURES].values
x = StandardScaler().fit_transform(x)    # should be (103, 4) shaped array

# Set up significance counter
significance_counter = 0
max_rho = 0
rho = np.array([])

# Run iterations
for n in range(ITERATIONS):
    ## shuffle array
    np.random.shuffle(x)    # shuffles alone first axis, the 103

    ## compute PCA
    principal_components = pca.fit_transform(x)    # should be (103, 4) shaped array

    ## add up if something is significant
    for i in range(4):
        # get correlation -- [0]
        # don't worry if it is a negative correlation -- abs()
        rho_i = abs(spearmanr(principal_components[:, i], data.loc[:, 'hr'].values)[0])
        rho = np.append(rho, rho_i)

        ### count significant correlations
        if rho_i > SPEARMAN_CUT:
            significance_counter += 1

        ### update max rho
        if rho_i > max_rho:
            max_rho = rho_i

    ### update progress
    if ITERATIONS > 300:
        if n % 50 == 0:
            # over print a percentage
            # need an extra \n in next print statement
            print(f'Shuffling --- {100*n/ITERATIONS:.2f}% done', end='\r')

# report
np.savetxt('resources/PCASignificance.csv', rho, delimiter=',')
# need an extra \n because of \r in previous print statement
print(f'\n\nWith {ITERATIONS:,} iterations, there are {significance_counter} "significant" correlations.')
print(f'That is {significance_counter/(ITERATIONS*4)}% of possible correlations.')
print(f'The maximum correlation is {max_rho:.5f}')

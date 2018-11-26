"""script_package_data_for_unity.py --

converts my results to a stan pickel file

Also I should save out the age Gaussian mixture numbers to a csv.
"""
from glob import glob
import pickle
import re

import numpy as np
import pandas as pd
from astropy.table import Table
from sklearn.mixture import GaussianMixture
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

FIT_AGES = True
# todo make this a list and loop over it
DATASETS = ['campbell', 'campbellG']    # campbell, campbellG, both
N_AGE_MIX = 1    # should be >0
N_DATASETS = len(DATASETS)

# Get the first set of files to understand the number of SN.
files = glob(f'resources/SN*_{DATASETS[0]}_chain.tsv')
N_SNE = len(files)
# end = -11 - len(DATASET)
# get the numbers after the SN.
snids = map(lambda x: re.search('(?<=SN)\d*', x).group(0), files)
snids = list(map(int, snids))

amplitudes = np.zeros((N_DATASETS, N_SNE, N_AGE_MIX))
means = np.zeros((N_DATASETS, N_SNE, N_AGE_MIX))
stds = np.zeros((N_DATASETS, N_SNE, N_AGE_MIX))


# todo loop over results
for i, DATASET in enumerate(DATASETS):
    files = glob(f'resources/SN*_{DATASET}_chain.tsv')
    # todo redundant, but this comment took longer to write then the code took to run multiple times.



    # Get Host Properties
    ######################
    ## get age representations
    print('Importing age representations.')
    data = np.loadtxt(f'resources/age_{DATASET}_{N_AGE_MIX}gaus_representation.csv', delimiter=',')
    amplitudes[i] = data[:, 1:N_AGE_MIX+1]
    means[i] = data[:, N_AGE_MIX+1:2*N_AGE_MIX+1]
    stds[i] = data[:, 2*N_AGE_MIX+1:3*N_AGE_MIX+1]



## mass & redshift (may have lost a digit)
print('Pulling in mass')
mass = pd.read_csv('resources/kcorrect_stellarmass.csv', index_col='GAL',
                    usecols=['GAL', 'stellarmass', 'redshift'])
mass['stellarmass'] = np.log10(mass['stellarmass'])
mass = mass.loc[snids]
# mass.values       # a (103, 1) numpy array
# print(mass.head())
# print(mass.loc[[762, 1032, 1371]]).   # Check that masses are correct.


# Get SN Properties
###################
## redshift via either the mass, or campbell values

## Campbell stuff, mb, x1, c, covariance
print('Pulling in LC data')
lc_data = Table.read('data/SDSS_Photometric_SNe_Ia.fits')
lc_data = lc_data.to_pandas().set_index('CID')
if DATASET == 'both':
    lc_data = lc_data.loc[snids[:N_SNE//2]]
else:
    lc_data = lc_data.loc[snids]
lc_data = lc_data[['Z', 'Z_ERR', 'X0', 'X0_ERR', 'X1', 'X1_ERR', 'COLOR',
                   'COLOR_ERR', 'C01', 'C00', 'C11', 'C22', 'C02', 'C12',
                   'MU_MB', 'MU_ERR']]
# somehow this seems like a different redshift.
# The values in mass seem to match BOSS's host values. So may be miss matched hosts?
# print(lc_data.head())
# print(mass.loc[[1794, 6057]])


# Make sure all SN properties have the same order of content
# print(snids[:5])
# print(mass.index[:5])
# print(lc_data.index[:5])
# print(ages[:5])


print('Putting it all together')

m_b = -2.5*np.log10(lc_data['X0'].values)

if N_AGE_MIX == 1:
    # just make it a Gaussian parameter
    
    # Make covariance matrix
    # if DATASET == 'both':
    #     # todo define cov size as 4 + length of list (if Gaussian).
    #     cov = np.zeros((N_SNE//2, 6, 6))
    #     cov[:, 3, 3] = np.array([0.3**2]*(N_SNE//2))
    #     cov[:, 4, 4] = stds[:N_SNE//2].reshape(N_SNE//2)
    #     cov[:, 5, 5] = stds[N_SNE//2:].reshape(N_SNE//2)
    #     obs_mBx1c = [[m_b[i], lc_data['X1'].iloc[i], lc_data['COLOR'].iloc[i], 
    #                  mass['stellarmass'].iloc[i], means[:N_SNE//2], means[N_SNE//2:]] for i in range(N_SNE//2)]

    # else:
    cov = np.zeros((N_SNE, 4+N_DATASETS, 4+N_DATASETS))
    cov[:, 3, 3] = np.array([0.3**2]*N_SNE)
    # todo what if I N_DATASETS>1
    for i in range(N_DATASETS):
        cov[:, 4+i, 4+i] = stds[i].reshape(N_SNE)
    # Saddly this needs to be hard coded
    # np.stack also does not help
    if N_DATASETS == 1:
        obs_mBx1c = [[m_b[j], lc_data['X1'].iloc[j], lc_data['COLOR'].iloc[j], 
                      mass['stellarmass'].iloc[j], means[0, j]] for j in range(N_SNE)]
    elif N_DATASETS == 2:
        # (N_DATASETS, N_SNE, N_AGE_MIX)
        obs_mBx1c = [[m_b[j], lc_data['X1'].iloc[j], lc_data['COLOR'].iloc[j], 
                      mass['stellarmass'].iloc[j], means[0, j], means[1, j]] for j in range(N_SNE)]
    else:
        raise RuntimeError('Adding ages to obs_mBx1c is hard coded and needs your attention.')

    # With N_AGE_MIX == 1, don't use the non-gaussian parameters
    n_props=4+N_DATASETS
    n_non_gaus_props = 0
    n_age_mix = 0
    age_gaus_mean = np.zeros((n_non_gaus_props, N_SNE, n_age_mix))
    age_gaus_std = np.zeros((n_non_gaus_props, N_SNE, n_age_mix))
    age_gaus_A = np.zeros((n_non_gaus_props, N_SNE, n_age_mix))

else:
    # Make covariance matrix
    # todo what if I have multiple datasets?
    cov = np.zeros((N_SNE, 4, 4))
    obs_mBx1c=[[m_b[i], lc_data['X1'].iloc[i], lc_data['COLOR'].iloc[i], 
                           mass['stellarmass'].iloc[i]] for i in range(N_SNE)]
    
    n_props=4
    n_non_gaus_props = N_DATASETS
    n_age_mix=N_AGE_MIX
    #these likely don't work.
    age_gaus_mean=np.expand_dims(means, 0)
    age_gaus_std=np.expand_dims(stds, 0)
    age_gaus_A=np.expand_dims(amplitudes, 0)



# update other parts of the covariance matrix.
cov[:, 0, 0] = (2.5/np.log(10))**2*(lc_data['C00'].values/lc_data['X0']**2)
cov[:, 0, 1] = (-2.5/np.log(10))*(lc_data['C01'].values/lc_data['X0'])
cov[:, 0, 2] = (-2.5/np.log(10))*(lc_data['C02'].values/lc_data['X0'])
cov[:, 1, 0] = (-2.5/np.log(10))*(lc_data['C01'].values/lc_data['X0'])
cov[:, 1, 1] = lc_data['C11'].values
cov[:, 1, 2] = lc_data['C12'].values
cov[:, 2, 0] = (-2.5/np.log(10))*(lc_data['C02'].values/lc_data['X0'])
cov[:, 2, 1] = lc_data['C12'].values
cov[:, 2, 2] = lc_data['C22'].values


# Get data set name
if N_DATASETS == 1:
    data_name = f'{DATASETS[0]}_{N_AGE_MIX}gaus'
else:
    data_name = f''
    for i in DATASETS:
        data_name += f'{i}'
    data_name += f'_{N_AGE_MIX}gaus'

print('Saving results')
pickle.dump(dict(     # general properties
                 n_sne=N_SNE, n_props=n_props, n_non_gaus_props=n_non_gaus_props, n_sn_set=1,
                 sn_set_inds=[0]*N_SNE,
                      # redshifts
                 z_helio=mass['redshift'], z_CMB=mass['redshift'],
                      # Gaussian defined properties
                 obs_mBx1c=obs_mBx1c, obs_mBx1c_cov=cov,
                      # Non-Gaussian properties, aka age
                 n_age_mix=n_age_mix, age_gaus_mean=age_gaus_mean,
                 age_gaus_std=age_gaus_std, age_gaus_A=age_gaus_A,
                      # Other stuff that does not really need to change
                 do_fullDint=0, outl_frac_prior_lnmean=-4.6, outl_frac_prior_lnwidth=1,
                 lognormal_intr_prior=0, allow_alpha_S_N=0),
            open(f'{data_name}_forUnity.pkl', 'wb'))

print('Done!')

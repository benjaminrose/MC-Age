"""script_package_results_for_unity.py

converts my results to a stan pickel file

Also I should save out the age Gaussian mixture numbers to a csv.
"""
# todo fit ages once and compile unity data separates.
from glob import glob
import re

import numpy as np
from sklearn.mixture import GaussianMixture
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt


def main(DATASET='campbell', N_AGE_MIX=1):
    """Fit and save the the an N_AGE_MIX number of Gaussians mixture to the given DATASET.
    
    DATASET: str
        campbell, campbellG, ...,
    N_AGE_MIX: int
        Number of Gaussian mixtures to use. Should be >0.
    """
    files = glob(f'resources/SN*_{DATASET}_chain.tsv')
    N_SNE = len(files)
    # end = -11 - len(DATASET)
    # get the numbers after the SN.
    snids = map(lambda x: re.search('(?<=SN)\d*', x).group(0), files)
    snids = list(map(int, snids))


    model = GaussianMixture(N_AGE_MIX)
    amplitudes = np.zeros((N_SNE, N_AGE_MIX))
    means = np.zeros((N_SNE, N_AGE_MIX))
    stds = np.zeros((N_SNE, N_AGE_MIX))

    print(f'Fitting ages to {N_AGE_MIX} Gaussians')
    pdf = PdfPages(f'resources/age_{DATASET}_{N_AGE_MIX}gaus_representation_preview.pdf')

    for i, f in enumerate(files):
        data = np.genfromtxt(f, delimiter='\t')
        data = data[:, 7]

        model.fit(np.expand_dims(data, 1))

        amplitudes[i] = model.weights_.reshape(N_AGE_MIX)
        means[i] = model.means_.reshape(N_AGE_MIX)
        stds[i] = np.sqrt(model.covariances_).reshape(N_AGE_MIX)

        plt.figure()
        plt.hist(data, bins=np.linspace(-5, 20, 200))
        plt.hist(model.sample(1020000)[0], alpha=0.5, bins=np.linspace(-5, 20, 200))
        plt.title(f)
        
        pdf.savefig()
        plt.close()

        if (i+1)%10 == 0:
            print(f'Finished with the {i+1}th age fit')

    pdf.close()

    # if DATASET != 'both':
    ages = np.column_stack((snids, amplitudes, means, stds))
    # todo update the header to match the number of Gaussians used.
    np.savetxt(f'resources/age_{DATASET}_{N_AGE_MIX}gaus_representation.csv', ages, delimiter=',',
               header='sn id, amp_1, amp_2, amp_3, mean_1, mean_2, mean_2, std_1, std_2, std_3')
    
    print(f'Done with {N_AGE_MIX} Gaussian mixture for {DATASET}.')



if __name__ == '__main__':
    N = [1, 2]
    D = ['campbell', 'campbellG']
    for i in N:
        for j in D:
            main(j, i)
    print('Done!')
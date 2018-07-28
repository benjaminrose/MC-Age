"""CombineAges.py -- Script to combine MCMC ages into one csv file."""

__author__ = '''* Benjamin Rose
* brose3@nd.edu
* benjamin.rose@me.com
* University of Notre Dame
* 2018-07-27
* Python 3.6'''

import glob
import re

import numpy as np
from astropy.table import Table

import util

def combine_ages(dataset: str) -> None:
    """Combine the ages into one file.

    Saves the median +/- a default interval.

    Parameters
    ----------
    dataset :
        Selects a dataset to combine. Works with ``'gupta'``, ``'messier'``,
        ``'circle'``, ``'campbell'``, ``'campbellG'``, ``'riess'``, or
        ``'riessL'``,

    Raises
    ------
    ValueError
        Raised if and invalid ``dataset`` is used.

    """
    # import & calculate median (& +/- intervals)
    if dataset == 'gupta':
        # find files with 4 or more characters after 'SN'
        # files = glob.glob('../resources/SN????*_chain.tsv')
        # will be 4 or more numbers and dataset name
        files = glob.glob('../resources/SN????*_gupta_chain.tsv')
    elif dataset == 'messier':
        # find files with 2 and 3 digits after 'SN'
        # files1 = glob.glob('../resources/SN'+'[0-9]'*2+'_chain.tsv')
        # files2 = glob.glob('../resources/SN'+'[0-9]'*3+'_chain.tsv')
        # files = files1 + files2  # Yeah python list concatenation!!
        # will be anything and dataset name
        files = glob.glob('../resources/SN*_messier_chain.tsv')
    elif dataset == 'circle':
        # find files with 1 digits after 'SN'
        files = glob.glob('../resources/SN[0-9]_chain.tsv')
        # will be 1 number and dataset name
        # files = glob.glob('../resources/SN[0-9]_circle_chain.tsv')
    elif dataset == 'campbell':
        # files = glob.glob('../resources/SN????*_chain.tsv')
        # will be 4 or more numbers and dataset name
        files = glob.glob('../resources/SN????*_campbell_chain.tsv')
    elif dataset == 'campbellG':
        # files = glob.glob('../resources/SN????*_chain.tsv')
        # will be 4 or more numbers and dataset name
        files = glob.glob('../resources/SN????*_campbellG_chain.tsv')
    elif dataset == 'riess':
        # The H-0 calibration sample global age
        files = glob.glob('../resources/SN*_riess_chain.tsv')
    elif dataset == 'riessL':
        # The H-0 calibration sample age at the location of SN
        files = glob.glob('../resources/SN*_riessL_chain.tsv')
    else:
        raise ValueError
    print("Collecting {} ages for the {} stellar populations".format(
          len(files), dataset))

    ages = np.array([])

    for f in files:
        print('reading {}'.format(f))
        #this is very slow
        data = Table.read(f, format='ascii.csv', delimiter='\t',
                          data_start=2)
        data = data.to_pandas()
        data.dropna(inplace=True)
        # calculate data an append it to ages variable
        # Is there was a way that does not calculate `len` each time?
        if len(ages) == 0:
            ages = np.array([util.median(data['age'])])
        else:
            ages = np.vstack((ages, np.array([util.median(data['age'])])))


    # get SN ID's
    # Don't get 'SN' because that will force it to be a string and mess up `to_save`
    print('getting SN IDs')
    p = re.compile(r'\d+')
    # add a holding value so I can append to
    names = np.array([1])
    for i in files:
        # save as float so `to_save` works
        names = np.append(names, float(p.findall(i)[0]))
    # remove holding value
    names = names[1:]
    # "reshape" so np.hstack works
    names = names.reshape(len(names), 1)

    # save names and ages to one variable
    to_save = np.hstack((names, ages))

    # Save data
    print('saving data')
    header = 'sn id\tage\tlower limit\tupper limit\n\tGyr\tGyr\tGyr'
    with open('../resources/ages_{}.tsv'.format(dataset), 'wb') as f:
        np.savetxt(f, to_save, delimiter='\t', header=header)


if __name__ == '__main__':
    combine_ages('riessL')

""" util.py -- holds utility functions for fspsage

Benjamin Rose
brose3@nd.edu
benjamin.rose@me.com
University of Notre Dame
2017-11-24
Python 3.5
"""

import glob
import re

import numpy as np
from astropy.table import Table
import fsps

def combine_ages(dataset):
    """combines the ages into one file

    Parameters
    ----------
    dataset :
        Selects a dataset to combine. Works with 'gupta', 'messier',
        'circle', or 'campbell'.

    Raise
    -----
    ValueError
        Raised if and invalid `dataset` is used.
    """

    # import & calculate median (& +/- intervals)
    if dataset == 'gupta':
        # find files with 4 or more characters after 'SN'
        # files = glob.glob('resources/SN????*_chain.tsv')
        # will be 4 or more numbers and dataset name
        files = glob.glob('resources/SN????*_gupta_chain.tsv')
    elif dataset == 'messier':
        # find files with 2 and 3 digits after 'SN'
        # files1 = glob.glob('resources/SN'+'[0-9]'*2+'_chain.tsv')
        # files2 = glob.glob('resources/SN'+'[0-9]'*3+'_chain.tsv')
        # files = files1 + files2  # Yeah python list concatenation!!
        # will be anything and dataset name
        files = glob.glob('resources/SN*_messier_chain.tsv')
    elif dataset == 'circle':
        # find files with 1 digits after 'SN'
        files = glob.glob('resources/SN[0-9]_chain.tsv')
        # will be 1 number and dataset name
        # files = glob.glob('resources/SN[0-9]_circle_chain.tsv')
    elif dataset == 'campbell':
        # files = glob.glob('resources/SN????*_chain.tsv')
        # will be 4 or more numbers and dataset name
        files = glob.glob('resources/SN????*_campbell_chain.tsv')
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
            ages = np.array([median(data['age'])])
        else:
            ages = np.vstack((ages, np.array([median(data['age'])])))
   

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
    with open('resources/ages_{}.tsv'.format(dataset), 'wb') as f:
        np.savetxt(f, to_save, delimiter='\t', header=header)

def median(data, interval=34):
    """calculates the median of a 1D array and ``interval`` size quartile uncertainties.

    Parameters
    ----------
    data: numpy.nparray
        Needs to be an 1D array.

    interval: int
        The percentile (to one side) you want for uncertainties. 

    Returns
    -------
    numpy.nparray
        Three values: median, upper and lower uncertainties.
    """
    low = 50 - interval     # 16 by default
    high = 50 + interval    # 84 by default
    
    # map(
    #     lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
    #     zip(*np.percentile(samples, [16, 50, 84], axis=0))
    # )

    # w/ 500,000 runs each, a few non-physical ages slip in and produce a nan
    v = np.nanpercentile(data, [low, 50, high])
    # return np.array([v[1], v[2]-v[1], v[1]-v[0]])
    return np.array([v[1], v[0], v[2]])

def mode(data, interval=34):
    """calculates and returns the "mode" of a continuous variable and also returns the interval percentile.

    Parameters
    ----------
    data: numpy.nparray
        Needs to be an 1D array.

    interval: int
        The percentile (to one side) you want for range 

    Returns
    -------
    numpy.nparray
        Three values: median, lower value and upper value. Note these are not
        uncertainties but the value of the parameter at these intervals.
    """
    low = 50 - interval     # 16 by default
    high = 50 + interval    # 84 by default
    
    
    hist = np.histogram(data, bins='fd')
    # hist[0] are frequencies
    # hist[1] are the bins

    # there has to be a better way to get the index of the max value of hist
    peak_val = np.max(hist[0])

    #What if two mosts?

    mode_index = np.where(np.equal(hist[0], peak_val))[0]
    if mode_index.size > 1:
        if mode_index[0] + 1 == mode_index[1]:
            # take the middle bin edge if mode is two consecutive bins
            mode = hist[1][mode_index[1]]
        else:
            # if there are two non-consecutive modes
            # raise RuntimeWarning('There exists two non-consecutive modes')
            print('\nTwo non-consecutive modes')
            print(mode_index)
            print(hist[1][mode_index])
            print('\n')
            mode = hist[1][mode_index[0]]
    else:
        # hist[1] is bin edges. I need the mid-point
        mode = (hist[1][mode_index] + hist[1][mode_index+1])/2
    # remove the excess numpy array wrapper
    if mode.size == 1:
        # most times mode is an ndarray of size 1. Onetime it was a float.
        # It looks like the two non-consecutive modes workaround causes this.
        if type(mode) is not np.float64:
            try:
                mode = mode[0]
            except:
                import pdb; pdb.set_trace()
    else:
        raise ValueError('somehow there are still two modes')

    v = np.nanpercentile(data, [16, 84])

    return np.array([mode, v[0], v[1]])

def get_sp():
    """returns the default stellar population.

    This is for constancy throughout the program.
    """
    return fsps.StellarPopulation(zcontinuous=2, cloudy_dust=True,
                                 add_neb_emission = True, sfh=5)


if __name__ == '__main__':
    # data = np.append(np.random.randn(600)*4, (np.random.randn(1000)+7))
    # data = np.random.randn(1000)
    # print('mode result: ', mode(data))#, interval=49))
    # import matplotlib.pyplot as plt
    # plt.hist(data, bins='auto')
    # plt.savefig('del.pdf')

    combine_ages('messier')
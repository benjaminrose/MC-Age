""" util.py -- holds utility functions for fspsage

* Benjamin Rose
* brose3@nd.edu
* benjamin.rose@me.com
* University of Notre Dame
* 2017-11-24
* Python 3.5
"""
import numpy as np
import pandas as pd
from sfdmap import SFDMap   # needs FITS dust maps installed
from extinction import fitzpatrick99

# ignore fsps on Read The Docs
try:
    from os import environ
    environ['READTHEDOCS']
except KeyError:
    import fsps

def median(data, interval=34):
    """Calculate the median and ``interval`` size quartile uncertainties.

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
    """Calculate the "mode" of a continuous variable.

    Also returns the interval percentile.

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
                                  add_neb_emission=True, sfh=5)

def correct_dust(data: pd.DataFrame, inplace: bool = False) -> pd.DataFrame:
    """Corrects a data set of SED for dust.

    Estimates the E(B-V) values from Schlegel, Finkbeiner & Davis (1998),
    but uses the scaling of 0.86 is to reflect the recalibration
    by Schlafly & Finkbeiner (2011). Also uses Fitzpatrick (1999) dust
    extinction function, with an R_V = 3.1

    Note this acts strange if passed part of a data frame. Either use:
    ``correct_dust(df)`` or ``correct_dust(df.loc[0:2].copy())`` to make
    sure the a copy/instance is passed rather than a view. **Caveat:** this
    was discovered and tested before ``inplace`` was set as an option.
    Using the default of ``inplace=False`` seems to fix this at first glance.

    Parameters
    ----------
    data: pandas.DataFrame
        The data that you want dust corrected. Need headers 'u', 'g',
        'r', 'i', 'z', 'ra', 'dec'. RA and Dec should be floats in degrees.

    inplace: bool (False)
        Defaults to making an internal copy of the DataFrame and passing
        that back. If set to `False` this will operate directly on in
        the DataFrame.

    Returns
    -------
    pandas.DataFrame
        It returns the same data set, but with the 'u', 'g', 'r', 'i',
        'z' dust corrected

    """
    # With out this, we operate directly in the DataFrame.
    if not inplace:
        data = data.copy()

    # read from file, not environment variable
    dust_map = SFDMap('sfddata-master')

    #set up useful numbers
    SED = np.array([3543, 4770, 6231, 7625, 9134])   # assumes ugriz in AA

    # for each deredden:
    # https://stackoverflow.com/questions/16476924/how-to-iterate-over-rows-in-a-dataframe-in-pandas#16476974
    for index, row in data.iterrows():
        # get E(B-V) at location, difference in extinction between the B and V bands
        # get RA and Dec and then expand it out to the two needed arguments
        ebv = dust_map.ebv(*row.loc[['ra', 'dec']])

        # deredden the DataFrame
        data.loc[index, ['u', 'g', 'r', 'i', 'z']] -= fitzpatrick99(SED, 3.1*ebv)


    #return data after deredening
    return data

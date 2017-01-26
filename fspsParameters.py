""" fspsParameters.py -- A file to help determine the correct FSPS parameters to use.

Benjamin Rose
brose3@nd.edu
benjamin.rose@me.com
University of Notre Dame
2017-01-20
Python 3
"""
import numpy as np
import fsps

def test_compute_vega_mags():
    """
    Does explicitly turning off `compute_vega_mags` matter in the output of fsps?

    # Parameters

    sp : fsps.StellarPopulation
        An instance of the FSPS StellarPopulation class. 
    """
    filters = ['sdss_u', 'sdss_g', 'sdss_r', 'sdss_i', 'sdss_z']
    age = 13.7

    #somewhere someone said to only have one instance of 
    #`fsps.StellarPopulation()` but to test the defaults I need a clean 
    #instance.
    sp = fsps.StellarPopulation()

    # test default (should be `False`)
    mags1 = sp.get_mags(tage=age, bands=filters)

    #test like what many do (set to `False`)
    #We have to reinitialize, because `compute_vega_mags` can only be changed 
    #at initiation.
    sp = fsps.StellarPopulation(compute_vega_mags = False)
    mags2 = sp.get_mags(tage=age, bands=filters)

    #compare to the milimag (ish)
    #could use `np.allclose` if I want to drop this down to a single boolean
    return np.isclose(mags1, mags2)

def test_add_neb_emission():
    """
    This function shows the difference of turning on the nebular emissions. We 
    expect that for young stars we can get a *** r-i slop with nebular 
    emissions on, but every other case (old-with/without, young/without) does 
    not.
    """
    # set up basics
    sdss_bands = ['sdss_u', 'sdss_g', 'sdss_r', 'sdss_i', 'sdss_z']
    sp = fsps.StellarPopulation(zcontinuous=2, logzsol=0, dust1=1.0, dust2=0.5,
        cloudy_dust=True, sfh=4)

    # create SED for old stellar population-without emission
    sp.params['tau'] = 0.5
    sp.params['sf_start'] = 1
    #set age to a redshift of 0.2
    oldwo = sp.get_mags(tage=11.4, redshift=0.2, bands=sdss_bands)

    # create SED for old stellar population-with emission
    sp.params['add_neb_emission'] = True
    oldw = sp.get_mags(tage=11.4, redshift=0.2, bands=sdss_bands)

    # create SED for young stellar population-without emission
    sp.params['tau'] = 5.0
    sp.params['sf_start'] = 9
    youngwo = sp.get_mags(tage=11.4, redshift=0.2, bands=sdss_bands)

    # create SED for young stellar population-with emission
    sp.params['add_neb_emission'] = True
    youngw = sp.get_mags(tage=11.4, redshift=0.2, bands=sdss_bands)

    #scale them so they all have the same r-band mag (of 45 "mags")
    #this simulates the best fit that determines the total stellar mass.
    print('pre correction:  ', oldwo)
    def correctSED(SED):
        delta = 45 - SED[2]
        return SED + delta
    oldwo = correctSED(oldwo)
    oldw = correctSED(oldw)
    youngwo = correctSED(youngwo)
    youngw = correctSED(youngw)
    print('post correction: ', oldwo)


    # plot 
    import matplotlib.pyplot as plt
    import seaborn as sns
    fig = plt.figure('neb emission test')
    ax = fig.add_subplot(111)
    # x = [0,1,2,3,4]
    x = [3551, 4686, 6166, 7480, 8932]
    plt.plot(x, oldwo, label='old - without')
    plt.plot(x, oldw, '-.', label='old - with')
    plt.plot(x, youngwo, ':', label='young - without')
    plt.plot(x, youngw, '--', label='young - with')
    plt.xticks(x)
    ax.set_xticklabels(['u','g','r','i','z'])
    plt.gca().invert_yaxis()
    plt.legend(loc=4)
    plt.show()

    print(oldwo)
    print(youngwo)
    print(oldw)
    print(youngw)
    return None
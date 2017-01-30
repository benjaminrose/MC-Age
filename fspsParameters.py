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
    sp.params['tau'] = 10.0
    sp.params['sf_start'] = 6.4
    #set age to a redshift of 0.2
    oldwo = sp.get_mags(tage=11.4, redshift=0.2, bands=sdss_bands)

    # create SED for old stellar population-with emission
    sp.params['add_neb_emission'] = True
    oldw = sp.get_mags(tage=11.4, redshift=0.2, bands=sdss_bands)

    # create SED for young stellar population-without emission
    sp.params['add_neb_emission'] = False
    # sp.params['tau'] = 10.0
    sp.params['sf_start'] = 11
    youngwo = sp.get_mags(tage=11.4, redshift=0.2, bands=sdss_bands)

    # create SED for young stellar population-with emission
    sp.params['add_neb_emission'] = True
    youngw = sp.get_mags(tage=11.4, redshift=0.2, bands=sdss_bands)

    # create SED for youngest stellar population-without emission
    sp.params['add_neb_emission'] = False
    # sp.params['tau'] = 10.0
    sp.params['sf_start'] = 11.35
    youngestwo = sp.get_mags(tage=11.4, redshift=0.2, bands=sdss_bands)

    # create SED for youngest stellar population-with emission
    sp.params['add_neb_emission'] = True
    youngestw = sp.get_mags(tage=11.4, redshift=0.2, bands=sdss_bands)

    #scale them so they all have the same r-band mag (of 45 "mags")
    #this simulates the best fit that determines the total stellar mass.
    def correctSED(SED):
        delta = 21 - SED[2]
        return SED + delta
    oldwo = correctSED(oldwo)
    oldw = correctSED(oldw)
    youngwo = correctSED(youngwo)
    youngw = correctSED(youngw)
    youngestwo = correctSED(youngestwo)
    youngestw = correctSED(youngestw)

    # plot 
    import matplotlib.pyplot as plt
    import seaborn as sns
    fig = plt.figure('neb emission test')
    ax = fig.add_subplot(111)
    # x = [0,1,2,3,4]
    x = [3551, 4686, 6166, 7480, 8932]
    plt.plot(x, oldwo, label='5 Gyr - no emission')
    plt.plot(x, oldw, '--', label='5 Gyr - emission')
    plt.plot(x, youngwo, '-', label='400 Myr - no emission')
    plt.plot(x, youngw, '-.', label='400 Myr - emission')
    plt.plot(x, youngestwo, '--', label='50 Myr - no emission')
    plt.plot(x, youngestw, ':', label='50 Myr - emission')
    plt.xticks(x)
    ax.set_xticklabels(['u','g','r','i','z'])
    # ax.set_yticklabels(['u','g','r','i','z'])
    plt.ylabel('theoretical magnitudes')
    plt.gca().invert_yaxis()
    plt.legend(loc=4)
    plt.savefig('figures/add_neb_emission.pdf')
    plt.show()

    print(oldwo)
    print(oldw)
    print(youngwo)
    print(youngw)
    print(youngestwo)
    print(youngestw)
    return None

def test_neg_ri_color(NebEmission=False):
    """
    Tried to see if the r-i color was age only or just nebular emission.
    Saves a figure to either `figures/oldVSyoung.pdf` or`figures/oldVSyoung_withNebEmission.pdf` depending on input parameter.

    Was helpful for determing if negative r-i color was from age alone or 
    needed `add_neb_emission`.

    Parameters 
    ----------
    NebEmission : bool
        Determines if it outputs a figure of varying ages with our without emission turned on. 
    """
    sdss_bands = ['sdss_u', 'sdss_g', 'sdss_r', 'sdss_i', 'sdss_z']
    sp = fsps.StellarPopulation(zcontinuous=2, logzsol=0, dust1=1.0, dust2=0.5,
        cloudy_dust=True, sfh=4)
    if NebEmission:
        sp.params['add_neb_emission'] = True

    # create SED for oldest stellar
    sp.params['tau'] = 0.5
    sp.params['sf_start'] = 1
    #set age to a redshift of 0.2
    oldest = sp.get_mags(tage=11.4, redshift=0.2, bands=sdss_bands)

    # create SED for old stellar
    sp.params['tau'] = 2
    sp.params['sf_start'] = 4
    #set age to a redshift of 0.2
    old = sp.get_mags(tage=11.4, redshift=0.2, bands=sdss_bands)

    # create SED for young stellar
    sp.params['tau'] = 5
    sp.params['sf_start'] = 7
    #set age to a redshift of 0.2
    young = sp.get_mags(tage=11.4, redshift=0.2, bands=sdss_bands)

    # create SED for younger stellar
    sp.params['tau'] = 10
    sp.params['sf_start'] = 11
    #set age to a redshift of 0.2
    younger = sp.get_mags(tage=11.4, redshift=0.2, bands=sdss_bands)

    # create SED for youngest stellar
    sp.params['tau'] = 10
    sp.params['sf_start'] = 11.35
    #set age to a redshift of 0.2
    youngest = sp.get_mags(tage=11.4, redshift=0.2, bands=sdss_bands)

    def correctSED(SED):
        delta = 45 - SED[2]
        return SED + delta
    oldest = correctSED(oldest)
    old = correctSED(old)
    young = correctSED(young)
    younger = correctSED(younger)
    youngest = correctSED(youngest)

    # plot 
    import matplotlib.pyplot as plt
    import seaborn as sns
    fig = plt.figure('neb emission test')
    ax = fig.add_subplot(111)
    # x = [0,1,2,3,4]
    x = [3551, 4686, 6166, 7480, 8932]
    plt.plot(x, oldest, label='oldest - 10.5 Gyr')
    plt.plot(x, old, '-.', label='old - 7.4 Gyr')
    plt.plot(x, young, ':', label='young - 4.4 Gyr')
    plt.plot(x, younger, ':', label='younger - 400 Myr')
    plt.plot(x, youngest, '--', label='youngest - 50 Myr')
    plt.xticks(x)
    ax.set_xticklabels(['u','g','r','i','z'])
    plt.gca().invert_yaxis()
    plt.legend(loc=4)
    if NebEmission:
        plt.savefig('figures/oldVSyoung_withNebEmission.pdf')
    else:
        plt.savefig('figures/oldVSyoung.pdf')
    plt.show()
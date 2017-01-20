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
    """

    # read in data

    # determine best candidate

    # get best fit model without nebular emission

    # get best fit model with nebular emission

    #compare two results
    return False
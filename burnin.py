"""burnin.py - a method to verify mcmc burnin length

Benjamin Rose
2017-06-22
Python 3.5
"""
import logging

import numpy as np

from calculateAge import calculateSFH

module_logger = logging.getLogger("fsps-age.burnin")

def burnin(SED, SEDerr, redshift, SNID=None):
    """
    """
    # Test inputs
    if not isinstance(SED, (list, np.ndarray)):
        raise TypeError('SED must be of type list or np.ndarray')
    if not isinstance(SEDerr, (list, np.ndarray)):
        raise TypeError('SEDerr must be of type list or np.ndarray')
    if not len(SED) == len(SEDerr) == 5:
        raise ValueError("The SED and it's errors need to represent the 5 SDSS filters")
    if redshift <= 0:
        raise ValueError("Redshifts should be greater than zero.")
    
    # Set up logger
    logger = logging.getLogger("fsps-age.burnin.burnin")

    # call `calculateSFH with burnin
    # size should be (50, 1000, 7)
    sampler = calculateSFH(SED, SEDerr, redshift, SNID, burnin=True)

    samples = sampler.chain             # size == (nwalker, nsteps, ndim)
    lnprob = sampler.lnprobability     # size == (nwalkers, nsteps)


    #save data
    # Note header should be:
    # logzsol, dust2, tau, tStart, sfTrans, sfSlope, c
    # dex, 1/Gyr, Gyr, Gyr, , mag
    np.save('resources/burnin/samples-{}'.format(SNID), samples)
    np.save('resources/burnin/lnprob-{}'.format(SNID), lnprob)
    logger.info('saved samples and ln probability')
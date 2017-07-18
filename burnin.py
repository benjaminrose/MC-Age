"""burnin.py - a method to verify mcmc burnin length

Benjamin Rose
2017-06-22
Python 3.5
"""
import logging

import numpy as np

from calculateAge import calculateSFH

module_logger = logging.getLogger("fsps-age.burnin")

def burnin(SED, SEDerr, redshift):
    """
    """
    # Set up logger
    logger = logging.getLogger("fsps-age.burnin.burnin")

    # call `calculateSFH with burnin
    # size should be (1000, 28, 7)
    sampler = calculateSFH(SED, SEDerr, redshift, burnin=True)

    samples = sampler.chain             # size == (nwalker, nsteps, ndim)
    samples = sampler.lnprobability     # size == (nwalkers, nsteps)


    #save data
    # Note header should be:
    # logzsol, dust2, tau, tStart, sfTrans, sfSlope, c
    # dex, 1/Gyr, Gyr, Gyr, , mag
    np.save('resources/burnin/samples', samples)
    np.save('resources/burnin/lnprob', arr)
    logger.info('saved samples and ln probability')
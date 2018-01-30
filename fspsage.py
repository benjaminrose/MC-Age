"""fsps-age.py -- Estimates the age of a Photometric SED using FSPS.

Science goal:
    Check local environment effects on SNIa by looking for correlations
    between HR and the age of the local environment calculated from SDSS Scene
    Modeling Photometry

Usage:
    fspsage.py burnin OBJID
    fspsage.py run DATASET JOBID JOBLENGTH [-d | --debug]
    fspsage.py (-h | --help)
    fspsage.py --version

Option:
    burnin          run a shorter run on specific objects only
    OBJID           the SN (or Messier) ID of the object to observe
    run             estimate age for a given data set
    DATASET         analyses circle, messier, gupta or campbell data sets
    JOBID           the ID for the piece of the data set to be analyzed
    JOBLENGTH       the total number of objects looked at
    -d --debug      run shorter and with more logs
    -h --help       show this screen
    --version       show version

Benjamin Rose
brose3@nd.edu
benjamin.rose@me.com
University of Notre Dame
2017-01-19
Python 3.5
"""
"""future options.
    main.py testFspsParameters [-d | --debug]
    main.py calculateAge (--sn=SNID) [-d | --debug]
    main.py age NAME U G R I Z U_ERR G_ERR R_ERR I_ERR Z_ERR
    main.py local [-d | --debug]
"""

__author__ = "Benjamin Rose"
__version__ = "0.2"
__license__ = "MIT"

import logging

from docopt import docopt
import numpy as np
import fsps

def formatLogging():
    """Setting up the logger parts that are the same for all commands
    """
    # update logging formating
    formatter = logging.Formatter(
        '%(asctime)s:%(name)s:%(levelname)s - %(message)s',
        datefmt='%Y-%m-%dT%I:%M:%S')
    fh.setFormatter(formatter)
    # add handler to logger object
    logger.addHandler(fh)
    logger.info('Starting')


def testFspsParameters():
    import fspsParameters

    # print(fspsParameters.test_compute_vega_mags())
    # fspsParameters.test_neg_ri_color()
    # fspsParameters.test_add_neb_emission()


def calculateAge(SNID):
    """
    # Get MCMC running and be able to calculate an age with uncertainties.
    """
    import calculateAge
    #####
    # unit tests
    #####
    # runFSPS()
    # - should return a list with 5 arguments
    # sp = fsps.StellarPopulation(zcontinuous=2, logzsol=0, dust1=1.0,
    #    dust2=0.5, cloudy_dust=True, sfh=4)
    # redshift = 0.084
    logzsol = 0.0
    dust2 = 0.1
    tau = 1
    tStart = 1
    sfTrans = 10
    sfSlope = 1
    # results = calculateAge.runFSPS(sp, redshift, logzsol, dust2, tau,
    #    tStart, sfTrans, sfSlope)
    # print('runFSPS(): ', results)
    # currently succeeds! We get a list with 5 parameters

    # lnlike()
    # - should return a float
    # [logzsol, dust2, tau, tStart, sfTrans, sfSlope, c]`.
    c = -20
    theta = [logzsol, dust2, tau, tStart, sfTrans, sfSlope, c]
    # from SN12781, or 2006er just using the values in the file
    # SED = np.array([24.41, 23.92, 23.08, 22.68, 22.01])
    # SEDerr = np.array([0.49, 0.10, 0.05, 0.05, 0.10])
    # redshift = 0.084
    #from SN10028, or ? just using the values in the file
    SED = np.array([21.22, 19.45, 18.64, 18.27, 17.98])
    SEDerr = np.array([0.041, 0.004, 0.019, 0.012, 0.004])
    redshift = 0.065

    # from SN15776 (global) that should be red and dead, but on 2017-05-11 was young and dusty
    SNID = 15776
    SED = np.array([23.14426, 21.00639, 19.41827, 18.82437, 18.46391])
    SEDerr = np.array([0.8009404, 0.04814964, 0.01899451, 0.01771959, 0.04554904])
    redshift = 0.305
    # results = calculateAge.lnlike(theta, SED, SEDerr, redshift, sp)
    # print('lnlike(): ')
    # print(type(results))
    # print(results)
    # print('where is this X < 0 coming from?')

    # calculateSFH()
    # sp = fsps.StellarPopulation(zcontinuous=2, logzsol=0, dust1=1.0, dust2=0.5,
            # cloudy_dust=True, sfh=4)
    # SED = np.array([24.41, 23.92, 23.08, 22.68, 22.01])
    # SEDerr = np.array([0.49, 0.10, 0.05, 0.05, 0.10])
    # redshift = 0.084
    # results = calculateAge.calculateSFH(SED, SEDerr, redshift, 10028)#, sp=sp)
    # print('calculateSFH(): ', results)
    # results = calculateAge.calculateAge(redshift, SED, SEDerr, SNID=10028)
    # results = calculateAge.calculateAge(redshift, SED, SEDerr, SNID=10028, debug=True)
    # print('calculateAge(): ', results)
    # Currently Fails!!


    # from SN12781, or 2006er just using the values in the file
    # SED = np.array([24.41, 23.92, 23.08, 22.68, 22.01])
    # SEDerr = np.array([0.49, 0.10, 0.05, 0.05, 0.10])
    # redshift = 0.084
    # results = calculateAge.calculateAge(redshift, SED, SEDerr)
    # print('calcualteAge()')


def redoGupta(cli):
    """
    # Test on global SED's

    We want to redo what Gupta did to make sure we can actually do something
    before we analyze on new data.

    :param cli:
        The dictionary of the CLI constructed by :docop:.
    """
    """ run redoGupta.py

    Parameters
    ----------
    cli : dictionary
        The result of docopt's parsing of the CLI.
    """
    import redoGupta

    if not cli['DATASET'] in ['gupta', 'messier', 'circle', 'campbell']:
        raise ValueError("Please use '--help' option for details on DATASET")
    
    redoGupta.redoGupta(int(cli['JOBID']), int(cli['JOBLENGTH']),
                        cli['--debug'], cli['DATASET'])


def burnin(cli):
    """
    Runs a smaller `emcee.sampler` to see how the sampling is progressing.
    Hopefully these values will not effect how it samples, but it will likely
    effect how well it sampled.
    """
    import calculateAge as age

    all_data = {'15776': {'SNID': 15776,
                          'SED':  np.array([23.14426, 21.00639, 19.41827,
                                            18.82437, 18.46391]),
                          'SEDerr': np.array([0.8009404, 0.04814964,
                                              0.01899451, 0.01771959,
                                              0.04554904]),
                          'redshift': 0.305
                         },
                '63': {'SNID': 63,
                       'SED': np.array([11.006, 11.076, 10.179, 8.925, 9.317]),
                       'SEDerr': np.array([0.002, 0.002, 0.001, 0.001, 0.002]),
                       'redshift': 0.001681
                      },
                '101': {'SNID': 101,
                        'SED':  np.array([13.671, 12.219, 11.529, 11.237,
                                          10.883]),
                        'SEDerr': np.array([0.004, 0.002, 0.002, 0.002,
                                            0.002]),
                        'redshift': 0.000804
                       },
               }
    try:
        data = all_data[cli['OBJID']]
    except KeyError as e:
        raise ValueError('The OBJID argument can only be 15776 or 101 at this'
                         ' time.')

    age.calculateSFH(data['SED'], data['SEDerr'], data['redshift'], data['SNID'],
                 burnin=True)


if __name__ == '__main__':
    # parse docopts
    cli = docopt(__doc__, version=__version__)

    # Setup logger
    ## initiate
    ## for all submodules to log with these settings they need to have the same start of a name.
    logger = logging.getLogger("fsps-age")
    ## set level
    if cli['--debug']:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)

    # create handler object to be able to change settings.
    ## choose saving location
    if cli['run']:
        fh = logging.FileHandler('logs/global/fsps-age_{}.log'.format(cli['JOBID']))
        formatLogging()
        redoGupta(cli)
    elif cli['burnin']:
        fh = logging.FileHandler('logs/burnin/fsps-age.log')
        formatLogging()
        burnin(cli)
    else:
        fh = logging.FileHandles('logs/fsps-age.log')
        formatLogging()

    logger.info('Done')

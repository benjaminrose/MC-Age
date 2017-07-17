"""burnin.py - 

Benjamin Rose
2017-06-22
python 3.5
"""
from calculateAge import 

module_logger = logging.getLogger("fsps-age.burnin")

def brunin(SED, SEDerr, redshift):
    """
    """
    # Set up logger
    logger = logging.getLogger("fsps-age.burnin.burnin")

    # call `calculateSFH with burnin
    # size should be (1000, 28, 7)
    samples = calculateSFH(SED, SEDerr, redshift, burin=False)

    #save data
    # Note header should be:
    # logzsol, dust2, tau, tStart, sfTrans, sfSlope, c
    # dex, 1/Gyr, Gyr, Gyr, , mag
    np.save('burnin', samples)
    logger.info('saved )
"""calculateAge.py -- This routine takes the *ugriz* and returns age PDF

Benjamin Rose
brose3@nd.edu
benjamin.rose@me.com
University of Notre Dame
2017-01-27
Python 3.6
"""
"""
Function outline, here is what each function calls or required data. Model parameters are **bold**. File is in inverse order of this.

- `calculateAge()` -- calculates mean mass weighted age
    - `calculateSFH()` -- calculates SFH parameters via MCMC
        - `lnprob()` -- sum of the ln(prior probability) and ln(likelihood)
            - `lnprior()` -- natural log of the prior probabilities (-inf or 0)
            - `lnlike()` -- calculates the likelihood of the model
                - `runFSPS()` -- the model
                    - redshift
                    - **star formation history model parameters**
                - `SED` in magnitudes
                - SED uncertainty (`SEDerr`) in magnitudes
                - **constant** that represents the observations seeing more then 1 solar mass of stars.
        - `lambda nll()` -- the negative of `lnlike()`, used for initial guess.
    - `starFormation()` -- the functional form of our SFH
"""
import logging
import platform
import warnings
import time

import numpy as np
from scipy import integrate
import fsps
import emcee

#used to write HDF5 tables, requires `h5py`
from astropy.table import Table
# from astropy.cosmology import WMAP9 as cosmo   # or make my own
from astropy.cosmology import FlatLambdaCDM

cosmo = FlatLambdaCDM(H0=70, Om0=0.27) ## Kessler 2009a but flat and 2-sig figs (Gupta 2011, §2.3)
# CampbellCosmo = FlatLambdaCDM(H0=73.8, Om0=0.24) 
module_logger = logging.getLogger("localEnvironments.calculateAge")

def runFSPS(sp, redshift, logzsol, dust2, tau, tStart, sfTrans, sfSlope):
    """Calculates expected SED for given stellar population (`sp`) and allows for variations in `redshift`, `logzsol`, `dust2`, `tau`, `tStart`, `sfTrans`, and `sfSlope`. It makes sense to use this function to minimize over metallicity, dust, and SFH. 
    
    Parameters
    ----------
    sp : fsps.fsps.StellarPopulation
        An FSPS StellarPopulation that we want to calculated the resulting SED.
    redshift : float
        The redshift of where this object is observed. 
    logzsol : float
        Parameter describing the metallicity, given in units of log(Z/Z⊙).
    dust2 : float
        Dust parameter describing the attenuation of old stellar light, i.e. where t > `dust_tesc`. `dust-test` is left at its default of 7.0 Gyr.
    tau : float
        The e-folding time for the SFH, in Gyr. The range is 0.1 < `tau` < 100.
    tStart : float
        Start time of the SFH, in Gyr.
    sfTrans : float
        Truncation time of the SFH, in Gyr. If set to 0.0, there is no truncation.
    sfSlope : float
        The slope of the SFR after time `sfTrans`.
    
    Returns
    -------
    list (or np.array ?)
        The output of `fsps.StellarPopulation.get_mags()` for the given redshift: the *ugriz* magnitudes of 1 solar mass. http://dan.iel.fm/python-fsps/current/stellarpop_api/#fsps.StellarPopulation.get_mags Currently *ugriz* are hard-coded in and not configurable. This could be changed. 

    """
    sdss_bands = ['sdss_u', 'sdss_g', 'sdss_r', 'sdss_i', 'sdss_z']

    #set metallicity
    sp.params['logzsol'] = logzsol

    #set Dust
    #young stars see dust1 & dust2. This produces young stars with 3x optical depth. 
    #using default `dust_tesc` of 7.0 Gyr
    # Conroy 2009 says Charlot 2000 recommends (dust1, dust2) = (1.0, 0.3)
    dust1 = 2.0*dust2
    sp.params['dust1'] = dust1
    sp.params['dust2'] = dust2

    #set SFH
    sp.params['tau'] = tau
    sp.params['sf_start'] = tStart
    sp.params['sf_trunc'] = sfTrans
    sp.params['sf_slope'] = sfSlope

    #calculate age of universe when observed light was emitted
    tage = cosmo.age(redshift).to('Gyr').value

    return sp.get_mags(tage=tage, redshift=redshift, bands=sdss_bands)


def lnlike(theta, magnitudes, magerr, redshift, sp):
    """
    log-likelyhood of the model to be minimized.: 
    ln p(y|x,σ,m,b,f) = −1/2 sum_n{(y_n−m*x_n−b)^2/s^2_n+ln(2π*s^2_n)}
    s^2_n = σ^2_n+f^2(m*x_n+b)^2

    Parameters
    ----------
    theta : list
        A list of the values of the current iteration of model parameters & a fitting constant to raise the mass from 1 solar mass (FSPS) to reality. Should be in order expected by `calculateAge.runFSPS()` with a constant on the end: `[logzsol, dust2, tau, tStart, sfTrans, sfSlope, c]`.
    magnitudes : list
        A list of the magnitudes of the stellar populations. Currently needs to be *ugriz* only.
    magerr : np.array
        A list of the associated uncertainties for `magnitudes`. Needs to be the same length.
    redshift : float
        The known redshift of the object
    sp : fsps.fsps.StellarPopulation
        The stellar population from which the age is derived

    See Also
    --------
    runFSPS : ? explains theta parameters. 
    """
    #calculate what the FSPS model is for the appropriate values in `theta`.
    model = runFSPS(sp, redshift, *theta[:-1])

    #test if `magnitudes` and `magerr` and `model` are all the same length
    if not len(magnitudes) == len(magerr) == len(model):
        raise ValueError('The inputs `magnitudes` and `magerr` need to be the same length as the result of `calculateAge.runFSPS()`. The lengths are {}, {}, and {} respectively.'.format(len(magnitudes), len(magerr), len(model)))

    #the inverse of the square of the uncertainty
    #aka the inverse of `magerr` squared. Just like in chi-square (before the summation)
    #section 4.2.6 of Statistics Data Mining ... (also look at 4.2.3)
    inv_sigma2 = 1.0/np.square(magerr)

    #todo(why does http://dan.iel.fm/emcee/current/user/line/#maximum-likelihood-estimation not have `np.log(2*np.pi*inv_sigma2)`)
    #todo(Do I need 1/2 at front or 5/2? https://www.statlect.com/fundamentals-of-statistics/normal-distribution-maximum-likelihood)
    #do these constants not matter because the `theta` that maximizes this function will be the same anyways and the value at maximum is not all that important? -- correct for 2*pi. Tested in `learningMCMC.py`.
    ln_like = -0.5*(np.sum((magnitudes-model-theta[-1])**2*inv_sigma2 
                   - np.log(2*np.pi*inv_sigma2))
                  )
    return ln_like


def lnprior(theta, redshift):
    """flat-top, log-priors for the parameter space to be searched: p(model)
    
    Parameters
    ----------
    theta : list
        Needs to be in order `[logzsol, dust2, tau, tStart, sfTrans, sfSlope]`. Variables are defined in runFSPS.
    
    Returns
    -------
    float
        Returns flat top ln-probability off by a constant: `0.0` if it is 
        within all prior bounds, negative infinity if at least one is outside 
        its prior's range.

    See Also
    --------
    runFSPS : ? explains theta parameters. 
    """
    # Conroy 2009 says Charlot 2000 recommends (dust1, dust2) = (1.0, 0.3)
    logzsol, dust2, tau, tStart, sfTrans, sfSlope, c = theta
    age = cosmo.age(redshift).to('Gyr').value
    # initial guess is `[0., 0., 1., 1., 10., 1., -20.]`
    # If `sfTrans` set to 0.0, there is no truncation.
    # should we allow sfTrans < tStart?
    if (-3.0  < logzsol < 0.5      and 
        0.0   < dust2   < 1.75     and 
        0.1   < tau     < 10.0     and 
        0.5   < tStart  < sfTrans  and 
        1.0   < sfTrans <= age     and
        -20.0 < sfSlope < 20.0     and 
        -35.0 < c       < -15.0):
        
        return 0.0
    
    return -np.inf


def lnprob(theta, magnitudes, magerr, redshift, sp):
    """
    log-probability of distributions (log-priors + log-likelyhood): 
    p(m,b,f|x,y,σ) ∝ p(m,b,f)*p(y|x,σ,m,b,f)

    Parameters
    ----------
    theta : array-like
        A list of the values of the current iteration of model parameters. 
        Should be in order expected by `calculateAge.runFSPS()`: [logzsol, 
        dust2, tau, tStart, sfTrans, sfSlope, c].
    magnitudes : array-like
        A list of the magnitudes of the stellar populations. Currently needs to be *ugriz* only.
    magerr : array-like
        A list of the associated uncertainties for `magnitudes`.
    redshift : float
        The known redshift of the object
    sp : fsps.fsps.StellarPopulation
        The stellar population from which ??
    
    Returns
    -------
    float, np.inf
        Returns negative infinity if theta is outside the prior range, otherwise returns the prior-potability times the likelihood-probability of the inputs. 

    See Also
    --------
    lnprior : 
    lnlike : 
    """
    lp = lnprior(theta, redshift)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, magnitudes, magerr, redshift, sp)


def calculateSFH(SED, SEDerr, redshift, SNID=None, sp=None, debug=False):
    """calculates the SHF. Save posterior distributions to ???.
    
    Parameters
    ----------
    SED : array-like
        The *ugriz* SED in magnitudes of the stellar population what we want 
        to fit a SFH to. 
    SEDerr : array-like
        Optional if `SED` is `False` and `x` is the SFH parameters.
    redshift : float
        The redshift of the region where the age is being calculated.
    SNID : int, optional
        This is the ID number for the SN. Used for 'unique' ID while saving 
        data. If nothing is given, MCMC chain is not saved.
    sp : fsps.fsps.StellarPopulation, optional
        An optional argument of an FSPS StellarPopulation if you don't want make a new one with the default settings.
    debug : bool
        Flag to have MCMC run incredibly short and in no way accurately. Also 
        does not save resulting chain. Should take around ~12 mins to get a 
        value of one SN.
    
    Returns
    -------
    ? : ? 
        PDF of age

    """
    #set up logger
    logger = logging.getLogger("localEnvironments.calculateAge.calculateSFH")
    logger.info('called calculateSFH')
    logger.debug('arguments are: {}, {}, {}, {}, {}'.format(SED, SEDerr, redshift, sp, SNID))

    #set up StellarPopulation if need be
    if sp is None:
        logger.debug('no stellar population argument passed')
        sp = fsps.StellarPopulation(zcontinuous=2, 
                  cloudy_dust=True, add_neb_emission = True,
                  sfh=5)
        logger.debug('no stellar population now created')

    #Setup MCMC
    logger.debug('initializing MCMC')
    if debug:
        ndim, nwalkers = 7, 14
        logger.debug('running a shorter MCMC run')
        nsteps = 20
        burnInSize = 10
        maxLikilhoodSize = 10
    else:
        ndim, nwalkers = 7, 100
        nsteps = 3000 #6000 for long run
        burnInSize = 400
        maxLikilhoodSize = 300
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(SED, SEDerr, redshift, sp))

    # "Maximize" likelihood for initial guess
    # start with uniformly (not 100% of space) distributed walkers
    # run MCMC because it has bounds
    # then select position with highest likelihood value.
    pos = np.zeros((nwalkers, ndim))
    '''PRIOR BOUNDS! -- with age = cosmo.age(redshift).to('Gyr').value
        -3.0  < logzsol < 0.5      and 
        0.0   < dust2   < 1.75     and 
        0.1   < tau     < 10.0     and 
        0.5   < tStart  < sfTrans  and 
        1.0   < sfTrans <= age     and
        -20.0 < sfSlope < 20.0     and 
        -35.0 < c       < -15.0'''
    pos[:,0] = np.random.uniform(-2.9, 0.4, size=nwalkers)   #logzsol
    pos[:,1] = np.random.uniform(0.1, 1.25, size=nwalkers)   #dust2
    pos[:,2] = np.random.uniform(0.5, 5.0, size=nwalkers)    #tau
    pos[:,3] = np.random.uniform(1.0, 5.0, size=nwalkers)    #tStart
    age = cosmo.age(redshift).to('Gyr').value
    pos[:,4] = np.random.uniform(1.2, age - 0.1, size=nwalkers)  #sfTrans
    pos[:,5] = np.random.uniform(-9.0, 9.0, size=nwalkers)   #sfSlope
    pos[:,6] = np.random.uniform(-30, -20, size=nwalkers)    #c

    print('Running MCMC for initial position')
    logger.info('Running MCMC for initial position with {} walkers, for {} steps'.format(nwalkers, maxLikilhoodSize))
    pos, prob, state = sampler.run_mcmc(pos, maxLikilhoodSize)
    
    #get position with "maximum" likelihood from limited run
    best_pos = sampler.flatchain[sampler.flatlnprobability.argmax()]
    sampler.reset()
    print('Best position from initial search: ', best_pos)
    logger.info('Best position from initial search: {}'.format(best_pos))
    
    #Set up new start position as a Gaussian ball around "max" likelihood
    pos = emcee.utils.sample_ball(best_pos, best_pos/1000., size=nwalkers)
    logger.debug('Sample ball returned: {}'.format(pos))

    #Run full MCMC analysis
    print('Running full MCMC fit')
    logger.info('Running with {} walkers, for {} steps, using a burn in cut after {} steps'.format(nwalkers, nsteps, burnInSize))
    logger.info('Running full MCMC fit')
    sampler.run_mcmc(pos, nsteps)
    logger.debug('Finished full MCMC fit')

    #save temp results in case something else fails soon.
    uuid = '{:.2f}'.format(time.time())
    header = 'logzsol\tdust2\ttau\ttStart\tsfTrans\tsfSlope\tc\ndex\t\t1/Gyr\tGyr\tGyr\t\tmag'
    np.savetxt('resources/temp/chain_'+uuid+'.tsv'.format(SNID), 
                sampler.flatchain, delimiter='\t', header=header)
    logger.info('saved resources/temp/chain_'+uuid+'.tsv')

    #save acceptance fraction
    #note() acceptace_fraction has len == nwalkers
    logger.info('Acceptance fraction: {}'.format(sampler.acceptance_fraction)) 
    print('Acceptance fraction: {}'.format(sampler.acceptance_fraction))

    ###################
    #plot burn in test
    # samples = sampler.chain  ##size == (nwalkers, nsteps, ndim)
    # lnprop_resutls = sampler.lnprobability  ##size == (nwalkers, nsteps)
    # import matplotlib.pyplot as plt
    # f, axarr = plt.subplots(8, sharex=True)
    # axarr[0].plot(samples[:,:,0].T)
    # axarr[0].set_title('Testing for Burn-in value')
    # axarr[0].set_ylabel('logzsol')
    # axarr[1].plot(samples[:,:,1].T)
    # axarr[1].set_ylabel('dust')
    # axarr[2].plot(samples[:,:,2].T)
    # axarr[2].set_ylabel('tau')
    # axarr[3].plot(samples[:,:,3].T)
    # axarr[3].set_ylabel('tStart')
    # axarr[4].plot(samples[:,:,4].T)
    # axarr[4].set_ylabel('sfTrans')
    # axarr[5].plot(samples[:,:,5].T)
    # axarr[5].set_ylabel('sfSlope')
    # axarr[6].plot(samples[:,:,6].T)
    # axarr[6].set_ylabel('c')
    # axarr[7].plot(lnprop_resutls.T)
    # axarr[7].set_ylabel('ln')
    # plt.savefig('figures/burnin.pdf')
    # logger.info('saved "figures/burnin.pdf"')
    # plt.show()
    ###################

    #Select only after "burn in" is complete
    #flatchain works if you run burn in specularly then full run
    # samples = sampler.flatchain      #size == (nsteps*nwalkers, ndim)
    #This method is a bit strange, but cuts the "burn in" section with ease
    samples = sampler.chain[:, burnInSize:, :].reshape((-1, ndim))
    
    #Save basic results to standard out & log
    logzso, dust2, tau, tStart, sfTrans, sfSlope, c = map(
                            lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                            zip(*np.percentile(samples, [16, 50, 84], axis=0))
                        )
    logger.info('MCMC values for logzsol, dust2, tau, tStart, sfTrans, sfSlope, c:')
    logger.info('{}, {}, {}, {}, {}, {}, {}'.format(logzso, dust2, tau, tStart,
                                                    sfTrans, sfSlope, c))
    print('MCMC: ', logzso, dust2, tau, tStart, sfTrans, sfSlope, c)

    #only save if SNID is given
    if SNID:
        logger.debug('SNID was given and now saving results')
        #only make a plot if on a Mac
        if platform.system() == 'Darwin':
            logger.debug('platform is Darwin, now making corner plot')
            import corner
            fig = corner.corner(samples, labels=["$logZsol$", "$dust_2$", 
                                        "$tau$", "$t_{start}$", "$t_{trans}$",
                                        '$sf slope$', 'c'])
            fig.savefig("figures/SF_triangle.pdf")
            logger.info('saved figures/SF_triangle.pdf')
        else:
            logger.debug('platform is {}, no corner plot made'.format(platform.system()))

        ## save clean results to disk
        header = 'logzsol\tdust2\ttau\ttStart\tsfTrans\tsfSlope\tc\ndex\t\t1/Gyr\tGyr\tGyr\t\tmag'
        np.savetxt('resources/SN{}_chain.tsv'.format(SNID), samples, 
                    delimiter='\t', header=header)
        logger.info('saved resources/SN{}_chain.tsv')
    else:
        logger.debug('SNID was {}, no data was saved to disk.'.format(SNID))
    
    logger.debug('done running calculateSFH')
    ## return flat clean chain
    return samples


def calculateAge(redshift, x, SEDerr=None, isSED=True, SNID=None, sp=None, debug=False):
    """calculateAge either from a given Star Formation History (SFH) or from a 
    *ugriz* SED. If a SED is given (the default) then it calculates the SFH by 
    calling `calculateSFH()`.
    
    Parameters
    ----------
    x : array-like
        Either an SED of *ugriz* magnitudes or SFH parameters (`tau`, 
        `tStart`, `sfTrans`, `sfSlope` of Shimha 2014 equation 6) if `SED` is 
        set to `False`.
    redshift : float
        The redshift of the region where the age is being calculated.
    SEDerr : array-like, optional
        Optional if `SED` is `False` and `x` is the SFH parameters. 
    isSED : boolean, optional
        A flag to determine if x is an SED (default value) or star formation 
        history parameters. 
    SNID : int, optional
        This is the ID number for the SN. Used for 'unique' ID while saving 
        data. Needed if `SED = True`.
    sp : fsps.fsps.StellarPopulation, optional
        An optional argument of an FSPS StellarPopulation if you don't want 
        make a new one with the default settings of `calculateSFH()`.
    debug : bool
        Flag to have MCMC run incredibly short and in no way accurately. Also 
        does not save resulting chain. Should take around ~12 mins to get a 
        value of one SN.
    
    Returns
    -------
    sfh : np.array
        The mean star formation history parameters (tau, tStart, sfTrans, 
        sfSlope). Either calculated from the SED given in `x` or the 
        parameters given in `x` if `SED = False`. A full posterior 
        distribution is save by `calculateSFH()`
    age : float
        The mean mass-weighted(?) age from the given or calculated SFH
        parameters. If the SFH was fit, a full posterior distribution is save 
        at `???`.
    
    Raises
    ------
    LinAlgException
        If the matrix is not numerically invertible.

    See Also
    --------
    calculateSFH : A way to calculate a Star Formation History from an SED
    
    Notes
    -----
    
    Examples
    --------
    
    """
    logger = logging.getLogger("localEnvironments.calculateAge.calculateAge")
    logger.info('called calculateAge')

    #if SED, calculate SFH
    if isSED:
        logger.info('Calculating SFH')
        # no need for SNID. We can save in this function.
        samples = calculateSFH(x, SEDerr, redshift, debug=debug)
        # logger.info('importing SFH to speed!!!!!')
        # samples = np.genfromtxt('resources/SN0_chain.tsv', delimiter='\t')
        #extract variables
        logzso, dust2, tau, tStart, sfTrans, sfSlope, c = np.hsplit(samples, 7)
        #reshape these to be 1D arrays
        logzso = logzso.reshape(-1)
        dust2 = dust2.reshape(-1)
        tau = tau.reshape(-1)
        tStart = tStart.reshape(-1)
        sfTrans = sfTrans.reshape(-1)
        sfSlope = sfSlope.reshape(-1)
        c = c.reshape(-1)
    else:
        logger.info('using SFH that was given')
        #unpack input parameters
        #todo(4 should not be hard coded)
        tau, tStart, sfTrans, sfSlope = np.hsplit(x, 4)

    #with SFH -> Calculate age
    '''
    How do I maintain the information of the posterior distribution?
    each column of from the MCMC sample is a "possible" SFH. So we should 
    just calculate the expected age for each index. This gives us the same 
    number of ages as samples and gives a us a distribution for the age.
    '''

    ## Calculate Age!
    ageOfUniverse = cosmo.age(redshift)
    # lengthOfSF = ageOfUniverse.to('Gyr').value - tStart
    numerator = np.array([])
    denominator = np.array([])
    for j, k, l, m in zip(tau, tStart, sfTrans, sfSlope):
        #only need the value of `integrate.quad` not the absolute error
        numerator = np.append(numerator, integrate.quad(tStarFormation, k, ageOfUniverse.to('Gyr').value, args=(j, k, l, m))[0])
        denominator = np.append(denominator, integrate.quad(starFormation, k, ageOfUniverse.to('Gyr').value, args=(j, k, l, m))[0])
        if denominator[-1] == 0:
            logger.info('because of scipy issue need to use sample integration')
            logger.warning('''SFH: {}, {}, {}, {} 
                (emitted at z={}, ageOfUniverse={}) 
                for SN{} produced a zero integrated SFH in the age calculation.'''.format(j, k, l, m, redshift, ageOfUniverse.to('Gyr').value, SNID))
            warnings.warn('Getting zero integrated SFH, check log.')
            time, dx = np.linspace(k, ageOfUniverse.to('Gyr').value, num=8193,
                                    retstep=True)
            num_y, den_y = np.array([]), np.array([])
            for n in time:
                num_y = np.append(num_y, tStarFormation(n, j, k, l, m))
                den_y = np.append(den_y, starFormation(n, j, k, l, m))
            numerator[-1] = integrate.romb(num_y, dx)
            denominator[-1] = integrate.romb(den_y, dx)

    # This is from Gupta 2011 Equation 3 but with a change in t_0. He used t_0
    # = start of star formation, I use t_0 = big bang. Explained in FIndings 
    # on 2017-05-10.
    age = ageOfUniverse.to('Gyr').value - numerator/denominator

    #Warn if any age calculations produce NaN
    if np.isnan(age).any():
        logger.warning("part of the age calculations for SN{} are not a number".format(SNID))
    # get median +/- 1 sigma like from MCMC example. 
    #Reshape to (1,3) so all 3 values are passed to lambda over 1 iteration.
    #ignore NaNs. They should be such a small percent and will be followed up
    age_precentiels = list(map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), 
                               np.nanpercentile(age, [16, 50, 84]).reshape(1,3)))
    logger.info('Age is: '+ str(age_precentiels))

    #only save data/figure if running full run.
    if debug:
        logger.info('skipping saving of data while debugging')
    else:
        # if MCMC was ran --> save full MCMC results with age
        if isSED:
            logger.debug('saving full MCMC & age data')
            samples = np.append(samples, age.reshape(age.size, 1), axis=-1)

            header = 'logzsol\tdust2\ttau\ttStart\tsfTrans\tsfSlope\tc\tage\ndex\t\t1/Gyr\tGyr\tGyr\t\tmag\tGyr'

            np.savetxt('resources/SN{}_chain.tsv'.format(SNID), samples, 
                    delimiter='\t', header=header)
            logger.info('saved resources/SN{}_chain.tsv'.format(SNID))

            #only on a Mac, make corner plot of MCMC parameters & Age
            if platform.system() == 'Darwin':
                import corner
                fig = corner.corner(samples, labels=["$logZsol$", "$dust_2$", 
                                    "$tau$", "$t_{start}$", "$t_{trans}$", 
                                    '$sf slope$', 'c', 'age'])
                fig.savefig("figures/SF_Age_triangle.pdf")

        print(np.nanmedian(age))
        print(np.median(age))
    
    return np.nanmedian(age)

##################
#helper functions for Star Formation History
#All these functions use t_0 = Big Band. Gupta 2011 Age calculation 
#assumes t_0 is the start of star formation. Details in Findings 2017-05-10
##################
#http://stackoverflow.com/questions/32877587/rampfunction-code
ramp = lambda x: x if x >= 0 else 0 #Simha14's Delta function (eqn 6)

heavyside = lambda x: 1 if x >=0 else 0

def starFormation(t, tau, tStart, sfTrans, sfSlope):
    '''
    Defined piecewise around `sfTrans` time. 
    no variables can be arrays! ALso must return a float!
    '''
    if t <= sfTrans:
        sf = heavyside(t-tStart)*(t-tStart)*np.e**(-(t-tStart)/tau)
    else:
        sf = ( starFormation(sfTrans, tau, tStart, sfTrans, sfSlope) +
             sfSlope*ramp(t-sfTrans) ) #\Delta \def sfSlope*ramp() in Simha

    #only return a positive star formation
    if sf <= 0:
        return 0
    else:
        return sf

# This needs to multiply star formation by t-tStart for integral in numerator of
# age calculation. Details on why it is different then Gupta 2011 equation 3 is
# in Findings "Can I redo Gupta's results with our improved SFH?" on 2017-05-10.
#Note: args[1] needs to be `tStart` to work when passed to `starFormation()`
tStarFormation = lambda t, *args: (t-args[1])*starFormation(t, *args)
##################
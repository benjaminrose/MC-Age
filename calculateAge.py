"""calculateAge.py -- This routine takes the SED photometry and computes an 
age probability distribution

Benjamin Rose
brose3@nd.edu
benjamin.rose@me.com
University of Notre Dame
2017-01-27
Python 3.5
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

import numpy as np
from scipy import integrate
import fsps
import emcee
import astropy.units as u
# from astropy.cosmology import WMAP9 as cosmo   # a generic cosmology
from astropy.cosmology import FlatLambdaCDM

# Kessler 2009a but flat and 2-sig figs (Gupta 2011, §2.3)
cosmo = FlatLambdaCDM(H0=70, Om0=0.27)
# CampbellCosmo = FlatLambdaCDM(H0=73.8, Om0=0.24)

# The part before the "." needs to be the same name as the initial logger in
# fsps-age.py if we want to use its settings.
module_logger = logging.getLogger("fsps-age.calculateAge")

def runFSPS(sp, redshift, logzsol, dust2, tau, tStart, sfTrans, sfSlope):
    """Calculates expected SED for given stellar population (`sp`) and allows
    for variations in `redshift`, `logzsol`, `dust2`, `tau`, `tStart`,
    `sfTrans`, and `sfSlope`. Returns the apparent magnitude in *ugriz*.

    Note that the model produces 1 solar mass of stars, so a scaling factor
    may need to be applied.
    
    Parameters
    ----------
    sp : fsps.fsps.StellarPopulation
        An FSPS StellarPopulation that we want to calculated the resulting SED.
        Assumes `zcontinuous` and `sfh` are set to $> 0$ and $5$ respectively.
    redshift : float
        The redshift of the object is observed.
    logzsol : float
        Parameter describing the metallicity, given in units of log(Z/Z⊙).
        python-FSPS's `logzsol`. Not sure what FSPS parameter it is.
    dust2 : float
        Dust parameter describing the attenuation of old stellar light, i.e.
        where t > `dust_tesc`. `dust-test` is left at its default of 7.0 Gyr.
        FSPS's `dust2`.
    tau : float
        The e-folding time for the SFH, in Gyr. The range is 0.1 < `tau` < 100.
        FSPS's `tau`.
    tStart : float
        Start time of the SFH, in Gyr. FSPS's `sf_start`.
    sfTrans : float
        Transition time of the SFH from exponential to linear, in Gyr.
        FSPS's `sf_trunc`.
    sfSlope : float
        The slope of the linear SFR after time `sfTrans`.
    
    Returns
    -------
    np.ndarray
        A 1D array of the output of `fsps.StellarPopulation.get_mags()` for
        the given redshift: the *ugriz* magnitudes of 1 solar mass.
        dan.iel.fm/python-fsps/current/stellarpop_api/#fsps.StellarPopulation.get_mags
    """
    sdss_bands = ['sdss_u', 'sdss_g', 'sdss_r', 'sdss_i', 'sdss_z']

    # set metallicity
    sp.params['logzsol'] = logzsol

    # set dust
    # using default `dust_tesc` of 7.0 Gyr
    # young stars see dust1 & dust2. 
    # Conroy 2009 says Charlot 2000 recommends (dust1, dust2) = (1.0, 0.3)
    # This produces young stars with 3x optical depth. 
    dust1 = 2.0*dust2
    sp.params['dust1'] = dust1
    sp.params['dust2'] = dust2

    # set SF parameters
    sp.params['tau'] = tau
    sp.params['sf_start'] = tStart
    sp.params['sf_trunc'] = sfTrans
    sp.params['sf_slope'] = sfSlope

    # calculate age of universe when observed light was emitted
    tage = cosmo.age(redshift).to('Gyr').value

    # * `zmet` has no effect, 2017-08-15 lab notes
    # * since we don't use `add_igm_absortion` there is no difference between
    # `sp.get_mags(redshfit) and `sp.params['redshift']`
    # * `tage` is the age of the stellar population. This uses the same t_0
    # reference as the `sf_start`, `st_trunc`, etc., aka since the bing bang
    return sp.get_mags(tage=tage, redshift=redshift, bands=sdss_bands)


def lnlike(theta, magnitudes, magerr, redshift, sp):
    """
    log-likelyhood of the model to be minimized.
    Example of fitting a line with undersized errors: 
    ln p(y|x,σ,m,b,f) = −1/2 sum_n{(y_n−m*x_n−b)^2/s^2_n+ln(2π*s^2_n)}
    s^2_n = σ^2_n+f^2(m*x_n+b)^2

    Parameters
    ----------
    theta : array-like
        A list (or an array) of the values of the current iteration of model
        parameters & a fitting constant to raise the mass from 1 solar mass
        (FSPS) to reality. Should be in order expected by
        `calculateAge.runFSPS()` with a constant on the end:
        `[logzsol, dust2, tau, tStart, sfTrans, phi, c]` except
        `sfSlope` = tan(`phi`).
    magnitudes : np.ndarray
        A list of the magnitudes of the stellar populations. Currently needs
        to be *ugriz* only.
    magerr : np.ndarray
        A list of the associated uncertainties for `magnitudes`. Needs to be
        the same length.
    redshift : float
        The known redshift of the object
    sp : fsps.StellarPopulation
        The stellar population from which the age is derived

    Returns
    -------
    float
        The likelihood of the given model parameters (`theta`, `sp`) relative
        to the given data (`magnitudes`, `magerr`, `redshift`). 

        If `-np.inf` is returned that means either the model produced at least
        one filter more then 0.5 mags away from the observed value or the
        likelihood was less then -10,000. This is because (as of 2017-08-15)
        MCMC was selecting terrible fits.

    Raises
    ------
    ValueError
        if `theta` does not contain the expected 7 model parameters
    ValueError
        if `magnitudes` and `magerr` are not length 5, one value per SDSS
        filter

    See Also
    --------
    runFSPS : explains `theta` parameters, except `sfSlope` = tan(`phi`)
    lnprior : the prior probability
    lnprob : the posterior probability 
    """
    # test inputs
    if not len(theta) == 7:
        raise ValueError('Likelihood expects 7 model variables')

    if not len(magnitudes) == len(magerr) == 5:
        raise ValueError('The inputs `magnitudes` and `magerr` need to be' +
                         ' arrays of length 5. The lengths are {} and {}'
                         .format(len(magnitudes), len(magerr)) +
                         ' respectively.')

    # calculate what the FSPS model is for the appropriate values in `theta`.
    model = runFSPS(sp, redshift, *theta[:-2], np.tan(theta[-2])) + theta[-1]

    # test if `magnitudes` and `magerr` and `model` are all the same length
    # this verifies that `runFSPS()` has not changed what filters it is using.
    if not len(magnitudes) == len(magerr) == len(model):
        raise ValueError('The inputs `magnitudes` and `magerr` need to be' +
                         ' the same length as the result of' +
                         ' `calculateAge.runFSPS()`. The lengths are {}, {},'
                         .format(len(magnitudes), len(magerr)) + 
                         ' and {} respectively.'.format(len(model)))

    # the inverse of the square of the uncertainty
    # aka the inverse of `magerr` squared. Just like in chi-square (before the summation)
    # section 4.2.6 of Statistics Data Mining ... (also look at 4.2.3)
    inv_sigma2 = 1.0/np.square(magerr)

    #todo(why does http://dan.iel.fm/emcee/current/user/line/#maximum-likelihood-estimation not have `np.log(2*np.pi*inv_sigma2)`)
    #todo(Do I need 1/2 at front or 5/2? https://www.statlect.com/fundamentals-of-statistics/normal-distribution-maximum-likelihood)
    #do these constants not matter because the `theta` that maximizes this function will be the same anyways and the value at maximum is not all that important? -- correct for 2*pi. Tested in `learningMCMC.py`.
    lnlike_model = np.square(magnitudes-model)*inv_sigma2
    lnlike_error = -np.log(2*np.pi*inv_sigma2)
    assert len(lnlike_model) == 5, 'should still be 5 filters'
    assert len(lnlike_error) == 5, 'should still be 5 filters'
    lnlike_sum = np.sum(lnlike_model + lnlike_error)
    ln_like = -0.5*lnlike_sum

    return ln_like


def lnprior(theta, redshift):
    """flat-top, log-priors for the parameter space to be searched: p(model)
    
    Parameters
    ----------
    theta : list
        Needs to be in order `[logzsol, dust2, tau, tStart, sfTrans, phi]`.
        Variables are defined in runFSPS, except `sfSlope` = tan(`phi`).
    
    Returns
    -------
    float
        Returns flat top ln-probability off by a constant: `0.0` if it is 
        within all prior bounds, negative infinity if at least one is outside 
        its prior's range.

    See Also
    --------
    runFSPS : explains theta parameters, except `sfSlope` = tan(`phi`).
    """
    # Conroy 2009 says Charlot 2000 recommends (dust1, dust2) = (1.0, 0.3)
    logzsol, dust2, tau, tStart, sfTrans, phi, c = theta
    age = cosmo.age(redshift).to('Gyr').value
    # initial guess is `[0., 0., 1., 1., 10., 1., -20.]`
    # If `sfTrans` set to 0.0, there is no truncation.
    # should we allow sfTrans < tStart?

    # Initially set flat priors, except for variables to follow if-statement
    if (-2.5  < logzsol < 0.5           and
        0.0   <= dust2 <= 0.9           and
        0.1   < tau     < 10.0          and
        0.5   < tStart  < sfTrans - 2.0 and   #force at least 2 Gyr of tau
        2.5   < sfTrans <= age          and
        -1.520838 < phi < 1.520838      and    # from -20 to 20 in slope
        -45.0 < c       < -5.0):
        
        # logzsol #
        # http://www.nature.com/nature/journal/v534/n7608/full/nature18322.html
        # Mean metallicity vs redshift (figure 6)
        # not much change over 0 < redshift < 0.5
        # Simply, this distribution can be characterized as a normal
        # distribution (between logzsol & redshift) centered at -0.5, and
        # sigma = 0.5 dex
        CENTER_Z = -0.5
        SIGMA_Z = 0.6   #cover -1.5 in 2sigma

        # dust2 #
        #high dust (dust2>1) is (starburst) is ~1% of galaxies
        #also Conroy 2009 2.6 says 0.3 is fine for most.
        sigma = 0.3
        center = 0.0
        #return ln-prior of dust2. Centered at 0 with sigma from above. 
        #Note that this return only takes place if `dust2`>0, so this is only
        #the right and side of the Gaussian. 

        return -1*((center-dust2)**2/(2*sigma**2) +
                   np.log(np.sqrt(2*np.pi)*sigma) +
                   (CENTER_Z - logzsol)**2/(2*SIGMA_Z**2) +
                   np.log(np.sqrt(2*np.pi)*SIGMA_Z))
    
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
        dust2, tau, tStart, sfTrans, phi, c], except `sfSlope` = tan(`phi`).
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

def setUpMCMC(ndim, nwalkers, maxLikilhoodSize, sampler, redshift):
    """
    mostly-maximize likelihood to find initial MCMC position

    Parameters
    ----------
    ndim, nwalkers, maxLikilhoodSize, sampler
    """
    logger = logging.getLogger("fsps-age.calculateAge.setUpMCMC")
    logger.info('called setupMCMC')

    ############
    # mostly-maximize likelihood to find initial MCMC position
    ############
    # Use MCMC because it can simply force boundary conditions
    # Set up initial array
    pos = np.zeros((nwalkers, ndim))
    # PRIOR BOUNDS! -- with age = cosmo.age(redshift).to('Gyr').value
    #     -2.5  < logzsol < 0.5           and
    #     0.0   < dust2                   and
    #     0.1   < tau     < 10.0          and
    #     0.5   < tStart  < sfTrans - 2.0 and
    #     2.5   < sfTrans <= age          and
    #     -1.520838 < phi < 1.520838         and
    #     -45.0 < c       < -5.0):
    # Keep
    # Assuming you have enough walkers, These random numbers will "uniformly
    # fill" the 7-d parameter space. More walkers (or a smaller initial search
    # space) will make this assumption more valid.
    # Metallicity should not bee too metal poor. Minimize search space by not
    # removing the some of the edge. Also z=0.1 objects unlikely to be more
    # metal rich then the Sun.
    pos[:,0] = np.random.uniform(-1.0, 0.0, size=nwalkers)   # logzsol to ~1σ
    pos[:,1] = np.random.uniform(0.1, 0.4, size=nwalkers)    # dust2 to ~1σ
    pos[:,2] = np.random.uniform(0.5, 8.0, size=nwalkers)    # tau
    pos[:,3] = np.random.uniform(1.0, 5.0, size=nwalkers)    # tStart
    age = cosmo.age(redshift).to('Gyr').value
    pos[:,4] = np.random.uniform(3.0, age - 0.1, size=nwalkers)   # sfTrans
    # Keep this one a bit tight because it should have a diminishing influence
    # as it approaches it bounds.
    # It does have a MAJOR influence if it is large. This is HOW we get young
    # populations.
    # from -8 to +19 in slope
    pos[:,5] = np.random.uniform(-1.4464, 1.5182, size=nwalkers)   # phi
    # Keep values tight if expected value is well understood
    # Also this variable is highly influential to the likelihood
    pos[:,6] = np.random.uniform(-28, -15.0, size=nwalkers)    # c

    print('Running MCMC for initial position')
    logger.info('Running MCMC for initial position with {} walkers, for {} steps'.format(nwalkers, maxLikilhoodSize))
    pos, prob, state = sampler.run_mcmc(pos, maxLikilhoodSize)
    
    # get position with "maximum" likelihood from limited run
    best_pos = sampler.flatchain[sampler.flatlnprobability.argmax()]
    how_likely = sampler.flatlnprobability.max()
    print('Best position from initial search: ', best_pos)
    print('with a ln-likelihood of: ', how_likely)
    logger.info('Best position from initial search: {}'.format(best_pos))
    logger.info('with a ln-likelihood of: ', how_likely)
    logger.debug("Mean ln-probability for each walker: {}".format(
                sampler.lnprobability.mean(-1)))
    logger.debug("Max ln-probability for each walker: {}".format(
                sampler.lnprobability.max(-1)))
    sampler.reset()
    logger.debug("Reset sampler's chain and lnprobability arrays")

    return best_pos

# def runMCMC(ndim, nwalkers = 7, 300
#         burnInSize = 200
#         nsteps = 1000
#         sampler):
#     return

def calculateSFH(SED, SEDerr, redshift, SNID=None, sp=None, debug=False, 
                 burnin=False):
    """Calculates the SFH. Returns the cleaned chain.
    
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
        An optional argument of an FSPS StellarPopulation if you don't want 
        make a new one with the default settings of this analysis.
    debug : bool
        Flag to have MCMC run incredibly short and in no way accurately. Also 
        does not save resulting chain. Should take around ~12 mins to get a 
        value of one SN.
    burnin : bool
        Should MCMC run with a smaller number of walkers, smaller search for initial Maximum Likelihood location, and not throw away any 'burn-in' steps? This is used to test and validate what burn-in size should be used. This argument does nothing if `debug` is `True`.
    
    Returns
    -------
    samples : numpy.array 
        The resulting chain of the MCMC analysis. The shape will be (`nsteps`*`nwalkers`, `ndim`)

    or

    None
        Returns nothing if 

    """
    #set up logger
    logger = logging.getLogger("fsps-age.calculateAge.calculateSFH")
    logger.info('called calculateSFH')
    logger.debug('arguments are: {}, {}, {}, {}, {}'.format(SED, SEDerr, redshift, sp, SNID))

    #set up StellarPopulation if need be
    if sp is None:
        logger.debug('no stellar population argument passed')
        sp = fsps.StellarPopulation(zcontinuous=2, 
                  cloudy_dust=True, add_neb_emission = True,
                  sfh=5)
        logger.debug('stellar population now created')


    ############
    # Setup MCMC
    ############
    logger.debug('initializing MCMC parameters')
    if debug:
        # This is minimum value to test if the code runs. Resulting values will
        # be meaningless.
        ndim, nwalkers = 7, 14
        maxLikilhoodSize = 3
        burnInSize = 3
        nsteps = 6
    elif burnin:
        # This is a smaller test to see how the sampling is progressing.
        # Hopefully these values will not effect how it samples, but it will
        # likely effect how well it sampled.
        ndim, nwalkers = 7, 100
        maxLikilhoodSize = 10
        burnInSize = 0
        nsteps = 300
    else:
        # Gupta's grid search was over 2592 points (in 4D, roughly 7 points
        # per dimension).
        # MCMC should walk at least that many or I do not improve on his work.
        # Initial course search (via "maximum-likelihood") could be an order 
        # of magnitude less. So `nwalkers` ~ 200 and the `maxLikeilhoodSize`
        # would say how well that area is searched.
        # To search 7 dimensions 200 walkers would sample each dimension 2.1
        # times (200**(1/7)). If we wanted 4 samples per dimension we would
        # need 16,384 walkers! 3 samples per dimension would be 2,187.
        # `emcee`'s "stretch move" will fill in the gaps by pulling all the low likelihood walkers harder to the high likelihood zones. We just need to run this section long enough for a low likelihood walker to get paired with a high likelihood walker, ~N the number of walkers.
        ndim, nwalkers = 7, 300
        maxLikilhoodSize = 500    # increase from 300 just cause?
        burnInSize = 500          # 500 seems to be good.
        nsteps = 1350         
        # We want ~250,000 accepted values: nwalkers*(nsteps - burninSize).
        # We want this many because it looks good?
        # 750,000 seemed to be too much. There was no change in the
        # distributions from the first, second, third quarter million or the
        # whole distribution.
    logger.info('running {} walkers,\n\t\t {} initial search steps,\n\t\t {} final steps,\n\t\twith {} burnin steps remove'.format(nwalkers, maxLikilhoodSize, nsteps, burnInSize))

    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(SED, SEDerr, redshift, sp))

    best_pos = setUpMCMC(ndim, nwalkers, maxLikilhoodSize, sampler, redshift)

    ############
    # Set up new start position as a Gaussian ball around "max" likelihood.
    ############
    # This starts the run in a high-probability zone, reducing the need to cut
    # any burnin.
    # Make ball a gaussian of "1/2" expected sigma for each parameter
    # needs to one for every dimension -- This starts all the walkers in a
    # position that resembles the expected distribution. Reducing the need for
    # burnin to be cut out.
    spread = [0.1, 0.05, 0.1, 0.05, 0.05, 1.0, 0.1]
    # spread = [0.001, 0.0005, 0.001, 0.0005, 0.0005, 0.01, 0.001]   # is the result very small?
    pos = emcee.utils.sample_ball(best_pos, spread, size=nwalkers)
    logger.debug('Sample ball returned: {}'.format(pos))


    ############
    # Run full MCMC analysis
    ############
    print('Running full MCMC fit')
    logger.info('Running full MCMC with {} walkers, for {} steps, using a burn in cut after {} steps'.format(nwalkers, nsteps, burnInSize))
    sampler.run_mcmc(pos, nsteps)
    logger.debug('Finished full MCMC fit')

    ############
    # save temp results in case something else fails soon.
    ############
    header = 'logzsol\tdust2\ttau\ttStart\tsfTrans\tsfAngle\tc\ndex\t\t1/Gyr\tGyr\tGyr\t\tmag'
    # convert from phi to slope
    to_save = sampler.flatchain
    np.savetxt('resources/SN{}_chain.tsv'.format(SNID), sampler.flatchain, 
                delimiter='\t', header=header)
    logger.info('saved full results resources/SN{}_chain.tsv'.format(SNID))

    #save acceptance fraction & ln-probability
    logger.info("Mean ln-probability for each walker: {}".format(
                sampler.lnprobability.mean(-1)))
    logger.info("Max ln-probability for each walker: {}".format(
                sampler.lnprobability.max(-1)))
    #note() acceptace_fraction has len == nwalkers
    logger.info('Acceptance fraction: {}'.format(sampler.acceptance_fraction)) 
    print('Acceptance fraction: {}'.format(sampler.acceptance_fraction))

    #Select only after "burn in" is complete
    #flatchain works if you run burn in specularly then full run
    # samples = sampler.flatchain      #size == (nsteps*nwalkers, ndim)
    allSamples = sampler.chain             # size == (nwalker, nsteps, ndim)
    allLnProb = sampler.lnprobability     # size == (nwalkers, nsteps)

    # Note header should be:
    # logzsol, dust2, tau, tStart, sfTrans, sfAngle, c
    # dex, 1/Gyr, Gyr, Gyr, , mag
    np.save('resources/burnin/SN{}_samples'.format(SNID), allSamples)
    np.save('resources/burnin/SN{}_lnprob'.format(SNID), allLnProb)
    logger.debug('saved samples and ln probability')

    if burnin:
        # burnin test now done
        return None
         
    #This method is a bit strange, but cuts the "burn in" section with ease
    # And makes it a flat chain
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

    ## save trimmed results to disk
    header = 'logzsol\tdust2\ttau\ttStart\tsfTrans\tsfAngle\tc\ndex\t\t1/Gyr\tGyr\tGyr\t\tmag'
    np.savetxt('resources/SN{}_chain.tsv'.format(SNID), samples, 
                delimiter='\t', header=header)
    logger.info('saved resources/SN{}_chain.tsv'.format(SNID))
    
    logger.debug('done running calculateSFH')

    ## return flat clean chain
    return samples


def integrate_age(tau, tStart, sfTrans, sfSlope, redshift):
    """Integrates over star formation history to get the age mass weighted age
    of a stellar history. This Integrates over a star formation history defined
    in Rose & Garnavich 2017, the same as FSPS with `sfh` set to 5.

    Parameters
    ----------
    tau : array_like, real numbers
        The e-folding scale ($e^(t/tau)$) of the initial star formation burst.
    tStart : array_like, real numbers
        The time, since the big bang, when star formation start.
    sfTrans : array_like, real numbers
        The time, since the bing bang, when star formation transfers from an
        exponential decay to a linear rate.
    sfSlope : array_like, real numbers
        The is the slope of the late time star formation linear rate.
    redshift : float
        This observed redshift. This is used to calculate when the light was
        emitted. To convert from redshift to age of the universe a cosmology
        is assumed. Currently that is Kessler 2009a (also Gupta 2011, §2.3),
        but it could change to Campbell 2013.

    Returns
    -------
    np.ndarray
        The age for those particular star formation parameters. Will be the
        same length as the array_like SFH parameters.

    Notes
    -----
    - This will fail if passed an np.ndarray of a float/int, such as
    np.array(1). It will pass the input tests, but the zip (for the for-loop)
    will fail.
    """
    logger = logging.getLogger("fsps-age.calculateAge.integrate_age")
    logger.info('called integrate_age')

    #######
    # Test inputs
    #######

    # Test type & convert to np.ndarray
    if isinstance(tau, (int, float)):
        tau = np.asarray([tau])
    elif isinstance(tau, (np.ndarray, list, tuple)):
        tau = np.asarray(tau)
    elif not isinstance(tau, (np.ndarray)):
        raise TypeError('tau must be array_like')

    if isinstance(tStart, (int, float)):
        tStart = np.asarray([tStart])
    elif isinstance(tStart, (np.ndarray, list, tuple)):
        tStart = np.asarray(tStart)
    elif not isinstance(tStart, (np.ndarray)):
        raise TypeError('tStart must be array_like')
    
    if isinstance(sfTrans, (int, float)):
        sfTrans = np.asarray([sfTrans])
    elif isinstance(sfTrans, (np.ndarray, list, tuple)):
        sfTrans = np.asarray(sfTrans)
    elif not isinstance(sfTrans, (np.ndarray)):
        raise TypeError('sfTrans must be array_like')
    
    if isinstance(sfSlope, (int, float)):
        sfSlope = np.asarray([sfSlope])
    elif isinstance(sfSlope, (np.ndarray, list, tuple)):
        sfSlope = np.asarray(sfSlope)
    elif not isinstance(sfSlope, (np.ndarray)):
        raise TypeError('sfSlope must be array_like')

    # make sure objects are 0 or 1 dimensional
    if len(tau.shape) > 1:
        raise TypeError('tau cannot be more then 1 dimensional.')
    if len(tStart.shape) > 1:
        raise TypeError('tStart cannot be more then 1 dimensional.')
    if len(sfTrans.shape) > 1:
        raise TypeError('sfTrans cannot be more then 1 dimensional.')
    if len(sfSlope.shape) > 1:
        raise TypeError('sfSlope cannot be more then 1 dimensional.')

    # use size over len() to work with arrays of floats
    # but now ints and floats become arrays of a list
    if not tau.size == tStart.size == sfTrans.size == sfSlope.size:
        raise ValueError("The star formation parameters need to be equal lengths")

    # test redshift
    if not isinstance(redshift, (float)):
        raise TypeError('redshift must be a float')
    if redshift <= 0:
        raise ValueError("Redshifts should be greater than zero.")

    # Get correct variables
    ageOfUniverse = cosmo.age(redshift)
    # this is now an np.ndarray
    lengthOfSF = ageOfUniverse - tStart*u.Gyr     # Gupta's A value,
    sfTrans_gupta = sfTrans - tStart       # transition time, since start of SF


    # Calculate numerator and denominator for each SFH
    numerator = np.array([])
    denominator = np.array([])
    for j, l, m, a in zip(tau, sfTrans_gupta, sfSlope, 
                         lengthOfSF.to('Gyr').value):
        #only need the value of `integrate.quad` not the absolute error
        numerator = np.append(numerator, integrate.quad(t_star_formation_gupta, 0, a, args=(j, l, m))[0])
        denominator = np.append(denominator, integrate.quad(star_formation_gupta, 0, a, args=(j, l, m))[0])

        # if SF is too much of a burst, compute by hand.
        if denominator[-1] == 0:
            logger.info('because of scipy issue need to use sample integration')
            logger.warning('''SFH: {}, {}, {} 
                (emitted at z={}, lengthOfSF={}) 
                produced a zero integrated SFH in the age calculation.'''.format(j, l, m, redshift, a))
            time_, dx = np.linspace(0, a, num=8193,
                                    retstep=True)
            num_y, den_y = np.array([]), np.array([])
            for n in time_:
                num_y = np.append(num_y, t_star_formation_gupta(n, j, l, m))
                den_y = np.append(den_y, star_formation_gupta(n, j, l, m))
            numerator[-1] = integrate.romb(num_y, dx)
            denominator[-1] = integrate.romb(den_y, dx)

    # This is from Gupta 2011 Equation 3 but with a change in t_0. He used t_0
    # = start of star formation, I use t_0 = big bang. Explained in FIndings 
    # on 2017-05-10 & updated on 2017-05-15
    # return ageOfUniverse.to('Gyr').value - tStart + numerator/denominator

    # return a function that integrates over Gupta's variables.
    return lengthOfSF.to('Gyr').value - numerator/denominator

def calculateAge(redshift, SED, SEDerr, SNID, sp=None, debug=False):
    """calculateAge either from a given Star Formation History (SFH) or from a 
    *ugriz* SED. If a SED is given (the default) then it calculates the SFH by 
    calling `calculateSFH()`.
    
    Parameters
    ----------
    SED : array-like
        Either an SED of *ugriz* magnitudes or SFH parameters (`tau`, 
        `tStart`, `sfTrans`, `sfSlope` of Shimha 2014 equation 6) if `SED` is 
        set to `False`.
    redshift : float
        The redshift of the region where the age is being calculated.
    SEDerr : array-like
        Optional if `SED` is `False` and `SED` is the SFH parameters. 
    SNID : int
        This is the ID number for the SN. Used for 'unique' ID while saving 
        data. Needed if `SED = True`.
    sp : fsps.fsps.StellarPopulation, optional
        An optional argument of an FSPS StellarPopulation if you don't want 
        make a new one with the default settings of `calculateSFH()` and this 
        analysis.
    debug : bool
        Flag to have MCMC run incredibly short and in no way accurately. Also 
        does not save resulting chain. Should take around ~12 mins to get a 
        value of one SN.
    
    Returns
    -------
    sfh : np.array
        The mean star formation history parameters (tau, tStart, sfTrans, 
        sfSlope). Either calculated from the SED given in `SED` or the 
        parameters given in `SED` if `SED = False`. A full posterior 
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
    logger = logging.getLogger("fsps-age.calculateAge.calculateAge")
    logger.info('called calculateAge')

    # calculate SFH
    logger.info('Calculating SFH')
    # no need for SNID. We can save in this function.
    samples = calculateSFH(SED, SEDerr, redshift, SNID, debug=debug)
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

    #with SFH -> Calculate age
    '''
    How do I maintain the information of the posterior distribution?
    each column of from the MCMC sample is a "possible" SFH. So we should 
    just calculate the expected age for each index. This gives us the same 
    number of ages as samples and gives a us a distribution for the age.
    '''

    age = integrate_age(tau, tStart, sfTrans, sfSlope, redshift)

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
        logger.debug('saving full MCMC & age data')
        samples = np.append(samples, age.reshape(age.size, 1), axis=-1)

        header = 'logzsol\tdust2\ttau\ttStart\tsfTrans\tsfAngle\tc\tage\ndex\t\t1/Gyr\tGyr\tGyr\t\tmag\tGyr'

        np.savetxt('resources/SN{}_chain.tsv'.format(SNID), samples, 
                delimiter='\t', header=header)
        logger.info('saved resources/SN{}_chain.tsv'.format(SNID))

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

def star_formation_gupta(t, tau, sfTrans, sfSlope):
    '''
    Defined piecewise around `sfTrans` time. 
    no variables can be arrays! ALso must return a float!
    '''
    if t <= sfTrans:
        sf = t*np.e**(-t/tau)
    else:
        sf = ( star_formation_gupta(sfTrans, tau, sfTrans, sfSlope) +
             sfSlope*ramp(t-sfTrans) ) #\Delta \def sfSlope*ramp() in Simha

    #only return a positive star formation
    if sf <= 0:
        return 0
    else:
        return sf

def t_star_formation_gupta(t, tau, sfTrans, sfSlope):
    return t*star_formation_gupta(t, tau, sfTrans, sfSlope)
##################
"""
Usage (local):
    pytest -v -rx --cov-config=.coveragerc --cov --runxfail

Usage-CRC (production):
    pytest -v -rx --runxfail

Usage-quick (no StellarPopulations)
    pytest -v -rx --cov-config=.coveragerc --cov --runxfail
"""
import pytest
import numpy as np
import fsps

import calculateAge

longrun = pytest.mark.skipif(False, reason="Runs too long for quick tests.")

class BaseTestCase():
    #only set up one stellar population
    sp = fsps.StellarPopulation(zcontinuous=2, 
                  cloudy_dust=True, add_neb_emission = True,
                  sfh=5)

    # only used for `lnlike()`, 'lnprior()' & 'lnprob()' so sfSlope is now phi
    theta_low = [0, -1, 0, 0, 0, np.arctan(-20.1), -35.1]   # outside bounds
    theta_high =[-3, -1, 10.1, 13, 13, np.arctan(20.1), -4.9]  # outside bounds
    theta_mean = [-0.5, 0.0, 1.0, 2.0, 10.0, 0.0, -25]     # at Gaussian means
    theta = [-0.3, 0.2, 1.0, 2.0, 10.0, np.arctan(5.0), -25]

    # from circle 1
    SED = [20.36, 18.76, 17.99, 17.67, 17.39]
    # is this the same as in code? -- What should the error be
    # Gupta is 0.05 to 0.5, lets say 0.1
    # Messier is 0.001
    # circle was currently 0.01 
    SED_err = [0.01, 0.01, 0.01, 0.01, 0.01]
    redshift = 0.05
    # details in 2017-07-28 lab notebook
    sf_parameters1 = [-0.5, 0.1, 0.5, 1.5, 9.0, -1.0, -25]  # old
    sf_fit_params1 = [-2.50, 0.01, 7.17, 7.94, 10.40, -5.25, -23.48]
    c = -25.0
    # from circle 3
    sf_parameters3 = [-0.5, 0.1, 7.0, 3.0, 10.0, 15.0, -25] # young


class TestRunFSPS(BaseTestCase):
    @pytest.mark.xfail(reason="Too tight of a tolerance for changes in newer FSPS.")
    @longrun
    def test_modelStillWorks_oldFSPS(self):
        """These have been calculated before (on crc) when creating new circle 
        test. Check out lab notes on 2017-07-28"""
        # is it the same to the 1/100th of a magnitude?
        assert np.allclose(calculateAge.runFSPS(self.sp, self.redshift,
                                                *self.sf_parameters1[:-1]),
                           np.array(self.SED) - self.sf_parameters1[-1],
                           atol=9e-03)

    @longrun
    def test_modelStillWorks_newFSPS(self):
        """These have been calculated before (on crc) when creating new circle 
        test. Check out lab notes on 2017-07-28"""
        # is it the same to the 1/100th of a magnitude?
        assert np.allclose(calculateAge.runFSPS(self.sp, self.redshift,
                                                *self.sf_parameters1[:-1]),
                           np.array(self.SED) - self.sf_parameters1[-1],
                           atol=2e-02)


class Test_lnlike(BaseTestCase):
    @longrun
    def test_rejectShortSED(self):
        """tests both the SED and the SED_err at the same time"""
        with pytest.raises(ValueError, match=r'need to be arrays of length 5'):
            calculateAge.lnlike(self.theta, [1], self.SED_err, self.redshift,
                                self.sp)
            calculateAge.lnlike(self.theta, self.SED, [1], self.redshift,
                                self.sp)
    
    @longrun
    def test_rejectLongSED(self):
        """tests both the SED and the SED_err at the same time"""
        with pytest.raises(ValueError, match=r'need to be arrays of length 5'):
            calculateAge.lnlike(self.theta, [1,2,3,4,5,6], self.SED_err,
                                self.redshift, self.sp)
            calculateAge.lnlike(self.theta, self.SED, [1,2,3,4,5,6],
                                self.redshift, self.sp)

    @longrun
    def test_rejectShortTheta(self):
        with pytest.raises(ValueError, match=r'^Likelihood expects 7'):
            calculateAge.lnlike(self.theta[:-1], self.SED, self.SED_err,
                                self.redshift, self.sp)

    @longrun
    def test_rejectLongTheta(self):
        with pytest.raises(ValueError, match=r'^Likelihood expects 7'):
            # some how I can't use self.theta.append(5)
            # also that changes the value of theta for future tests
            hold = [-0.3, 0.2, 1.0, 2.0, 10.0, 5.0, -25, 5]
            calculateAge.lnlike(hold, self.SED, self.SED_err, self.redshift,
                                self.sp)

    @longrun
    def test_likeChanges(self):
        """likelihood should be higher for the correct star formation parameters and lower when changing even 1 parameter"""
        # self.sf_parameters1 is circle 1 sf parameters
        sf_params1 = self.sf_parameters1[:]
        sf_params1[5] = np.arctan(sf_params1[5])    # convert from slope to phi
        sf_params1[2] = 7.0     # change tau
        sf_params2 = self.sf_parameters1[:]
        sf_params2[5] = np.arctan(sf_params2[5])    # convert from slope to phi
        sf_params2[3] = 3.0    # change t_start

        # need to some positive sfSlope to see change in t_trans
        # some how it needs a big change in transition. 
        # this has to do with forcing 1 solar mass of production total. 
        sf_params_trans = self.sf_parameters1[:]
        sf_params_trans[5] = np.arctan(2.0)    # and convert from slope to phi
        sf_params_trans_changed = sf_params_trans[:]
        sf_params_trans_changed[4] = 12.0    # change t_trans

        sf_params4 = self.sf_parameters1[:]
        sf_params4[5] = np.arctan(15.0)    # change sf_slope and convert to phi

        assert calculateAge.lnlike(self.sf_parameters1, self.SED, self.SED_err, self.redshift, self.sp) >calculateAge.lnlike(sf_params1, self.SED, self.SED_err, self.redshift, self.sp), "Changing tau should lower likelihood"
        assert calculateAge.lnlike(self.sf_parameters1, self.SED, self.SED_err, self.redshift, self.sp) > calculateAge.lnlike(sf_params2, self.SED, self.SED_err, self.redshift, self.sp), "Changing t_start should lower likelihood"
        # need another value of sf_slope to test sf_trans, so using
        # sf_params_trans rather then self.sf_parameters1
        n1 = calculateAge.lnlike(sf_params_trans, self.SED, self.SED_err, self.redshift, self.sp)
        n2 = calculateAge.lnlike(sf_params_trans_changed, self.SED, self.SED_err, self.redshift, self.sp)
        print(n1, n2)
        assert n1 > n2 , "Changing t_trans should lower likelihood"
        assert calculateAge.lnlike(self.sf_parameters1, self.SED, self.SED_err, self.redshift, self.sp) > calculateAge.lnlike(sf_params4, self.SED, self.SED_err, self.redshift, self.sp), "Changing sf_slope should lower likelihood"
        assert calculateAge.lnlike(self.sf_parameters1, self.SED, self.SED_err, self.redshift, self.sp) > calculateAge.lnlike(self.sf_parameters3, self.SED, self.SED_err, self.redshift, self.sp),  "Correct SF parameters should be more likely than another set"
    
    @longrun
    def test_idk(self):
        """ make sure correct parameters have a higher likelihood then median (or modal) parameters from fit. 

        see also similar test of `lnprob`.
        """
        # need to update sf_parameter `sfSlope` to `phi`.
        input1 = self.sf_parameters1[:]
        input1[5] = np.arctan(input1[5])
        input2 = self.sf_fit_params1[:]
        input2[5] = np.arctan(input2[5])

        assert calculateAge.lnlike(input1, self.SED, self.SED_err, self.redshift, self.sp) > calculateAge.lnlike(input2, self.SED, self.SED_err, self.redshift, self.sp),  "Input SF parameters (from circle test) should be more likely than median results of bad fits."


class Test_lnprior(BaseTestCase):
    def test_lowValues(self):
        assert calculateAge.lnprior(self.theta_low, self.redshift) == -np.inf, "Value should be outside prior range but ln-probability is not minus infinity"

    def test_highValues(self):
        assert calculateAge.lnprior(self.theta_high, self.redshift) == -np.inf, "Value should be outside prior range but ln-probability is not minus infinity"

    def test_transitionBeforeStart(self):
        # start with passing values
        theta = self.theta[:]
        # Change start to be too close to transition
        theta[3] = theta[4] - 1.9
        assert calculateAge.lnprior(theta, self.redshift) == -np.inf, "Value should be outside prior range but ln-probability is not minus infinity"

    def test_passingValues(self):
        assert calculateAge.lnprior(self.theta, self.redshift) > -np.inf, "Prior failing unexpectedly"

    def test_correctResult_mean(self):
        """Takes a calculated a prior probability from `self.theta_mean`, from http://www.wolframalpha.com/input/?i=-1*((0.0-0.0)**2%2F(2*0.3**2)+%2B+log(sqrt(2*pi)*0.3)+%2B+(-0.5+-+-0.5)**2%2F(2*0.6**2)+%2B+log(sqrt(2*pi)*0.6)), then compares to see if calcuateAge.lnprior calculated the same value.

        Probability comes from `self.theta_mean`: 
        CENTER_Z = -0.5
        SIGMA_Z = 0.5
        logzsol = CENTER_Z
        center = 0.0
        sigma = 0.3
        dust2 = center
        sfSlope = 0.0
        theta = [logzsol, dust2, 1.0, 2.0, 10.0, sfSlope, -25]
        expected = -1*((center-dust2)**2/(2*sigma**2) +
                   np.log(np.sqrt(2*np.pi)*sigma) +
                   (CENTER_Z - logzsol)**2/(2*SIGMA_Z**2) +
                   np.log(np.sqrt(2*np.pi)*SIGMA_Z))
        """
        # It looks like the return value has 16 decimals
        expected = -0.12307863831741855
        assert np.isclose(calculateAge.lnprior(self.theta_mean, self.redshift),
                          expected), "Prior not returning expected result."

    def test_correctResult_other(self):
        """Takes a calculated a prior probability from `self.theta`, from http://www.wolframalpha.com/input/?i=-1*((0.0-0.2)**2%2F(2*0.3**2)+%2B+log(sqrt(2*pi)*0.3)+%2B+(-0.5+-+-0.3)**2%2F(2*0.6**2)+%2B+log(sqrt(2*pi)*0.6)), then compares to see if calcuateAge.lnprior calculated the same value.

        Probability comes from `self.theta_mean`: 
        CENTER_Z = -0.5
        SIGMA_Z = 0.5
        logzsol = -0.3
        center = 0.0
        sigma = 0.3
        dust2 = 0.2
        sfSlope = 5.0
        theta = [logzsol, dust2, 1.0, 2.0, 10.0, sfSlope, -25]
        expected = -1*((center-dust2)**2/(2*sigma**2) +
                   np.log(np.sqrt(2*np.pi)*sigma) +
                   (CENTER_Z - logzsol)**2/(2*SIGMA_Z**2) +
                   np.log(np.sqrt(2*np.pi)*SIGMA_Z))
        """
        # It looks like the return value has 16 decimals
        expected = -0.4008564160951964
        calculated = calculateAge.lnprior(self.theta, self.redshift)
        assert np.isclose(calculated, expected), "Prior not returning " \
               "expected result. {} {}".format(calculated, expected)

    def test_priorChanges(self):
        assert calculateAge.lnprior(self.theta_mean, self.redshift) > calculateAge.lnprior(self.theta, self.redshift) , "Mean values should have a higher prior probability than other accepted values."

    def test_logzsol25(self):
        """why is logzsol = -2.5 so popular?"""
        theta = self.theta_mean[:]
        theta[0] = -2.5
        assert calculateAge.lnprior(self.theta_mean, self.redshift) > calculateAge.lnprior(theta, self.redshift) , "Mean values should have a higher prior than logzsol == -2.5."

        theta = self.theta[:]
        theta[0] = -2.5
        assert calculateAge.lnprior(self.theta_mean, self.redshift) > calculateAge.lnprior(theta, self.redshift) , "'Normal' values should have a higher prior than logzsol == -2.5."


class Test_lnprob(BaseTestCase):
    @longrun
    def test_posteriorFail(self):
        theta = [0, -1, 0, 0, 0, np.arctan(-20.1), -35.1]
        assert calculateAge.lnprob(theta, self.SED, self.SED_err, self.redshift, self.sp) == -np.inf, "Value should be outside prior range but ln-probability is not minus infinity"

    @longrun
    def test_posteriorPass(self):
        theta = [-0.3, 0.2, 1.0, 2.0, 10.0, np.arctan(5.0), -25]
        assert calculateAge.lnprob(theta, self.SED, self.SED_err, self.redshift, self.sp) == calculateAge.lnprior(theta, self.redshift) + calculateAge.lnlike(theta, self.SED, self.SED_err, self.redshift, self.sp)

    @longrun
    def test_posteriorChange(self):
        """The correct values, the ones used to construct `self.SED`, of
        circle test 1 should have a higher posterior probability the MCMC
        results found on 2017-08-01 from CRC: 153454.1 ("globalCircle-07-31")
        """
        # sf_parameters1 needs to be converted from `sfSlope` to `phi`
        input_ = self.sf_parameters1[:]
        input_[5] = np.arctan(input_[5])
        assert calculateAge.lnprob(input_, self.SED, self.SED_err, self.redshift, self.sp) > calculateAge.lnprob([-2.5, 0.01, 7.17, 7.94, 10.40, np.arctan(-5.24), -23.48], self.SED, self.SED_err, self.redshift, self.sp)

    @longrun
    def test_idk(self):
        """ make sure correct parameters have a higher likelihood then median (or modal) parameters from fit. 

        see also similar test of `lnlike`.
        """
        # sf_parameters1 & sf_fit_params1& need to be converted from `sfSlope`
        # to `phi`
        input1 = self.sf_parameters1[:]
        input1[5] = np.arctan(input1[5])
        input2 = self.sf_fit_params1[:]
        input2[5] = np.arctan(input2[5])

        assert calculateAge.lnprob(input1, self.SED, self.SED_err, self.redshift, self.sp) > calculateAge.lnprob(input2, self.SED, self.SED_err, self.redshift, self.sp),  "Input SF parameters (from circle test) should be more likely than median results of bad fits."

    @longrun
    def test_same_overtime(self):
        """
        Burn-in plots see the likelihood getting better over time. Final
        results of MCMC return a very bad fit as the "best-fit". But other
        tests prove that the very bad "best-fit" is truly (via likelihood and
        proprietor calculations) very bad and worse then the expected fit.
        """
        # sf_parameters1 needs to be converted from `sfSlope` to `phi`
        input_ = self.sf_parameters1[:]
        input_[5] = np.arctan(input_[5])

        run1 = calculateAge.lnprob(input_, self.SED, self.SED_err,
                                  self.redshift, self.sp)
        run2 = calculateAge.lnprob(input_, self.SED, self.SED_err,
                                  self.redshift, self.sp)
        # skip over many
        for i in range(100):
            _ = calculateAge.lnprob(input_, self.SED,
                                    self.SED_err, self.redshift, self.sp)
        run3 = calculateAge.lnprob(input_, self.SED, self.SED_err,
                                  self.redshift, self.sp)
        assert run1 == run2 == run3, "Results shouldn't change over runs"


class TestIntegrateAge(BaseTestCase):

    def test_integrate_age_argument_amount(self):
        with pytest.raises(TypeError, match=r'missing \d required positional arguments'):
            calculateAge.integrate_age([1])
        with pytest.raises(TypeError, match=r'positional arguments but \d were given$'):
            calculateAge.integrate_age([1], [2], 3, 4, 5, 6)

    def test_integrate_age_argument_type(self):
        with pytest.raises(TypeError, match=r'must be array_like'):
            calculateAge.integrate_age({'a': 1}, {'b': 1}, {'c': 1}, {'d': 1},
                                       {'e': 0.05})

    def test_ints_floats_inputs(self):
        """
        should calculate an age and return a array of floats. Just check the first element.
        """
        assert isinstance(calculateAge.integrate_age(0.1, 1.0, 13.0, 0.0, 0.05)[0], float)

    def test_lists_tuples_inputs(self):
        '''no need to test list of length 1. Python often treats that like a
        float.
        '''
        assert isinstance(calculateAge.integrate_age([0.1, 0.1], [1.0, 1.0], [13.0, 13.0], [0.0, 0.0], 0.05)[0], float)

    def test_np_ndarrays_inputs(self):
        '''no need to test array of length 1. Python often treats that like a
        float.
        '''
        assert isinstance(calculateAge.integrate_age(np.array([0.1, 0.1]), np.array([1.0, 1.0]), np.array([13.0, 13.0]), [0.0, 0.0], 0.05)[0], float)

    def test_failure_with_2d_array(self):
        """
        Test tau, tStart, sfTrans, & sfSlope to fail on two dimensional inputs.
        Redshift is not array-like so this test is not needed.
        """
        with pytest.raises(TypeError, match=r'cannot be more then 1 dimensional.$'):
            calculateAge.integrate_age(np.random.randn(3,3), 2, 3, 4, 0.5)
        with pytest.raises(TypeError, match=r'cannot be more then 1 dimensional.$'):
            calculateAge.integrate_age(1, np.random.randn(3,3), 3, 4, 0.5)
        with pytest.raises(TypeError, match=r'cannot be more then 1 dimensional.$'):
            calculateAge.integrate_age(1, 2, np.random.randn(3,3), 4, 0.5)
        with pytest.raises(TypeError, match=r'cannot be more then 1 dimensional.$'):
            calculateAge.integrate_age(1, 2, 3, np.random.randn(3,3), 0.5)

    def test_failure_with_unequal_arrays(self):
        with pytest.raises(ValueError, match=r'need to be equal lengths$'):
            calculateAge.integrate_age(1, [2,2], 3, 4, 0.05)

    def test_integrate_age_argument_redshift(self):
        with pytest.raises(ValueError, match=r'greater than zero'):
            calculateAge.integrate_age(np.array([1, 2]), np.array([0.1, 0.2]), np.array([0.1, 0.2]), np.array([1, 2]), -0.4)

    # also acts as a test to change in tStart
    def test_correct_young_age(self):
        assert np.isclose(12, calculateAge.integrate_age(0.1, 1.0, 13.0,
                                                         0.0, 0.05), 
                          atol = 0.05)

    def correct_old_age(self):
        assert np.isclose(3, calculateAge.integrate_age(0.1, 10.0, 13.0,
                                                        0.0, 0.05),
                          atol = 0.05)

    def test_change_tau(self):
        """With a larger tau, we should have a younger population
        having a large sfSlope will mess up this test.
        """
        # this could possibly use part of `self.sf_parameters1`
        tau_high = np.array([9.0])    # should be closer to a flat sfh
        tau_low = np.array([1.0])    # should be a shorter/sharper burst
        tStart = np.array([3.0])
        sfTrans = np.array([12.0])
        sfSlope = np.array([-1.0])
        redshift = 0.05

        assert calculateAge.integrate_age(tau_high, tStart, sfTrans, sfSlope, redshift) < calculateAge.integrate_age(tau_low, tStart, sfTrans, sfSlope, redshift), "A larger tau, we should have a younger population"

    def test_change_sfTrans(self):
        """
        This needs some sfSlope to have any effect
        """
        tau = np.array([1.0])
        tStart = np.array([3.0])
        sfTrans_early = np.array([5.0])
        sfTrans_late = np.array([12.0])
        sfSlope = np.array([15.0])
        redshift = 0.05

        assert calculateAge.integrate_age(tau, tStart, sfTrans_late, sfSlope, redshift) < calculateAge.integrate_age(tau, tStart, sfTrans_early, sfSlope, redshift), "A later sfTrans should have a younger population."

    def test_chagne_sfSlope(self):
        """
        Should have an early sfTrans so that sfSlope has a magnified effect
        """
        tau = np.array([1.0])
        tStart = np.array([3.0])
        sfTrans = np.array([5.0])
        sfSlope_large = np.array([15.0])
        sfSlope_small = np.array([-5.0])
        redshift = 0.05

        assert calculateAge.integrate_age(tau, tStart, sfTrans, sfSlope_large, redshift) < calculateAge.integrate_age(tau, tStart, sfTrans, sfSlope_small, redshift), "A larger sfSlope should have a younger population."

    def test_romb_method(self):
        """test the method that requires integration via scipy.integrate.romb()
        Sadly this current SFH should not be allowed since sfTrans is not 2Gyrs
        past tStart.
        """
        assert isinstance(calculateAge.integrate_age(np.array([2.72, 2.72]), np.array([6.701, 6.701]), np.array([6.7067, 6.7067]), np.array([-19.33, -19.3]), 0.0043)[0], float)


class Test_SFH(BaseTestCase):
    def test_ramp_lowbounds(self):
        assert calculateAge.ramp(-0.01) == 0
        assert calculateAge.ramp(-5) == 0
        assert calculateAge.ramp(-500) == 0

    def test_ramp_inbounds(self):
        assert calculateAge.ramp(0) == 0
        assert calculateAge.ramp(1) == 1
        assert calculateAge.ramp(100) == 100

    def test_t_sfr(self):
        assert calculateAge.t_star_formation_gupta(5, 4, 8, 0) == 5*calculateAge.star_formation_gupta(5, 4, 8, 0)

    # def test_a(self):
    #     """ Circle 1 should be older then Circle 3"""
    #     assert calculateAge.starFormation(10, *self.sf_parameters1[:-1]) > calculateAge.starFormation(10, *self.sf_parameters3[:-1])
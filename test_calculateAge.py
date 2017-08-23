"""
Usage:
    pytest -v -rx --cov-config=.coveragerc --cov --runxfail
"""
import pytest
import numpy as np
import fsps

import calculateAge

class BaseTestCase():
    #only set up one stellar population
    sp = fsps.StellarPopulation(zcontinuous=2, 
                  cloudy_dust=True, add_neb_emission = True,
                  sfh=5)

    theta_low = [0, -1, 0, 0, 0, -20.1, -35.1]   # outside bounds
    theta_high =[-3, -1, 10.1, 13, 13, 20.1, -4.9]    # outside bounds
    theta_mean = [-0.5, 0.0, 1.0, 2.0, 10.0, 0.0, -25]     # at Gaussian means
    theta = [-0.3, 0.2, 1.0, 2.0, 10.0, 5.0, -25]

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
    def test_modelStillWorks_oldFSPS(self):
        """These have been calculated before (on crc) when creating new circle 
        test. Check out lab notes on 2017-07-28"""
        # is it the same to the 1/100th of a magnitude?
        assert np.allclose(calculateAge.runFSPS(self.sp, self.redshift,
                                                *self.sf_parameters1[:-1]),
                           np.array(self.SED) - self.sf_parameters1[-1],
                           atol=9e-03)

    def test_modelStillWorks_newFSPS(self):
        """These have been calculated before (on crc) when creating new circle 
        test. Check out lab notes on 2017-07-28"""
        # is it the same to the 1/100th of a magnitude?
        assert np.allclose(calculateAge.runFSPS(self.sp, self.redshift,
                                                *self.sf_parameters1[:-1]),
                           np.array(self.SED) - self.sf_parameters1[-1],
                           atol=2e-02)


class Test_lnlike(BaseTestCase):
    def test_rejectShortSED(self):
        """tests both the SED and the SED_err at the same time"""
        with pytest.raises(ValueError, match=r'need to be arrays of length 5'):
            calculateAge.lnlike(self.theta, [1], self.SED_err, self.redshift,
                                self.sp)
            calculateAge.lnlike(self.theta, self.SED, [1], self.redshift,
                                self.sp)

    def test_rejectLongSED(self):
        """tests both the SED and the SED_err at the same time"""
        with pytest.raises(ValueError, match=r'need to be arrays of length 5'):
            calculateAge.lnlike(self.theta, [1,2,3,4,5,6], self.SED_err,
                                self.redshift, self.sp)
            calculateAge.lnlike(self.theta, self.SED, [1,2,3,4,5,6],
                                self.redshift, self.sp)

    def test_rejectShortTheta(self):
        with pytest.raises(ValueError, match=r'^Likelihood expects 7'):
            calculateAge.lnlike(self.theta[:-1], self.SED, self.SED_err,
                                self.redshift, self.sp)

    def test_rejectLongTheta(self):
        with pytest.raises(ValueError, match=r'^Likelihood expects 7'):
            # some how I can't use self.theta.append(5)
            # also that changes the value of theta for future tests
            hold = [-0.3, 0.2, 1.0, 2.0, 10.0, 5.0, -25, 5]
            calculateAge.lnlike(hold, self.SED, self.SED_err, self.redshift,
                                self.sp)

    def test_likeChanges(self):
        """likelihood should be higher for the correct star formation parameters and lower when changing even 1 parameter"""
        # self.sf_parameters1 is circle 1 sf parameters
        sf_params1 = self.sf_parameters1[:]
        sf_params1[2] = 7.0     # change tau
        sf_params2 = self.sf_parameters1[:]
        sf_params2[3] = 3.0    # change t_start
        sf_params3 = self.sf_parameters1[:]
        sf_params3[4] = 10.0    # change t_trans
        sf_params4 = self.sf_parameters1[:]
        sf_params4[5] = 15.0    # change sf_slope

        assert calculateAge.lnlike(self.sf_parameters1, self.SED, self.SED_err, self.redshift, self.sp) >calculateAge.lnlike(sf_params1, self.SED, self.SED_err, self.redshift, self.sp), "Changing tau should lower likelihood"
        assert calculateAge.lnlike(self.sf_parameters1, self.SED, self.SED_err, self.redshift, self.sp) > calculateAge.lnlike(sf_params2, self.SED, self.SED_err, self.redshift, self.sp), "Changing t_start should lower likelihood"
        # assert calculateAge.lnlike(self.sf_parameters1, self.SED, self.SED_err, self.redshift, self.sp) > calculateAge.lnlike(sf_params3, self.SED, self.SED_err, self.redshift, self.sp), "Changing t_trans should lower likelihood"
        assert calculateAge.lnlike(self.sf_parameters1, self.SED, self.SED_err, self.redshift, self.sp) > calculateAge.lnlike(sf_params4, self.SED, self.SED_err, self.redshift, self.sp), "Changing sf_slope should lower likelihood"
        assert calculateAge.lnlike(self.sf_parameters1, self.SED, self.SED_err, self.redshift, self.sp) > calculateAge.lnlike(self.sf_parameters3, self.SED, self.SED_err, self.redshift, self.sp),  "Correct SF parameters should be more likely than another set"
    def test_idk(self):
        """ make sure correct parameters have a higher likelihood then median (or modal) parameters from fit. 

        see also similar test of `lnprob`.
        """
        assert calculateAge.lnlike(self.sf_parameters1, self.SED, self.SED_err, self.redshift, self.sp) > calculateAge.lnlike(self.sf_fit_params1, self.SED, self.SED_err, self.redshift, self.sp),  "Input SF parameters (from circle test) should be more likely than median results of bad fits."


class Test_lnprior(BaseTestCase):
    def test_lowValues(self):
        assert calculateAge.lnprior(self.theta_low, self.redshift) == -np.inf, "Value should be outside prior range but ln-probability is not minus infinity"

    def test_highValues(self):
        assert calculateAge.lnprior(self.theta_high, self.redshift) == -np.inf, "Value should be outside prior range but ln-probability is not minus infinity"

    def test_passingValues(self):
        assert calculateAge.lnprior(self.theta, self.redshift) > -np.inf, "Prior failing unexpectedly"

    def test_correctResult_mean(self):
        """Takes a calculated a prior probability from `self.theta_mean`, from http://www.wolframalpha.com/input/?i=-1*(1.5*log(1+%2B+0.0**2)+%2B+(0.0-0.0)**2%2F(2*0.3**2)+%2B+log(sqrt(2*pi)*0.3)+%2B+(-0.5+%2B+0.5)**2%2F(2*0.5**2)+%2B+log(sqrt(2*pi)*0.5)), then compares to see if calcuateAge.lnprior calculated the same value.

        Probability comes from `self.theta_mean`: 
        CENTER_Z = -0.5
        SIGMA_Z = 0.5
        logzsol = CENTER_Z
        center = 0.0
        sigma = 0.3
        dust2 = center
        sfSlope = 0.0
        theta = [logzsol, dust2, 1.0, 2.0, 10.0, sfSlope, -25]
        expected = -1*(1.5*np.log(1 + sfSlope**2) + (center-dust2)**2/(2*sigma**2) + np.log(np.sqrt(2*np.pi)*sigma) + (CENTER_Z - logzsol)**2/(2*SIGMA_Z**2) + np.log(np.sqrt(2*np.pi)*SIGMA_Z))
        """
        # It looks like the return value has 16 decimals
        expected = 0.0592429184765358
        assert np.isclose(calculateAge.lnprior(self.theta_mean, self.redshift),
                          expected), "Prior not returning expected result."

    def test_correctResult_other(self):
        """Takes a calculated a prior probability from `self.theta`, from http://www.wolframalpha.com/input/?i=-1*(1.5*log(1+%2B+5.0**2)+%2B+(0.0-0.2)**2%2F(2*0.3**2)+%2B+log(sqrt(2*pi)*0.3)+%2B+(-0.5+%2B+0.3)**2%2F(2*0.5**2)+%2B+log(sqrt(2*pi)*0.5)), then compares to see if calcuateAge.lnprior calculated the same value.

        Probability comes from `self.theta_mean`: 
        CENTER_Z = -0.5
        SIGMA_Z = 0.5
        logzsol = -0.3
        center = 0.0
        sigma = 0.3
        dust2 = 0.2
        sfSlope = 5.0
        theta = [logzsol, dust2, 1.0, 2.0, 10.0, sfSlope, -25]
        expected = -1*(1.5*np.log(1 + sfSlope**2) + 
                    (center-dust2)**2/(2*sigma**2) + 
                    np.log(np.sqrt(2*np.pi)*sigma) + 
                    (CENTER_Z - logzsol)**2/(2*SIGMA_Z**2) + 
                    np.log(np.sqrt(2*np.pi)*SIGMA_Z))
        """
        # It looks like the return value has 16 decimals
        expected = -5.1301241107779094
        assert np.isclose(calculateAge.lnprior(self.theta, self.redshift),
                         expected), "Prior not returning expected result."

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
    # pass
    def test_posteriorFail(self):
        theta = [0, -1, 0, 0, 0, -20.1, -35.1]
        assert calculateAge.lnprob(theta, self.SED, self.SED_err, self.redshift, self.sp) == -np.inf, "Value should be outside prior range but ln-probability is not minus infinity"

    def test_posteriorPass(self):
        theta = [-0.3, 0.2, 1.0, 2.0, 10.0, 5.0, -25]
        assert calculateAge.lnprob(theta, self.SED, self.SED_err, self.redshift, self.sp) == calculateAge.lnprior(theta, self.redshift) + calculateAge.lnlike(theta, self.SED, self.SED_err, self.redshift, self.sp)

    def test_posteriorChange(self):
        """The correct values, the ones used to construct `self.SED`, of
        circle test 1 should have a higher posterior probability the MCMC
        results found on 2017-08-01 from CRC: 153454.1 ("globalCircle-07-31")
        """
        assert calculateAge.lnprob(self.sf_parameters1, self.SED, self.SED_err, self.redshift, self.sp) > calculateAge.lnprob([-2.5, 0.01, 7.17, 7.94, 10.40, -5.24, -23.48], self.SED, self.SED_err, self.redshift, self.sp)

    def test_idk(self):
        """ make sure correct parameters have a higher likelihood then median (or modal) parameters from fit. 

        see also similar test of `lnlike`.
        """
        assert calculateAge.lnprob(self.sf_parameters1, self.SED, self.SED_err, self.redshift, self.sp) > calculateAge.lnprob(self.sf_fit_params1, self.SED, self.SED_err, self.redshift, self.sp),  "Input SF parameters (from circle test) should be more likely than median results of bad fits."

    def test_same_overtime(self):
        """
        Burn-in plots see the likelihood getting better over time. Final
        results of MCMC return a very bad fit as the "best-fit". But other
        tests prove that the very bad "best-fit" is truly (via likelihood and
        proprietor calculations) very bad and worse then the expected fit.
        """
        run1 = calculateAge.lnprob(self.sf_parameters1, self.SED, self.SED_err,
                                  self.redshift, self.sp)
        run2 = calculateAge.lnprob(self.sf_parameters1, self.SED, self.SED_err,
                                  self.redshift, self.sp)
        # skip over many
        for i in range(100):
            _ = calculateAge.lnprob(self.sf_parameters1, self.SED,
                                    self.SED_err, self.redshift, self.sp)
        run3 = calculateAge.lnprob(self.sf_parameters1, self.SED, self.SED_err,
                                  self.redshift, self.sp)
        assert run1 == run2 == run3, "Results shouldn't change over runs"


class Test_SFH(BaseTestCase):
    def test_ramp_lowbounds(self):
        assert calculateAge.ramp(-0.01) == 0
        assert calculateAge.ramp(-5) == 0
        assert calculateAge.ramp(-500) == 0

    def test_ramp_inbounds(self):
        assert calculateAge.ramp(0) == 0
        assert calculateAge.ramp(1) == 1
        assert calculateAge.ramp(100) == 100

    def test_heavyside_lowbounds(self):
        assert calculateAge.heavyside(-0.01) == 0
        assert calculateAge.heavyside(-5) == 0
        assert calculateAge.heavyside(-500) == 0

    def test_heavyside_inbounds(self):
        assert calculateAge.heavyside(0) == 1
        assert calculateAge.heavyside(1) == 1
        assert calculateAge.heavyside(100) == 1

class TestIntegrateAge(BaseTestCase):

    def test_integrate_age_argument_amount(self):
        with pytest.raises(TypeError, match=r'missing \d required positional arguments'):
            calculateAge.integrate_age([1])
        with pytest.raises(TypeError, match=r'positional arguments but \d were given$'):
            calculateAge.integrate_age([1], [2], 3, 4, 5, 6)

    def test_integrate_age_argument_type(self):
        with pytest.raises(TypeError, match=r'must be of type'):
            calculateAge.integrate_age(1, 1, 1, 1, 0.05)

    # def test if all have the same length:

    def test_integrate_age_argument_redshift(self, match=r'greater than zero'):
        with pytest.raises(ValueError):
            calculateAge.integrate_age(np.array([1, 2]), np.array([0.1, 0.2]), np.array([0.1, 0.2]), np.array([1, 2]), -0.4)

    def test_change_tau(self):
        """With a larger tau, we should have a younger population"""
        # this could possibly use part of `self.sf_parameters1`
        tau_high = np.array([7.0])    # should be closer to a flat sfh
        tau_low = np.array([1.0])    # should be a shorter/sharper burst
        tStart = np.array([3.0])
        sfTrans = np.array([10.0])
        sfSlope = np.array([15.0])
        redshift = 0.05

        assert calculateAge.integrate_age(tau_high, tStart, sfTrans, sfSlope, redshift) < calculateAge.integrate_age(tau_low, tStart, sfTrans, sfSlope, redshift), "A larger tau, we should have a younger population"

    # def test_a(self):
    #     """ Circle 1 should be older then Circle 3"""
    #     assert calculateAge.starFormation(10, *self.sf_parameters1[:-1]) > calculateAge.starFormation(10, *self.sf_parameters3[:-1])
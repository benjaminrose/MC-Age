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
    SED_err = [0.1, 0.1, 0.1, 0.1, 0.1]
    redshift = 0.05
    # details in 2017-07-28 lab notebook
    sf_parameters1 = [-0.5, 0.1, 0.5, 1.5, 9.0, -1.0, -25]  # old
    c = -25.0
    # from circle 3
    sf_parameters3 = [-0.5, 0.1, 7.0, 3.0, 10, 15.0, -25] # young


class Test_runFSPS(BaseTestCase):
    @pytest.mark.xfail(reason="Too tight of a tolerance. Newest FSPS adds too much variability.")
    def test_modelStillWorks_oldFSPS(self):
        """These have been calculated before (on crc) when creating new circle 
        test. Check out lab notes on 2017-07-28"""
        # is it the same to the 1/100th of a magnitude?
        assert np.allclose(calculateAge.runFSPS(self.sp, self.redshift,
                                                *self.sf_parameters1[:-1]),
                           np.array(self.SED) - self.sf_parameters1[-1],
                           atol=1e-02)

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
        """likelihood should be higher for the correct star formation parameters"""
        assert calculateAge.lnlike(self.sf_parameters1, self.SED, self.SED_err, self.redshift, self.sp) > calculateAge.lnlike(self.sf_parameters3, self.SED, self.SED_err, self.redshift, self.sp)


class Test_lnprior(BaseTestCase):
    def test_lowValues(self):
        assert calculateAge.lnprior(self.theta_low, self.redshift) == -np.inf, "Value should be outside prior range but ln-probability is not minus infinity"

    def test_highValues(self):
        assert calculateAge.lnprior(self.theta_high, self.redshift) == -np.inf, "Value should be outside prior range but ln-probability is not minus infinity"

    def test_passingValues(self):
        assert calculateAge.lnprior(self.theta, self.redshift) > -np.inf, "Prior failing unexpectedly"

    def test_correctResult(self):
        # should be the same as self.theta_mean
        CENTER_Z = -0.5
        SIGMA_Z = 0.5
        logzsol = CENTER_Z
        center = 0.0
        sigma = 0.3
        dust2 = center
        sfSlope = 0.0
        theta = [logzsol, dust2, 1.0, 2.0, 10.0, sfSlope, -25]
        expected = -1*(1.5*np.log(1 + sfSlope**2) + (center-dust2)**2/(2*sigma**2) + np.log(np.sqrt(2*np.pi)*sigma) + (CENTER_Z - logzsol)**2/(2*SIGMA_Z**2) + np.log(np.sqrt(2*np.pi)*SIGMA_Z))
        assert calculateAge.lnprior(theta, self.redshift) == expected, "Prior not returning expected result."

        # should be the same at self.theta
        logzsol = -0.3
        dust2 = 0.2
        sfSlope = 5.0
        theta = [logzsol, dust2, 1.0, 2.0, 10.0, sfSlope, -25]
        expected = -1*(1.5*np.log(1 + sfSlope**2) +
                   (center-dust2)**2/(2*sigma**2) +
                   np.log(np.sqrt(2*np.pi)*sigma) +
                   (CENTER_Z - logzsol)**2/(2*SIGMA_Z**2) +
                   np.log(np.sqrt(2*np.pi)*SIGMA_Z))
        assert calculateAge.lnprior(theta, self.redshift) == expected, "Prior not returning expected result."

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

    # def test_a(self):
    #     """ Circle 1 should be older then Circle 3"""
    #     assert calculateAge.starFormation(10, *self.sf_parameters1[:-1]) > calculateAge.starFormation(10, *self.sf_parameters3[:-1])
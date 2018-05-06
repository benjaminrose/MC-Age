import pytest
import numpy as np
import fsps

import util

class TestSP():
    def test_get_sp(self):
        """
        """
        sp = util.get_sp()
        assert type(sp) is fsps.fsps.StellarPopulation
        # zcontinuous can't be tested this way
        # sp.params['zcontinuous'] == 2
        assert sp.params['cloudy_dust'] == True
        assert sp.params['add_neb_emission'] == True
        assert sp.params['sfh'] == 5

    # def metals_dust_effects_magnitudes(self):

    def test_dust_effects_magnitudes(self):
        """
        Varying the dust parameters, like what the MCMC does, should
        change the output SED.
        """
        #getting a new stellar population is not a great idea for speed.
        sp = util.get_sp()
        sdss_bands = ['sdss_u', 'sdss_g', 'sdss_r', 'sdss_i', 'sdss_z']

        # mean of prior
        sp.params['dust1'] = 0.6
        sp.params['dust2'] = 0.3
        one = sp.get_mags(tage=13, bands=sdss_bands)

        # one sigma of prior
        sp.params['dust1'] = 1.0
        sp.params['dust2'] = 0.5
        two = sp.get_mags(tage=13, bands=sdss_bands)

        assert not np.allclose(one, two)

    # def sf_start_dust_effects_magnitudes(self):

    # def sf_tau_dust_effects_magnitudes(self):

    # def sf_transition_dust_effects_magnitudes(self):

    # def sf_slope_dust_effects_magnitudes(self):


class TestMode():
    mode = util.mode(np.random.randn(1000))
    def test_return_length(self):
        """
        """
        assert len(self.mode) == 3
    def test_return_value(self):
        assert abs(self.mode[0]) < 0.5, "mode is not close (<0.5) to zero."

class TestMedian():
    median = util.median(np.random.randn(1000))
    def test_return_length(self):
        """
        """
        assert len(self.median) == 3
    def test_return_value(self):
        assert abs(self.median[0]) < 0.1, "median is not close (<0.1) to zero."
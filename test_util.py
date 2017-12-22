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
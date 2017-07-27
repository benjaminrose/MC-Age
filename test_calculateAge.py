# import sys
# sys.path.append('../')

# content of test_sysexit.py
import pytest
import numpy as np

import calculateAge

class Test_lnprior():
    def test_lowValues(self):
        redshift = 0.1
        theta = [0, -1, 0, 0, 0, -20.1, -35.1]
        assert calculateAge.lnprior(theta, redshift) == -np.inf, "Value should be outside prior range but ln-probability is not minus infinity"
        assert 2 == 2, "what says it fails?"

    def test_highValues(self):
        redshift = 0.1
        theta = [-3, -1, 10.1, 13, 13, 20.1, -4.9]
        assert calculateAge.lnprior(theta, redshift) == -np.inf, "Value should be outside prior range but ln-probability is not minus infinity"

    def test_passingValues(self):
        redshift = 0.1
        theta = [0, 1, 1, 2, 5, 0, -25]
        assert calculateAge.lnprior(theta, redshift) > -np.inf, "Should this fail?"

class Test_SFH():
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
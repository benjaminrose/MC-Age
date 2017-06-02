# import sys
# sys.path.append('../')

# content of test_sysexit.py
import pytest
import numpy as np

import calculateAge

class Test_calcualteAge_failures:
    def test_lnprior(self):
        redshift = 0.1
        theta = [0, 0, 0, 0, 0, 0, 0]
        assert calculateAge.lnprior(theta, redshift) == -np.inf, "Value should be outside prior range but ln-probability is not minus infinity"

class Test_calcualteAge_passes:
    def test_middle_of_prior_space(self):
        redshift = 0.1
        theta = [0, 1, 1, 2, 5, 0, -25]
        assert calculateAge.lnprior(theta, redshift) > -np.inf, "Should this fail?"
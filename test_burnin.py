import pytest
import numpy as np

import burnin

class Test_burnin():
    def test_burnin_argument_amount(self):
        with pytest.raises(TypeError, match=r'missing \d required positional arguments'):
            burnin.burnin([1])
        with pytest.raises(TypeError, match=r'positional arguments but \d were given$'):
            burnin.burnin([1],[2],3,4,5)

    def test_burnin_argument_type(self):
        with pytest.raises(TypeError, match=r'must be of type'):
            burnin.burnin(1, 1, 0.05)

    def test_burnin_argument_lists(self):
        with pytest.raises(ValueError, match=r'5 SDSS filters'):
            burnin.burnin([1], [2], 0.05, 1415)

    def test_burnin_argument_redshift(self, match=r'greater than zero'):
        with pytest.raises(ValueError):
            burnin.burnin([1, 2, 3, 4, 5], [0.1, 0.2, 0.3, 0.4, 0.5], -0.4)
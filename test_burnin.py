import pytest
import numpy as np

import burnin

class Test_burnin():
    def test_burnin_argument_amount(self):
        with pytest.raises(TypeError):
            burnin.burnin(1)
        with pytest.raises(TypeError):
            burnin.burnin(1,2,3,4,5)

    def test_burnin_argument_lists(self):
        with pytest.raises(TypeError):
            burnin.burnin(1, 2, 0.05, 1415)

    def test_burnin_argument_redshift(self):
        with pytest.raises(ValueError):
            burnin.burnin([1, 2, 3, 4, 5], [0.1, 0.2, 0.3, 0.4, 0.5], -0.4)
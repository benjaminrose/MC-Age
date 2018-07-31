import pytest
import numpy as np
import fsps

import redoGupta

class TestRedoGupta():
    """ Test running data

    This function should be tested. The 4 month import bug was here. But I
    combined too much. I can't test that it works, because calling this run my
    analyses for days! Here are a few tests that I can do.
    """
    def testBadDataSet(self):
        """Test failure with several invalid data set/path names"""
        with pytest.raises(ValueError, match=r'Invalid dataset'):
            redoGupta.redoGupta(1, dataset='not good')
        with pytest.raises(ValueError, match=r'Invalid dataset'):
            redoGupta.redoGupta(1, dataset='bad')
        with pytest.raises(ValueError, match=r'Invalid dataset'):
            redoGupta.redoGupta(1, dataset='one more once')

        # use these tests if can start to use arbitrary file names
        # with pytest.raises(FileNotFoundError, match=r'not good'):
        #     redoGupta.redoGupta(1, dataset='not good')
        # with pytest.raises(FileNotFoundError, match=r'bad'):
        #     redoGupta.redoGupta(1, dataset='bad')
        # with pytest.raises(FileNotFoundError, match=r'more'):
        #     redoGupta.redoGupta(1, dataset='one more once')
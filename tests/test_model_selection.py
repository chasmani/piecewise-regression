from piecewise_regression.main import Fit
from piecewise_regression.model_selection import ModelSelection

import numpy as np
import unittest
from importlib.machinery import SourceFileLoader

import os
import sys
sys.path.insert(1, os.path.join(sys.path[0], '..'))

DATA_SOURCE = "tests/data/data.txt"


class TestFit(unittest.TestCase):

    def test_it_works_with_differnet_options(self):
        """
        Chekc the the n_breakpoints converged True/False boolean gives the
        same answer with model selection and fit
        """
        # This is influenced by random seeds, so might not pass
        # with different random seeds (although it doensnt mean it isn't
        # working)
        np.random.seed(2)

        data = SourceFileLoader('data', DATA_SOURCE).load_module()

        xx = np.array(data.MUGGEO_1_XX)
        yy = np.array(data.MUGGEO_1_YY)

        ms = ModelSelection(xx, yy, n_boot=20)

        # For each n_breakpoints, chekc the ModelSelection vs fit results
        for n_breakpoints in range(1,10):
            fit = Fit(xx, yy, n_breakpoints=n_breakpoints, verbose=False, n_boot=20)

            fit_converged = fit.get_results()["converged"]
            ms_converged = ms.model_summaries[n_breakpoints]["converged"]
            print(n_breakpoints, fit_converged, ms_converged)
            

            self.assertEqual(fit_converged, ms_converged)

if __name__ == '__main__':
    unittest.main()

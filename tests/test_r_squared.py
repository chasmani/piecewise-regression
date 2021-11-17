import piecewise_regression.r_squared_calc as r_squared_calc
import numpy as np
import unittest
from importlib.machinery import SourceFileLoader

import os
import sys
sys.path.insert(1, os.path.join(sys.path[0], '..'))

DATA_SOURCE = "tests/data/data.txt"


class TestRSquared(unittest.TestCase):

    def test_some_data(self):

        ff = np.linspace(0, 10)
        yy = np.linspace(0, 10)

        rss, tss, r_2, adjusted_r_2 = r_squared_calc.get_r_squared(yy, ff, 1)

        self.assertEqual(r_2, 1)
        self.assertEqual(adjusted_r_2, 1)

        ff = [0.1, 0.8, 0.4, -4, 6, 12, 14, 1]
        yy = [1, 2, 3, 4, 5, 6, 6, 6]

        rss, tss, r_2, adjusted_r_2 = r_squared_calc.get_r_squared(yy, ff, 1)

        # Value calculated from sklearn's r2_score function
        r_2_from_sklearn = -6.405023255813953
        self.assertEqual(r_2_from_sklearn, r_2)

        data = SourceFileLoader('data', DATA_SOURCE).load_module()

        ff = np.array(data.BP_1_FF)
        yy = np.array(data.BP_1_YY)

        rss, tss, r_2, adjusted_r_2 = r_squared_calc.get_r_squared(yy, ff, 6)
        # Value calculated from sklearn's r2_score function
        r_2_from_sklearn = 0.9990626123719015
        self.assertEqual(r_2_from_sklearn, r_2)


if __name__ == '__main__':
    unittest.main()

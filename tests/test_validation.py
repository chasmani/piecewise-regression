import copy

import pandas as pd
import numpy as np
import unittest

from breakpoint_regression import Fit

class TestValidation(unittest.TestCase):

    def test_with_invalid_data_types(self):

        xx = np.linspace(0,10)
        yy = np.linspace(0,10)

        KWARGS = {
                "xx":xx,
                "yy":yy,
                "n_breakpoints":2,
                "start_values":[1,2]
            }

        # Lots of invalid data types
        for test_variable in ["xx", "yy", "n_breakpoints", "start_values"]:
            for invalid_value in [None, "hi", 12.1, 12, [1,1,1,1,1], 0, []]:
                
                new_kwargs = copy.deepcopy(KWARGS)
                new_kwargs[test_variable] = invalid_value

                self.assertRaises(ValueError, Fit, **new_kwargs)

        # Initial guesses outside data range
        bp_fit = Fit(yy, xx, n_breakpoints=2, start_values=[-100, 10000])

        

    def test_without_enough_args(self):

        xx = np.linspace(0,10)
        yy = np.linspace(0,10)

        self.assertRaises(TypeError, Fit, xx, yy, n_breakpoints=3)
        self.assertRaises(TypeError, Fit, xx, yy, start_values=[3])
        self.assertRaises(TypeError, Fit, xx, yy)
        self.assertRaises(TypeError, Fit, xx, n_breakpoints=3, start_values=[3])
        self.assertRaises(TypeError, Fit, yy, n_breakpoints=3, start_values=[3])


        

if __name__ == '__main__':
    unittest.main()
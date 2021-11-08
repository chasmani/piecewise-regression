from piecewise_regression import Fit
import copy

import numpy as np
import unittest

import os
import sys
sys.path.insert(1, os.path.join(sys.path[0], '..'))


class TestValidation(unittest.TestCase):

    def test_with_invalid_data_types(self):

        xx = np.linspace(0, 10)
        yy = np.linspace(0, 10)

        KWARGS = {
            "xx": xx,
            "yy": yy,
            "start_values": [1, 2]
        }

        # Lots of invalid data types
        for test_variable in ["xx", "yy", "start_values"]:
            for invalid_value in [None, "hi", 12.1, 12, 0, []]:

                new_kwargs = copy.deepcopy(KWARGS)
                new_kwargs[test_variable] = invalid_value

                self.assertRaises(ValueError, Fit, **new_kwargs)

        for test_variable in [
                "max_iterations", "tolerance",
                "min_distance_between_breakpoints",
                "min_distance_to_edge"]:
            for invalid_value in [None, "hi", [1, 1, 1, 1, 1], [], 0, -3]:

                new_kwargs = copy.deepcopy(KWARGS)
                new_kwargs[test_variable] = invalid_value

                self.assertRaises(ValueError, Fit, **new_kwargs)

    def test_without_enough_args(self):

        xx = np.linspace(0, 10)
        yy = np.linspace(0, 10)

        self.assertRaises(ValueError, Fit, xx, yy)
        self.assertRaises(TypeError, Fit, xx, start_values=[3])
        self.assertRaises(TypeError, Fit, yy, start_values=[3])


if __name__ == '__main__':
    unittest.main()

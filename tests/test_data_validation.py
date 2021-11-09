from piecewise_regression.data_validation import (
        validate_positive_number,
        validate_boolean,
        validate_positive_integer,
        validate_list_of_numbers,
        validate_non_negative_integer
    )
import unittest

import os
import sys
sys.path.insert(1, os.path.join(sys.path[0], '..'))


class TestDataValidation(unittest.TestCase):

    def test_positive_number(self):

        valid_inputs = [0.111, 3, 6770, 23.1, 7, 0.00001]

        invalid_inputs = [-1, 0, "hi", [1, 1], True, False, None]

        for valid in valid_inputs:
            var_return = validate_positive_number(valid, "test")
            self.assertEqual(var_return, valid)

        for invalid in invalid_inputs:
            test_kwargs = {"var": invalid, "var_name": "test"}
            self.assertRaises(
                ValueError,
                validate_positive_number,
                **test_kwargs)

    def test_boolean(self):

        valid_inputs = [True, False]

        invalid_inputs = [-1, 0, "hi", [1, 1], 22, 1.11]

        for valid in valid_inputs:
            var_return = validate_boolean(valid, "test")
            self.assertEqual(var_return, valid)

        for invalid in invalid_inputs:
            test_kwargs = {"var": invalid, "var_name": "test"}
            self.assertRaises(
                ValueError,
                validate_boolean,
                **test_kwargs)

    def test_positive_integer(self):

        valid_inputs = [1, 2, 31, 6543]

        invalid_inputs = [-1, 0, "hi", [1, 1], 1.11, True, False, None]

        for valid in valid_inputs:
            var_return = validate_positive_integer(valid, "test")
            self.assertEqual(var_return, valid)

        for invalid in invalid_inputs:
            test_kwargs = {"var": invalid, "var_name": "test"}
            self.assertRaises(
                ValueError,
                validate_positive_integer,
                **test_kwargs)

    def test_non_negative_integer(self):

        valid_inputs = [1, 2, 0, 31, 6543]

        invalid_inputs = [-1, "hi", [1, 1], 1.11, True, False, None]

        for valid in valid_inputs:
            var_return = validate_non_negative_integer(valid, "test")
            self.assertEqual(var_return, valid)

        for invalid in invalid_inputs:
            test_kwargs = {"var": invalid, "var_name": "test"}
            self.assertRaises(
                ValueError,
                validate_non_negative_integer,
                **test_kwargs)

    def test_list_of_numbers(self):

        valid_inputs = [
            [1, 1, 1],
            [0.1, 43, 12],
            [2, 1]
            ]

        invalid_inputs = [
            -1, "hi", [1, "h"], 1.11, True, False, None,
            [], [1]]

        for valid in valid_inputs:
            var_return = validate_list_of_numbers(valid, "test", 2)
            self.assertListEqual(list(var_return), list(valid))

        for invalid in invalid_inputs:
            test_kwargs = {
                "var": invalid, "var_name": "test",
                "min_length": 2
                }
            self.assertRaises(
                ValueError,
                validate_list_of_numbers,
                **test_kwargs)


if __name__ == '__main__':
    unittest.main()

import pandas as pd
import unittest

from breakpoint_regression.breakpoint_regression_1_bp import breakpoint_iterate

class TestSegmented(unittest.TestCase):

    def test_on_simple_data_1_bp(self):

    	filename = "data/test_data_simple_1_bp.csv"

    	# Results from running r segmented pacakge on the same data
    	bp_result = 39
    	intercept = 3.19471 
    	alpha_1 = -0.08223 
    	alpha_2 = 0.08410

    	df = pd.read_csv(filename, header=0)
    	print(df)
    	xx = df["x"].to_numpy()
    	yy = df["y"].to_numpy()
    	print(xx, yy)
    	result = breakpoint_iterate(xx, yy, 1)
    	print(result)






#        self.assertEqual(sum([1, 2, 3]), 6, "Should be 6")


if __name__ == '__main__':
    unittest.main()
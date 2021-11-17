from piecewise_regression.main import NextBreakpoints

import numpy as np
import unittest
from importlib.machinery import SourceFileLoader

import os
import sys
sys.path.insert(1, os.path.join(sys.path[0], '..'))

DATA_SOURCE = "tests/data/data.txt"


class TestNextBreakpoint(unittest.TestCase):

    def test_against_muggeo_r_package_data_1(self):
        """
        Muggeo uses slightly different packages and methods etc, so just check
        values are very close, not exact
        Starting from Muggeo's converged breakpoint values, I am iterating
        once, slight change
        The NextBreakpoint class is very vanilla, so in this example is getting
        in some local minima
        """

        data = SourceFileLoader('data', DATA_SOURCE).load_module()

        xx = np.array(data.MUGGEO_1_XX)
        yy = np.array(data.MUGGEO_1_YY)

        # Choose some bps values from Muggeo converged values
        bps = np.array([7.152, 12.161])

        next_fit = NextBreakpoints(xx, yy, bps)

        # Check statistics from breakpoints etc found by Muggeo
        # 1. MUggeo rss for these brreakpoints are
        muggeo_rss = 190.02

        self.assertAlmostEqual(
            muggeo_rss, next_fit.residual_sum_squares, places=1)

        muggeo_c = 100.64655
        muggeo_c_se = 0.23081
        muggeo_c_t = 436.07

        muggeo_alpha = -4.18526
        muggeo_alpha_se = 0.05583
        muggeo_alpha_t = -74.97

        muggeo_beta1 = -3.65462
        muggeo_beta1_se = 0.11405
        muggeo_beta1_t = -32.05

        muggeo_beta2 = 3.81336
        muggeo_beta2_se = 0.11405
        muggeo_beta2_t = 34.45

        muggeo_bp1 = 7.152
        muggeo_bp1_se = 0.101

        muggeo_bp2 = 12.161
        muggeo_bp2_se = 0.095

        estimates = next_fit.estimates

        self.assertAlmostEqual(
            muggeo_c, estimates["const"]["estimate"], places=1)
        self.assertAlmostEqual(
            muggeo_alpha, estimates["alpha1"]["estimate"], places=1)
        self.assertAlmostEqual(
            muggeo_beta1, estimates["beta1"]["estimate"], places=1)
        self.assertAlmostEqual(
            muggeo_beta2, estimates["beta2"]["estimate"], places=1)
        self.assertAlmostEqual(
            muggeo_bp1, estimates["breakpoint1"]["estimate"], places=1)
        self.assertAlmostEqual(
            muggeo_bp2, estimates["breakpoint2"]["estimate"], places=1)

        self.assertAlmostEqual(muggeo_c_se, estimates["const"]["se"], places=1)
        self.assertAlmostEqual(
            muggeo_alpha_se, estimates["alpha1"]["se"], places=1)
        self.assertAlmostEqual(
            muggeo_beta1_se, estimates["beta1"]["se"], places=1)
        self.assertAlmostEqual(
            muggeo_beta2_se, estimates["beta2"]["se"], places=1)
        self.assertAlmostEqual(
            muggeo_bp1_se, estimates["breakpoint1"]["se"], places=1)
        self.assertAlmostEqual(
            muggeo_bp2_se, estimates["breakpoint2"]["se"], places=1)

        self.assertAlmostEqual(
            muggeo_c_t, estimates["const"]["t_stat"], delta=1)
        self.assertAlmostEqual(
            muggeo_alpha_t, estimates["alpha1"]["t_stat"], delta=1)
        self.assertAlmostEqual(
            muggeo_beta1_t, estimates["beta1"]["t_stat"], delta=1)
        self.assertAlmostEqual(
            muggeo_beta2_t, estimates["beta2"]["t_stat"], delta=1)

        muggeo_r_squared = 0.9991
        muggeo_adj_r_squared = 0.999

        self.assertAlmostEqual(muggeo_r_squared, next_fit.r_squared, places=2)
        self.assertAlmostEqual(muggeo_adj_r_squared,
                               next_fit.adjusted_r_squared, places=2)

    def test_against_muggeo_r_package_data_2(self):
        """
        Muggeo uses slightly different packages and methods etc, so just check
        values are very close, not exact
        Starting from Muggeo's converged breakpoint values, I am iterating
        once, slight change
        The NextBreakpoint class is very vanilla, so in this example is
        getting in some local minima
        """

        data = SourceFileLoader('data', DATA_SOURCE).load_module()

        xx = np.array(data.MUGGEO_2_XX)
        yy = np.array(data.MUGGEO_2_YY)

        # Choose some bps values from Muggeo converged values
        bps = np.array([1.608])

        next_fit = NextBreakpoints(xx, yy, bps)

        # Check statistics from breakpoints etc found by Muggeo
        # 1. MUggeo rss for these brreakpoints are
        muggeo_rss = 34.70127

        self.assertAlmostEqual(
            muggeo_rss, next_fit.residual_sum_squares, places=1)

        muggeo_c = 4.86206
        muggeo_c_se = 0.31020
        muggeo_c_t = 15.67

        muggeo_alpha = 0.94472
        muggeo_alpha_se = 0.06744
        muggeo_alpha_t = 14.01

        muggeo_beta1 = 1.95719
        muggeo_beta1_se = 0.11008
        muggeo_beta1_t = 17.78

        muggeo_bp1 = 1.608
        muggeo_bp1_se = 0.291

        estimates = next_fit.estimates

        self.assertAlmostEqual(
            muggeo_c, estimates["const"]["estimate"], places=1)
        self.assertAlmostEqual(
            muggeo_alpha, estimates["alpha1"]["estimate"], places=1)
        self.assertAlmostEqual(
            muggeo_beta1, estimates["beta1"]["estimate"], places=1)
        self.assertAlmostEqual(
            muggeo_bp1, estimates["breakpoint1"]["estimate"], places=1)

        self.assertAlmostEqual(muggeo_c_se, estimates["const"]["se"], places=1)
        self.assertAlmostEqual(
            muggeo_alpha_se, estimates["alpha1"]["se"], places=1)
        self.assertAlmostEqual(
            muggeo_beta1_se, estimates["beta1"]["se"], places=1)
        self.assertAlmostEqual(
            muggeo_bp1_se, estimates["breakpoint1"]["se"], places=1)

        self.assertAlmostEqual(
            muggeo_c_t, estimates["const"]["t_stat"], delta=1)
        self.assertAlmostEqual(
            muggeo_alpha_t, estimates["alpha1"]["t_stat"], delta=1)
        self.assertAlmostEqual(
            muggeo_beta1_t, estimates["beta1"]["t_stat"], delta=1)

        muggeo_r_squared = 0.9911
        muggeo_adj_r_squared = 0.9903

        self.assertAlmostEqual(muggeo_r_squared, next_fit.r_squared, places=2)
        self.assertAlmostEqual(muggeo_adj_r_squared,
                               next_fit.adjusted_r_squared, places=2)


if __name__ == '__main__':
    unittest.main()

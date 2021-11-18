from piecewise_regression.main import Fit, Muggeo

import numpy as np
import unittest
import matplotlib.pyplot as plt
from importlib.machinery import SourceFileLoader

import os
import sys
sys.path.insert(1, os.path.join(sys.path[0], '..'))

DATA_SOURCE = "tests/data/data.txt"


class TestFit(unittest.TestCase):

    def test_it_works_with_differnet_options(self):

        xx = np.linspace(0, 1, 100)
        yy = np.linspace(0, 1, 100)

        Fit(xx, yy, n_breakpoints=2, verbose=False, n_boot=3)
        Fit(xx, yy, n_breakpoints=1, verbose=False, n_boot=3)
        Fit(xx, yy, n_breakpoints=4, verbose=False, n_boot=3)
        Fit(xx, yy, n_breakpoints=7, verbose=False, n_boot=3)
        Fit(xx, yy, start_values=[0.1, 0.6], verbose=False)
        Fit(xx, yy, start_values=[2, 3], verbose=False)
        Fit(xx, yy, start_values=[2], verbose=True)

    def test_against_muggeo_r_package_data_1(self):
        """
        Muggeo uses slightly different packages and methods etc, so just check
        values are very close, not exact
        Starting from Muggeo's converged breakpoint values, I am iterating
        once, slight change
        The NextBreakpoint class is very vanilla, so in this example is
        getting in some local minima
        """

        data = SourceFileLoader('data', DATA_SOURCE).load_module()

        xx = np.array(data.MUGGEO_1_XX)
        yy = np.array(data.MUGGEO_1_YY)

        # Choose some bps values from Muggeo converged values
        bps = np.array([7, 13])

        np.random.seed(1)
        fit = Fit(xx, yy, start_values=bps, verbose=False, n_boot=10)

        best_fit = fit.best_muggeo.best_fit

        # Check statistics from breakpoints etc found by Muggeo
        # 1. MUggeo rss for these brreakpoints are
        muggeo_rss = 190.02

        self.assertAlmostEqual(
            muggeo_rss, best_fit.residual_sum_squares, places=0)

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

        estimates = best_fit.estimates

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

        self.assertAlmostEqual(muggeo_r_squared, best_fit.r_squared, places=2)
        self.assertAlmostEqual(muggeo_adj_r_squared,
                               best_fit.adjusted_r_squared, places=2)

        # Estimates from get_results
        estimates = fit.get_results()["estimates"]

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

        fit = Fit(xx, yy, start_values=bps, verbose=False)

        best_fit = fit.best_muggeo.best_fit

        # Check statistics from breakpoints etc found by Muggeo
        # 1. MUggeo rss for these brreakpoints are
        muggeo_rss = 34.70127

        self.assertAlmostEqual(
            muggeo_rss, best_fit.residual_sum_squares, places=1)

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

        estimates = best_fit.estimates

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

        self.assertAlmostEqual(muggeo_r_squared, best_fit.r_squared, places=2)
        self.assertAlmostEqual(muggeo_adj_r_squared,
                               best_fit.adjusted_r_squared, places=2)

    def test_against_muggeo_r_package_data_3(self):
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

        # Set just number of bps
        fit = Fit(xx, yy, n_breakpoints=1, verbose=False)

        best_fit = fit.best_muggeo.best_fit

        # Check statistics from breakpoints etc found by Muggeo
        # 1. MUggeo rss for these brreakpoints are
        muggeo_rss = 34.70127

        self.assertAlmostEqual(
            muggeo_rss, best_fit.residual_sum_squares, places=1)

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

        estimates = best_fit.estimates

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

        self.assertAlmostEqual(muggeo_r_squared, best_fit.r_squared, places=2)
        self.assertAlmostEqual(muggeo_adj_r_squared,
                               best_fit.adjusted_r_squared, places=2)

    def test_n_boot_zero_does_muggeo(self):
        """
        """
        np.random.seed(2)

        data = SourceFileLoader('data', DATA_SOURCE).load_module()

        xx = np.array(data.MUGGEO_1_XX)
        yy = np.array(data.MUGGEO_1_YY)

        # Choose some bps values from Muggeo converged values
        bps = np.array([7, 13])

        fit = Fit(xx, yy, start_values=bps, verbose=False, n_boot=0)

        best_fit = fit.best_muggeo.best_fit

        fit_muggeo = Muggeo(xx, yy, start_values=bps,
                            n_breakpoints=2, verbose=False)
        best_fit_2 = fit_muggeo.best_fit

        self.assertDictEqual(best_fit.estimates, best_fit_2.estimates)

    def test_summary_created(self):

        np.random.seed(2)

        data = SourceFileLoader('data', DATA_SOURCE).load_module()

        xx = np.array(data.MUGGEO_1_XX)
        yy = np.array(data.MUGGEO_1_YY)

        # Choose some bps values from Muggeo converged values
        bps = np.array([7, 13])

        fit = Fit(xx, yy, start_values=bps, verbose=False, n_boot=0)

        summary = fit.summary()
        self.assertIn('Davies', summary)
        self.assertIn('Model Parameters', summary)
        self.assertIn('Degrees of Freedom', summary)
        self.assertIn('Sum of Squares', summary)
        self.assertIn('Converged', summary)
        self.assertIn('Estimate', summary)
        self.assertIn('alphas', summary)

    def test_non_converging(self):

        np.random.seed(1)

        xx = np.linspace(0, 10, 10)
        yy = np.linspace(0, 10, 10)

        fit = Fit(xx, yy, start_values=[1], n_boot=0)

        summary = fit.summary()

        self.assertIn('Algorithm did not converge', summary)

        ax = plt.gca()
        initial_lines_count = len(ax.lines)

        fit.plot_bootstrap_restarting_rss_history()

        ax = plt.gca()

        self.assertGreater(len(ax.lines), initial_lines_count)

        ax = plt.gca()
        initial_lines_count = len(ax.lines)

        fit.plot_bootstrap_restarting_history()

        ax = plt.gca()

        self.assertGreater(len(ax.lines), initial_lines_count)

        # Estimates from get_results
        results = fit.get_results()

        self.assertEqual(results["converged"], False)
        self.assertEqual(results["estimates"], None)
        self.assertEqual(results["bic"], None)
        self.assertEqual(results["rss"], None)


class TestPlots(unittest.TestCase):

    def test_main_plot_happens(self):

        np.random.seed(2)

        data = SourceFileLoader('data', DATA_SOURCE).load_module()

        xx = np.array(data.MUGGEO_1_XX)
        yy = np.array(data.MUGGEO_1_YY)

        # Choose some bps values from Muggeo converged values
        bps = np.array([7, 13])

        print("Here: ", plt.gcf().number)

        fit = Fit(xx, yy, start_values=bps, verbose=False, n_boot=0)

        ax = plt.gca()
        initial_lines_count = len(ax.lines)

        fit.plot()

        ax = plt.gca()

        self.assertEqual(len(ax.lines), 3 + initial_lines_count)

    def test_main_breakpoint_history(self):

        np.random.seed(2)

        data = SourceFileLoader('data', DATA_SOURCE).load_module()

        xx = np.array(data.MUGGEO_1_XX)
        yy = np.array(data.MUGGEO_1_YY)

        # Choose some bps values from Muggeo converged values
        bps = np.array([7, 13])

        print("Here: ", plt.gcf().number)

        fit = Fit(xx, yy, start_values=bps, verbose=False, n_boot=0)

        ax = plt.gca()
        initial_lines_count = len(ax.lines)

        fit.plot_bootstrap_restarting_history()

        ax = plt.gca()

        print(initial_lines_count, len(ax.lines))

        self.assertGreater(len(ax.lines), initial_lines_count)

    def test_main_rss_breakpoint_history(self):

        np.random.seed(2)

        data = SourceFileLoader('data', DATA_SOURCE).load_module()

        xx = np.array(data.MUGGEO_1_XX)
        yy = np.array(data.MUGGEO_1_YY)

        # Choose some bps values from Muggeo converged values
        bps = np.array([7, 13])

        fit = Fit(xx, yy, start_values=bps, verbose=False, n_boot=0)

        ax = plt.gca()
        initial_lines_count = len(ax.lines)

        fit.plot_bootstrap_restarting_rss_history()

        ax = plt.gca()

        self.assertEqual(len(ax.lines), 1 + initial_lines_count)

    def test_main_best_muggeo_breakpoint_history(self):

        np.random.seed(2)

        data = SourceFileLoader('data', DATA_SOURCE).load_module()

        xx = np.array(data.MUGGEO_1_XX)
        yy = np.array(data.MUGGEO_1_YY)

        # Choose some bps values from Muggeo converged values
        bps = np.array([7, 13])

        fit = Fit(xx, yy, start_values=bps, verbose=False, n_boot=0)

        ax = plt.gca()
        initial_lines_count = len(ax.lines)

        fit.plot_best_muggeo_breakpoint_history()

        ax = plt.gca()

        self.assertGreater(len(ax.lines), initial_lines_count)


if __name__ == '__main__':
    unittest.main()

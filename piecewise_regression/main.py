
import warnings

import numpy as np
import scipy.stats
import matplotlib.pyplot as plt
import statsmodels.api as sm

try:
    import piecewise_regression.davies as davies
    import piecewise_regression.r_squared_calc as r_squared_calc
    from piecewise_regression.data_validation import (
        validate_positive_number,
        validate_boolean,
        validate_positive_integer,
        validate_list_of_numbers,
        validate_non_negative_integer
    )
except ImportError:
    import davies
    import r_squared_calc
    from data_validation import (
        validate_positive_number,
        validate_boolean,
        validate_positive_integer,
        validate_list_of_numbers,
        validate_non_negative_integer
    )


class NextBreakpoints:
    """
    One iteration of Muggeo's segmented regression algorithm. Gets
    the next breakpoints.
    Also calculates interesting statistics. This expects data validation
    and error
    handling are done at a higher level.

    :param xx: Data series in x-axis for fitting (same axis as the breakpoints)
    :type xx: list

    :param yy: Data series in y-axis for fitting
    :type yy: list

    :param current_breakpoints: The starting breakpoints for this iteration
    :type current_breakpoints: list

    """

    def __init__(self,
                 # list(float) or numpy(float). REQUIRED. Data series in x-axis
                 xx,
                 # list(float) or numpy(float). REQUIRED. Data series in y-axis
                 yy,
                 # list(float) or numpy(float). REQUIRED. Current values of
                 # breakpoint positions
                 current_breakpoints
                 ):

        # Data validation done at a higher level
        self.xx = xx
        self.yy = yy
        self.current_breakpoints = current_breakpoints
        self.n_breakpoints = len(current_breakpoints)

        self.next_breakpoints = None
        self.raw_params = None
        self.covariance_matrix = None

        self.breakpoint_fit()

        # All estimate data saved in dictionary
        self.estimates = {}

        # Don't really need to do this at this point. But it is very quick and
        # nice to have the record of these as we go
        self.calculate_all_estimates()
        self.calculate_all_standard_errors()
        self.calculate_all_confidence_intervals()
        self.calculate_all_t_stats()

        # R squared etc
        self.residual_sum_squares = None
        self.total_sum_squares = None
        self.r_squared = None
        self.adjusted_r_squared = None
        self.bic = None

        self.calculate_r_squared()
        self.calculate_bayesian_information_criterion()

    def breakpoint_fit(self):
        """
        Fit the linear approximation given the current breakpoint guesses.
        Sets the next breakpoints
        and the params from the fit. The params are of the form [c, a,
        beta_hats, gamma_hats]
        """
        Z = np.array([self.xx])
        # Convert data based on breakpoints

        UU = [(self.xx - bp) * np.heaviside(self.xx - bp, 1)
              for bp in self.current_breakpoints]
        VV = [np.heaviside(self.xx - bp, 1) for bp in self.current_breakpoints]

        Z = np.concatenate((Z, UU, VV))
        Z = Z.T
        Z = sm.add_constant(Z, has_constant='add')

        results = sm.OLS(endog=self.yy, exog=Z).fit()

        self.raw_params = results.params
        self.covariance_matrix = results.cov_params()

        # First two params are a and c in the line equation
        # Beta hats are the next group of params, same length as breakpoints
        beta_hats = results.params[2:2 + len(self.current_breakpoints)]
        # Gamma hats are the last group of params, same length as breakpoints
        gamma_hats = results.params[2 + len(self.current_breakpoints):]
        # The next breakpoints are calculated iteratively
        self.next_breakpoints = self.current_breakpoints \
            - gamma_hats / beta_hats

    def calculate_all_estimates(self):
        """
        Extract estiamtes from the params and saves in self.estimates
        """
        params = self.raw_params

        # Extract the exstimates form the correct locations in the params
        const_estimate = params[0]
        beta_estimates = params[2:self.n_breakpoints + 2]
        breakpoint_estimates = self.next_breakpoints

        self.estimates["const"] = {"estimate": const_estimate}
        for bp_i in range(self.n_breakpoints):
            self.estimates["beta{}".format(
                bp_i + 1)] = {"estimate": beta_estimates[bp_i]}
            self.estimates["breakpoint{}".format(
                bp_i + 1)] = {"estimate": breakpoint_estimates[bp_i]}
        # Also calculate alphas
        for alpha_i in range(self.n_breakpoints + 1):
            alpha_estimate = np.sum(params[1:alpha_i + 2])
            self.estimates["alpha{}".format(
                alpha_i + 1)] = {"estimate": alpha_estimate}

    def get_alpha_standard_errors(self):
        """
        Get the standard errors for the alphas (gradients of segments)
        """
        cov_matrix = self.covariance_matrix

        # Alphas are calculated as the sum of alpha_0 and betas up that part
        # of the regression line
        # The var of each alpha is the sum of the covariance matrix up to it
        # Removing the intercept column and row.
        # var(alpha_k) = var(alpha_1) + sum var(beta_j) +
        #    2*sum_{i=1,j=2}^k *cov(alpha, betas))
        # var(alpha_k) = sum_{i, j} cov(alpha and betas)
        alpha_vars = []
        for alpha_n in range(self.n_breakpoints + 1):
            alpha_cov_matrix = cov_matrix[1:alpha_n + 2, 1:alpha_n + 2]
            alpha_vars.append(np.sum(alpha_cov_matrix))

        alpha_ses = np.sqrt(alpha_vars)
        return alpha_ses

    def get_bp_standard_errors(self):
        """
        Get the standard errors of the breakpoints.
        Considering bp = gamma/beta + bp_0, the
        standard error of the breakpoint estaimtes can be found using the
        ratio/delta method.
        See e.g. Muggeo (2003) for clarification
        """
        cov_matrix = self.covariance_matrix
        params = self.raw_params

        bp_vars = []

        # For each breakpoint, calcaulte the variance of the estimator
        for bp_n in range(self.n_breakpoints):
            beta_index = 2 + bp_n
            gamma_index = 2 + self.n_breakpoints + bp_n

            beta = params[beta_index]
            gamma = params[gamma_index]
            gamma_var = cov_matrix[gamma_index, gamma_index]
            beta_var = cov_matrix[beta_index, beta_index]
            gamma_beta_covar = cov_matrix[beta_index, gamma_index]

            # From Muggeo (2003). The sign before the covariance term is
            # opposite to Muggeo, this is because the gamma is defined with
            # the opposite sign
            # The calculation is equivalent to Muggeos.
            bp_var = (gamma_var + beta_var * (gamma / beta)**2 -
                      2 * (gamma / beta) * gamma_beta_covar) / (beta**2)
            bp_vars.append(bp_var)

        bp_ses = np.sqrt(bp_vars)
        return bp_ses

    def get_const_standard_error(self):
        """
        Get the constant standard error from the covariance matrix
        """
        # Covariance matrix is [c, alpha, betas, gammas]
        # Constant variance is just the top left cell in covariance matrix
        cov_matrix = self.covariance_matrix
        c_var = cov_matrix[0][0]
        return np.sqrt(c_var)

    def get_beta_standard_errors(self):
        """
        Get the beta estimates standard errors from the covariance matrix
        """
        # Covariance matrix is [c, alpha, betas, gammas]
        # Beta variances are along the diagonal of the covariance matrix
        cov_matrix = self.covariance_matrix
        cov_diagonal = np.diagonal(cov_matrix)
        beta_vars = cov_diagonal[2:self.n_breakpoints + 2]
        return np.sqrt(beta_vars)

    def calculate_all_standard_errors(self):
        """
        Calculate standard errors for all the variables of interest.
        Save to the self.estimates dictionary
        """
        const_ses = self.get_const_standard_error()
        self.estimates["const"]["se"] = const_ses

        beta_ses = self.get_beta_standard_errors()
        bp_ses = self.get_bp_standard_errors()
        for bp_i in range(self.n_breakpoints):
            self.estimates["beta{}".format(bp_i + 1)]["se"] = beta_ses[bp_i]
            self.estimates["breakpoint{}".format(
                bp_i + 1)]["se"] = bp_ses[bp_i]
        alpha_ses = self.get_alpha_standard_errors()
        for alpha_i in range(self.n_breakpoints + 1):
            self.estimates["alpha{}".format(
                alpha_i + 1)]["se"] = alpha_ses[alpha_i]

    def calculate_all_confidence_intervals(self):
        """
        Calculate all confidence intervals, based on t-distribution and
        standard errors.
        """
        # Estimates
        dof = len(self.xx) - 2 - 2 * self.n_breakpoints
        t_const = scipy.stats.t.ppf(0.975, dof)

        # Iterate over the estimate dictionary, add confidence intervals
        # to all estimators
        for estimator_name, details in self.estimates.items():
            confidence_interval = (
                details["estimate"] - t_const * details["se"],
                details["estimate"] + t_const * details["se"])
            details["confidence_interval"] = confidence_interval

    def calculate_all_t_stats(self):
        """
        Get t stats for all the estimators
        """
        dof = len(self.xx) - 2 - 2 * self.n_breakpoints
        for estimator_name, details in self.estimates.items():
            # Breakpoint t stats don't make sense
            # Don't exist in the null model - nuisance parameter
            # H_0 isn't bp=0, it's that bp doesn't exist
            if "breakpoint" in estimator_name:
                details["t_stat"] = "-"
                details["p_t"] = "-"
            else:
                t_stat = details["estimate"] / details["se"]
                p_t = scipy.stats.t.sf(np.abs(t_stat), dof) * 2
                details["t_stat"] = t_stat
                if "beta" in estimator_name:
                    details["p_t"] = "-"
                else:
                    details["p_t"] = p_t

    def get_predicted_yy(self):
        """
        Get the model predictions for each of the xx data points
        """
        params = self.raw_params
        breakpoints = self.next_breakpoints
        # Extract what we need from params etc
        intercept_hat = params[0]
        alpha_hat = params[1]
        beta_hats = params[2:2 + len(breakpoints)]

        yy_predicted = []

        yy_predicted = intercept_hat + alpha_hat * self.xx
        for bp_count in range(len(breakpoints)):
            yy_predicted += beta_hats[bp_count] * \
                np.maximum(self.xx - breakpoints[bp_count], 0)

        return yy_predicted

    def calculate_r_squared(self):
        """
        Calculate R squared from the fitted model.
        """
        yy_predicted = self.get_predicted_yy()
        n_params = 2 * self.n_breakpoints + 2

        rss, tss, r2, adjr2 = r_squared_calc.get_r_squared(
            self.yy,
            yy_predicted,
            n_params)

        self.residual_sum_squares = rss
        self.total_sum_squares = tss
        self.r_squared = r2
        self.adjusted_r_squared = adjr2

    def calculate_bayesian_information_criterion(self):
        """
        Calculates the Bayesian Information Criterion of the fitted model.
        Assuming normal noise, uses the standard version for OLS models.
        This should hold for breakpoint regression models, because the BIC
        is based on the likelihood of the data
        given the model. That likelihood function doesn't involve the
        breakpoint values - it just depends on distances
        of the data to the fitted model predictions. Also depends on the error
        in the noise term being constant.
        """
        n = len(self.xx)  # No. data points
        k = 2 + 2 * self.n_breakpoints  # No. model parameters
        rss = self.residual_sum_squares
        self.bic = n * np.log(rss / n) + k * np.log(n)


class Muggeo:
    """
    Muggeo's iterative segmented regression method. This is a simple version.
    Errors are handled at a higher level in the Fit object. See Muggeo (2003).
    If the breakpoints get too close together, or get outside
    (or near the edge) of the data range,
    the algorithm is stopped because this is likely to be a local
    minimum that is difficult to escape, as well
    as possibly generating errors in the iterative procedue.

    :param xx: Data series in x-axis for fitting (same axis as the breakpoints)
    :type xx: list of floats

    :param yy: Data series in y-axis for fitting
    :type yy: list of floats

    :param n_breakpoints: The number of breakpoints to fit
    :type n_breakpoints: positive int

    :param start_values: A list of initial guesses for the breakpoints
    :type start_values: list floats

    :param verbose: If True, prints out updates to the terminal
    :type verbose: bool

    :param max_iterations: How many iterations before stopping if not converged
    :type max_iterations: positive int

    :param tolerance: How close breakpoints from pervious iterations must be
        to consider converged.
    :type tolerance: positive float

    :param min_distance_between_breakpoints: The minimum allowed distance
        between breakpoints, as a proportion of the data range.
    :type min_distance_between_breakpoints: positive float

    :param min_distance_between_breakpoints: The minimum allowed distance from
        the edge of data to a breakpoint, as a proportion of the data range.
    :type min_distance_between_breakpoints: positive float

    """

    def __init__(self,
                 # list(float) or numpy(float). REQUIRED. Data series in x-axis
                 xx,
                 # list(float) or numpy(float). REQUIRED. Data series in y-axis
                 yy,
                 n_breakpoints,  # int. REQUIRED. Number of breakpoints
                 # list(float) or numpy(float). REQUIRED. Initial guesses for
                 # breakpoint positions
                 start_values=None,
                 # Boolean. whether to print progress to terminal.
                 verbose=False,
                 max_iterations=30,  # Positive int. Maximum iterations of
                 # Muggeo algorithm if not converged
                 tolerance=10**-5,  # Positive float. If breakpoints change
                 # less than the tolerance then the algorithm has converged
                 # Positive float. The minimum required distance between
                 # breakpoints, as a proportion of the data range.
                 min_distance_between_breakpoints=0.01,
                 # Positive float. Minimum distance from edge of data to a
                 # breakpoint, as a proportion of the data range.
                 min_distance_to_edge=0.02,
                 ):

        self.verbose = verbose
        if self.verbose:
            print("\nInstantiating Muggeo . . . with start_values = {}".format(
                start_values))

        # validation is done at a higher level
        self.xx = xx
        self.yy = yy

        self.n_breakpoints = n_breakpoints
        self.min_distance_between_breakpoints = \
            min_distance_between_breakpoints
        self.min_distance_to_edge = min_distance_to_edge

        if start_values is None:
            start_values = self._generate_breakpoints()

        self.start_values = self._validate_start_values(start_values)

        self.max_iterations = max_iterations
        self.tolerance = tolerance
        self.verbose = verbose

        self.stop = False

        if self.n_breakpoints != len(self.start_values):
            raise ValueError(
                "n_breakpoints is not the same as the length of start_values")

        self.fit_history = []
        self.best_fit = None
        self.converged = False

        # Records the reason why the algorithm stopped
        self.stop_reason = None

        self.fit()

    def fit(self):
        """
        Runs the breakpoint iterative procedure
        """
        if self.verbose:
            print("Running Muggeo's iterative algorithm . . . ")
        while not self.stop:
            # Do the fit
            # Get the current breakpoints. If first fit then use start_vakues,
            # otherwise use fit history
            if len(self.fit_history) == 0:
                current_breakpoints = self.start_values
            else:
                current_breakpoints = self.fit_history[-1].next_breakpoints

            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    IteratedFit = NextBreakpoints(
                        self.xx, self.yy, current_breakpoints)
                    self.fit_history.append(IteratedFit)
                    self.stop_or_not()
                    if self.verbose:
                        print("Next breakpoints are ",
                              IteratedFit.next_breakpoints)
            except Exception as e:
                self.stop = True
                self.stop_reason = "Error encountered: " + str(e)

        # Select the best fit from the fit history by finding the smallest rss
        if len(self.fit_history) > 0:
            self.best_fit = min(
                self.fit_history, key=lambda x: x.residual_sum_squares)

    def stop_or_not(self):
        """
        Test to see if the iterative procedure should stop.
        Stop if it's converged, if max_iterations reached, of the breakpoints
        fall into values that
        are too close together or outside of the allowed rance.
        """
        # Stop if maximum iterations reached
        if len(self.fit_history) > self.max_iterations:
            self.stop_reason = "Algorithm stopped as max iterations reached"
            self.stop = True

        if not self._are_breakpoint_values_far_apart(
                self.fit_history[-1].next_breakpoints):
            self.stop_reason = "Breakpoint values too close together"
            self.stop = True

        if not self._are_breakpoint_values_within_range(
                self.fit_history[-1].next_breakpoints):
            self.stop_reason = "Breakpoint values outside range"
            self.stop = True

        # Stop if tolerance reached - small change between this and last
        # breakpoints
        if len(self.fit_history) > 1:
            breakpoint_differences = self.fit_history[-2].next_breakpoints - \
                self.fit_history[-1].next_breakpoints
            if np.max(np.abs(breakpoint_differences)) <= self.tolerance:
                self.stop_reason = "Algorithm converged on breakpoint values"
                self.stop = True
                self.converged = True

        # Stop if the algorithm is iterating back to previous values, within
        # tolerance
        if len(self.fit_history) > 2:
            breakpoint_two_step_differences = \
                self.fit_history[-3].next_breakpoints - \
                self.fit_history[-1].next_breakpoints
            if np.max(np.abs(breakpoint_two_step_differences)) \
                    <= self.tolerance:
                self.stop_reason = "Algorithm converged on breakpoint values"
                self.stop = True
                self.converged = True

    def _validate_start_values(self, start_values):
        """
        Validate the breakpoint start_values.
        Should be a list of numbers, or a numpy list of numbers.
        They should also be not too close to the edge of the data, and not too
        close together - to avoid the Muggeo algorithm diverging

        :param start_values: A list of initial guesses for the breakpoints
        :type start_values: list floats

        """

        start_values = validate_list_of_numbers(
            start_values, "start_values", min_length=1)

        if not self._are_breakpoint_values_within_range(start_values):
            self.stop = True
            self.stop_reason = ("start_values are not within allowed range."
                                "Try changing min_distance_to_edge")

        if not self._are_breakpoint_values_far_apart(start_values):
            self.stop = True
            self.stop_reason = ("start_values are too close together."
                                "Try changing "
                                "min_distance_between_breakpoints")
        return start_values

    def _are_breakpoint_values_within_range(self, breakpoints):
        """
        Check that breakpoints are self.min_distance_to_edge * range away for
        the edge of the data in x

        :param breakpoints: A list of breakpoints
        :type breakpoints: list floats

        """
        min_allowed_bp = np.quantile(self.xx, self.min_distance_to_edge)
        max_allowed_bp = np.quantile(self.xx, 1 - self.min_distance_to_edge)

        for bp in breakpoints:
            if bp <= min_allowed_bp or bp >= max_allowed_bp:
                return False
        return True

    def _are_breakpoint_values_far_apart(self, breakpoints):
        """
        Check if breakpoint values are
        self.min_distance_between_breakpoints*range away from each other

        :param breakpoints: A list of breakpoints
        :type breakpoints: list floats

        """
        min_distance = np.diff(np.sort(breakpoints))

        # numpy ptp gives the range of the data, closeness relative to that
        min_distance_allowed = self.min_distance_between_breakpoints * \
            np.ptp(self.xx)

        if (min_distance <= min_distance_allowed).any():
            return False
        return True

    def _generate_breakpoints(self):
        """
        Randomly generate some breakpoint values
        """
        # Get breakpoints within allowed range
        min_allowed_bp = np.quantile(self.xx, self.min_distance_to_edge)
        max_allowed_bp = np.quantile(self.xx, 1 - self.min_distance_to_edge)
        start_values = np.random.uniform(
            low=min_allowed_bp, high=max_allowed_bp, size=self.n_breakpoints)
        if self.verbose:
            print("Generating some random breakpoints: ", start_values)
        return start_values


class Fit:
    """
    Fit a peicewise (segmented) regression model to data.
    Uses bootstrap restarting to avoid local minima.
    Requires either n_breakpoints of start_values.
    if no start_vaues are given, they are instead uniformly randomly
    generated across range of data.
    This is the main user facing object and input data is validated mainly
    at this level.

    :param xx: Data series in x-axis for fitting (same axis as the breakpoints)
    :type xx: list of floats

    :param yy: Data series in y-axis for fitting.
    :type yy: list of floats

    :param n_breakpoints: The number of breakpoints to fit.
    :type n_breakpoints: positive int

    :param start_values: A list of initial guesses for the breakpoints.
    :type start_values: list floats

    :param n_boot: How many times to run the bootstrap restarting procedure.
        Set to zero for no bootstrap restarting.
    :type n_boot: non-negative int

    :param verbose: If True, prints out updates to the terminal.
    :type verbose: bool

    :param max_iterations: How many iterations before stopping if not
        converged, in the Muggeo iterative procedure.
    :type max_iterations: positive int

    :param tolerance: How close breakpoints from pervious iterations must
        be to consider converged.
    :type tolerance: positive float

    :param min_distance_between_breakpoints: The minimum allowed distance
        between breakpoints, as a proportion of the data range.
    :type min_distance_between_breakpoints: positive float

    :param min_distance_between_breakpoints: The minimum allowed distance from
        the edge of data to a breakpoint, as a proportion of the data range.
    :type min_distance_between_breakpoints: positive float

    """

    def __init__(self,
                 # list(float) or numpy(float). REQUIRED Data series in x-axis
                 xx,
                 # list(float) or numpy(float). REQUIRED Data series in y-axis
                 yy,
                 # list(float) or numpy(float). Initial guesses for breakpoint
                 # positions
                 start_values=None,
                 # int. If not start_values, the number of breakpoints to fit.
                 # REQUIRED if no start_values
                 n_breakpoints=None,
                 n_boot=100,  # Positive int. The number of times to run the
                 # bootstrap restarting. n_boot=0 runs the Muggeo algorithm
                 # with no bootstrap
                 # Boolean. whether to print progress to terminal.
                 verbose=False,
                 max_iterations=30,  # Positive int. Maximum iterations of
                 # Muggeo algorithm if not converged
                 tolerance=10**-5,  # Positive float. If breakpoints change
                 # less than the tolerance then the algorithm has converged
                 # Positive float. The minimum required distance between
                 # breakpoints, as a proportion of the data range.
                 min_distance_between_breakpoints=0.01,
                 # Positive float. Minimum distance from edge of data to a
                 # breakpoint, as a proportion of the data range.
                 min_distance_to_edge=0.02,
                 ):

        self.verbose = validate_boolean(verbose, "verbose")

        if self.verbose:
            print("\nInstantiating Fit . . . ")

        # Validate all input data
        self.xx = validate_list_of_numbers(xx, "xx", min_length=3)
        self.yy = validate_list_of_numbers(yy, "yy", min_length=3)

        if len(self.yy) != len(self.xx):
            raise ValueError("x and y data series must be the same size")

        self.n_boot = validate_non_negative_integer(n_boot, "n_boot")
        self.max_iterations = validate_positive_integer(
            max_iterations, "max_iterations")
        self.tolerance = validate_positive_number(tolerance, "tolerance")
        self.min_distance_between_breakpoints = validate_positive_number(
            min_distance_between_breakpoints,
            "min_distance_between_breakpoints")
        self.min_distance_to_edge = validate_positive_number(
            min_distance_to_edge, "min_distance_to_edge")

        # We need either start_values or n_breakpoints
        if n_breakpoints is None:
            self.n_breakpoints = None
        else:
            self.n_breakpoints = validate_positive_integer(
                n_breakpoints, "n_breakpoints")

        if start_values is None:
            self.start_values = None
        else:
            self.start_values = validate_list_of_numbers(
                start_values, "start_values", min_length=1)
            self.n_breakpoints = len(self.start_values)

        if start_values is None and n_breakpoints is None:
            raise ValueError(
                "Fit algorithm requires either start_values or n_breakpoints")

        self.bootstrap_history = []
        self.best_muggeo = None
        self.stop = False

        self.bootstrap_restarting()

        self.davies = davies.davies_test(self.xx, self.yy)

    def get_results(self):
        """
        Return a small dictionary with key results form the fit.
        Useful for using this code in a larger analysis. E.g. ModelSelection
        """
        results = {
            "davies": self.davies,
        }

        if self.best_muggeo:
            results["estimates"] = self.best_muggeo.best_fit.estimates
            results["bic"] = self.best_muggeo.best_fit.bic
            results["rss"] = self.best_muggeo.best_fit.residual_sum_squares
            results["converged"] = True
        else:
            results["converged"] = False
            results["estimates"] = None
            results["bic"] = None
            results["rss"] = None
        return results

    def bootstrap_restarting(self):
        """
        The main fitting algorithm. Begins by doing a fit based on
        Muggeo's algorithm.
        if n_boot = 0 we stop there. Otherwise we do some bootstrap restarting.
        Bootstrap Restarting escapes local minima.
        Each bootstrap restart:
        - We take the best current breakpoints, and get new data by running a
            non-parametric bootstrap by resampling data.
        - Then run a Muggeo fit on the new data and best current breakpoints.
            This gives new breakpoint values.
        - Then run a Muggeo fit again with the original data and these new
            breakpoint values.
        - Throughout, keep track of the history of fits and the best_muggeo
            fit that converged - defined as the lowest residual sum of squares.
        """
        min_d_between_bps = self.min_distance_between_breakpoints
        muggeo_fit = Muggeo(
            xx=self.xx,
            yy=self.yy,
            start_values=self.start_values,
            n_breakpoints=self.n_breakpoints,
            max_iterations=self.max_iterations,
            tolerance=self.tolerance,
            verbose=self.verbose,
            min_distance_between_breakpoints=min_d_between_bps,
            min_distance_to_edge=self.min_distance_to_edge)

        self.bootstrap_history.append(muggeo_fit)

        # best_muggeo is the best converged muggeo
        if muggeo_fit.converged:
            self.best_muggeo = muggeo_fit

        # Iterate bootstraps
        for i in range(self.n_boot):

            # Best breakpoints are either from best converged muggeo so far,
            # or start values, or randomly generated
            if self.best_muggeo and np.random.uniform() < 0.5:
                best_bps = self.best_muggeo.best_fit.next_breakpoints
            else:
                best_bps = self.start_values

            # Get some new breakpoint values from a bootstrapped fit
            # Non parametric bootstrap by resampling from data
            xx_boot, yy_boot = self.bootstrap_data(self.xx, self.yy)
            bootstrap_fit = Muggeo(
                xx_boot, yy_boot,
                start_values=best_bps,
                n_breakpoints=self.n_breakpoints,
                max_iterations=self.max_iterations,
                tolerance=self.tolerance,
                verbose=self.verbose,
                min_distance_between_breakpoints=min_d_between_bps,
                min_distance_to_edge=self.min_distance_to_edge)
            if bootstrap_fit.converged:
                bootstrap_bps = bootstrap_fit.best_fit.next_breakpoints
            else:
                # Give it something - even though these breakpoints are
                # already run
                # prefer this to using breaks in the for loop
                bootstrap_bps = best_bps

            # Do a new fit with the new breakpoint values
            next_muggeo = Muggeo(
                self.xx,
                self.yy,
                start_values=bootstrap_bps,
                n_breakpoints=self.n_breakpoints,
                max_iterations=self.max_iterations,
                tolerance=self.tolerance,
                verbose=self.verbose,
                min_distance_between_breakpoints=min_d_between_bps,
                min_distance_to_edge=self.min_distance_to_edge)
            self.bootstrap_history.append(next_muggeo)

            # If we get a converged answer, see if this new fit is the best
            if next_muggeo.converged:
                # If there is already a converged best_muggeo, see if this one
                # is better
                if self.best_muggeo:
                    if next_muggeo.best_fit.residual_sum_squares \
                            < self.best_muggeo.best_fit.residual_sum_squares:
                        self.best_muggeo = next_muggeo
                # If there is not already a converged best_muggeo, use this
                # fit instead
                else:
                    self.best_muggeo = next_muggeo

    def bootstrap_data(self, xx, yy):
        """
        Non parametric bootstrap, randomly sample data points with replacement.
        Return bootstrapped data of same length as oriignal data.

        :param xx: Data series in x-axis.
        :type xx: list of floats

        :param yy: Data series in y-axis.
        :type yy: list of floats

        """
        n = len(xx)
        # Get bootstrap samples as array index locations
        boot_indices = np.random.choice(n, size=n, replace=True)

        xx_boot = xx[boot_indices]
        yy_boot = yy[boot_indices]
        return xx_boot, yy_boot

    def plot_data(self, **kwargs):
        """
        Plot the data as a scatter plot.
        Passes any kwargs to the matplotlib scatter function, e.g. color="red".

        """
        plt.scatter(self.xx, self.yy, **kwargs)

    def plot_fit(self, **kwargs):
        """
        Plot the fitted model as a series of straight lines.
        Passes any kwargs to the matplotlib plot function, e.g. color="red".

        """
        if not self.best_muggeo:
            print("Algorithm didn't converge. No fit to plot.")
        else:
            # Get the final results from the fitted model variables
            # Params are in terms of [intercept, alpha, betas, gammas]
            final_params = self.best_muggeo.best_fit.raw_params
            breakpoints = self.best_muggeo.best_fit.next_breakpoints

            # Extract what we need from params etc
            intercept_hat = final_params[0]
            alpha_hat = final_params[1]
            beta_hats = final_params[2:2 + len(breakpoints)]

            xx_plot = np.linspace(min(self.xx), max(self.xx), 100)

            # Build the fit plot segment by segment. Betas are defined as
            # difference in gradient from previous section
            yy_plot = intercept_hat + alpha_hat * xx_plot
            for bp_count in range(len(breakpoints)):
                yy_plot += beta_hats[bp_count] * \
                    np.maximum(xx_plot - breakpoints[bp_count], 0)

            plt.plot(xx_plot, yy_plot, **kwargs)

    def plot_breakpoints(self, **kwargs):
        """
        Plot the breakpoint locations as vertical lines.
        Passes kwargs to the matplotlib function, e.g. color="red".

        """
        if not self.best_muggeo:
            print("Algorithm didn't converge. No breakpoints to plot")
        else:
            breakpoints = self.best_muggeo.best_fit.next_breakpoints

            for bp in breakpoints:
                plt.axvline(bp, **kwargs)

    def plot_breakpoint_confidence_intervals(self, **kwargs):
        """
        Plot the breakpoint confidence intervals as vertical shaded regions.
        Passes kwargs to the matplotlib function, e.g. color="red".

        """

        if not self.best_muggeo:
            print("Algorithm didn't converge. No breakpoint estimates to plot")
        else:
            estimates = self.best_muggeo.best_fit.estimates

            for bp_i in range(self.best_muggeo.n_breakpoints):
                bp_ci = estimates["breakpoint{}".format(
                    bp_i + 1)]["confidence_interval"]
                plt.axvspan(bp_ci[0], bp_ci[1], alpha=0.1)

    def plot_best_muggeo_breakpoint_history(self, **kwargs):
        """
        Plot the history of the breakpoints as they iterate.
        History of the best_muggeo fit.

        """
        if not self.best_muggeo:
            print("Algorithm didn't converge. No meaningful history to plot")
        else:
            # Get the data from the best_muggeo in a form for plotting
            breakpoint_history = [self.best_muggeo.start_values]
            for fit_details in self.best_muggeo.fit_history:
                breakpoint_history.append(fit_details.next_breakpoints)
            breakpoint_history_series = zip(*breakpoint_history)

            # Plot a history for each breakpoint
            count = 0
            for bh in breakpoint_history_series:
                count += 1
                plt.plot(range(0, len(bh)), bh,
                         label="Breakpoint {}".format(count), **kwargs)
                plt.xlabel("Muggeo Iteration")
                plt.ylabel("Breakpoint")

    def plot_bootstrap_restarting_history(self, **kwargs):
        """
        Plot the history of the breakpoint values as they iterate.
        History of the bootstrap restarting procedure.

        """
        if not self.best_muggeo:
            print("Algorithm didn't converge. Plotting breakpoint history")

        # Get the data from the best_muggeo in a form for plotting

        breakpoint_history = [self.start_values]

        for muggeo_fit in self.bootstrap_history:
            if muggeo_fit.best_fit:
                breakpoint_history.append(muggeo_fit.best_fit.next_breakpoints)
            else:
                breakpoint_history.append(None)

        breakpoint_history_series = zip(*breakpoint_history)

        # Plot a history for each breakpoint
        count = 0
        for bh in breakpoint_history_series:
            count += 1
            plt.plot(range(0, len(bh)), bh,
                     label="Breakpoint {}".format(count), **kwargs)
            plt.xlabel("Bootstrap Iteration")
            plt.ylabel("Breakpoint")

    def plot_bootstrap_restarting_rss_history(self, **kwargs):
        """
        Plot the history of the residual sum of squares.
        History of the bootstrap restarting algorithm.

        """
        if not self.best_muggeo:
            print("Algorithm didn't converge. Plotting rss history anyway")

        # Get the data from the best_muggeo in a form for plotting

        rss_history = []

        for muggeo_fit in self.bootstrap_history:
            if muggeo_fit.best_fit:
                rss_history.append(muggeo_fit.best_fit.residual_sum_squares)
            else:
                rss_history.append(None)

        plt.plot(range(1, len(rss_history) + 1), rss_history, **kwargs)
        plt.xlabel("Bootstrap Iteration")
        plt.ylabel("Residual Sum of Squares")

    def plot(self):
        """
        Plot the full fit including the data, fitted model, breakpoint
        positions and breakpoint confidence intervals.
        Doesn't allow control of matplotlib kwargs for style changes.
        """
        self.plot_data()
        self.plot_fit()
        self.plot_breakpoints()
        self.plot_breakpoint_confidence_intervals()

    def summary(self):
        """
        Print a summary of the best fit, along the lines of the summary
        given by python's statsmodels OLS fit.
        """

        if not self.best_muggeo:
            summary = (
                "Algorithm did not converge. Try different n_breakpoints, "
                "different start_values, or start_values=None\n")

            summary += "Summary of why the algorithm did not converge:\n"
            run_count = 1
            for muggeo_fit in self.bootstrap_history:
                summary += "Run {}: {} \n".format(
                    run_count,
                    muggeo_fit.stop_reason)
                run_count += 1
            print(summary)

            return summary

        else:
            header = "\n{:^70}\n".format("Breakpoint Regression Results")

            line_length = 100
            double_line = "=" * line_length + "\n"
            single_line = "-" * line_length + "\n"

            # Overview
            n_obs = len(self.xx)
            n_model_params = 2 + 2 * self.n_breakpoints
            dof = n_obs - n_model_params
            no_obs_text = "{:<20} {:>20}\n".format("No. Observations", n_obs)
            no_model_parameters_text = "{:<20} {:>20}\n".format(
                "No. Model Parameters", n_model_params)
            dof_text = "{:<20} {:>20}\n".format("Degrees of Freedom", dof)
            rss_text = "{:<20} {:>20.6}\n".format(
                "Res. Sum of Squares",
                self.best_muggeo.best_fit.residual_sum_squares)
            tss_text = "{:<20} {:>20.6}\n".format(
                "Total Sum of Squares",
                self.best_muggeo.best_fit.total_sum_squares)
            r_2_text = "{:<20} {:>20.6f}\n".format(
                "R Squared", self.best_muggeo.best_fit.r_squared)
            adj_r_2_text = "{:<20} {:>20.6f}\n".format(
                "Adjusted R Squared",
                self.best_muggeo.best_fit.adjusted_r_squared)
            converged_text = "{:<20} {:>20s}\n".format(
                "Converged: ", str(self.best_muggeo.converged))

            overview = double_line + no_obs_text + no_model_parameters_text + \
                dof_text + rss_text + tss_text + r_2_text + adj_r_2_text + \
                converged_text + double_line

            # Table of results
            table_header_template = ("{:<15} {:>12} {:>12} {:>12} "
                                     "{:>12} {:>12} {:>12}\n")

            table_header = table_header_template.format(
                "", "Estimate", "Std Err", "t", "P>|t|", "[0.025", "0.975]")

            table_row_template = ("{:<15} {:>12.6} {:>12.3} {:>12.5} "
                                  " {:>12.3} {:>12.5} {:>12.5}\n")

            table_contents = ""

            beta_names = ["beta{}".format(i + 1)
                          for i in range(self.n_breakpoints)]
            bp_names = ["breakpoint{}".format(i + 1)
                        for i in range(self.n_breakpoints)]

            model_estimator_names = ["const", "alpha1"] + beta_names + bp_names

            estimates = self.best_muggeo.best_fit.estimates

            for est_name in model_estimator_names:
                estimator_row = table_row_template.format(
                    est_name,
                    estimates[est_name]["estimate"],
                    estimates[est_name]["se"],
                    estimates[est_name]["t_stat"],
                    estimates[est_name]["p_t"],
                    estimates[est_name]["confidence_interval"][0],
                    estimates[est_name]["confidence_interval"][1])
                table_contents += estimator_row

            table_contents += single_line

            table_contents += (
                "These alphas(gradients of segments) are estimated"
                "from betas(change in gradient)\n")

            alpha_names = ["alpha{}".format(
                alpha_i + 1) for alpha_i in range(1, self.n_breakpoints + 1)]

            table_contents += single_line

            for est_name in alpha_names:
                estimator_row = table_row_template.format(
                    est_name, estimates[est_name]["estimate"],
                    estimates[est_name]["se"],
                    estimates[est_name]["t_stat"], estimates[est_name]["p_t"],
                    estimates[est_name]["confidence_interval"][0],
                    estimates[est_name]["confidence_interval"][1])
                table_contents += estimator_row

            table_contents += double_line

            table = double_line + table_header + single_line + table_contents

            davies_result = (
                "Davies test for existence of at least "
                "1 breakpoint: p={:.6} (e.g. p<0.05 means reject null "
                "hypothesis of no breakpoints "
                " at 5% significance)".format(self.davies))

            summary = header + overview + table + davies_result + "\n\n"

            print(summary)

            return summary


if __name__ == "__main__":
    pass

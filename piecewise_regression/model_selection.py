

import numpy as np
import statsmodels.api as sm

try:
    import piecewise_regression.r_squared_calc as r_squared_calc
    from piecewise_regression.main import Fit
except ImportError:
    import r_squared_calc
    from main import Fit


class ModelSelection:
    """
    Experimental - uses simple BIC based on simple linear model.
    """

    def __init__(
            self, xx, yy, max_breakpoints=10, n_boot=100,
            max_iterations=30, tolerance=10**-5,
            min_distance_between_breakpoints=0.01, min_distance_to_edge=0.02,
            verbose=True):

        # The actual fit model objects
        self.models = []
        # The model summary data
        self.model_summaries = []

        self.stop = False

        if verbose:
            print("Running fit with n_breakpoint = 0 . . ")

        self.no_breakpoint_fit(xx, yy)

        min_d_between_bps = min_distance_between_breakpoints
        for k in range(1, max_breakpoints + 1):
            if verbose:
                print("Running fit with n_breakpoint = {} . . ".format(k))
            bootstrapped_fit = Fit(
                xx, yy, n_breakpoints=k, verbose=False,
                n_boot=n_boot, max_iterations=max_iterations,
                tolerance=tolerance,
                min_distance_between_breakpoints=min_d_between_bps,
                min_distance_to_edge=min_distance_to_edge)
            fit_summary = bootstrapped_fit.get_results()
            fit_summary["n_breakpoints"] = k
            self.model_summaries.append(fit_summary)
            self.models.append(bootstrapped_fit)

        self.summary()

    def summary(self):

        header = "\n{:^70}\n".format("Breakpoint Model Comparision Results")

        line_length = 100
        double_line = "=" * line_length + "\n"
        single_line = "-" * line_length + "\n"

        table_header_template = "{:<15} {:>12} {:>12} {:>12} \n"
        table_header = table_header_template.format(
            "n_breakpoints", "BIC", "converged", "RSS")
        table_row_template = "{:<15} {:>12.5} {:>12} {:>12.5} \n"

        table_contents = header
        table_contents += double_line

        table_contents += table_header
        table_contents += single_line

        for model_summary in self.model_summaries:

            if model_summary["converged"]:
                model_row = table_row_template.format(
                    model_summary["n_breakpoints"], 
                    model_summary["bic"],
                    str(model_summary["converged"]), 
                    model_summary["rss"])
            else:
                model_row = table_row_template.format(
                    model_summary["n_breakpoints"], "", 
                    str(model_summary["converged"]), "")

            table_contents += model_row

        print(table_contents)

        print("Min BIC (Bayesian Information Criterion) suggests best model")

    def no_breakpoint_fit(self, xx, yy):

        Z = np.array([xx])
        Z = Z.T
        Z = sm.add_constant(Z, has_constant='add')
        # Basic OLS fit
        results = sm.OLS(endog=np.array(yy), exog=Z).fit()

        # get the predicted values
        ff = [(results.params[0] + results.params[1] * x) for x in xx]

        # Get Rss
        rss, tss, r_2, adjusted_r_2 = r_squared_calc.get_r_squared(
            yy, ff, n_params=2)

        # Calcualte BIC
        n = len(xx)  # No. data points
        k = 2  # No. parameters
        bic = n * np.log(rss / n) + k * np.log(n)

        fit_data = {
            "bic": bic,
            "n_breakpoints": 0,
            "estimates": {},
            "converged": True,
            "rss": rss
        }

        fit_data["estimates"]["const"] = results.params[0]
        fit_data["estimates"]["alpha1"] = results.params[1]

        self.model_summaries.append(fit_data)


if __name__ == "__main__":
    
    pass



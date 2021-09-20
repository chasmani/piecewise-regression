
API
================================================

Main
----------
The main module includes the Fit function, which runs the bootstrap restarting algorithm.

.. automodule:: piecewise_regression.main
	:members:

Model selection
---------------------
The model selection module is experimental. It compares models with different `n_breakpoints` using the Bayesian Information Criterion.

.. automodule:: piecewise_regression.model_selection
	:members:

Davies test
----------------------
Implements the Davies hypothesis test for existence of at least one breakpoint.

.. automodule:: piecewise_regression.davies
	:members:


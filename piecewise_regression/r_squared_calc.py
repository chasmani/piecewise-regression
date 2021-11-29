
import numpy as np


def get_r_squared(yy, ff, n_params):

    n_data = len(yy)
    yy_mean = np.mean(yy)

    # Calculate residual and total sum of squares
    residual_sum_squares = 0
    total_sum_squares = 0
    for i in range(n_data):
        residual_sum_squares += (yy[i] - ff[i])**2
        total_sum_squares += (yy[i] - yy_mean)**2

    # R Squares
    r_squared = 1 - residual_sum_squares / total_sum_squares

    # Adjusted R squared
    adj_r_squared = 1 - (1 - r_squared) * \
        (n_data - 1) / (n_data - n_params - 1)

    return residual_sum_squares, total_sum_squares, r_squared, adj_r_squared

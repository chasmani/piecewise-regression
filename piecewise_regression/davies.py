
import numpy as np
import scipy
import math


def get_test_statistic_wald(xx, yy, theta):


    import statsmodels.api as sm

    Z = np.array([xx])

    UU = [(xx - theta ) * np.heaviside(xx-theta, 1)]
    VV = [np.heaviside(xx - theta, 1)]

    Z = np.concatenate((Z, UU, VV))
    Z = Z.T
    Z = sm.add_constant(Z, has_constant='add')

    results = sm.OLS(endog=yy, exog=Z).fit()

    beta_hat = results.params[2]
    se_beta_hat = results.bse[2]
    return beta_hat/se_beta_hat


def davies_test(xx, yy, k=10, alternative="two_sided"):
    """
    Significance test for the existence of a breakpoint
    Null hypothesis is that there is no breakpoint, or that the change in
    gradient is zero.
    Alternative hypothesis is that there is a breakpoint, with a non-zero
    change in gradient.
    The change is gradient is a function of the breakpoint position.
    The breakpoint posiition is a nuisannce parameter that only exists in the
    alternative hypothesis.
    Based on Davies (1987), "Hypothesis Testing when a nuisance parameter is
    present only under the alternative".

    :param xx: Data series in x-axis (same axis as the breakpoints).
    :type xx: list of floats

    :param yy: Data series in y-axis.
    :type yy: list of floats

    :param k: A control parameter that determines the number of points to
        consider within the xx range.
    :type k: int

    :param alternative: Whether to consider a two-sided hypothesis test,
        or a one sided test with change of gradient greater or less than zero.
        For existence of a breakpoint, use "two-sided".
    :type alternative: str. One of "two_sided", "less", "greater"

    """
    # Centre the x values - makes no difference to existence of a breakpoint
    # The Davies test has this as an assumption
    xx_davies = xx - np.mean(xx)
    yy_davies = yy

    # As in Muggeo's R package "segmented", cut from second to second to last
    # data point
    # Need more data in the xx than in [L,U] for the test to work
    # Take off points form edge to be conservative 
    L = xx_davies[2]
    U = xx_davies[-3]

    # More thetas is better
    thetas = np.linspace(L, U, k)
    # For each value of theta, compute a test statistic
    test_stats = []
    for theta in thetas:
        test_stat = get_test_statistic_wald(xx_davies, yy_davies, theta)
        test_stats.append(test_stat)
    if alternative == "two_sided":
        # Two sided test, M as defined by Davies
        M = np.max(np.abs(test_stats))
    elif alternative == "less":
        M = np.abs(np.min(test_stats))
    elif alternative == "greater":
        M = np.max(test_stats)

    # Use formulas from Davies
    V = 0
    for i in range(len(thetas) - 1):
        V += np.abs(test_stats[i + 1] - test_stats[i])

    p = scipy.stats.norm.cdf(-M) + V * np.exp(-.5 * M **
                                              2) * 1 / (np.sqrt(8 * math.pi))

    if alternative == "two_sided":
        return p * 2
    else:
        return p


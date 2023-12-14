
import os
import sys
sys.path.insert(1, os.path.join(sys.path[0], '..'))


from piecewise_regression.davies import davies_test
import numpy as np



def check_p_values(alternative="two_sided", breakpoint_1=0, beta_1=0, breakpoint_2=0, beta_2=0, noise_scale=1):

    p_count_20 = 0
    p_count_10 = 0
    p_count_5 = 0
    p_count_2 = 0
    p_count_1 = 0

    p_values = []

    sample_size = 1000

    for seed in range(sample_size):
        np.random.seed(seed)

        xx_bp, yy_bp = generate_data(breakpoint_1, beta_1, breakpoint_2, beta_2, noise_scale)

        p = davies_test(xx_bp, yy_bp, alternative=alternative)

        if p < 0.2:
            p_count_20 += 1
        if p < 0.05:
            p_count_5 += 1
        if p < 0.02:
            p_count_2 += 1
        if p < 0.01:
            p_count_1 += 1

        p_values.append(p)

    # plt.hist(p_values)
    # plt.show()

    print("P values from empirical testing, along with expected p-values")
    print("p-values \t\t\t0.2  \t0.05\t0.02\t0.01")
    print("fraction less than \t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}".format(
        p_count_20/sample_size, p_count_5/sample_size,
        p_count_2/sample_size, p_count_1/sample_size))


def generate_data(breakpoint_1=0, beta_1=0, breakpoint_2=0, beta_2=0, noise_scale=1):

    intercept = 5
    alpha = 1
    n_points = 100

    xx_bp = np.linspace(-9.5, 9.5, n_points)

    yy_bp = intercept + alpha*xx_bp + \
        beta_1 * np.maximum(xx_bp - breakpoint_1, 0) + \
        beta_2 * np.maximum(xx_bp - breakpoint_2, 0) + \
        np.random.normal(size=len(xx_bp), scale=noise_scale)

    return xx_bp, yy_bp


def test_under_null_hypothesis_no_breakpoints():

    print("For all of these checks, data is generated without breakpoints.")
    print("We repeat this many times, and record how many times gave a signifincat result.")
    print("We report the significance level and fraction passing the test.")
    print("The fraction testing should be below the significance level")
    # Assume no breakpoint
    check_p_values(alternative="two_sided")
    check_p_values(alternative="less")
    check_p_values(alternative="greater")

    # Low noise
    check_p_values(alternative="two_sided", noise_scale=0.1)
    check_p_values(alternative="less", noise_scale=0.1)
    check_p_values(alternative="greater", noise_scale=0.1)

    # High noise
    check_p_values(alternative="two_sided", noise_scale=10)
    check_p_values(alternative="less", noise_scale=10)
    check_p_values(alternative="greater", noise_scale=10)    


def test_with_breakpoints():
    print("\n\nFor all of these checks, data is generated with breakpoints.")
    print("We repeat this many times, and record how many times gave a signifincat result.")
    print("We report the significance level and fraction passing the test.")
    print("The fraction testing should be much higher than the significance level")
    print("This will depend on the change in gradient being a greater scale than the noise")


    breakpoints = [-1,0,5]
    beta_1s = [0.1,-0.2, -1, 3, 5]
    noises = [0.1, 1, 10]
    for breakpoint_1 in breakpoints:
        for beta_1 in beta_1s:
            for noise in noises:
                for alternative in ["two_sided"]:
                    print("\nbeta_1 is {}. noise is {}".format(beta_1, noise))
                    check_p_values(alternative=alternative, breakpoint_1=breakpoint_1, beta_1=beta_1, noise_scale=noise)


def test_with_two_breakpoints():
    print("\n\nFor all of these checks, data is generated with breakpoints.")
    print("We repeat this many times, and record how many times gave a signifincat result.")
    print("We report the significance level and fraction passing the test.")
    print("The fraction testing should be much higher than the significance level")
    print("This will depend on the change in gradient being a greater scale than the noise")


    breakpoints = [-1,0,5]
    beta_1s = [0.1,-0.2, -1, 3,5]
    noises = [0.1, 1, 10]
    for breakpoint_1 in breakpoints:
        for beta_1 in beta_1s:
            for noise in noises:
                for alternative in ["two_sided"]:
                    print("\nbeta_1 is {}. noise is {}".format(beta_1, noise))
                    breakpoint_2 = breakpoint_1 + 3
                    beta_2 = - beta_1
                    check_p_values(alternative=alternative, breakpoint_1=breakpoint_1, beta_1=beta_1, breakpoint_2=breakpoint_2, beta_2=beta_2, noise_scale=noise)


def test_one_sided_with_breakpoints():

    print("\n\nFor all of these checks, data is generated with breakpoints.")
    print("Here we look at one-side results.")
    print("We should see alternating high and low pass rates")
    breakpoints = [-3,-1,0,2,4,5]
    beta_1s = [3,5,-1,1,10]
    for breakpoint_1 in breakpoints:
        for beta_1 in beta_1s:
            for alternative in ["greater", "less"]:
                print("\nActual beta_1 is {}, looking at {} than 0".format(beta_1, alternative))
                check_p_values(alternative=alternative, breakpoint_1=breakpoint_1, beta_1=beta_1)


if __name__ == "__main__":
    
    #test_under_null_hypothesis_no_breakpoints()
    test_with_two_breakpoints()
    #test_one_sided_with_breakpoints()




import numpy as np


def validate_boolean(var, var_name):
    if isinstance(var, bool):
        return var
    else:
        raise ValueError(
            "{} must be a Boolean: True or False".format(var_name))


def validate_positive_integer(var, var_name):
    if isinstance(var, bool):
        raise ValueError("{} must be a positive Integer".format(var_name))
    if isinstance(var, int) and var > 0:
        return var
    else:
        raise ValueError("{} must be a positive Integer".format(var_name))


def validate_non_negative_integer(var, var_name):

    if isinstance(var, bool):
        raise ValueError("{} must be a non-negative Integer".format(var_name))

    if isinstance(var, int) and var >= 0:
        return var
    else:
        raise ValueError("{} must be a non-negative Integer".format(var_name))


def validate_positive_number(var, var_name):
    if isinstance(var, bool):
        raise ValueError("{} must be a Float".format(var_name))
    if (isinstance(var, float) or isinstance(var, int)) and var > 0:
        return var
    else:
        raise ValueError("{} must be a Float".format(var_name))


def validate_list_of_numbers(var, var_name, min_length):
    """
    Allowed types:
            List of integers of floats
            Numpy array of integers or floats
    """
    value_error_text = "{} must be a list of numbers with minimum length {}"
    value_error_text = value_error_text.format(var_name, min_length)

    # If its a list, convert it to a numpy array
    if isinstance(var, list):
        var = np.array(var)

    # If its not a numpy array at this point, raise a value error
    if not isinstance(var, np.ndarray):
        raise ValueError(value_error_text)

    # Check the array has numebrs in it
    if not np.issubdtype(var.dtype, np.number):
        raise ValueError(value_error_text)

    if len(var) < min_length:
        raise ValueError(value_error_text)

    return var

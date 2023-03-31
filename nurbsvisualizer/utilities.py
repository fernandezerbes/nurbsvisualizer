"""
Common algorithms used for the computation and visualization of NURBS.
"""

import numpy as np


def cox_de_boor(knot_span, polynomial_degree, parameter, knot_vector):
    """
    Evaluate the B-spline basis function.

    Arguments
    ---------
    knot_span : int
        The knot_span where the basis is evaluated.
    polynomial_degree : int
        The polynomial degree of the basis.
    parameter : float
        The parameter for which the basis is evaluated.
    knot_vector : array_like
            The array of knots to construct the basis, with shape (1, ).

    Returns
    -------
    float
        The basis function evaluation.

    Raises
    ------
    ValueError
        If parameter is not between the first and last knot values.
    """
    if parameter < knot_vector[0] or parameter > knot_vector[-1]:
        raise ValueError(
            "parameter value should be between the first and last knot values"
        )

    # Base case
    if polynomial_degree == 0:
        is_end_of_parameter_space = np.isclose(
            parameter, knot_vector[knot_span + 1]
        ) and np.isclose(parameter, knot_vector[-1])
        if (
            knot_vector[knot_span] <= parameter < knot_vector[knot_span + 1]
            or is_end_of_parameter_space
        ):
            return 1.0
        else:
            return 0.0
    # Recursion
    else:
        # Numerator on the left hand side
        n1 = parameter - knot_vector[knot_span]
        # Denominator on the left hand side
        d1 = knot_vector[knot_span + polynomial_degree] - knot_vector[knot_span]
        # Quotient on the left hand side
        q1 = cox_de_boor_division(n1, d1)

        # Numerator on the right hand side
        n2 = knot_vector[knot_span + polynomial_degree + 1] - parameter
        # Denominator on the right hand side
        d2 = knot_vector[knot_span + polynomial_degree + 1] - knot_vector[knot_span + 1]
        # Quotient on the right hand side
        q2 = cox_de_boor_division(n2, d2)

        # Cox-De Boor's recursion formula
        return q1 * cox_de_boor(
            knot_span, polynomial_degree - 1, parameter, knot_vector
        ) + q2 * cox_de_boor(
            knot_span + 1, polynomial_degree - 1, parameter, knot_vector
        )


def cox_de_boor_division(numerator, denominator):
    """
    Handle division by zero defined for the Cox-De Boor's algorithm.

    Arguments
    ---------
    numerator : float
        The numerator of the division.
    denominator : float
        The denominator of the division.

    Returns
    -------
    float
        The division of numerator by denominator.
    """
    if np.isclose(denominator, 0.0):
        return 0.0

    else:
        return numerator / denominator


def evaluate_bspline_basis(
    polynomial_degree, knot_vector, control_points_count, samples
):
    """
    Evaluate the B-spline basis function for the given samples.

    Arguments
    ---------
    polynomial_degree : int
        The polynomial degree of the basis.
    knot_vector : array_like
            The array of knots to construct the basis, with shape (1, ).
    control_points_count: int
        The number of control points.
    samples : array_like
        The array of samples coordinates in the parameter space, with shape (1, samples_count).

    Returns
    -------
    ndarray
        The array of basis functions evaluations, with shape (control_points_count, len(samples)).
    """
    # Initialize the array of basis evaluations
    evaluations = np.zeros((control_points_count, len(samples)))

    # Evaluate the shape functions for the given samples using the Cox-De Boor's algorithm
    for i in range(control_points_count):
        first_nonzero_sample = knot_vector[i]
        last_nonzero_sample = knot_vector[i + polynomial_degree + 1]

        # Perform the computation only where the shape function i has support
        for j, sample in enumerate(samples):
            if first_nonzero_sample <= sample <= last_nonzero_sample:
                evaluations[i, j] = cox_de_boor(
                    i, polynomial_degree, sample, knot_vector
                )

    return evaluations


def get_index_to_nearest(array, value):
    """
    Get index of the nearest array entry to a given value.

    Arguments
    ---------
    array : array_like
        The array to search for the nearest value.
    value : float
        The value to look for.

    Returns
    -------
    ndarray
        The index of the nearest array entry.
    """
    return (np.abs(array - value)).argmin()


def replace_nearest(array, value):
    """
    Replace the nearest array entry to a given value for the given value.

    Arguments
    ---------
    array : array_like
        The array to replace the nearest value.
    value : float
        The value to replace.
    """
    # Get index of the nearest array entry to the given value
    index = get_index_to_nearest(array, value)

    # Replace the value in index
    array[index] = value


def is_open(knot_vector, polynomial_degree):
    """
    Check whether a knot vector is open.

    Arguments
    ---------
    knot_vector : array_like
        The knot vector.
    polynomial_degree : int
        The polynomial degree of the curve generated by knot_vector.

    Returns
    -------
    bool
        True if the knot vector is open or False otherwise.
    """
    # Check whether the first and last (polynomial_degree + 1) knots are the same
    for i in range(1, polynomial_degree + 1):
        if knot_vector[0] != knot_vector[i] or knot_vector[-1] != knot_vector[-1 - i]:
            return False
    if (
        knot_vector[polynomial_degree + 1] == knot_vector[0]
        or knot_vector[-1 - (polynomial_degree + 1)] == knot_vector[-1]
    ):
        return False
    return True


def is_non_descending(knot_vector):
    """
    Check whether a knot vector is in non-descending order.

    Arguments
    ---------
    knot_vector : array_like
        The knot vector.

    Returns
    -------
    bool
        True if the knot vector is in non-descending order or False otherwise.
    """
    for i in range(len(knot_vector) - 1):
        if knot_vector[i + 1] < knot_vector[i]:
            return False
    return True


def generate_open_knot(control_points_count, polynomial_degree):
    """
    Generate a uniform open knot vector in the normalized range.

    Arguments
    ---------
    control_points_count : int
        The number of control points.
    polynomial_degree : int
        The polynomial degree.

    Returns
    -------
    np.ndarray
        True uniform open knot vector in the normalized range.

    Raises
    ------
    ValueError
        If control_points_count is smaller than polynomial_degree + 1.
    """
    if control_points_count < polynomial_degree + 1:
        raise ValueError(
            "control_points_count should be at least polynomial_degree + 1"
        )

    knot_vector = np.zeros(control_points_count + polynomial_degree + 1)

    # Set the last (polynomial_degree + 1) to 1
    knot_vector[control_points_count:] += 1

    # Calculate the delta between inner knots
    delta = 1 / (control_points_count - polynomial_degree)

    # Compute inner knots
    for i in range(polynomial_degree + 1, control_points_count):
        knot_vector[i] = knot_vector[i - 1] + delta

    return knot_vector

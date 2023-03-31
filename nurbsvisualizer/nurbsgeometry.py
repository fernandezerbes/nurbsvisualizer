import numpy as np

from nurbsvisualizer.bsplinegeometry import BsplineCurve, BsplineSurface
from nurbsvisualizer.utilities import evaluate_bspline_basis, is_open


class NurbsCurve(BsplineCurve):
    """
    A class used for the definition of NURBS curves.

    Attributes
    ----------
    control_points : array_like
        The array of spatial coordinates of the control points, with shape (dimensions, control_points_count).
    weights: array_like
        The array of weights, with shape (1, control_points_count).
    knot_vector : array_like
        The array of knots to construct the basis, with shape (1, ).
    samples_count : int, optional, default=101
        The number of samples used to evaluate the curve.
    """

    def __init__(self, control_points, weights, knot_vector, samples_count=101):
        """
        Construct NURBS curve objects.

        Arguments
        ---------
        control_points : array_like
            The array of spatial coordinates of the control points, with shape (dimensions, control_points_count).
        weights: array_like
            The array of weights, with shape (1, control_points_count).
        knot_vector : array_like
            The array of knots to construct the basis, with shape (1, ).
        samples_count : int, optional, default=101
            The number of samples used to evaluate the curve.
        """
        BsplineCurve.__init__(self, control_points, knot_vector, samples_count)
        self.weights = weights

    @property
    def weights(self):
        """
        Get the weights.

        Returns
        -------
        ndarray
            The array of weights.
        """
        return self._weights

    @weights.setter
    def weights(self, value):
        """
        Set the weights.

        Raises
        ------
        ValueError
            If value doesn't have the right dimensions.
        TypeError
            If value isn't a np.ndarray or a tuple or a list.
        """
        if np.shape(value)[0] != self.control_points_count:
            raise ValueError(
                "the number of weights should be the same as the number of control points"
            )
        elif isinstance(value, np.ndarray):
            self._weights = value
        elif isinstance(value, (tuple, list)):
            self._weights = np.array(value)
        else:
            raise TypeError("weights should be a tuple, list, or np.ndarray")

    @BsplineCurve.knot_vector.setter
    def knot_vector(self, value):
        """
        Set the knot vector.

        This setter overrides @BsplineCurve.knot_vector.setter.

        Raises
        ------
        ValueError
            If value doesn't have the right dimensions.
            If value is not in non-descending order.
        TypeError
            If value isn't a np.ndarray or a tuple or a list.
        NotImplementedError
            If value is a non-open knot vector.
        """
        if not is_open(value, np.shape(value)[0] - self.control_points_count - 1):
            raise NotImplementedError("non-open knot vectors aren't supported")

        # Call the B-spline knot_vector setter
        BsplineCurve.knot_vector.fset(self, value)

    def evaluate_basis(self):
        """
        Evaluate the NURBS basis functions.

        This method overrides BsplineCurve.evaluate_basis.
        """
        # Evaluate the BSpline basis with the Cox-De Boor's algorithm
        evaluations = evaluate_bspline_basis(
            self.polynomial_degree,
            self.knot_vector,
            self.control_points_count,
            self.samples,
        )

        # Multiply basis evaluations with the corresponding weights (numerator of rational basis formula)
        weighted_evaluations = np.reshape(self.weights, (-1, 1)) * evaluations

        # Sum of basis evaluations at the same sample (denominator of rational basis formula).
        weighted_evaluations_sum = weighted_evaluations.sum(axis=0)

        # Compute the NURBS basis evaluations and return
        self._basis_evals = weighted_evaluations / weighted_evaluations_sum


class NurbsSurface(BsplineSurface):
    """
    A class used for the definition of NURBS surfaces.

    Attributes
    ----------
    control_points : array_like
        The array of spatial coordinates of the control points, with shape (3, control_points_count).
    weights: array_like
        The array of weights, with shape (control_points_count).
    knot_vector : array_like
        The array of knots to construct the basis, with shape (2, ).
    samples_count : array_like, optional, default=(51, 51)
        The number of samples used to evaluate the surface in each parameter dimension.
    """

    def __init__(self, control_points, weights, knot_vector, samples_count=(51, 51)):
        """
        Construct NURBS surface objects.

        Arguments
        ---------
        control_points : array_like
            The array of spatial coordinates of the control points, with shape (3, control_points_count).
        weights: array_like
            The array of weights, with shape (control_points_count).
        knot_vector : array_like
            The array of knots to construct the basis, with shape (2, ).
        samples_count : array_like, optional, default=(51, 51)
            The number of samples used to evaluate the surface in each parameter dimension.
        """
        BsplineSurface.__init__(self, control_points, knot_vector, samples_count)
        self.weights = weights

    @property
    def weights(self):
        """
        Get the weights.

        Returns
        -------
        ndarray
            The array of weights.
        """
        return self._weights

    @weights.setter
    def weights(self, value):
        """
        Set the weights.

        Raises
        ------
        ValueError
            If value doesn't have the right dimensions.
        TypeError
            If value isn't a np.ndarray or a tuple or a list.
        """
        if np.shape(value) != self.control_points[0].shape:
            raise ValueError(
                "the shape of the weights grid doesn't match the shape of the control points grid"
            )
        elif isinstance(value, np.ndarray):
            self._weights = value
        elif isinstance(value, (tuple, list)):
            self._weights = np.array(value)
        else:
            raise TypeError("weights should be a tuple, list, or np.ndarray")

    @BsplineSurface.knot_vector.setter
    def knot_vector(self, value):
        """
        Set the knot vector.

        This setter overrides @BsplineSurface.knot_vector.setter.

        Raises
        ------
        ValueError
            If value doesn't have the right dimensions.
            If value is not in non-descending order.
        TypeError
            If value isn't a np.ndarray or a tuple or a list.
        NotImplementedError
            If value is a non-open knot vector.
        """
        if not is_open(
            value[0], np.shape(value[0])[0] - self.control_points_count[0] - 1
        ) or not is_open(
            value[1], np.shape(value[1])[0] - self.control_points_count[1] - 1
        ):
            raise NotImplementedError("non-open knot vectors aren't supported")

        # Call the B-spline knot_vector setter
        BsplineSurface.knot_vector.fset(self, value)

    def evaluate_basis(self):
        """
        Evaluate the NURBS basis functions.

        This method overrides BsplineSurface.evaluate_basis.
        """
        # Initialize the array of basis evaluations
        self._basis_evals = np.zeros(
            (
                self.control_points_count[0],
                self.control_points_count[1],
                self.samples_count[0],
                self.samples_count[1],
            )
        )

        # Evaluate the BSpline basis for each parameter space dimension with the Cox-De Boor's algorithm
        evaluations_xi = evaluate_bspline_basis(
            self.polynomial_degree[0],
            self.knot_vector[0],
            self.control_points_count[0],
            self.samples[0],
        )
        evaluations_eta = evaluate_bspline_basis(
            self.polynomial_degree[1],
            self.knot_vector[1],
            self.control_points_count[1],
            self.samples[1],
        )

        # Compute the tensor product of basis evaluations in each parameter dimension and multiply each entry of the
        # result with the corresponding weight (numerator of rational basis formula).
        for n, Ni_evals in enumerate(evaluations_xi):
            for m, Nj_evals in enumerate(evaluations_eta):
                a = np.reshape(Nj_evals, (-1, 1)) * Ni_evals * self.weights[n][m]
                self._basis_evals[n][m] = a

        # Compute the sum of basis functions at the same sample in each dimension (denominator of rational basis
        # formula). The sum is done over the first and second axis of the evaluations array.
        basis_sums = np.sum(self._basis_evals, (0, 1))

        # Compute the NURBS basis evaluations
        for n in range(self.control_points_count[0]):
            for m in range(self.control_points_count[1]):
                self._basis_evals[n][m] = np.divide(self._basis_evals[n][m], basis_sums)

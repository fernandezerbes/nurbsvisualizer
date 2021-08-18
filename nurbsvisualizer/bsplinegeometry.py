import numpy as np
from nurbsvisualizer.utilities import replace_nearest, evaluate_bspline_basis, is_non_descending


class BsplineCurve:
    """
    A class used for the definition of B-spline curves.

    Attributes
    ----------
    control_points : array_like
        The array of spatial coordinates of the control points, with shape (dimensions, control_points_count).
    knot_vector : array_like
        The array of knots to construct the basis, with shape (1, ).
    samples_count : int, optional, default=101
        The number of samples used to evaluate the curve.
    """

    def __init__(self, control_points, knot_vector, samples_count=101):
        """
        Construct B-spline curve objects.

        Arguments
        ---------
        control_points : array_like
            The array of spatial coordinates of the control points, with shape (dimensions, control_points_count).
        knot_vector : array_like
            The array of knots to construct the basis, with shape (1, ).
        samples_count : int, optional, default=101
            The number of samples used to evaluate the curve.
        """

        self.control_points = control_points
        self.knot_vector = knot_vector
        self.samples_count = samples_count
        self._samples = None
        self._basis_evals = None
        self._curve_evals = None

    @property
    def control_points(self):
        """
        Get the control points.

        Returns
        -------
        ndarray
            The array of control points.
        """
        return self._control_points

    @control_points.setter
    def control_points(self, value):
        """
        Set the control points.

        Raises
        ------
        ValueError
            If value doesn't have the right dimensions.
        TypeError
            If value isn't a np.ndarray or a tuple or a list.
        """
        if np.shape(value)[0] != 2 and np.shape(value)[0] != 3:
            raise ValueError("control_points should have 2 or 3 dimensions")
        elif isinstance(value, np.ndarray):
            self._control_points = value
        elif isinstance(value, (tuple, list)):
            self._control_points = np.array(value)
        else:
            raise TypeError("control_points should be a tuple, list, or np.ndarray")

    @property
    def knot_vector(self):
        """
        Get the knot vector.

        Returns
        -------
        ndarray
            The array of knot vector.
        """
        return self._knot_vector

    @knot_vector.setter
    def knot_vector(self, value):
        """
        Set the knot vector.

        Raises
        ------
        ValueError
            If value doesn't have the right dimensions.
            If value is not in non-descending order.
        TypeError
            If value isn't a np.ndarray or a tuple or a list.
        """
        if np.ndim(value) != 1:
            raise ValueError("knot_vector should have dimension 1")
        elif np.shape(value)[0] - self.control_points_count < 2:
            raise ValueError("knot_vector should have at least {} knots".format(self.control_points_count + 2))
        elif not is_non_descending(value):
            raise ValueError("knot_vector should be in non-descending order")
        elif isinstance(value, np.ndarray):
            self._knot_vector = value
        elif isinstance(value, (tuple, list)):
            self._knot_vector = np.array(value)
        else:
            raise TypeError("knot_vector should be a tuple, list, or np.ndarray")

    @property
    def polynomial_degree(self):
        """
        Get the polynomial degree.

        Returns
        -------
        int
            The polynomial degree.
        """
        return self.knot_vector.shape[0] - self.control_points_count - 1

    @property
    def dimensions(self):
        """
        Get the geometric dimension of the curve.

        Returns
        -------
        int
            The geometric dimension of the curve.
        """
        return self.control_points.shape[0]

    @property
    def control_points_count(self):
        """
        Get control points count.

        Returns
        -------
        int
            The control points count.
        """
        return self.control_points.shape[1]

    @property
    def samples(self):
        """
        Get samples for basis and curve evaluations.

        Returns
        -------
        ndarray
            The array of samples coordinates in the parameter space, with shape (1, samples_count).
        """
        if self._samples is None:
            self.generate_samples()
        return self._samples

    def generate_samples(self):
        """
        Generate the samples coordinates in the parameter space.

        The generated samples are evenly spaced in the parameter space, except at the inner knot positions. In these
        positions, the samples closest to the inner knots are replaced by the latest. This is done because samples at
        the knot positions are needed in the case of discontinuous basis functions.
        """
        # Generate an evenly spaced array of samples
        self._samples = np.linspace(self.knot_vector[0], self.knot_vector[-1], self.samples_count)

        # Get the array of distinct knots
        knots_to_replace = self.get_unique_knots()

        # Replace samples closest to inner knots
        for knot in knots_to_replace[1:-1]:
            replace_nearest(self._samples, knot)

    def get_unique_knots(self):
        """
        Get the array of distinct knots.

        For example, if the knots are [0, 0, 0, 0.5, 1, 1, 1], then the distinct knots are [0, 0.5, 1].

        Returns
        -------
        ndarray
            The distinct knots.
        """
        return np.unique(self.knot_vector)

    @property
    def basis_evals(self):
        """
        Get the basis evaluations.

        Returns
        -------
        ndarray
            The array of basis functions evaluations, with shape (control_points_count, samples_count).
        """
        if self._basis_evals is None:
            self.evaluate_basis()
        return self._basis_evals

    def evaluate_basis(self):
        """
        Evaluate the B-spline basis functions.

        Returns
        -------
        ndarray
            The array of basis function evaluations.
        """
        # Evaluate the BSpline basis with the Cox-De Boor's algorithm
        self._basis_evals = evaluate_bspline_basis(self.polynomial_degree, self.knot_vector,
                                                             self.control_points_count, self.samples)

    @property
    def curve_evals(self):
        """
        Get the curve evaluations.

        Returns
        -------
        ndarray
            The array of curve evaluations, with shape (dimensions, samples_count).
        """
        if self._curve_evals is None:
            self.evaluate_curve()
        return self._curve_evals

    def evaluate_curve(self):
        """Evaluate the curve."""
        # Multiply each control point by the corresponding basis function
        self._curve_evals = np.dot(self.control_points, self.basis_evals)


class BsplineSurface:
    """
    A class used for the definition of B-spline surfaces.

    Attributes
    ----------
    control_points : array_like
        The array of spatial coordinates of the control points, with shape (3, control_points_count).
    knot_vector : array_like
        The array of knots to construct the basis, with shape (2, ).
    samples_count : array_like, optional, default=(51, 51)
        The number of samples used to evaluate the surface in each parameter dimension.
    """

    def __init__(self, control_points, knot_vector, samples_count=(51, 51)):
        """
        Construct B-spline surface objects.

        Arguments
        ---------
        control_points : array_like
            The array of spatial coordinates of the control points, with shape (3, control_points_count).
        knot_vector : array_like
            The array of knots to construct the basis, with shape (2, ).
        samples_count : array_like, optional, default=(51, 51)
            The number of samples used to evaluate the surface in each parameter dimension.
        """

        self.control_points = control_points
        self.knot_vector = knot_vector
        self.samples_count = samples_count
        self._samples = None
        self._basis_evals = None
        self._surface_evals = None

    @property
    def control_points(self):
        """
        Get the control points.

        Returns
        -------
        ndarray
            The array of control points.
        """
        return self._control_points

    @control_points.setter
    def control_points(self, value):
        """
        Set the control points.

        Raises
        ------
        ValueError
            If value doesn't have the right dimensions.
        TypeError
            If value isn't a np.ndarray or a tuple or a list.
        """
        if np.shape(value)[0] != 3:
            raise ValueError("control_points should have 3 dimensions")
        elif isinstance(value, np.ndarray):
            self._control_points = value
        elif isinstance(value, (tuple, list)):
            self._control_points = np.array(value)
        else:
            raise TypeError("control_points should be a tuple, list, or np.ndarray")

    @property
    def knot_vector(self):
        """
        Get the knot vector.

        Returns
        -------
        ndarray
            The array of knot vector.
        """
        return self._knot_vector

    @knot_vector.setter
    def knot_vector(self, value):
        """
        Set the knot vector.

        Raises
        ------
        ValueError
            If value doesn't have the right dimensions.
            If value is not in non-descending order.
        TypeError
            If value isn't a np.ndarray or a tuple or a list.
        """
        if np.shape(value)[0] != 2:
            raise ValueError("knot_vector should have one knot vector for each parametric direction")
        elif np.ndim(value[0]) != 1 or np.ndim(value[1]) != 1:
            raise ValueError("each knot vector in knot_vector should have dimension 1")
        elif np.shape(value[0])[0] - self.control_points_count[0] < 2:
            raise ValueError("knot_vector should have at least {} knots in the first direction".format(
                             self.control_points_count[0] + 2))
        elif np.shape(value[1])[0] - self.control_points_count[1] < 2:
            raise ValueError("knot_vector should have at least {} knots in the second direction".format(
                             self.control_points_count[1] + 2))
        elif not is_non_descending(value[0]) or not is_non_descending(value[1]):
            raise ValueError("knot_vector should be in non-descending order")
        elif isinstance(value, np.ndarray):
            self._knot_vector = value
        elif isinstance(value, (tuple, list)):
            self._knot_vector = np.array(value)
        else:
            raise TypeError("knot_vector should be a tuple, list, or np.ndarray")

    @property
    def polynomial_degree(self):
        """
        Get the polynomial degree in each parametric direction.

        Returns
        -------
        (int, int)
            The polynomial degree in each parametric direction.
        """
        return np.shape(self.knot_vector[0])[0] - self.control_points_count[0] - 1,\
               np.shape(self.knot_vector[1])[0] - self.control_points_count[1] - 1

    @property
    def dimensions(self):
        """
        Get the geometric dimension of the surface.

        Returns
        -------
        int
            The geometric dimension of the surface.
        """
        return self.control_points.shape[0]

    @property
    def control_points_count(self):
        """
        Get control points count in each parametric direction.

        Returns
        -------
        (int, int)
            The control points count in each parametric direction.
        """
        return self.control_points[0].shape

    @property
    def samples(self):
        """
        Get samples for basis and surface evaluations.

        Returns
        -------
        samples : tuple of ndarray
            The arrays of samples coordinates in each parameter space dimension.
        """
        if self._samples is None:
            self.generate_samples()
        return self._samples

    def generate_samples(self):
        """
        Generate the samples coordinates in the parameter space.

        The generated samples are evenly spaced in the parameter space, except at the inner knot positions. In these
        positions, the samples closest to the inner knots are replaced by the latest. This is done because samples at
        the knot positions are needed in the case of discontinuous basis functions.
        """
        # Generate an evenly spaced array of samples for parameter space dimension
        samples_xi = np.linspace(self.knot_vector[0][0], self.knot_vector[0][-1], self.samples_count[0])
        samples_eta = np.linspace(self.knot_vector[1][0], self.knot_vector[1][-1], self.samples_count[1])

        # Get the array of distinct knots for each parameter space dimension
        knots_to_replace_xi, knots_to_replace_eta = self.get_unique_knots()

        # Replace samples closest to inner knots in each parameter space dimension
        for knot_xi in knots_to_replace_xi[1:-1]:
            replace_nearest(samples_xi, knot_xi)
        for knot_eta in knots_to_replace_eta[1:-1]:
            replace_nearest(samples_eta, knot_eta)

        self._samples = (samples_xi, samples_eta)

    def get_unique_knots(self):
        """
        Get an array of distinct knots for each parameter space dimension.

        For example, if the knots are ([0, 0, 0, 0.5, 1, 1, 1], [1, 1, 1, 1.5, 2, 2, 2]), then the distinct knots are
        ([0, 0.5, 1], [1, 1.5, 2]).

        Returns
        -------
        tuple of ndarray
            The distinct knots for each parameter space dimension.
        """
        return np.unique(self.knot_vector[0]), np.unique(self.knot_vector[1])

    @property
    def basis_evals(self):
        """
        Get the basis evaluations.

        Returns
        -------
        basis_evals : ndarray
            The array of basis functions evaluations, with shape (control_points_count, samples_count).
        """
        if self._basis_evals is None:
            self.evaluate_basis()
        return self._basis_evals

    def evaluate_basis(self):
        """Evaluate the B-spline basis functions."""
        # Initialize the array of basis evaluations
        self._basis_evals = np.zeros((self.control_points_count[0], self.control_points_count[1],
                                      self.samples_count[0], self.samples_count[1]))

        # Evaluate the BSpline basis for each parameter space dimension with the Cox-De Boor's algorithm
        evaluations_xi = evaluate_bspline_basis(self.polynomial_degree[0], self.knot_vector[0],
                                                          self.control_points_count[0], self.samples[0])
        evaluations_eta = evaluate_bspline_basis(self.polynomial_degree[1], self.knot_vector[1],
                                                           self.control_points_count[1], self.samples[1])

        # Compute the tensor product of basis evaluations in each parameter dimension.
        for n, Ni_evals in enumerate(evaluations_xi):
            for m, Nj_evals in enumerate(evaluations_eta):
                self._basis_evals[n][m] = np.reshape(Nj_evals, (-1, 1)) * Ni_evals

    @property
    def surface_evals(self):
        """
        Get the surface evaluations.

        Returns
        -------
        ndarray
            The array of surface evaluations, with shape (3, samples_count).
        """
        if self._surface_evals is None:
            self.evaluate_surface()
        return self._surface_evals

    def evaluate_surface(self):
        """Evaluate the surface."""
        # Multiply each control point by the corresponding basis function

        # Initialize the array of surface evaluations
        self._surface_evals = np.zeros((self.dimensions, self.samples_count[0], self.samples_count[1]))

        # Loop over each dimension (x, y, z) and multiply each control point coordinate with the corresponding basis
        for dim, coordinate_grid in enumerate(self.control_points):
            for i in range(self.control_points_count[0]):
                for j in range(self.control_points_count[1]):
                    self._surface_evals[dim] += self.control_points[dim][i][j] * self.basis_evals[i][j]

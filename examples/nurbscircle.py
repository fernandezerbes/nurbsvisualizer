import numpy as np

from nurbsvisualizer.nurbsgeometry import NurbsCurve
from nurbsvisualizer.visualizer import CurveVisualizer

#################
# 2D NURBS CIRCLE
#################

# Define knot vector
knot_vector = [
    0,
    0,
    0,
    np.pi / 2,
    np.pi / 2,
    np.pi,
    np.pi,
    3 * np.pi / 2,
    3 * np.pi / 2,
    2 * np.pi,
    2 * np.pi,
    2 * np.pi,
]

# Define control points
control_points = [[-1, -1, 0, 1, 1, 1, 0, -1, -1], [0, 1, 1, 1, 0, -1, -1, -1, 0]]

# Define weights
weights = [
    1,
    1 / np.sqrt(2),
    1,
    1 / np.sqrt(2),
    1,
    1 / np.sqrt(2),
    1,
    1 / np.sqrt(2),
    1,
]

# Create a nurbs curve object
nurbs_circle = NurbsCurve(control_points, weights, knot_vector)

# Visualize the curve and basis
CurveVisualizer(nurbs_circle)

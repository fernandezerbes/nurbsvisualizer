from nurbsvisualizer.bsplinegeometry import BsplineCurve
from nurbsvisualizer.visualizer import CurveVisualizer


##################
# 3D BSPLINE CURVE
##################

# Define 3D control points
control_points = [[0, 1, 4, 7, 10, 15, 17, 20],   # x values
                  [0, 5, 0, 8, -5, 5, 1, -3],     # y values
                  [2, 3, -5, 8, -1, 4, 9, 0]]     # z values

# Define knot vector
knot_vector = [0, 0, 0, 0, 0.3, 0.5, 0.65, 0.8, 1, 1, 1, 1]

# Create a B-spline curve object
bspline_curve_3d = BsplineCurve(control_points, knot_vector)

# Visualize curve and basis
CurveVisualizer(bspline_curve_3d)

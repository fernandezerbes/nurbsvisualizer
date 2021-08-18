from nurbsvisualizer.nurbsgeometry import NurbsSurface
from nurbsvisualizer.visualizer import SurfaceVisualizer


###############
# NURBS SURFACE
###############

# Define control point coordinates grids
x = [[0, 1, 2, 3],
     [0, 1, 2, 3],
     [0, 1, 2, 3],
     [0, 1, 2, 3]]

y = [[0, 0, 0, 0],
     [1, 1, 1, 1],
     [2, 2, 2, 2],
     [3, 3, 3, 3]]

z = [[0, 1, 1, 0],
     [1, 2, 2, 1],
     [1, 2, 2, 1],
     [0, 1, 1, 0]]

# Collect coordinate grids in an array
control_points = [x, y, z]

# Define weights
weights = [[0.5, 1, 1, 0.8],
           [0.5, 1, 2, 0.8],
           [0.5, 1, 2, 0.8],
           [0.5, 1, 2, 0.8]]

# Define knot vector
knot_vector = [[0, 0, 0, 0, 1, 1, 1, 1],     # xi values
               [0, 0, 0, 1, 2, 2, 2]]        # eta values

# Create a NURBS surface object
nurbs_surface = NurbsSurface(control_points, weights, knot_vector)

# Visualize the surface and basis
SurfaceVisualizer(nurbs_surface)

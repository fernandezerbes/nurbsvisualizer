from nurbsvisualizer.bsplinegeometry import BsplineSurface
from nurbsvisualizer.visualizer import SurfaceVisualizer


#################
# BSPLINE SURFACE
#################

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

# Define knot vector
knot_vector = [[0, 0, 0, 0, 1, 1, 1, 1],     # xi values
               [0, 0, 0, 1, 2, 2, 2]]        # eta values

# Create a B-Spline surface object
bspline_surface = BsplineSurface(control_points, knot_vector)

# Visualize the surface and basis
SurfaceVisualizer(bspline_surface)

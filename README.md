# NURBSvisualizer

NURBSvisualizer is a tool for the generation and interactive visualization of
B-spline and NURBS basis functions, curves, and surfaces.

## Usage

### Install the nurbsvisualizer package

1. Clone the repository

```shell
git clone https://github.com/FernandezErbes/nurbsvisualizer.git
```

2. Go to the repository folder and install the package with [pip](https://pypi.org/project/pip/)

```shell
pip install .
```

3. Create your own drivers by following the [examples](https://github.com/FernandezErbes/nurbsvisualizer/tree/main/examples)


### Curve visualization

#### Object creation and visualization

1.  Import ```NurbsCurve``` from ```nurbsgeometry``` (alternatively ```BsplineCurve``` from ```bsplinegeometry```):

```python
from nurbsvisualizer.nurbsgeometry import NurbsCurve
```


2. Import ```CurveVisualizer``` from ```visualizer```:

```python
from nurbsvisualizer.visualizer import CurveVisualizer
```


3.  Set up the control points:

```python
control_points = [[0, 1, 4, 7, 10, 15, 17, 20],   # x values
                  [0, 5, 0, 8, -5, 5, 1, -3],     # y values
                  [2, 3, -5, 8, -1, 4, 9, 0]]     # z values
```


4. Define the weights (for NURBS only): 

```python
weights = [0.5, 1, 1, 2, 0.8, 1, 2.5, 0.2]
```


5. Define the knot vector: 

```python
knot_vector = [0, 0, 0, 0, 0.3, 0.5, 0.65, 0.8, 1, 1, 1, 1]
```


6. Create the curve object: 

```python
nurbs_curve_3D = NurbsCurve(control_points, weights, knot_vector)
```


7. Visualize the curve and basis: 

```python
CurveVisualizer(nurbs_curve_2d)
```

![curve_plot](https://user-images.githubusercontent.com/62465061/130354532-21a82925-5ba3-41d8-abce-280fceb64000.gif)

*Note - Alternative automatic generation of uniform open knot vectors (here for
8 control points and polynomial degree 2):*

```python
import utilities
knot_vector = utilities.generate_open_knot(8, 2)
```

#### Interaction


*  Move the mouse in the basis plot to see the corresponding point in the curve.
*  Select/deselect a basis curve by clicking to highlight it and label/unlabel the corresponding control point.
*  Pick and drag the curve plot to rotate it (only for 3D curves).
*  Press the *C* key to toggle on and off the curve.
*  Press the *P* key to toggle on and off the control polygon.
*  Press the *K* key to toggle on and off the points in the curve corresponding to knots.


---
### Surface visualization

#### Object creation and visualization

1.  Import ```NurbsSurface``` from ```nurbsgeometry``` (alternatively ```BsplineSurface``` from ```bsplinegeometry```):

```python
from nurbsvisualizer.nurbsgeometry import NurbsSurface
```


2. Import ```SurfaceVisualizer``` from ```visualizer```:

```python
from nurbsvisualizer.visualizer import SurfaceVisualizer
```


3.  Set up the control points:

```python
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
```


4. Define the weights (for NURBS only): 

```python
weights = [[0.5, 1, 1, 0.8],
           [0.5, 1, 2, 0.8],
           [0.5, 1, 2, 0.8],
           [0.5, 1, 2, 0.8]]
```


5. Define the knot vectors: 

```python
knot_vector = [[0, 0, 0, 0, 1, 1, 1, 1],     # xi values
               [0, 0, 0, 1, 2, 2, 2]]        # eta values
```


6. Create the surface object: 

```python
nurbs_surface = nurbsgeometry.NurbsSurface(control_points, weights, knot_vector)
```


7. Visualize the surface and basis: 

```python
SurfaceVisualizer(nurbs_surface)
```

![surface_plot](https://user-images.githubusercontent.com/62465061/130354540-3fde5a8e-d71b-49cc-888a-6226b1f16854.gif)

#### Interaction

*  Scroll up/down in the basis plot to change the shown basis. Alternatively press *Up/Down* keys.
*  Pick and drag the plots to rotate them.
*  Press the *S* key to toggle on and off the surface.
*  Press the *P* key to toggle on and off the control polygon.

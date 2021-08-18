import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import proj3d
from matplotlib.text import OffsetFrom
import numpy as np
from matplotlib import cm
from matplotlib.ticker import AutoLocator
from matplotlib.widgets import Cursor
from nurbsvisualizer.utilities import get_index_to_nearest
from nurbsvisualizer.bsplinegeometry import BsplineCurve, BsplineSurface
from nurbsvisualizer.nurbsgeometry import NurbsCurve, NurbsSurface

# Set up matplotlib global parameters
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['font.size'] = 12
matplotlib.rcParams['axes.labelsize'] = 14
matplotlib.rcParams["figure.figsize"] = [15, 8]
matplotlib.rcParams["axes.titlepad"] = 45

# Disable default 's' key shortcut for saving the plot. Use 'Ctrl + s' instead.
plt.rcParams['keymap.save'].remove('s')
# Disable default 'k' key shortcut for changing x scale of the plot.
plt.rcParams['keymap.xscale'].remove('k')
# Disable default 'p' key shortcut for changing panning.
plt.rcParams['keymap.pan'].remove('p')


class CurveVisualizer:
    """
    A class used for visualization of B-spline and NURBS curves.

    Attributes
    ----------
    curve : BsplineCurve or NurbsCurve
        The curve object to visualize.
    """

    def __init__(self, curve):
        """
        Construct B-spline or NURBS curves.

        Arguments
        ---------
        curve : BsplineCurve or NurbsCurve
            The curve object to visualize.
        """
        self.curve = curve
        self.is_3D = self.curve.dimensions == 3

        # Create figure for the plots
        self.fig = plt.figure()
        self.set_title()

        # Initialize attributes to store the matplotlib artists and states for future manipulation
        self.control_point_label = None
        self.curve_art = None
        self.control_polygon_art = None
        self.control_polygon_marker_art = None
        self.basis_function_art = []
        self.picked_basis_index = None
        self.knot_points_art = None
        self.knot_lines_art = []

        # Create basis plot
        self.ax_basis = self.fig.add_subplot(122)
        self.basis_plot_setup()
        self.plot_basis()

        # Create vertical cursor line for basis plot
        self.cursor = Cursor(self.ax_basis, horizOn=False, vertOn=True, linewidth=0.5, color="grey")

        # Create curve plot
        if self.is_3D:
            self.ax_curve = self.fig.add_subplot(121, projection='3d')
            self.ax_curve.view_init(elev=30, azim=45)
        else:
            self.ax_curve = self.fig.add_subplot(121)
        self.curve_plot_setup()
        self.plot_curve()
        self.plot_control_polygon()
        self.plot_knot_points()
        self.plot_knot_lines()

        # Define events for interactive plots
        self.key_press_event = self.fig.canvas.mpl_connect('key_press_event', self.on_key_press)
        self.motion_notify_event_basis = self.fig.canvas.mpl_connect('motion_notify_event', self.on_motion)
        self.pick_event = self.fig.canvas.mpl_connect('pick_event', self.on_pick)
        if self.is_3D:
            self.motion_notify_event_curve = self.fig.canvas.mpl_connect('motion_notify_event',
                                                                         self.update_control_point_label_position)

        # Show the plot
        plt.tight_layout(pad=2, w_pad=2.5)
        plt.show()

    def set_title(self):
        """Set the window title depending on the type of curve and dimensions."""
        if type(self.curve) is BsplineCurve:
            if self.is_3D:
                title = '3D B-spline Curve Visualizer'
            else:
                title = '2D B-spline Curve Visualizer'
        elif type(self.curve) is NurbsCurve:
            if self.is_3D:
                title = '3D NURBS Curve Visualizer'
            else:
                title = '2D NURBS Curve Visualizer'
        self.fig.canvas.set_window_title(title)

    def basis_plot_setup(self):
        """Set up the basis plot format."""
        # Set up plot title depending on the type of curve
        if type(self.curve) is BsplineCurve:
            self.ax_basis.set_title('B-spline Basis Functions')
        elif type(self.curve) is NurbsCurve:
            self.ax_basis.set_title('NURBS Basis Functions')

        # Set up plot axis labels
        self.ax_basis.set_xlabel('$\\xi$')
        self.ax_basis.set_ylabel('$R_{i,p}(\\xi)$')

        # Set up plot plot ticks
        unique_knots = self.curve.get_unique_knots()
        self.ax_basis.set_xticks(unique_knots)
        self.ax_basis.tick_params(axis='x', which='major', length=10, color='k', width=2, direction='inout')

        # Show minor ticks only if there are no inner knots
        if len(unique_knots) == 2:
            self.ax_basis.xaxis.set_minor_locator(AutoLocator())
            self.ax_basis.tick_params(which='minor', length=4, color='k')

    def plot_basis(self):
        """Plot the basis functions."""
        for i, basis in enumerate(self.curve.basis_evals):
            # Create current basis label
            basis_label = "$R_{" + str(i) + "," + str(self.curve.polynomial_degree) + "}(\\xi)$"

            # Plot the current basis
            art, = self.ax_basis.plot(self.curve.samples, basis, label=basis_label, picker=5)

            # Append artist and pick state to the corresponding attributes
            self.basis_function_art.append(art)

    def curve_plot_setup(self):
        """Set up the curve plot format."""
        # Set up plot title depending on the type of curve
        if type(self.curve) is BsplineCurve:
            self.ax_curve.set_title('B-spline Curve')
        elif type(self.curve) is NurbsCurve:
            self.ax_curve.set_title('NURBS Curve')

        # Set up plot aspect (not available on 3D Axes)
        if not self.is_3D:
            self.ax_curve.set_aspect('equal', adjustable='datalim')

        # Set up plot axis labels
        self.ax_curve.set_xlabel('$x$')
        self.ax_curve.set_ylabel('$y$')
        if self.is_3D:
            self.ax_curve.set_zlabel('$z$')

    def plot_curve(self):
        """Plot the curve."""
        if self.is_3D:
            self.curve_art, = self.ax_curve.plot(self.curve.curve_evals[0],
                                                 self.curve.curve_evals[1],
                                                 self.curve.curve_evals[2],
                                                 visible=True, label='Curve', marker='o', markevery=[])

        else:
            self.curve_art, = self.ax_curve.plot(self.curve.curve_evals[0],
                                                 self.curve.curve_evals[1],
                                                 visible=True, label='Curve', marker='o', markevery=[])

    def plot_control_polygon(self):
        """Plot the control polygon."""
        if self.is_3D:
            self.control_polygon_art, = self.ax_curve.plot(self.curve.control_points[0],
                                                           self.curve.control_points[1],
                                                           self.curve.control_points[2],
                                                           '-.', color='k', linewidth=1, visible=True,
                                                           label='Control polygon',
                                                           marker='s', mec='k', mfc='w', ms=5, zorder=1)

            # Plot markers for interaction
            self.control_polygon_marker_art, = self.ax_curve.plot(self.curve.control_points[0],
                                                                  self.curve.control_points[1],
                                                                  self.curve.control_points[2],
                                                                  ls='none', marker='s', mec='k', mfc='r', ms=5,
                                                                  markevery=[], zorder=3)

        else:
            self.control_polygon_art, = self.ax_curve.plot(self.curve.control_points[0],
                                                           self.curve.control_points[1],
                                                           '-.', color='k', linewidth=1, visible=True,
                                                           label='Control polygon',
                                                           marker='s', mec='k', mfc='w', ms=5, zorder=1)

            # Plot markers for interaction
            self.control_polygon_marker_art, = self.ax_curve.plot(self.curve.control_points[0],
                                                                  self.curve.control_points[1],
                                                                  ls='none', marker='s', mec='k', mfc='r', ms=5,
                                                                  markevery=[], zorder=3)

    def plot_knot_points(self):
        """Plot the curve points corresponding to inner knots."""
        # Get unique knots
        knots = self.curve.get_unique_knots()

        # Get indices of the samples array where the inner knots are placed
        indices = [get_index_to_nearest(self.curve.samples, knot) for knot in knots[1: -1]]

        # Get curve coordinates corresponding to inner knots
        x = self.curve.curve_evals[0][indices]
        y = self.curve.curve_evals[1][indices]

        # Plot the points
        if self.is_3D:
            z = self.curve.curve_evals[2][indices]  # Get the z coordinate of the curve
            self.knot_points_art, = self.ax_curve.plot(x, y, z, ls="none", marker='s', mec='k', mfc='w', ms=5)

        else:
            self.knot_points_art, = self.ax_curve.plot(x, y, ls="none", marker='s', mec='k', mfc='w', ms=5)

    def plot_knot_lines(self):
        """Plot vertical lines corresponding to knots."""
        # Get unique knots
        knots = self.curve.get_unique_knots()

        # Plot a vertical line for every knot
        for knot in knots:
            x = [knot, knot]
            y = [0, 1]
            art, = self.ax_basis.plot(x, y, linestyle='--', color='gray', linewidth=0.5)
            self.knot_lines_art.append(art)

    def on_key_press(self, event):
        """Handle the key press event."""
        if event.key == 'c':
            # Change visibility of curve
            self.curve_art.set_visible(not self.curve_art.get_visible())

        elif event.key == 'p':
            # Change visibility of control polygon
            self.control_polygon_art.set_visible(not self.control_polygon_art.get_visible())

        elif event.key == 'k':
            # Change visibility of curve points corresponding to inner knots
            self.knot_points_art.set_visible(not self.knot_points_art.get_visible())

            # Change visibility of vertical lines corresponding to knots
            for line in self.knot_lines_art:
                line.set_visible(not line.get_visible())

        # Draw plot to show changes
        plt.draw()

    def on_motion(self, event):
        """Handle the motion notify event."""
        if self.ax_basis == event.inaxes:
            # If mouse is over parameter domain
            if self.curve.knot_vector[0] <= event.xdata <= self.curve.knot_vector[-1]:
                # Get index of the sample corresponding to the mouse position
                nearest_point_idx = get_index_to_nearest(self.curve.samples, event.xdata)

                # Visualize marker in nearest_point_idx
                self.curve_art.set_markevery([nearest_point_idx])

            else:
                self.curve_art.set_markevery([])

            # Draw plot to show changes
            plt.draw()

    def on_pick(self, event):
        """Handle the pick event."""
        self.update_picked_basis_index(event)
        self.update_control_point_facecolor()
        self.update_control_point_label()
        self.update_basis_linewidth()
        self.update_basis_label()
        self.update_basis_opacity()

        # Draw plot to show changes
        plt.draw()

    def update_picked_basis_index(self, event):
        """
        Update the picked basis index.

        Arguments
        ---------
        event : mouseevent
            The mouse event.
        """
        # Get index of the new picked basis in the basis_function_art list
        new_picked_basis_index = self.basis_function_art.index(event.artist)

        # If there was no basis picked before or the new picked basis is different the current picked basis
        if self.picked_basis_index is None or self.picked_basis_index != new_picked_basis_index:
            self.picked_basis_index = new_picked_basis_index
        else:
            self.picked_basis_index = None

    def update_control_point_facecolor(self):
        """Update the label of the control point based on the picked basis curve."""
        if self.picked_basis_index is None:
            self.control_polygon_marker_art.set_markevery([])
        else:
            self.control_polygon_marker_art.set_markevery([self.picked_basis_index])

    def update_control_point_label(self):
        """Update the label of the control point based on the picked basis curve."""

        if self.control_point_label is not None:
            self.control_point_label.remove()

        if self.picked_basis_index is not None:
            # Set up the control point label
            label = self.get_current_control_point_label()
            bbox_props = dict(boxstyle='round,pad=0.25', fc='white', alpha=0.5)
            arrow_props = dict(arrowstyle='->', connectionstyle='arc3,rad=0')
            # Get the coordinates and update label
            if self.is_3D:
                x, y, z = self.get_current_control_point_coords()

                # Get the projection of the actual view
                x2, y2, _ = proj3d.proj_transform(x, y, z, self.ax_curve.get_proj())

                # Annotate in the position of the label
                self.control_point_label = self.ax_curve.annotate(label, xy=(x2, y2), xytext=(-20, 20),
                                                                  textcoords='offset points', ha='right', va='bottom',
                                                                  bbox=bbox_props, arrowprops=arrow_props)
            else:
                x, y = self.get_current_control_point_coords()
                self.control_point_label = self.ax_curve.annotate(label, xy=(x, y), xytext=(-20, 20),
                                                                  textcoords='offset points', ha='right', va='bottom',
                                                                  bbox=bbox_props, arrowprops=arrow_props)
        else:
            self.control_point_label = None

    def get_current_control_point_coords(self):
        """
        Get the coordinates of the current (highlighted) control point.

        Returns
        -------
        numpy.ndarray
            The x, y, z coordinates of the current control point.
        """
        return self.curve.control_points[:, self.picked_basis_index]

    def get_current_control_point_label(self):
        """Get the label of the current control point."""
        return '$\\mathbf{P}_{' + str(self.picked_basis_index + 1) + '}$'

    def update_control_point_label_position(self, event):
        """Update control point label position when rotating the surface plot."""
        if self.ax_curve == event.inaxes and self.picked_basis_index is not None:
            # Get the coordinates of the current control point
            x, y, z = self.get_current_control_point_coords()

            # Get the projection of the actual view
            x2, y2, _ = proj3d.proj_transform(x, y, z, self.ax_curve.get_proj())

            # Set the coordinates of the label according to projection
            self.control_point_label.xy = x2, y2
            self.control_point_label.update_positions(self.fig.canvas.renderer)

            # Draw plot to show changes
            plt.draw()

    def update_basis_linewidth(self):
        """Update the linewidth of the picked basis."""
        for i, art in enumerate(self.basis_function_art):
            if i == self.picked_basis_index:
                art.set_linewidth(3)
            elif art.get_linewidth() != 1.5:
                art.set_linewidth(1.5)

    def update_basis_opacity(self):
        """Update the opacity of the basis."""
        for i, art in enumerate(self.basis_function_art):
            if i == self.picked_basis_index or self.picked_basis_index is None:
                if art.get_alpha() != 1:
                    art.set_alpha(1)
            else:
                art.set_alpha(0.3)

    def update_basis_label(self):
        """Update the label of the picked basis."""
        if self.picked_basis_index is None:
            self.ax_basis.set_ylabel('$R_{i,p}(\\xi)$')
        else:
            self.ax_basis.set_ylabel(self.get_current_basis_label())

    def get_current_basis_label(self):
        """Get the label of the current basis function."""
        return '$R_{' + str(self.picked_basis_index + 1) + ',' + str(self.curve.polynomial_degree) + '}(\\xi)$'


class SurfaceVisualizer:
    """
    A class used for visualization of B-spline and NURBS surfaces.

    Attributes
    ----------
    surface : BsplineSurface or NurbsSurface
        The surface object to visualize.
    """

    def __init__(self, surface):
        """
        Construct B-spline or NURBS surfaces.

        Arguments
        ---------
        surface : BsplineSurface or NurbsSurface
            The surface object to visualize.
        """
        self.surface = surface

        # Create figure for the plots
        self.fig = plt.figure()
        self.set_title()

        # Initialize attributes to store the matplotlib artists and states for future manipulation
        self.control_point_label = None
        self.surface_art = None
        self.control_polygon_art = [[], []]
        self.basis_function_art = None
        self.current_basis_i_idx = 0
        self.current_basis_j_idx = 0

        # Create basis plot
        self.ax_basis = self.fig.add_subplot(122, projection='3d')
        self.ax_basis.view_init(elev=30, azim=45)
        self.basis_plot_setup()
        self.plot_basis()

        # Create surface plot
        self.ax_surface = self.fig.add_subplot(121, projection='3d')
        self.ax_surface.view_init(elev=30, azim=45)
        self.surface_plot_setup()
        self.plot_surface()
        self.plot_control_polygon()

        # Define events for interactive plots
        self.scroll_event = self.fig.canvas.mpl_connect('scroll_event', self.on_scroll)
        self.key_press_event = self.fig.canvas.mpl_connect('key_press_event', self.on_key_press)
        self.motion_notify_event = self.fig.canvas.mpl_connect('motion_notify_event',
                                                               self.update_control_point_label_position)

        # Show the plot
        plt.tight_layout(pad=3, w_pad=2.5)
        plt.show()

    def set_title(self):
        """Set the window title depending on the type of surface."""
        if type(self.surface) is BsplineSurface:
            self.fig.canvas.set_window_title('B-spline Surface Visualizer')
        if type(self.surface) is NurbsSurface:
            self.fig.canvas.set_window_title('NURBS Surface Visualizer')

    def basis_plot_setup(self):
        """Set up the basis plot format."""
        # Set up plot title depending on the type of surface
        if type(self.surface) is BsplineSurface:
            self.ax_basis.set_title('B-spline Basis Functions')
        if type(self.surface) is NurbsSurface:
            self.ax_basis.set_title('NURBS Basis Functions')

        # Set up plot axis labels
        self.ax_basis.set_xlabel('$\\xi$')
        self.ax_basis.set_ylabel('$\\eta$')
        self.ax_basis.xaxis.set_rotate_label(False)
        self.ax_basis.yaxis.set_rotate_label(False)
        self.ax_basis.zaxis.set_rotate_label(False)

        # Set up z axis limits
        self.ax_basis.set_zlim(0, 1)

        # Set up plot plot ticks
        self.ax_basis.set_xticks(self.surface.get_unique_knots()[0])
        self.ax_basis.set_yticks(self.surface.get_unique_knots()[1])

    def plot_basis(self):
        """Plot the current basis function."""
        # Get the current basis label
        label = self.get_current_basis_label()

        # Set the z axis label
        self.ax_basis.set_zlabel(label)

        # Plot the current basis function
        xi_grid, eta_grid = np.meshgrid(self.surface.samples[0], self.surface.samples[1])
        self.basis_function_art = self.ax_basis.plot_surface(xi_grid, eta_grid,
                                                             self.surface.basis_evals[self.current_basis_i_idx][
                                                                 self.current_basis_j_idx], antialiased=False)

    def get_current_basis_label(self):
        """Get the label of the current basis function."""
        return '$R_{' + str(self.current_basis_i_idx + 1) + str(self.current_basis_j_idx + 1) + '}(\\xi, \\eta)$'

    def surface_plot_setup(self):
        """Set up the surface plot format."""
        # Set up plot title depending on the type of surface
        if type(self.surface) is BsplineSurface:
            self.ax_surface.set_title('B-spline Surface')
        if type(self.surface) is NurbsSurface:
            self.ax_surface.set_title('NURBS Surface')

        # Set up plot labels
        self.ax_surface.set_xlabel('$x$')
        self.ax_surface.set_ylabel('$y$')
        self.ax_surface.set_zlabel('$z$')

    def plot_surface(self):
        """Plot the surface."""
        self.surface_art = self.ax_surface.plot_surface(self.surface.surface_evals[0],
                                                        self.surface.surface_evals[1],
                                                        self.surface.surface_evals[2],
                                                        cmap=cm.jet, alpha=0.8, antialiased=False)
        self.update_control_point_label()

    def update_control_point_label(self, direction=None):
        """Update the label of the current control point."""

        # Remove the previous control point label
        if self.control_point_label is not None:
            self.control_point_label.remove()

        # Get the coordinates of the current control point
        x, y, z = self.get_current_control_point_coords()

        # Update control point label
        label = self.get_current_control_point_label()

        # Get the projection of the actual view
        x2, y2, _ = proj3d.proj_transform(x, y, z, self.ax_surface.get_proj())

        # Annotate in the position of the label
        bbox_props = dict(boxstyle='round,pad=0.25', fc='white', alpha=0.5)
        arrow_props = dict(arrowstyle='->', connectionstyle='arc3,rad=0')

        self.control_point_label = self.ax_surface.annotate(label, xy=(x2, y2), xytext=(-20, 20),
                                                            textcoords='offset points', ha='right', va='bottom',
                                                            bbox=bbox_props, arrowprops=arrow_props)

    def get_current_control_point_coords(self):
        """
        Get the coordinates of the current (highlighted) control point.

        Returns
        -------
        numpy.ndarray
            The x, y, z coordinates of the current control point.
        """
        return self.surface.control_points[:, self.current_basis_i_idx, self.current_basis_j_idx]

    def update_control_point_label_position(self, event):
        """Update control point label position when rotating the surface plot."""

        # Get the coordinates of the current control point
        x, y, z = self.get_current_control_point_coords()

        # Get the projection of the actual view
        x2, y2, _ = proj3d.proj_transform(x, y, z, self.ax_surface.get_proj())

        # Set the coordinates of the label according to projection
        self.control_point_label.xy = x2, y2
        self.control_point_label.update_positions(self.fig.canvas.renderer)

        # Draw plot to show changes
        plt.draw()

    def get_current_control_point_label(self):
        """Get the label of the current control point."""
        return '$\\mathbf{P}_{' + str(self.current_basis_i_idx + 1) + str(self.current_basis_j_idx + 1) + '}$'

    def plot_control_polygon(self):
        """Plot the control polygon."""
        # Plot lines connecting control points in the first direction
        for row in range(self.surface.control_points_count[0]):
            art, = self.ax_surface.plot(self.surface.control_points[0, row, :],
                                        self.surface.control_points[1, row, :],
                                        self.surface.control_points[2, row, :],
                                        "-.", color="k", linewidth=1,
                                        marker='s', mec='k', mfc='r', ms=5, markevery=[], zorder=2)
            self.control_polygon_art[0].append(art)

        # Plot lines connecting control points in the second direction
        for column in range(self.surface.control_points_count[1]):
            art, = self.ax_surface.plot(self.surface.control_points[0, :, column],
                                        self.surface.control_points[1, :, column],
                                        self.surface.control_points[2, :, column],
                                        "-.", color="k", linewidth=1,
                                        marker='s', mec='k', mfc='w', ms=5, zorder=1)
            self.control_polygon_art[1].append(art)

        # Highlight the current control point
        self.control_polygon_art[0][self.current_basis_i_idx].set_markevery([self.current_basis_j_idx])

    def on_scroll(self, event):
        """Handle the scroll event."""
        # Get direction of scroll (up or down)
        direction = event.button  # Scroll can be with direction up or down

        # Update plots based on current basis and direction
        self.update_plots_new_basis(direction)

    def update_plots_new_basis(self, direction=None):
        """
        Update the plots based on a direction.

        Arguments
        ---------
        direction : str, 'up' or 'down', default=None
            The direction to update the plots.
        """
        # Remove current basis from plot
        self.basis_function_art.remove()

        # Update plots for the new basis
        self.update_basis_idx(direction)
        self.plot_basis()
        self.update_control_polygon_marker(direction)
        self.update_control_point_label()

        # Draw plot to show changes
        plt.draw()

    def update_basis_idx(self, direction=None):
        """
        Update the basis index based on a direction.

        Arguments
        ---------
        direction : str, 'up' or 'down', default=None
            The direction to update the basis index.
        """
        # Get number of control points in each parameter space direction
        points_count_xi, points_count_eta = self.surface.control_points_count

        # Compute maximum index in each parameter space direction
        i_max = points_count_xi - 1
        j_max = points_count_eta - 1

        # Update index based on current index and direction
        if direction == 'up':
            if self.current_basis_j_idx == j_max:
                self.current_basis_j_idx = 0
                if self.current_basis_i_idx == i_max:
                    self.current_basis_i_idx = 0
                else:
                    self.current_basis_i_idx += 1
            else:
                self.current_basis_j_idx += 1
        elif direction == 'down':
            if self.current_basis_j_idx == 0:
                self.current_basis_j_idx = j_max
                if self.current_basis_i_idx == 0:
                    self.current_basis_i_idx = i_max
                else:
                    self.current_basis_i_idx -= 1
            else:
                self.current_basis_j_idx -= 1

    def update_control_polygon_marker(self, direction=None):
        """
        Update the control polygon marker.

        Arguments
        ---------
        direction : str, 'up' or 'down', default=None
            The direction to update the control polygon marker.
        """
        # Get number of control points in each parameter space direction
        points_count_xi, points_count_eta = self.surface.control_points_count

        # Compute maximum index in each parameter space direction
        i_max = points_count_xi - 1
        j_max = points_count_eta - 1

        # Highlight control polygon marker based on current index and direction
        if direction == 'up':
            if self.current_basis_j_idx == 0:
                if self.current_basis_i_idx == 0:
                    self.control_polygon_art[0][i_max].set_markevery([])
                else:
                    self.control_polygon_art[0][self.current_basis_i_idx - 1].set_markevery([])
            else:
                self.control_polygon_art[0][self.current_basis_i_idx].set_markevery([])
        elif direction == 'down':
            if self.current_basis_j_idx == j_max:
                if self.current_basis_i_idx == i_max:
                    self.control_polygon_art[0][0].set_markevery([])
                else:
                    self.control_polygon_art[0][self.current_basis_i_idx + 1].set_markevery([])
            else:
                self.control_polygon_art[0][self.current_basis_i_idx].set_markevery([])
        self.control_polygon_art[0][self.current_basis_i_idx].set_markevery([self.current_basis_j_idx])

    def on_key_press(self, event):
        """Handle the key press event."""
        if event.key == 's':
            # Change visibility of surface
            self.surface_art.set_visible(not self.surface_art.get_visible())

        elif event.key == 'p':
            # Change visibility of control polygon
            for direction in self.control_polygon_art:
                for polygon in direction:
                    polygon.set_visible(not polygon.get_visible())

        elif event.key == 'up' or 'down':
            # Update plots based on current basis and direction
            self.update_plots_new_basis(event.key)

        # Draw plot to show changes
        plt.draw()

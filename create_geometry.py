"""
Class object which discretises the rotor blade(s) into bound and trailing vortex
filaments
"""
import numpy as np


class BladeGeometry:

    def __init__(self, radius, tsr, v_inf, n_blades, n_span, n_theta):
        # todo: check whether non-dim of span_arr is needed
        self.radius = radius
        self.tsr = tsr
        self.v_inf = v_inf
        self.n_blades = n_blades
        self.n_span = n_span
        self.n_theta = n_theta
        self.span_arr = np.linspace(0.2, 1.0, n_span)
        self.theta_arr = np.linspace(0, 2 * np.pi, n_theta)

        self.cp = np.zeros((3, n_span, 10))  # coord; normal; tangential
        self.rings = {"x1": np.zeros((n_blades * n_span, 2 * n_theta + 2)),
                      "x2": np.zeros((n_blades * n_span, 2 * n_theta + 2)),
                      "y1": np.zeros((n_blades * n_span, 2 * n_theta + 2)),
                      "y2": np.zeros((n_blades * n_span, 2 * n_theta + 2)),
                      "z1": np.zeros((n_blades * n_span, 2 * n_theta + 2)),
                      "z2": np.zeros((n_blades * n_span, 2 * n_theta + 2)),
                      "gamma": np.zeros((n_blades, n_span, 2 * n_theta + 2))
                      }
        self.bladepanels = np.zeros((n_blades, n_span, 4 * 3))  # empty dict to
        # store
        # blade geometry

    def discretize_blade(self):
        # TODO: solve control points

        # TODO: solve rings

        # TODO: solve blade panels

        raise NotImplementedError

    def _compute_ring(self):
        # TODO: compute bound vortex filaments

        # TODO: compute trailing vortex filaments

        # TODO: rotate all filaments

        # TODO: redefine in self.rings

        raise NotImplementedError

    def _compute_cp(self):
        # TODO: compute coordinates

        # TODO: compute tangential unit vector

        # TODO: compute normal unit vector

        raise NotImplementedError
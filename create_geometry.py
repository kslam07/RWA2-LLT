import numpy as np


"""
Class object which discretises the rotor blade(s) into bound and trailing vortex
filaments
"""

class BladeGeometry:

    def __init__(self, radius, tsr, v_inf, n_blades, n_span, n_theta):

        # todo: check whether non-dim of span_arr is needed
        self.radius = radius
        self.tsr = tsr
        self.v_inf = v_inf
        self.n_blades = n_blades
        self.span_arr = np.linspace(0.2, 1.0, n_span)
        self.theta_arr = np.linspace(0, 2*np.pi, n_theta)

        self.cp = np.zeros((3, n_span, 3))  # coord; normal; tangential
        self.rings = {"x1": np.zeros((n_theta, n_theta, n_blades)),
                      "x2": np.zeros((n_theta, n_theta, n_blades)),
                      "y1": np.zeros((n_theta, n_theta, n_blades)),
                      "y2": np.zeros((n_theta, n_theta, n_blades)),
                      "z1": np.zeros((n_theta, n_theta, n_blades)),
                      "z2": np.zeros((n_theta, n_theta, n_blades))
                      }
        self.bladepanels = {}  # empty dict to store blade geometry

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
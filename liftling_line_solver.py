"""
Runs frozen vortex wake solver // lifting line theory
"""

import numpy as np


class LiftingLineSolver:

    def __init__(self, geo, u_rot, r_rotor, weight=0.3, tol=1e-6, n_iter=1000):
        """

        :param geo: BladeGeometry class
        :param u_rot: rotational velocity [rad/s]
        :param r_rotor: rotor radius [m]
        :param weight: weighting factor for next step
        :param tol: stopping criteria
        :param n_iter: number of iterations
        """
        self.geo = geo
        self.u_rot = u_rot
        self.r_rotor = r_rotor
        self.weight = weight
        self.tol = tol
        self.n_iter = n_iter

        self.gamma = np.zeros((geo.n_blades, geo.n_theta, 1))

    def _compute_circ(self):
        # TODO: compute circulation of each filament
        raise NotImplementedError

    def _compute_induced_velocity(self):
        # TODO: compute induced velocity of rings on control points
        raise NotImplementedError

    def _compute_loads_blade(self):
        # TODO: compute normal, tangential force and circulation of the
        #  control points
        raise NotImplementedError

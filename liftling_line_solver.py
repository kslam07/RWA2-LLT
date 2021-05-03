"""
Runs frozen vortex wake solver // lifting line theory
"""

import numpy as np
import pandas as pd

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

        self.geo.rings["gamma"] = self.geo.rings["gamma"]*(1-self.weight) + (self.weight*self.gamma)

        return self.geo.rings

    def _compute_induced_velocity(self):

        # TODO: compute induced velocity of rings on control points

        raise NotImplementedError

    def _geo_blade(self):

        radial = [0, 0.3, .5, .8, 1]

        # chord_dist = [.05, .04, .03, .02, .015]
        # twist_dist = [-12, -8, -5, -4, 0.]

        pitch = 2

        chord = 3 * (1 - self.r_rotor) + 1

        twist = -14 * (1 - self.r_rotor)

        result = [chord, twist + pitch]

        return result

    def _polarAirfoil(self):

        data = pd.read_excel('polar_DU95W180.xlsx').to_numpy()
        data = data[3:]
        data = np.array(data, dtype=float)

        alphaRad = np.radians(data[:, 0])
        fcl = np.polyfit(alphaRad, data[:, 1], 3)
        fcd = np.polyfit(alphaRad, data[:, 2], 3)
        amax = max(alphaRad)
        amin = min(alphaRad)

        return [alphaRad, fcl, fcd, amax, amin]

    def _compute_loads_blade(self, V_norm=20, V_tan=10):

        V_mag2 = (V_norm**2 + V_tan**2)         # Velocity magnitude squared
        phi = np.atan(V_norm/V_tan)            # Inflow angle

        # Get chord and twist
        [chord, twist] = self._geo_blade()

        alpha = twist + phi*180/np.pi

        [alphaRad, fcl, fcd, amax, amin] = self._polarAirfoil()

        cl = fcl[0] + fcl[1] * alpha + fcl[2] * (alpha**2) + fcl[3] * (alpha**3)
        cd = fcd[0] + fcd[1] * alpha + fcd[2] * (alpha**2) + fcd[3] * (alpha**3)

        L = 0.5 * V_mag2 * cl * chord
        D = 0.5 * V_mag2 * cd * chord

        F_norm = L * np.cos(phi) + D * np.sin(phi)
        F_tan = L * np.sin(phi) - D * np.cos(phi)

        Gamma = 0.5 * np.sqrt(V_mag2) * cl * chord

        return [F_norm, F_tan, Gamma]

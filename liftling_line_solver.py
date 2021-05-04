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

        # rotor discretization
        self.geo = geo
        # rotor properties
        self.u_rot = u_rot
        self.r_rotor = r_rotor

        # solver settings
        self.weight = weight
        self.tol = tol
        self.n_iter = n_iter

        # airfoil variables
        self._polarAirfoil()  # computes minimum/maximum alpha, polynomials for CL and CD alpha approximations

    def _compute_circ(self, gamma, weight):
        self.geo.rings[-1] = self.geo.rings[-1] * (1 - weight) + (weight * gamma)

    def _velocity_3D_from_vortex_filament(self, core):
        r_gamma = self.geo.rings[-1]

        x1 = self.geo.rings[0]
        y1 = self.geo.rings[2]
        z1 = self.geo.rings[4]
        x2 = self.geo.rings[1]
        y2 = self.geo.rings[3]
        z2 = self.geo.rings[5]
        xc = self.geo.cp[0]
        yc = self.geo.cp[1]
        zc = self.geo.cp[2]

        R1 = np.sqrt((xc - x1) ** 2 + (yc - y1) ** 2 + (zc - z1) ** 2)
        R2 = np.sqrt((xc - x2) ** 2 + (yc - y2) ** 2 + (zc - z2) ** 2)

        R12_xx = ((yc - y1) * (zc - z2)) - ((zc - z1) * (yc - y2))
        R12_xy = -((xc - x1) * (zc - z2)) + ((zc - z1) * (xc - x2))
        R12_xz = ((xc - x1) * (yc - y2)) - ((yc - y1) * (xc - x2))

        R12_sq = (R12_xx ** 2) + (R12_xy ** 2) + (R12_xz ** 2)

        R01 = ((x2 - x1) * (xc - x1)) + ((y2 - y1) * (yc - y1)) + ((z2 - z1) * (zc - z1))
        R02 = ((x2 - x1) * (xc - x2)) + ((y2 - y1) * (yc - y2)) + ((z2 - z1) * (zc - z2))

        # check if target point is in the vortex filament core,
        # and modify to solid body rotation
        R12_sq = np.where(R12_sq < core ** 2, core ** 2, R12_sq)
        R1 = np.where(R1 < core, core, R1)
        R2 = np.where(R1 < core, core, R2)

        K = r_gamma / 4 / np.pi / R12_sq * (R01 / R1 - R02 / R2)

        U = K * R12_xx
        V = K * R12_xy
        W = K * R12_xz

        return [U, V, W]

    def _compute_induced_velocity(self):
        V_ind = np.zeros(3)
        core = 0.00001
        temp_v = self._velocity_3D_from_vortex_filament(core)
        V_ind += temp_v

        return V_ind

    def _geo_blade(self):
        # radial = [0, 0.3, .5, .8, 1]
        # chord_dist = [.05, .04, .03, .02, .015]
        # twist_dist = [-12, -8, -5, -4, 0.]

        pitch = 2  # in deg

        chord = 3 * (1 - self.r_rotor) + 1

        twist = -14 * (1 - self.r_rotor)  # in deg

        result = [chord, np.radians(twist + pitch)]

        return result

    def _polarAirfoil(self):
        data = pd.read_excel('polar_DU95W180.xlsx').to_numpy()
        data = data[3:]
        data = np.array(data, dtype=float)

        alphaRad = np.radians(data[:, 0])
        self.amax = max(alphaRad)
        self.amin = min(alphaRad)
        self.fcl = np.polyfit(alphaRad, data[:, 1], 3)
        self.fcd = np.polyfit(alphaRad, data[:, 2], 3)

    def _compute_loads_blade(self, v_norm=20, v_tan=10):
        V_mag2 = (v_norm ** 2 + v_tan ** 2)  # Velocity magnitude squared
        phi = np.arctan(v_norm / v_tan)  # Inflow angle

        # Get chord and twist
        [chord, twist] = self._geo_blade()

        alpha = twist + phi * 180 / np.pi

        # apply alpha constraints
        alpha[alpha < self.amin] = self.amin
        alpha[alpha > self.amax] = self.amax

        cl = np.polyval(self.fcl, alpha)
        cd = np.polyval(self.fcd, alpha)

        L = 0.5 * V_mag2 * cl * chord
        D = 0.5 * V_mag2 * cd * chord

        F_norm = L * np.cos(phi) + D * np.sin(phi)
        F_tan = L * np.sin(phi) - D * np.cos(phi)

        Gamma = 0.5 * np.sqrt(V_mag2) * cl * chord

        return [F_norm, F_tan, Gamma]

    def _initialize_solver(self):
        uvw_mat = np.zeros((3, self.geo.n_span, self.geo.n_span))
        gamma_new = np.shape(self.geo.rings[-1].shape)  # initialize gamma array to store new circulation

        # update Gamma given Gamma matrix, weight, and new Gamma
        self._compute_circ(gamma_new, self.weight)  # updates self.geo itself

        # compute [ui, vi, wi] based on vortex strength and distance
        # between control point and vortex
        v_induced = self._compute_induced_velocity()
        uvw_mat[0] = v_induced[0]  # u
        uvw_mat[1] = v_induced[1]  # v
        uvw_mat[2] = v_induced[2]  # w

        return uvw_mat

    def run_solver(self):
        for i in range(self.n_iter):
            # todo: update circulation

            # todo: calculate velocity, circulation, control points
            pass

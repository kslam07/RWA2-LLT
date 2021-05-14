"""
Class object which discretises the rotor blade(s) into bound and trailing vortex
filaments
"""
import numpy as np
import matplotlib.pyplot as plt


class BladeGeometry:

    def __init__(self, radius, tsr, v_inf, n_blades, n_span, n_theta, spacing,
                 a, xshift=0, yshift=100, zshift=0):
        # todo: check whether non-dim of span_arr is needed
        self.radius = radius
        self.tsr = tsr
        self.v_inf = v_inf
        self.n_blades = n_blades
        self.n_span = n_span
        self.n_theta = n_theta
        self.a = a
        self.span_arr = np.linspace(0.2, 1.0, n_span)
        self.theta_arr = np.linspace(0, 100 * np.pi, n_theta)
        self.cp = np.zeros((n_blades * (n_span - 1), 9))  # coord; normal; tangential
        self.bladepanels = np.zeros((n_blades * (n_span-1), 4 * 3))  # empty dict to
        
        self.xshift = xshift
        self.yshift = yshift
        self.zshift = zshift

        self.rRootRatio = 0.2
        self.rTipRatio = 1.0
        self.spacing = spacing
        self.filaments = np.zeros((7, n_blades * (n_span - 1), 2 * n_theta + 1))
        self.discretize_spanwise()
        self.discretize_blade()
        self.compute_ring()
        self._compute_cp()

    def discretize_spanwise(self):
        if self.spacing == 'equal':
            rsegment = np.linspace(self.rRootRatio, 1, self.n_span)

        elif self.spacing == 'cosine':
            midpoint = (1 - self.rRootRatio) / 2
            angleInc = np.pi / (self.n_span - 1)
            curAngle = angleInc
            rsegment = np.zeros(self.n_span)
            rsegment[0] = self.rRootRatio
            rsegment[-1] = 1
            for idx in range(1, self.n_span):
                rsegment[idx] = self.rRootRatio + midpoint * (1 - np.cos(curAngle))
                curAngle += angleInc

        else:
            print("Choose equal or cosine spacing")

        self.span_arr = rsegment
        return

    def discretize_blade(self):
        for blade in range(self.n_blades):
            bladeRot = 2 * np.pi / self.n_blades * blade
            angle1 = np.deg2rad(- 14 * (1 - self.span_arr[:-1]) + 2)
            angle2 = np.deg2rad(- 14 * (1 - self.span_arr[1:]) + 2)
            chord1 = (3 * (1 - self.span_arr[:-1]) + 1)/self.radius
            chord2 = (3 * (1 - self.span_arr[1:]) + 1)/self.radius

            # define the 4 corners
            p1 = [-0.25 * chord1 * np.sin(-angle1), self.span_arr[:-1], 0.25 * chord1 * np.cos(angle1)]
            p2 = [-0.25 * chord2 * np.sin(-angle2), self.span_arr[1:], 0.25 * chord2 * np.cos(angle2)]
            p3 = [0.75 * chord2 * np.sin(-angle2), self.span_arr[1:], -0.75 * chord2 * np.cos(angle2)]
            p4 = [0.75 * chord1 * np.sin(-angle1), self.span_arr[:-1], -0.75 * chord1 * np.cos(angle1)]

            # rotate coordinates
            p1Rot = np.column_stack([np.zeros(self.n_span-1), p1[1] * np.cos(bladeRot) - p1[2] * np.sin(bladeRot),
                                     p1[1] * np.sin(bladeRot) + p1[2] * np.cos(bladeRot)])
            p2Rot = np.column_stack([np.zeros(self.n_span-1), p2[1] * np.cos(bladeRot) - p2[2] * np.sin(bladeRot),
                                     p2[1] * np.sin(bladeRot) + p2[2] * np.cos(bladeRot)])
            p3Rot = np.column_stack([np.zeros(self.n_span-1), p3[1] * np.cos(bladeRot) - p3[2] * np.sin(bladeRot),
                                     p3[1] * np.sin(bladeRot) + p3[2] * np.cos(bladeRot)])
            p4Rot = np.column_stack([np.zeros(self.n_span-1), p4[1] * np.cos(bladeRot) - p4[2] * np.sin(bladeRot),
                                     p4[1] * np.sin(bladeRot) + p4[2] * np.cos(bladeRot)])

            # write to bladepanels
            self.bladepanels[blade * (self.n_span-1):blade * (self.n_span-1) + (self.n_span-1), :] = \
                np.column_stack((p1Rot, p2Rot, p3Rot, p4Rot))*self.radius
                
        return

    def compute_ring(self):
        # loop over different blades
        for blade in range(self.n_blades):
            r = self.span_arr
            bladeRot = 2 * np.pi / self.n_blades * blade

            # loop             
            for idx, span in enumerate(self.span_arr[:-1]):
                data_arr = np.empty((0, 7))
                chord1 = (3 * (1 - span) + 1)/self.radius
                angle1 = np.deg2rad(- 14 * (1 - span) + 2)

                # Define bound filament
                # [x1, y1, z1, x2 ,y2, z2, gamma]
                data_arr = np.vstack((data_arr, [0, r[idx], 0, 0, r[idx + 1], 0, 0]))

                # Define trailing filaments at x1
                data = [chord1 * np.sin(-angle1), r[idx], -chord1 * np.cos(angle1), 0, r[idx], 0, 0]
                data_arr = np.vstack((data_arr, data))

                # Loop over trailing filaments at x1
                for index, theta in enumerate(self.theta_arr[:-1]):
                    xt = data_arr[-1, 0]
                    yt = data_arr[-1, 1]
                    zt = data_arr[-1, 2]
                    dx = (self.theta_arr[index + 1] - theta) / (self.tsr*self.v_inf/self.v_inf/(1-self.a))
                    dy = (np.cos(-self.theta_arr[index + 1]) - np.cos(-theta)) * r[idx]
                    dz = (np.sin(-self.theta_arr[index + 1]) - np.sin(-theta)) * r[idx]
                    data_arr = np.vstack((data_arr, [xt + dx, (yt + dy), (zt + dz), xt, yt, zt, 0]))

                # trailing filaments at x2
                chord2 = (3 * (1 - r[idx + 1]) + 1)/self.radius
                angle2 = np.deg2rad(- 14 * (1 - r[idx + 1]) + 2)
                data = [0, r[idx + 1], 0, chord2 * np.sin(-angle2), r[idx + 1], -chord2 * np.cos(angle2), 0]
                data_arr = np.vstack((data_arr, data))

                # Loop over trailing filaments at x2
                for index, theta in enumerate(self.theta_arr[:-1]):
                    xt = data_arr[-1, 3]
                    yt = data_arr[-1, 4]
                    zt = data_arr[-1, 5]
                    dx = (self.theta_arr[index + 1] - theta) / (self.tsr*self.v_inf/self.v_inf/(1-self.a))
                    dy = (np.cos(-self.theta_arr[index + 1]) - np.cos(-theta)) * r[idx + 1]
                    dz = (np.sin(-self.theta_arr[index + 1]) - np.sin(-theta)) * r[idx + 1]
                    data_arr = np.vstack((data_arr, [xt, yt, zt, xt + dx, (yt + dy), (zt + dz), 0]))

                # rotate the filaments to correspond with blade orientation
                for filament in data_arr:
                    temp1 = filament.copy()  # store non-rotated filament for subsequent rotation
                    filament[1] = temp1[1] * np.cos(bladeRot) - temp1[2] * np.sin(bladeRot)  # y1
                    filament[2] = temp1[1] * np.sin(bladeRot) + temp1[2] * np.cos(bladeRot)  # z1
                    filament[4] = temp1[4] * np.cos(bladeRot) - temp1[5] * np.sin(bladeRot)  # y2
                    filament[5] = temp1[4] * np.sin(bladeRot) + temp1[5] * np.cos(bladeRot)  # z2

                self.filaments[:, blade * (self.n_span - 1) + idx, :] = data_arr.T * self.radius

        return
    

    def _compute_cp(self):
        segmentCenter = self.radius*(self.span_arr + (np.roll(self.span_arr, -1) - self.span_arr) / 2)[:-1]
        self.centerPoints = segmentCenter
        for blade in range(self.n_blades):
            bladeRot = 2 * np.pi / self.n_blades * blade
            angle = np.deg2rad(2 - 14 * (1 - segmentCenter))

            # bound edge
            boundEdge = np.column_stack(
                (np.zeros(self.n_span - 1), segmentCenter * np.cos(bladeRot), segmentCenter * np.sin(bladeRot)))
            tangVect = np.column_stack(
                (np.cos(angle), np.sin(angle) * np.sin(bladeRot), -np.sin(angle) * np.cos(bladeRot)))
            normVect = np.column_stack(
                (-np.sin(angle), np.cos(angle) * np.sin(bladeRot), -np.cos(angle) * np.cos(bladeRot)))
            # Assign to cp
            cp = np.column_stack((boundEdge, tangVect, normVect))
            # return [coord, norm, tang] x,y,z
            self.cp[blade * (self.n_span - 1):blade * (self.n_span - 1) + self.n_span - 1, :] = cp
        return
    
    def doubleRotor(self):
        # shift control points
        rotor2=self.cp.copy()
        rotor2[:,0]+=self.xshift
        rotor2[:,1]+=self.yshift
        rotor2[:,2]+=self.zshift
        self.cp=np.vstack((self.cp,rotor2))
        
        # shift blade 
        rotor2=self.bladepanels.copy()
        rotor2[:,(0,3,6,9)]+=self.xshift
        rotor2[:,(1,4,7,10)]+=self.yshift
        rotor2[:,(2,5,8,11)]+=self.zshift
        self.bladepanels=np.vstack((self.bladepanels,rotor2))
        
        # shift filaments
        f=self.filaments
        rotor2=self.filaments.copy()
        rotor2[(0,3),:,:]+=self.xshift
        rotor2[(1,4),:,:]+=self.yshift
        rotor2[(2,5),:,:]+=self.zshift
        
        x1=self.bladepanels=np.vstack((f[0],rotor2[0])).T
        y1=self.bladepanels=np.vstack((f[1],rotor2[1])).T
        z1=self.bladepanels=np.vstack((f[2],rotor2[2])).T
        x2=self.bladepanels=np.vstack((f[3],rotor2[3])).T
        y2=self.bladepanels=np.vstack((f[4],rotor2[4])).T
        z2=self.bladepanels=np.vstack((f[5],rotor2[5])).T
        g=self.bladepanels=np.vstack((f[6],rotor2[6])).T
        self.filaments=np.dstack((x1,y1,z1,x2,y2,z2,g)).T
    
    def doubleRotorUpdate(self):
        # shift filaments
        idxMid = int(np.shape(self.filaments)[1]/2-1)
        self.filaments[(0,3),idxMid:,:][:idxMid]+=self.xshift
        self.filaments[(1,4),idxMid:,:][:idxMid]+=self.yshift
        self.filaments[(2,5),idxMid:,:][:idxMid]+=self.zshift
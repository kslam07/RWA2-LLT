"""
Class object which discretises the rotor blade(s) into bound and trailing vortex
filaments
"""
import numpy as np
import matplotlib.pyplot as plt


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
        self.cp = np.zeros((n_blades*(n_span-1),9)) # coord; normal; tangential
        self.rings = {"x1": np.zeros((n_blades * (n_span-1), n_theta-1)),
                      "x2": np.zeros((n_blades * (n_span-1), n_theta-1)),
                      "y1": np.zeros((n_blades * (n_span-1), n_theta-1)),
                      "y2": np.zeros((n_blades * (n_span-1), n_theta-1)),
                      "z1": np.zeros((n_blades * (n_span-1), n_theta-1)),
                      "z2": np.zeros((n_blades * (n_span-1), n_theta-1)),
                      "gamma": np.zeros((n_blades, (n_span-1), n_theta-1))
                      }
        self.bladepanels = np.zeros((n_blades*n_span, 4 * 3))  # empty dict to
        

        self.rRootRatio = 0.2
        self.rTipRatio = 1.0
        self.spacing = 'cosine'
        self.discretize_spanwise()
        self.filaments=[]
        self.test=np.zeros((7,n_blades*(n_span-1),2*n_theta+1))
        # self.bladepanels = {}  # empty dict to store blade geometry
        
    def discretize_spanwise(self):
        if self.spacing == 'equal':
            rsegment = np.linspace(self.rRootRatio,1,self.n_span)
            
        elif self.spacing == 'cosine':
            midpoint = (1-self.rRootRatio)/2
            angleInc = np.pi/(self.n_span-1)
            curAngle = angleInc;
            rsegment = np.zeros(self.n_span)
            rsegment[0] = self.rRootRatio
            rsegment[-1] = 1
            for idx in range(1,self.n_span):
                rsegment[idx] = self.rRootRatio+midpoint*(1-np.cos(curAngle))
                curAngle+=angleInc;

        else:
            print ("Choose equal or cosine spacing")
            
        self.span_arr=rsegment
        return 

    def discretize_blade(self):
        # TODO: solve control points
        data_arr=np.empty((0,4*3))
        segmentCenter=(self.span_arr+(np.roll(self.span_arr,-1)-self.span_arr)/2)[:-1]
        self._compute_cp()
        for blade in range(self.n_blades):
            bladeRot = 2 * np.pi / self.n_blades*blade
            angle = np.deg2rad(3 * (1 - segmentCenter) + 1 - 14 * (1 - segmentCenter))
        
        # TODO: solve rings
        # TODO: solve blade panels
            angle1 = np.deg2rad(- 14 * (1 - self.span_arr) + 2)
            angle2 = np.deg2rad(- 14 * (1-np.roll(self.span_arr,-1)) + 2)
            chord1 = 3 * (1 - self.span_arr) + 1
            chord2 = 3 * (1 - np.roll(self.span_arr,-1)) + 1

            # define the 4 corners
            p1=[-0.25*chord1*np.sin(-angle1),self.span_arr,0.25*chord1*np.cos(angle1)]
            p2=[-0.25*chord2*np.sin(-angle2),np.roll(self.span_arr,-1),0.25*chord2*np.cos(angle2)]
            p3=[0.75*chord2*np.sin(-angle2),np.roll(self.span_arr,-1),-0.75*chord2*np.cos(angle2)]
            p4=[0.75*chord1*np.sin(-angle1),self.span_arr,-0.75*chord1*np.cos(angle1)]

            # rotate coordinates
            p1Rot=np.column_stack([np.zeros(self.n_span),p1[1]*np.cos(bladeRot)-p1[2]*np.sin(bladeRot),p1[1]*np.sin(bladeRot)+p1[2]*np.cos(bladeRot)])
            p2Rot=np.column_stack([np.zeros(self.n_span),p2[1]*np.cos(bladeRot)-p2[2]*np.sin(bladeRot),p2[1]*np.sin(bladeRot)+p2[2]*np.cos(bladeRot)])
            p3Rot=np.column_stack([np.zeros(self.n_span),p3[1]*np.cos(bladeRot)-p3[2]*np.sin(bladeRot),p3[1]*np.sin(bladeRot)+p3[2]*np.cos(bladeRot)])
            p4Rot=np.column_stack([np.zeros(self.n_span),p4[1]*np.cos(bladeRot)-p4[2]*np.sin(bladeRot),p4[1]*np.sin(bladeRot)+p4[2]*np.cos(bladeRot)])

            # write to bladepanels
            self.bladepanels[blade*self.n_span:blade*self.n_span+self.n_span,:]=\
                np.column_stack((p1Rot,p2Rot,p3Rot,p4Rot))

            data_arr=np.vstack((data_arr,np.column_stack((p1Rot,p2Rot,p3Rot,p4Rot))))

        return
    
    
    def comp_ring(self):
        self.discretize_spanwise()
        
        # loop over different blades
        for blade in range(self.n_blades):
            r = self.span_arr
            bladeRot = 2 * np.pi / self.n_blades*blade
            
            # loop             
            for idx,span in enumerate(self.span_arr[:-1]):
                data_arr = np.empty((0,7))
                chord1 = 3*(1-span)+1
                angle1 = np.deg2rad(- 14 * (1 - span) + 2)

                # Define bound filament
                # [x1, y1, z1, x2 ,y2, z2, gamma]
                data_arr = np.vstack((data_arr,[0, r[idx], 0, 0, r[idx+1], 0, 0]))
                
                # Define trailing filaments at x1
                data = [chord1*np.sin(-angle1), span, -chord1*np.cos(angle1), 0, span, 0, 0]
                data_arr = np.vstack((data_arr,data))
                
                # Loop over trailing filaments at x1
                for index, theta in enumerate(self.theta_arr[:-1]):
                    if index>0:
                        xt=data_arr[-1,0]
                        yt=data_arr[-1,1]
                        zt=data_arr[-1,2]
                        dx=(self.theta_arr[index+1]-theta)/self.tsr
                        dy=(np.cos(-self.theta_arr[index+1])-np.cos(-theta))*span
                        dz=(np.sin(-self.theta_arr[index+1])-np.sin(-theta))*span
                        # data_arr = np.vstack((data_arr,
                        #                       [xt+dx, yt+dy, zt+dz, xt, yt, zt, 0]))
                        data_arr = np.vstack((data_arr, 
                            [xt+dx, (yt+dy)*np.cos(bladeRot)-(zt+dz)*np.sin(bladeRot),
                            (yt+dy)*np.sin(bladeRot)+(zt+dz)*np.cos(bladeRot),
                            xt, yt*np.cos(bladeRot)-zt*np.sin(bladeRot),
                        yt*np.sin(bladeRot)+zt*np.cos(bladeRot),0]))
                    else:
                        xt=chord1*np.sin(-angle1)
                        yt=span
                        zt=-chord1*np.cos(angle1)
                        dx=(self.theta_arr[index+1]-theta)/self.tsr
                        dy=(np.cos(-self.theta_arr[index+1])-np.cos(-theta))*span
                        dz=(np.sin(-self.theta_arr[index+1])-np.sin(-theta))*span
                        # data_arr = np.vstack((data_arr,
                        #                       [xt+dx, yt+dy, zt+dz, xt, yt, zt, 0]))
                        data_arr = np.vstack((data_arr, 
                            [xt+dx, (yt+dy)*np.cos(bladeRot)-(zt+dz)*np.sin(bladeRot),
                            (yt+dy)*np.sin(bladeRot)+(zt+dz)*np.cos(bladeRot),
                            xt, yt*np.cos(bladeRot)-zt*np.sin(bladeRot),
                        yt*np.sin(bladeRot)+zt*np.cos(bladeRot),0]))
                    
                    
                # trailing filaments at x2
                chord2=3*(1-r[idx+1])+1
                angle2=np.deg2rad(- 14 * (1 - r[idx+1]) + 2)
                data=[0, r[idx+1], 0, chord2*np.sin(-angle2), r[idx+1], -chord2*np.cos(angle2), 0]
                data_arr = np.vstack((data_arr, data))
                
                # Loop over trailing filaments at x2
                for index, theta in enumerate(self.theta_arr[:-1]):
                    if index > 0:
                        xt=data_arr[-1,3]
                        yt=data_arr[-1,4]
                        zt=data_arr[-1,5]
                        dx=(self.theta_arr[index+1]-theta)/self.tsr
                        dy=(np.cos(-self.theta_arr[index+1])-np.cos(-theta))*r[idx+1]
                        dz=(np.sin(-self.theta_arr[index+1])-np.sin(-theta))*r[idx+1]
                        # data_arr = np.vstack((data_arr,
                                              # [xt, yt, zt, xt+dx, yt+dy, zt+dz, 0]))
                        data_arr = np.vstack((data_arr, 
                            [xt, yt*np.cos(bladeRot)-zt*np.sin(bladeRot),
                            yt*np.sin(bladeRot)+zt*np.cos(bladeRot),
                            xt+dx, (yt+dy)*np.cos(bladeRot)-(zt+dz)*np.sin(bladeRot),
                        (yt+dy)*np.sin(bladeRot)+(zt+dz)*np.cos(bladeRot),0]))
                    else:
                        xt=chord2*np.sin(-angle2)
                        yt=r[idx+1]
                        zt=-chord2*np.cos(angle2)
                        dx=(self.theta_arr[index+1]-theta)/self.tsr
                        dy=(np.cos(-self.theta_arr[index+1])-np.cos(-theta))*r[idx+1]
                        dz=(np.sin(-self.theta_arr[index+1])-np.sin(-theta))*r[idx+1]
                        # data_arr = np.vstack((data_arr,
                                              # [xt, yt, zt, xt+dx, yt+dy, zt+dz, 0]))
                        data_arr = np.vstack((data_arr, 
                            [xt, yt*np.cos(bladeRot)-zt*np.sin(bladeRot),
                            yt*np.sin(bladeRot)+zt*np.cos(bladeRot),
                            xt+dx, (yt+dy)*np.cos(bladeRot)-(zt+dz)*np.sin(bladeRot),
                        (yt+dy)*np.sin(bladeRot)+(zt+dz)*np.cos(bladeRot),0]))          

                self.test[:, blade*(self.n_span-1)+idx, :] = data_arr.T

        return
    
    def _compute_ring(self):
        self.discretize_spanwise()
        for blade in range(self.n_blades):
            return
        # TODO: compute bound vortex filaments
        # sets initial bound vortex filament
       
        
        # compute trailing filaments
        
       
            # TODO: compute trailing vortex filaments
            
            #starting point of x, y and z
               
            # compute increments dy, dz and dx
            

            # TODO: rotate all filaments
    
            # TODO: redefine in self.rings
    
        return

    def _compute_cp(self):
        # TODO: compute coordinates
        cp = np.zeros((self.n_blades, self.n_span-1, 9))
        self.discretize_spanwise()
        segmentCenter=(self.span_arr+(np.roll(self.span_arr,-1)-self.span_arr)/2)[:-1]
        for blade in range(self.n_blades):
            bladeRot = 2 * np.pi / self.n_blades*blade
            angle = np.deg2rad(2 - 14 * (1 - segmentCenter))
          
            #line 135-139 in carlos code
            # bound edge
            boundEdge=np.column_stack((np.zeros(self.n_span-1),segmentCenter*np.cos(bladeRot),segmentCenter*np.sin(bladeRot)))
            # TODO: compute tangential unit vector
            tangVect=np.column_stack((np.cos(angle),np.sin(angle)*np.sin(bladeRot),-np.sin(angle)*np.cos(bladeRot)))
            # TODO: compute normal unit vector 
            normVect=np.column_stack((-np.sin(angle),np.cos(angle)*np.sin(bladeRot),-np.cos(angle)*np.cos(bladeRot)))
            # Assign to cp
            cp=np.column_stack((boundEdge,tangVect,normVect))
            self.cp[blade*(self.n_span-1):blade*(self.n_span-1)+self.n_span-1,:]=cp
        return
        
solver = BladeGeometry(1, 6, 10, 3, 5, 5)
# solver._compute_cp()
# solver._compute_ring()
# solver.discretize_blade()
# solver.comp_ring()
# solver._compute_cp() # works
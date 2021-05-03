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
        self.theta_arr = np.linspace(0, 2 * np.pi, n_theta)[:-1]
        self.cp = np.zeros((n_blades*(n_span-1),9)) # coord; normal; tangential
        self.rings = {"x1": np.zeros((n_blades * n_span, n_theta-1)),
                      "x2": np.zeros((n_blades * n_span, n_theta-1)),
                      "y1": np.zeros((n_blades * n_span, n_theta-1)),
                      "y2": np.zeros((n_blades * n_span, n_theta-1)),
                      "z1": np.zeros((n_blades * n_span, n_theta-1)),
                      "z2": np.zeros((n_blades * n_span, n_theta-1)),
                      "gamma": np.zeros((n_blades, n_span, n_theta-1))
                      }
        self.bladepanels = np.zeros((n_blades, n_span, 4 * 3))  # empty dict to
        

        self.rRootRatio = 0.2
        self.rTipRatio = 1.0
        self.spacing = 'cosine'
        
        self.bladepanels = {}  # empty dict to store blade geometry
        
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
            print ("zet goeie spacing")
            
        self.spacing_arr=rsegment        
         
        return 

    def discretize_blade(self):
        # TODO: solve control points
        segmentCenter=(self.span_arr+(np.roll(self.span_arr,-1)-self.span_arr)/2)[:-1]
        self._compute_cp()
        for blade in range(self.n_blades):
            bladeRot = 2 * np.pi / self.n_blades*blade
            angle = np.deg2rad(3 * (1 - segmentCenter) + 1 - 14 * (1 - segmentCenter))
        # print(rsegment)
        # plt.scatter(self.bladesegments,np.ones(self.n_span))
        # plt.scatter(segmentCenter,np.ones(self.n_span-1))
        
        # TODO: solve rings
        

        # TODO: solve blade panels

        return

    def _compute_ring(self):
        
        # TODO: compute bound vortex filaments
        # sets initial bound vortex filament
        self.rings['x1'][0][0]=0
        self.rings['x2'][0][0]=0
        self.rings['y1'][0][0]=self.spacing_arr[0]
        self.rings['y2'][0][0]=self.spacing_arr[1]
        self.rings['z1'][0][0]=0
        self.rings['z2'][0][0]=0
        self.rings['gamma'][0][0]=0
        
        # compute trailing filaments
        
        for blade in range(self.n_blades):
            # TODO: compute trailing vortex filaments
            angle = np.deg2rad(3 * (1 - self.span_arr) + 1 - 14 * (1 - self.span_arr))
            chord = 3*(1-self.span_arr)+1
            # [x1, x2, y1, y2, z1, z2]
            coordList = [3*(1-self.span_arr)*np.sin(-angle),np.zeros(self.n_span),self.span_arr,self.span_arr,-3*(1-self.span_arr)*np.cos(angle),np.zeros(self.n_span)]
            
            #starting point of x, y and z
            xt=[]
            yt=[]
            zt=[]
            print(self.span_arr)
            print(self.theta_arr)
            print(np.cos(-np.roll(self.theta_arr,-1))-np.cos(-self.theta_arr))
            print(np.transpose(np.tile(self.span_arr,(self.n_theta,1))))
            
            # compute increments dy, dz and dx
            
            # np.transpose(np.tile(   ---   ),(self.n_theta,1))*self.span_arr
            dy=(np.cos(-np.roll(self.theta_arr,-1))-np.cos(-self.theta_arr))*np.transpose(np.tile(self.span_arr,(self.n_theta-1,1)))
            print(dy)
            # dy=np.cos(-np.roll(self.theta_arr,-1))-np.cos(-self.theta_arr)*self.span_arr
            print(dy)
            dz=(np.sin(-np.roll(self.theta_arr,-1))-np.sin(-self.theta_arr))*np.transpose(np.tile(self.span_arr,(self.n_theta-1,1)))
            dx=np.tile((np.roll(self.theta_arr,-1)-self.theta_arr)/self.tsr*self.radius,(self.n_span,1))
            print(coordList[0]+np.cumsum(dx,axis=0),'ok')
            self.rings['x1'][blade*(self.n_span-1)+1:blade*(self.n_span-1)+self.n_span+1,:]=coordList[0]+np.cumsum(dx,axis=0)
            
            # dz=np.sin(-np.roll(self.theta_arr,-1))-np.sin(-self.theta_arr)*self.span_arr
            # dx=(np.roll(self.theta_arr,-1)-self.theta_arr)/self.tsr*self.radius
            
            # use cumsum to compute total change
            # self.rings['x1'] = 
            
            # TODO: rotate all filaments
    
            # TODO: redefine in self.rings
    
            print (self.rings['x1'])
        return

    def _compute_cp(self):
        # TODO: compute coordinates
        cp = np.zeros((self.n_blades, self.n_span-1, 9))
        self.discretize_spanwise()
        segmentCenter=(self.spacing_arr+(np.roll(self.spacing_arr,-1)-self.spacing_arr)/2)[:-1]
        for blade in range(self.n_blades):
            bladeRot = 2 * np.pi / self.n_blades*blade
            angle = np.deg2rad(3 * (1 - segmentCenter) + 1 - 14 * (1 - segmentCenter))
          
            #line 135-139 in carlos code
            # bound edge
            boundEdge=np.column_stack((np.zeros(self.n_span-1),segmentCenter*np.cos(bladeRot),segmentCenter*np.sin(bladeRot)))
            tangVect=np.column_stack((np.cos(angle),np.sin(angle)*np.sin(bladeRot),-np.sin(angle)*np.cos(bladeRot)))
            normVect=np.column_stack((-np.sin(angle),np.cos(angle)*np.sin(bladeRot),-np.cos(angle)*np.cos(bladeRot)))
            cp=np.column_stack((boundEdge,tangVect,normVect))
            # cp[blade][0][1,2,3] = np.column_stack((np.zeros(self.n_span-1),segmentCenter*np.cos(bladeRot),segmentCenter*np.sin(bladeRot)))
            
            # TODO: compute tangential unit vector
            # cp[blade][1] = np.column_stack((np.cos(angle),np.sin(angle)*np.sin(bladeRot),-np.sin(angle)*np.cos(bladeRot)))
        
            # TODO: compute normal unit vector 
            # cp[blade][2] = np.column_stack((-np.sin(angle),np.cos(angle)*np.sin(bladeRot),-np.cos(angle)*np.cos(bladeRot)))
            # print(blade*self.n_span, blade*self.n_span+self.n_span-1)
            self.cp[blade*(self.n_span-1):blade*(self.n_span-1)+self.n_span-1,:]=cp
        # print(cp)
        # self.cp=np.reshape(cp,(self.n_blades*(self.n_span-1),9))
        return
        
solver = BladeGeometry(1, 6, 10, 3, 5, 6)
solver._compute_cp()
solver._compute_ring()

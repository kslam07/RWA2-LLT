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
        data_arr=np.empty((0,np.zeros(self.n_blades * (self.n_span-1))))
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
            angle1 = np.deg2rad(- 14 * (1 - self.span_arr) + 2)
            angle2 = np.deg2rad(- 14 * (1-np.roll(self.span_arr,-1)) + 2)
            chord1 = 3 * (1 - self.span_arr) + 1
            chord2 = 3 * (1 - np.roll(self.span_arr,-1)) + 1

            # define the 4 corners
            p1=[-0.25*chord1*np.sin(-angle1),self.span_arr,0.25*chord1*np.cos(angle1)]
            p2=[-0.25*chord1*np.sin(-angle1),np.roll(self.span_arr,-1),0.25*chord2*np.cos(angle2)]
            p3=[0.75*chord2*np.sin(-angle2),np.roll(self.span_arr,-1),-0.75*chord2*np.cos(angle2)]
            p4=[0.75*chord1*np.sin(-angle1),self.span_arr,-0.75*chord1*np.cos(angle1)]
            
            p1Rot=np.column_stack([np.zeros(self.n_span),p1[1]*np.cos(bladeRot)-p1[2]*np.sin(bladeRot),p1[1]*np.sin(bladeRot)+p1[2]*np.cos(bladeRot)])
            p2Rot=np.column_stack([np.zeros(self.n_span),p2[1]*np.cos(bladeRot)-p2[2]*np.sin(bladeRot),p2[1]*np.sin(bladeRot)+p2[2]*np.cos(bladeRot)])
            p3Rot=np.column_stack([np.zeros(self.n_span),p3[1]*np.cos(bladeRot)-p3[2]*np.sin(bladeRot),p3[1]*np.sin(bladeRot)+p3[2]*np.cos(bladeRot)])
            p4Rot=np.column_stack([np.zeros(self.n_span),p4[1]*np.cos(bladeRot)-p4[2]*np.sin(bladeRot),p4[1]*np.sin(bladeRot)+p4[2]*np.cos(bladeRot)])
            # print(self.bladepanels[blade*self.n_span:blade*self.n_span+self.n_span,:])
            # print(p1,p2,p3,p4)
            self.bladepanels[blade*self.n_span:blade*self.n_span+self.n_span,:]=\
                np.column_stack((p1Rot,p2Rot,p3Rot,p4Rot))
            # print(p1Rot,p2Rot)
            data_arr=np.vstack((data_arr,np.column_stack((p1Rot,p2Rot,p3Rot,p4Rot))))
            # print(p1Rot)
        self.rings=data_arr
        return
    
    
    def comp_ring(self):
        self.discretize_spanwise()
        data_arr = np.empty((0,7))
        # loop over different blades
        for blade in range(self.n_blades):
            r = self.span_arr
            
            
            for idx,span in enumerate(self.span_arr[:-1]):
                chord = 3*(1-span)+1
                angle = np.deg2rad(- 14 * (1 - span) + 2)

                # Define bound filament
                # [x1, y1, z1, x2 ,y2, z2, gamma]
                data_arr = np.vstack((data_arr,[0, r[idx], 0, 0, r[idx+1], 0, 0]))
                
                # Define trailing filaments at x1
                data = [chord*np.sin(-angle), span, -chord*np.cos(angle), 0, span, 0, 0]
                # print(data_arr,'1')
                data_arr = np.vstack((data_arr,data))
                # print(data_arr,'2\n')
                for index, theta in enumerate(self.theta_arr[:-1]):
                    xt=chord*np.sin(-angle)
                    yt=span
                    zt=-chord*np.cos(angle)
                    dx=(self.theta_arr[index+1]-theta)/self.tsr
                    dy=(np.cos(-self.theta_arr[index+1])-np.cos(-theta))*span
                    dz=(np.sin(-self.theta_arr[index+1])-np.sin(-theta))*span
                    data_arr = np.vstack((data_arr,
                                          [xt+dx, yt+dy, zt+dz, xt, yt, zt, 0]))
                    
                    
                # trailing filaments at x2
                data =[0, r[idx+1], 0, chord]
                # for theta in 
                # if blade==0:
                #     print(z1)
        
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
            angle = np.deg2rad(2 - 14 * (1 - self.span_arr))
            chord = 3*(1-self.span_arr)+1
            # [x1, x2, y1, y2, z1, z2]
            coordList = [3*(1-self.span_arr)*np.sin(-angle),np.zeros(self.n_span),self.span_arr,self.span_arr,-3*(1-self.span_arr)*np.cos(angle),np.zeros(self.n_span)]
            
            #starting point of x, y and z
            xt=[]
            yt=[]
            zt=[]
            # print(self.span_arr)
            # print(self.theta_arr)
            # print(np.cos(-np.roll(self.theta_arr,-1))-np.cos(-self.theta_arr))
            # print(np.transpose(np.tile(self.span_arr,(self.n_theta,1))))
            
            # compute increments dy, dz and dx
            
            # np.transpose(np.tile(   ---   ),(self.n_theta,1))*self.span_arr
            
            
            
            # check dyyyy
            dy=(np.cos(-np.roll(self.theta_arr,-1))-np.cos(-self.theta_arr))*np.transpose(np.tile(self.span_arr,(self.n_theta-1,1)))
            # print(dy)
            # dy=np.cos(-np.roll(self.theta_arr,-1))-np.cos(-self.theta_arr)*self.span_arr
            # print(dy)
            dz=(np.sin(-np.roll(self.theta_arr,-1))-np.sin(-self.theta_arr))*np.transpose(np.tile(self.span_arr,(self.n_theta-1,1)))
            dx=np.tile((np.roll(self.theta_arr,-1)-self.theta_arr)/self.tsr*self.radius,(self.n_span,1))
            # print(coordList[0]+np.cumsum(dx,axis=0),'ok')
            self.rings['x1'][blade*(self.n_span-1)+1:blade*(self.n_span-1)+self.n_span+1,:]=coordList[0]+np.cumsum(dx,axis=0)
            
            dz=np.sin(-np.roll(self.theta_arr,-1))-np.sin(-self.theta_arr)*self.span_arr
            dx=(np.roll(self.theta_arr,-1)-self.theta_arr)/self.tsr*self.radius
            
            # use cumsum to compute total change
            # self.rings['x1'] = 
            
            # TODO: rotate all filaments
    
            # TODO: redefine in self.rings
    
            # print (self.rings['x1'])
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
        

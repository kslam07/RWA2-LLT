"""
Plotting function
"""
from create_geometry import BladeGeometry
solver=BladeGeometry(50,8,10,3,25,25)
solver._compute_cp()
cp=solver.cp
solver.compute_ring()
solver.discretize_blade()
blade=solver.bladepanels
rings=solver.filaments

from matplotlib import pyplot
from mpl_toolkits.mplot3d import Axes3D
import random


fig = pyplot.figure()
ax = Axes3D(fig)


ax.scatter(rings[0,:24,1:26],rings[1,:24,1:26],rings[2,:24,1:26])
ax.scatter(rings[3,:24,1:26],rings[4,:24,1:26],rings[5,:24,1:26])
# ax.scatter(rings[0,24:48],rings[1,:24],rings[2,:24])
# ax.scatter(rings[0,24:48],rings[1,:24],rings[2,:24])
ax.set_xlabel('x')
ax.set_ylabel('y')
# ax.set_zlaebl('z')
pyplot.show()
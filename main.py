from create_geometry import BladeGeometry
from lifting_line_solver import LiftingLineSolver

prop_geo = BladeGeometry(radius=50.0, tsr=8, v_inf=10.0, n_blades=3, n_span=50, n_theta=50,
                spacing="equal",a=0,xshift=0,yshift=100,zshift=0)
solver = LiftingLineSolver(geo=prop_geo, r_rotor=50, tol=1e-3, n_iter=100,double_rotor=False)

solver.run_solver()
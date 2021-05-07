from create_geometry import BladeGeometry, doubleRotor
from lifting_line_solver import LiftingLineSolver

prop_geo = BladeGeometry(radius=1.0, tsr=8, v_inf=10.0, n_blades=3, n_span=5, n_theta=5)
solver = LiftingLineSolver(geo=prop_geo, u_rot=10/50, r_rotor=50, tol=1e-3, n_iter=100)

solver.run_solver()
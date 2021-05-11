from create_geometry import BladeGeometry, doubleRotor
from lifting_line_solver import LiftingLineSolver

prop_geo = BladeGeometry(radius=50.0, tsr=8, v_inf=10.0, n_blades=3, n_span=15, n_theta=15, spacing="equal")
solver = LiftingLineSolver(geo=prop_geo, r_rotor=50, tol=1e-3, n_iter=100)

solver.run_solver()
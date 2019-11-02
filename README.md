# vortex_panel_method
Calculating Vortex Strengths of Each Panel for the Vortex Panel Method with the Example of a Circular Cylinder

# Language and Environment
Written in Fortran90 compatible with the Linux GNU gfortran compiler. [LAPACK](https://github.com/Reference-LAPACK/lapack) is required for matrix inversion

![Sample Geometry](https://www.overleaf.com/project/5dba00bbb541450001db30df/file/5dba422ab541450001db5714)

# Input (input.in)
  - NPTS: # of total node points
  - FRANG: Angle of attack of the free stream flow
  - ISSYM: Symmetric flow (if 1) or asymmetric flow where the stagnation point occurs at θ = -150˚ (if 0)
  - X, Y: Nodes for vortex panels. Starts from the stangation point, indexed clockwise

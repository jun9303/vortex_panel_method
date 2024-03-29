# Vortex Panel Method
Calculating Vortex Strengths of Each Panel for the Vortex Panel Method with the Example of a Circular Cylinder

# Language and Environment
Written in Fortran90 compatible with the Linux GNU gfortran compiler. [LAPACK](https://github.com/Reference-LAPACK/lapack) is required for matrix inversion

# Sample Geometry
![Sample Geometry](https://drive.google.com/uc?export=view&id=1s_A3oSA1bRufd-UQqJIdCOeXWa3vUilE)

# Input (input.in)
  - NPTS: # of total node points
  - FRANG: Angle of attack of the free stream flow
  - ISSYM: Symmetric flow (if 1) or asymmetric flow where the stagnation point occurs at θ = -150˚ (if 0)
  - X, Y: Nodes for vortex panels. Starts from the stangation point, indexed clockwise

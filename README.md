# vortex_panel_method
Calculating Vortex Strengths of Each Panel for the Vortex Panel Method with the Example of a Circular Cylinder

# Language and Environment
Written in Fortran90 compatible with the Linux GNU gfortran compiler. [LAPACK](https://github.com/Reference-LAPACK/lapack) is required for matrix inversion

![Sample Geometry](https://doc-04-a8-docs.googleusercontent.com/docs/securesc/ffjn0859b4rebh0hh7sgijetk8hhlida/cmd6iomhnjokvmq8k6eaaitsepv1hp0m/1572681600000/00034033724614023160/08046271387926787997/1s_A3oSA1bRufd-UQqJIdCOeXWa3vUilE?e=view&nonce=jfmg01sh5ju9u&user=08046271387926787997&hash=8v5tk0fhovh4hhk202pruiea6sua6ad2)

# Input (input.in)
  - NPTS: # of total node points
  - FRANG: Angle of attack of the free stream flow
  - ISSYM: Symmetric flow (if 1) or asymmetric flow where the stagnation point occurs at θ = -150˚ (if 0)
  - X, Y: Nodes for vortex panels. Starts from the stangation point, indexed clockwise

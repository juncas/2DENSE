# 2DENSE
A 2 Dimensional Euler/Naviers-Stokes Equations solver. 

**2DENSE is still under developing and will be updated instantaneously.**

This is the original code for our paper ["An Adaptive Characteristic-wise Reconstruction WENO scheme for Gas Dynamic Euler
Equations"](https://arxiv.org/abs/1711.11288). 

## Time integral method
Third order TVD Runge-Kutta method.
## Riemann solvers
- Local Lax-Friedrichs splitting
- Global Lax-Friedrichs splitting
- Steger-Warming
- Characteristic-wise reconstruction with Roe-solver and global Lax-Friedrichs splitting
## Reconstruction methods
- 5th order upwind scheme
- 5th order WENO-JS scheme
- 5th order WENO-Z scheme
- 5th order AdaWENO scheme
## Pre-defined test problem setups
- Isentropic vortex convection problem
- Sedov problem
- Rayleigh-Taylor instability problem
- Richtmyer–Meshkov instability problem
- Double mach reflection problem
- Shock/Shear layer interaction problem
- Shock/Vortex Interaction problem

 

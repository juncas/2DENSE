import math
import scipy
import numpy as np
gamma = 1.4
Mach_s = 1.1

rho_u = 1.0
u_u = Mach_s*np.sqrt(gamma)
p_u = 1.0

rho = ((gamma+1)*Mach_s**2/(2.0+(gamma-1)*Mach_s**2))*rho_u
p = 1+2*gamma/(gamma+1)*(Mach_s**2-1)
u = (2.0+(gamma-1)*Mach_s**2)/((gamma+1)*Mach_s**2)*u_u

print("rho=%f\nu=%f\np=%f"%(rho,u,p))

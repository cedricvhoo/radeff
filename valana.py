# -*- coding: utf-8 -*-
"""
Created on Fri Jul  7 10:21:19 2023

@author: Cédric Van hoorickx
"""

# Import libraries
import math
import numpy as np
import matplotlib.pyplot as plt
from radfun import radeffYu,radeffpoly

## GIVEN ##
# Frequency
f = np.logspace(math.log10(10),math.log10(1e4),num=250)     # Frequency array
omega = 2*math.pi*f                         # Angular frequency

# Plot details
ymin = 1e-4                                 # Minimum radiation efficiency plotted
ymax = 1e1                                  # Maximum radiation efficiency plotted

# Air properties
c0 = 343                                    # Sound speed in air
k0 = omega/c0                               # WAvenumber air

# Plate properties
h = 6e-3                                    # Thickness
cL = 5200                                   # Longitudinal wave speed
rho = 2500                                  # Mass density
nu = 0.24                                   # Poisson ratio

# Structural properties
D = rho*cL**2*h**3/12                       # Bending stiffnes
m = rho*h                                   # Surface density
kf = omega**(1/2)*(m/D)**(1/4)              # Bending wavenumber

## EXAMPLE 1 ##
# Geometry
Matm = np.array([[1,1,1]])                  # Geometry matrix
Lxim = np.ones((1,3))                       # Vector with lengths of subplates in x direction
Lyim = np.ones((1,1))                       # Vector with lengths of subplates in y direction

# Radiation efficiency
s = radeffYu(k0,kf,sum(Lxim.T),sum(Lyim))
s1,s2 = radeffpoly(k0,kf,Matm,Lxim,Lyim)

# Plot
plt.figure(1)
plt.loglog(f,s,'k')                         # Plot the radiation efficiency of the total plate
plt.loglog(f,s1,'b')                        # Plot the first order approximation of the radiation efficiency
plt.loglog(f,s2,'r')                        # Plot the second order approximation of the radiation efficiency
plt.xlabel("Frequency [Hz]")
plt.ylabel("Radiation efficiency [-]")
plt.xlim([f[0],f[-1]])
plt.ylim([ymin,ymax])

## EXAMPLE 2 ##
# Geometry
Lx = 2                                      # Length of the plate in the x direction
Ly = 2                                      # Length of the plate in the y direction

Lxi = 1                                     # Length of the subplate in the x direction
Lyi = 1                                     # Lenght of the subplate in the y direction

# Radiation efficiency
s = radeffYu(k0,kf,Lx,Ly)                   # Compute the radiation efficiency of the total plate
si = radeffYu(k0,kf,Lxi,Lyi)                # Compute the radiation efficiency of the subplate
si2x = radeffYu(k0,kf,2*Lxi,Lyi)            # Compute the radiation efficiency of two adjacent subplates in the x direction
si2y = radeffYu(k0,kf,2*Lyi,Lxi)            # Compute the radiation efficiency of two adjacent subplates in the y direction
s1 = 4*si*1/4                               # Compute the first order approximation of the radiation efficiency
s2 = 2*si2x*2/4+2*si2y*2/4-4*si*1/4         # Compute the second order approximation of the radiation efficiency

# Geometry
Matm = np.array([[1,1],[1,1]])              # Geometry matrix
Lxim = np.ones((1,2))                       # Vector with lengths of subplates in x direction
Lyim = np.ones((2,1))                       # Vector with lengths of subplates in y direction

# Radiation efficiency
s = radeffYu(k0,kf,sum(Lxim.T),sum(Lyim))
s1,s2 = radeffpoly(k0,kf,Matm,Lxim,Lyim)

# Plot
plt.figure(2)
plt.loglog(f,s,'k')                         # Plot the radiation efficiency of the total plate
plt.loglog(f,s1,'b')                        # Plot the first order approximation of the radiation efficiency
plt.loglog(f,s2,'r')                        # Plot the second order approximation of the radiation efficiency
plt.xlabel("Frequency [Hz]")
plt.ylabel("Radiation efficiency [-]")
plt.xlim([f[0],f[-1]])
plt.ylim([ymin,ymax])

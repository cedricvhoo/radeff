# -*- coding: utf-8 -*-
"""
Created on Fri Jul  7 10:21:19 2023

@author: Cédric Van hoorickx
"""

# Import libraries
import math
import numpy as np
import matplotlib.pyplot as plt
from radfun import radeffpoly,progressbar
from scipy.special import j0


## GIVEN ##
# Frequency
f = np.logspace(math.log10(10),math.log10(1e4),num=250)     # Frequency array
omega = 2*math.pi*f                         # Angular frequency

# Plot details
ymin = 1e-4                                 # Minimum radiation efficiency plotted
ymax = 1e1                                  # Maximum radiation efficiency plotted

# Air properties
c0 = 343                                    # Sound speed in air
k0 = omega/c0                               # Wavenumber air

# Plate properties
h = 6e-3                                    # Thickness
cL = 5200                                   # Longitudinal wave speed
rho = 2500                                  # Mass density
nu = 0.24                                   # Poisson ratio

# Structural properties
D = rho*cL**2*h**3/12                       # Bending stiffnes
m = rho*h                                   # Surface density
kf = omega**(1/2)*(m/D)**(1/4)              # Bending wavenumber

# Geometry
Matm = np.array([[1,1],[1,0]])              # Geometry matrix 
Lxim = np.ones((1,2))                       # Vector with lengths of subplates in x direction
Lyim = np.ones((2,1))                       # Vector with lengths of subplates in y direction

# Mesh resolution ! Make this fine enough for computations at high frequencies!
dx = 0.05                                   # Mesh size in the x direction
dy = dx                                     # Mesh size in the y direction
fnmax = c0/(2*max(dx,dy))                   # Maximum frequency to be computed with mesh size
fn = f[f<=fnmax]                            # Frequency array for numerical computation

## NUMERICAL APPROACH ##
# Geometry
Lx = sum(Lxim.T)                            # (Maximum) length of the plate in the x direction
Ly = sum(Lyim)                              # (Maximum) length of the plate in the y direction

Lxi = Lxim[0,0]                             # Length of the subplate in the x direction
Lyi = Lyim[0,0]                             # Lenght of the subplate in the y direction

S = Lx*Ly-Lxi*Lyi                           # Total surface area

# Make X1 and X2 which contain all coordinates of the mesh
xA = np.arange(0,Lx+dx,dx)                  # x coordinates of part A (Lx * Lyi)
yA = np.arange(0,Lyi+dy,dy)                 # y coordinates of part A (Lx * Lyi)
xB = np.arange(0,Lxi+dx,dx)                 # x coordinates of part B (Lxi * Lyi)
yB = np.arange(Lyi+dy,Ly+dy,dy)             # y coordinates of part B (Lxi * Lyi)

XA = np.transpose(np.tile(xA,(1,len(yA))))                      # Repeat xA for all elements of yA
YA = np.reshape(np.tile(yA[:,np.newaxis],(1,len(xA))),(-1,1))   # Repeat yA for all elements of xA
XB = np.transpose(np.tile(xB,(1,len(yB))))                      # Repeat xB for all elements of yB
YB = np.reshape(np.tile(yB[:,np.newaxis],(1,len(xB))),(-1,1))   # Repeat yB for all elements of xB

X1 = np.vstack((XA,XB))                     # Combine x coordinates of parts A and B
Y1 = np.vstack((YA,YB))                     # Combine y coordinates of parts A and B

R = np.sqrt((X1-X1.T)**2+(Y1-Y1.T)**2)      # Compute the distance between the different mesh points
R[R==0] = 1e-8                              # Avoid numerical issues by replacing 0 by 1e-8

sn = np.zeros_like(fn)                      # Initialize the numerically computed radiation efficiency
for parf in range(len(fn)):
    progressbar((parf+1)/len(fn))           # Display a progress bar
    G = np.exp(-1j*k0[parf]*R)/(2*np.pi*R)  # Define the Green's function
    sn[parf] = np.real(1j*k0[parf]/S*np.sum(np.sum(G*j0(kf[parf]*R),axis=1),axis=0)*dx**2*dy**2)       # Perform numerical integration

## ANALYTICAL APPROACH ##
s1,s2 = radeffpoly(k0,kf,Matm,Lxim,Lyim)

## PLOT
plt.figure(3)
plt.loglog(fn,sn,'k')                       # Plot the radiation efficiency of the total plate [numerical solution]
plt.loglog(f,s1,'b')                        # Plot the first order approximation of the radiation efficiency [analytical solution]
plt.loglog(f,s2,'r')                        # Plot the second order approximation of the radiation efficiency [analytical solution]
plt.xlabel("Frequency [Hz]")
plt.ylabel("Radiation efficiency [-]")
plt.xlim([f[0],f[-1]])
plt.ylim([ymin,ymax])

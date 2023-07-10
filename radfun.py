# -*- coding: utf-8 -*-
"""
Created on Fri Jul  7 15:30:02 2023

@author: Cédric Van hoorickx
"""

# Import libraries
import math
import numpy as np
import scipy.integrate as integrate
import scipy.special as special
import sys

def radeffYu(k0,kf,Lxi,Lyi):
    """Computes the radiation efficiency of a rectangular plate according to:
        [1] Y. Yu and C. Hopkins. Reduced order integration for the radiation efficiency 
        of a rectangular plate. JASA Express Letters, 1(6):062801, 2021.
        
    Parameters
    ----------
    k0 : Array (Nf x 1)
        Wavenumber air [^{-1}]
    kf : Array (Nf x 1)
        Bending wavenumber of the plate [m^{-1}]
    Lxi : scalar
        Length of the plate in the x direction [m]
    Lyi : scalar
        Length of the plate in the y direction [m]

    Returns
    -------
    sigmai : Array (Nf x 1)
        Radiation efficiency of the plate [-]

    """
    sigmai = np.zeros(len(k0))        # initialize vector with radiation efficiencies
    for parom in range(len(k0)):
        A = -2*k0[parom]/(math.pi*Lxi*Lyi)         # constant, to be multiplied by the following integrals
        (I1,err) = integrate.quad(lambda x: np.imag((Lxi*Lyi*math.pi/2-Lyi*x-Lxi*x+x**2/2)*special.j0(kf[parom]*x)*np.exp(-1j*k0[parom]*x)),0,Lyi,limit=5000)
        (I2,err) = integrate.quad(lambda x: np.imag((Lxi*Lyi*math.asin(Lyi/x)-Lyi**2/2+Lxi*(x**2-Lyi**2)**(1/2)-Lxi*x)*special.j0(kf[parom]*x)*np.exp(-1j*k0[parom]*x)),Lyi,Lxi,limit=5000)
        (I3,err) = integrate.quad(lambda x: np.imag((Lxi*Lyi*(math.asin(Lyi/x)-math.acos(Lxi/x))+Lyi*(x**2-Lxi**2)**(1/2)+Lxi*math.sqrt(x**2-Lyi**2)-(Lxi**2+Lyi**2+x**2)/2)*special.j0(kf[parom]*x)*np.exp(-1j*k0[parom]*x)),Lxi,math.sqrt(Lxi**2+Lyi**2),limit=5000)
        sigmai[parom] = A*(I1+I2+I3)        # radiation efficiency
    return sigmai

def radeffpoly(k0,kf,Matm,Lxim,Lyim):
    """Computes the radiation efficiency of a polygonal plate
        
    Parameters
    ----------
    k0 : Array (Nf x 1)
        Wavenumber air [^{-1}]
    kf : Array (Nf x 1)
        Bending wavenumber of the plate [m^{-1}]
    Matm: Matrix (Ny x Nx)
        Matrix representing the geometry: equal to 1 if material is present, equal to 0 if not
    Lxim : Array (1 x Nx)
        Length of the subplates in the x direction [m]
    Lyim : Array (Ny x 1)
        Length of the subplates in the y direction [m]

    Returns
    -------
    sigma0 : Array (Nf x 1)
        Zero order approximation of the plate's radiation efficiency [-]
    sigma1 : Array (Nf x 1)
        First order approximation of the plate's radiation efficiency [-]

    """
    # Definitions
    Nex = Matm.shape[1]
    Ney = Matm.shape[0]
    Net = Nex*Ney

    Idx = np.tile(np.arange(0,Nex),(Ney,1))
    Idy = np.tile(np.arange(0,Ney).reshape(-1,1),(1,Nex))
    Idxv = np.reshape(Idx.T,-1)
    Idyv = np.reshape(Idy.T,-1)

    Matv = np.reshape(Matm,-1)

    # Surface areas
    Si = Lxim*Lyim
    S = np.sum(np.reshape(Si*Matm,-1))

    # Neighbouring elements
    nbe = np.empty((Net, 4))
    nbe[:] = np.nan

    indr1 = np.reshape(np.arange(0,Nex)*Ney+np.arange(1,Ney)[:,np.newaxis],-1)
    indr2 = np.reshape(np.arange(0,Nex)*Ney+np.arange(0,Ney-1)[:,np.newaxis],-1)
    indr3 = np.reshape(np.arange(1,Nex)*Ney+np.arange(0,Ney)[:,np.newaxis],-1)
    indr4 = np.reshape(np.arange(0,Nex-1)*Ney+np.arange(0,Ney)[:,np.newaxis],-1)

    nbe[indr1,0] = indr2
    nbe[indr2,1] = indr1
    nbe[indr3,2] = indr4
    nbe[indr4,3] = indr3

    for pare in range(Net):
        if Matv[pare] == 0:
            nbe[nbe == pare] = np.nan
            nbe[pare, :] = np.nan

    Nnb = np.sum(~np.isnan(nbe), axis=1)

    # Radiation efficiencies
    si = np.empty((Net,len(k0)))
    si[:] = np.nan
    for pare in range(Net):
        if Matv[pare] != 0:
            Lxe = max(Lxim.T[Idxv[pare]],Lyim[Idyv[pare]])
            Lye = min(Lxim.T[Idxv[pare]],Lyim[Idyv[pare]])
            si[pare,:] = radeffYu(k0,kf,Lxe,Lye)

    Sij = np.empty((Net, 4))
    Sij[:] = np.nan
    sij = np.empty((Net, 4, len(k0)))
    sij[:] = np.nan
    for pare in range(Net):
        for parn in range(4):
            if not np.isnan(nbe[pare, parn]):
                if parn <= 2:
                    Lxij = Lxim.T[Idxv[pare]]
                    Lyij = Lyim[Idyv[pare]]+Lyim[Idyv[int(nbe[pare,parn])]]
                else:
                    Lxij = Lxim.T[Idxv[pare]]+Lxim.T[Idxv[int(nbe[pare,parn])]]
                    Lyij = Lyim[Idyv[pare]]
                
                Sij[pare,parn] = Lxij*Lyij
                sij[pare,parn,:] = radeffYu(k0,kf,max(Lxij,Lyij),min(Lxij,Lyij))

    Sim = np.tile(np.reshape(Si,-1),(len(k0),1)).T
    sigma0 = np.sum(Sim/S*si,axis=0,where=~np.isnan(si))
    Sijm = np.tile(Sij.T,(len(k0),1,1)).T
    sig1a = 1/2*np.sum(np.sum(Sijm/S*sij,axis=0,where=~np.isnan(sij)),axis=0)
    Nnbm = np.tile(np.reshape(Nnb,-1),(len(k0),1)).T
    sig1b = np.sum((Nnbm-1)*Sim/S*si,axis=0,where=~np.isnan(si))
    sigma1 = sig1a-sig1b
    
    return sigma0,sigma1

def progressbar(progress):
    """Show a progress bar
    
    Parameters
    ----------
    progress : scalar between 0 and 1
        Progress ratio (0 at start, 1 at end)
    """
    barLength = 20 # Modify this to change the length of the progress bar
    status = ""
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
        status = "error: progress var must be float\r\n"
    if progress < 0:
        progress = 0
        status = "Halt...\r\n"
    if progress >= 1:
        progress = 1
        status = "Done...\r\n"
    block = int(math.floor(barLength*progress))
    text = "\rPercent: [{0}] {1:.1f}% {2}".format( "#"*block + "-"*(barLength-block), progress*100, status)
    sys.stdout.write(text)
    sys.stdout.flush()
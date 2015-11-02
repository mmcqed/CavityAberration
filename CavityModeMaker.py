from __future__ import division
import scipy
import scipy.misc
import scipy.fftpack
import cmath
import matplotlib.pyplot as plt
import scipy.linalg
from scipy.constants import pi,c
from CavityDefinitions import *



phiModeLib = (0.+0.j)*scipy.zeros((len(x),len(y),NumModes))
Full_C_phiModeLib = (0.+0.j)*scipy.zeros((len(x),len(y),NumModes))
Ideal_C_phiModeLib = (0.+0.j)*scipy.zeros((len(x),len(y),NumModes))


for i in scipy.linspace(0,NumModes-1,NumModes):
    l = 2*scipy.floor(scipy.sqrt(i))-(i-scipy.floor(scipy.sqrt(i))**2)
    m = (i-scipy.floor(scipy.sqrt(i))**2)

    phi_temp = philmAstig(l,m,x,y,w0x,w0y)
    
    phiModeLib[:,:,i] = phi_temp
    Full_C_phiModeLib[:,:,i] = iCFFT(Ghalf_Full*CFFT(rTotalAstig*iCFFT(Ghalf_Full**2*CFFT(rTotalAstig*iCFFT(Ghalf_Full*CFFT(phi_temp))))))
    Ideal_C_phiModeLib[:,:,i] = iCFFT(Ghalf_Helmholtz*CFFT(rOpAstig*iCFFT(Ghalf_Helmholtz**2*CFFT(rOpAstig*iCFFT(Ghalf_Helmholtz*CFFT(phi_temp))))))
    print i

scipy.save('CavityModeFunctions/phiModeLib.npy',phiModeLib)
scipy.save('CavityModeFunctions/Full_C_phiModeLib.npy',Full_C_phiModeLib)
scipy.save('CavityModeFunctions/Ideal_C_phiModeLib.npy',Ideal_C_phiModeLib)
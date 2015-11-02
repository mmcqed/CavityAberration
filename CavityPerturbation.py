from __future__ import division
import scipy
import scipy.misc
import scipy.fftpack
import cmath
import matplotlib.pyplot as plt
import scipy.linalg
from scipy.constants import pi,c
from CavityDefinitions import *
import sys


phiModeLib = scipy.load('/afs/ir.stanford.edu/users/p/a/papag/Documents/LevLab/CavityAberrationSimulation/CavityModeFunctions/phiModeLib.npy')
Full_C_phiModeLib = scipy.load('/afs/ir.stanford.edu/users/p/a/papag/Documents/LevLab/CavityAberrationSimulation/CavityModeFunctions/Full_C_phiModeLib.npy')
Ideal_C_phiModeLib = scipy.load('/afs/ir.stanford.edu/users/p/a/papag/Documents/LevLab/CavityAberrationSimulation/CavityModeFunctions/Ideal_C_phiModeLib.npy')



C_Lambda = (0.+0.j)*scipy.ones((NumModes))    #unperturbed eigenvalues 
Full_C_Lambda = (0.+0.j)*scipy.ones((NumModes))          #final eigenvalues

def CavTotalLambda(i):

    for j in scipy.linspace(0,i,i+1):
        
        Full_C_Lambda[j] = scipy.sum(phiModeLib[:,:,j]*Full_C_phiModeLib[:,:,i])*dx**2
        C_Lambda[j] = scipy.sum(phiModeLib[:,:,j]*Ideal_C_phiModeLib[:,:,i])*dx**2
        
        print i,j,C_Lambda[j]
    return C_Lambda,TotalC_Lambda
    

if __name__ == '__main__':
    i = int(sys.argv[1])-1
    scipy.save('/afs/ir.stanford.edu/users/p/a/papag/Documents/LevLab/CavityAberrationSimulation/CavDataDir/CLam'+str(i)+'.npy',CavTotalLambda(i))

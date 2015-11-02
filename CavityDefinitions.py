from __future__ import division
import scipy
import scipy.misc
import scipy.fftpack as ft
import cmath
import matplotlib.pyplot as plt
import scipy.linalg
from scipy.constants import pi,c
import numpy.polynomial.hermite


def CFFT(Axy):
    Apq = ft.fftshift(ft.fft2(ft.ifftshift(Axy)))
    return Apq
    
def iCFFT(Apq):
    Axy = ft.ifftshift(ft.ifft2(ft.fftshift(Apq)))
    return Axy
    
def philm(l,m,x,y,w0):
    
    normalize_factor = scipy.sqrt(2/2**(l+m)/scipy.misc.factorial(l)/scipy.misc.factorial(m)/pi)/w0
    hermite_x = numpy.polynomial.hermite.hermval(scipy.sqrt(2)*x/w0,scipy.concatenate((scipy.zeros(l),scipy.ones(1))))
    hermite_y = numpy.polynomial.hermite.hermval(scipy.sqrt(2)*y/w0,scipy.concatenate((scipy.zeros(m),scipy.ones(1))))
    phi_lm = normalize_factor*hermite_x*hermite_y*scipy.exp(-(x**2+y**2)/w0**2)
    return phi_lm
    
def philmAstig(l,m,x,y,w0x,w0y):
    normalize_factor = scipy.sqrt(2/2**(l+m)/scipy.misc.factorial(l)/scipy.misc.factorial(m)/pi)/scipy.sqrt(w0x*w0y)    
    hermite_x = numpy.polynomial.hermite.hermval(scipy.sqrt(2)*x/w0x,scipy.concatenate((scipy.zeros(l),scipy.ones(1))))
    hermite_y = numpy.polynomial.hermite.hermval(scipy.sqrt(2)*y/w0y,scipy.concatenate((scipy.zeros(m),scipy.ones(1))))
    phiAstig_lm = normalize_factor*hermite_x*hermite_y*scipy.exp(-(x**2/w0x**2+y**2/w0y**2))
    return phiAstig_lm

maxIndex = 20 # Considering all modes up to l+m = maxIndex
NumModes = (maxIndex/2+1)**2 # The number of modes that are included


# Cavity properties

ref = 0.999965 #reflectivity
ref = 0.999948

R = .01
Rx = R + 4.3e-6
Ry = Rx - 8.6e-6

L = R - 1.75e-6

#wavelength = 2*L/(25550+scipy.arccos(1-L/R)/pi)
wavelength = 2*L/(25550 + scipy.arccos(1-L/Rx)/(2*pi) + scipy.arccos(1-L/Ry)/(2*pi))
k = 2*pi/wavelength

w0 = scipy.sqrt(wavelength/pi*scipy.sqrt(L*R/2-L**2/4)) #waist size at the center radius
w0x = scipy.sqrt(wavelength/pi*scipy.sqrt(L*Rx/2-L**2/4))
w0y = scipy.sqrt(wavelength/pi*scipy.sqrt(L*Ry/2-L**2/4))



N = 512 #sample points


# Real space domain

maxLen = 5*scipy.sqrt(2)**(scipy.log2(N)-5)*w0
dx = 2.*maxLen/N # Discretization in realspace
xyspace = scipy.linspace(-maxLen,maxLen,N,endpoint=False)
x,y = scipy.meshgrid(xyspace,xyspace) #x,y coordinates



# Define reflection operators

rOp = scipy.exp(2*1j*k*(x**2+y**2)/2./R) # ideal reflection operator
rTotal = ref*scipy.exp(2*1j*k*(R-scipy.sqrt(R**2-x**2-y**2))) # actual reflection operator
delta_rOp = (rTotal-rOp) #first perturbation in rOp

rOpAstig = scipy.exp(2*1j*k*(x**2/Rx+y**2/Ry)/2.)
rTotalAstig = ref*scipy.exp(2*1j*k*(R-scipy.sqrt(R**2-x**2*(R/Rx)-y**2*(R/Ry))))


# k space domain

dk = 1/2./maxLen # discretization in kspace
PQspace = scipy.linspace(-N/4./maxLen,N/4./maxLen,N,endpoint=False)
P,Q = scipy.meshgrid(PQspace,PQspace)

# Define k space operators

Ghalf_Helmholtz = scipy.exp(1j*(L/2)*(4*pi**2*(P**2+Q**2)/2/k))  # Helmholtz k space propagator
Ghalf_Full = scipy.exp(1j*(L/2)*(k-scipy.sqrt(k**2-4*pi**2*(P**2+Q**2)))) # Full k space propagator
deltaGhalf = (Ghalf_Full-Ghalf_Helmholtz) # half length perturbation
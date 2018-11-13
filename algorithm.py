import numpy as np
import cmath
import math
#ieee document 7813017 used

def wavefront_func(x, y, R, C, A, S):
    c0 = ((R ** 2) * C * np.sin(2*A))/(4 * np.sqrt(6))
    z0 = (2 * np.sqrt(6) * x * y)
    c1 = ((R ** 2 * (S + C/2))/(4 * np.sqrt(3)))
    z1 = np.sqrt(3) * (2 * x ** 2 + 2 * y ** 2 - 1)
    c2 = ((R ** 2) * C * np.cos(2*A))/(4 * np.sqrt(6))
    z2 = np.sqrt(6) * (x ** 2 - y ** 2)
    
    W = c0*z0 + c1*z1 + c2*z2   
    return W
    

def S_calc(Sm, d):
    if 1/d < abs(Sm):
        S = Sm + 1/d
    elif 8 - 1/d < abs(Sm):
        S = Sm - (8 - 1/d)
    else:
        S = 0
    return S

def PSF(x, y, R, C, A, Sm, d):
    #still confused aqbout this pupil function malarkey
    P = 1
    #or it equals 0 outside of the projected aperture 
    wavelength = 550 #nm used for monochromatic simulation others are used for RGB need to perform seperately
    W = wavefront_func(x, y, R, C, A, S_calc(Sm, d))
    pupil_func = P * cmath.exp(cmath.sqrt(-1) * (2*cmath.pi)/wavelength) * W
    pupil_func_arr = []
    pupil_func_arr.append(pupil_func)
    psf_func = (np.fft.fft(pupil_func_arr) ** 2)

    return psf_func

print(PSF(1,2,3,4,5,6, 5))

#look at deconvolution processes probably implement a few 

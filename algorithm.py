import numpy as np
import cv2 


#ieee document 7813017 used

def wavefront_func(x, y, R, C, A, S_calc):
    c0 = ((R ** 2) * C * np.sin(2A))/(4 * np.sqrt(6))
    z0 = (2 * np.sqrt(6) * x * y)
    c1 = ((R ** 2 * (S + C/2))/(4 * np.sqrt(3)))
    z1 = np.sqrt(3) * (2 * x ** 2 + 2 * y ** 2 - 1)
    c0 = ((R ** 2) * C * np.cos(2A))/(4 * np.sqrt(6))
    z2 = np.sqrt(6) * (x ** 2 - y ** 2)
    
    W = c0*z0 + c1*z1 + c2*z2   
    return W
    

def S_calc(Sm, d):
    if 1/d < abs(Sm):
        S = Sm + 1/d
    else if 8 - 1/d < abs(Sm):
        S = Sm - (8 - 1/d)
    else:
        S = 0
    return S

def PSF(x, y, R, C, A, S):
    #still confused aqbout this pupil function malarkey
    P = 1
    #or it equals 0 outside of the projected aperture 
    wavelength = 550 #nm used for monochromatic simulation others are used for RGB need to perform seperately
    W = wavefront_func(x, y, R, C, A, S)
    psf_func = np.fft(P * exp(-1j * (2*Math.pi)/wavelength) * W) 

    return psf_func


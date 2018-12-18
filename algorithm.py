import numpy as np
import cmath
import math
import cv2
#ieee document 7813017 used

def image():
    img = cv2.imread('lena.png')
    height, width = img.shape[:2]
    img2 = np.empty(shape=(height,width))
    normalizedImg = np.empty(shape=(height,width))
    # print(img)
    for x in range(0, width):
        for y in range(0, height):
            img2[x][y] = ideal(x, y, 14, 1.25, 10, 0.5, 100, img)[0]/1 * 255
            
    cv2.imshow('image', img2)
    cv2.waitKey(0)
def ideal(x, y, R, C, A, Sm, d, img):
    weight = 2 #used to balance the two terms
    K = PSF(x, y, R, C, A, Sm, d)
    # top = np.multiply(weight * np.conj(np.fft.fft(K)), np.fft.fft(img))
    f1 = np.array([1, -1])
    f2 = (np.array([-1, 1]))
    # s1 = np.multiply(np.conj(np.fft.fft(f1)), np.fft.fft(f1))
    # s2 = np.multiply(np.conj(np.fft.fft(f2)), np.fft.fft(f2))
    # s3 = np.multiply(weight * np.conj(np.fft.fft(K)), np.fft.fft(K))
    # print(s1)
    # print(s2)
    # print(s3)
    # bottom = s1 + s2 + s3
    bottom = np.multiply(np.conj(np.fft.fft(f1)), np.fft.fft(f1)) + np.multiply(np.conj(np.fft.fft(f2)), np.fft.fft(f2)) + np.multiply(weight * np.conj(np.fft.fft(K)), np.fft.fft(K))
    
    top = np.multiply(weight, np.conj(np.fft.fft(K)))
    # bottom = s1 + s2 + s3
      
    
    Ip = np.fft.ifft(np.divide(top, bottom))
    # print(Ip)
    return Ip



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

# print(PSF(1,2,3,4,5,6, 5))
image()
#look at deconvolution processes probably implement a few 

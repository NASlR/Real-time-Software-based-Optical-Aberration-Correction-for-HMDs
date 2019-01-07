import numpy as np
import cmath
import math
import cv2
#ieee document 7813017 used

def image():
    img = cv2.imread('lena.png')
    img = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
    height, width = img.shape[:2]
    img2 = np.empty(shape=(height,width))
    img3 = np.empty(shape=(height,width))

    for x in range(0, width):
        for y in range(0, height):
            img2[x][y] = generalised_pupil_function(x, y, 1.5, -2, 45, 0.5, 0.0010, 550)
    
    K = abs(np.fft.fft(img2)) ** 2
    f1 = np.array([1, -1])
    f2 = (np.array([-1, 1]))
    weight = 2
    top = np.multiply(weight * np.conj(np.fft.fft(K)), np.fft.fft(img))
    s1 = np.multiply(np.conj(np.fft.fft(f1)), np.fft.fft(f1))
    s2 = np.multiply(np.conj(np.fft.fft(f2)), np.fft.fft(f2))
    s3 = np.multiply(weight * np.conj(np.fft.fft(K)), np.fft.fft(K))

    bottom = s3 # can't broadcast

    Ip = np.fft.ifft(np.divide(top, bottom))
    for x in range(0, width):
        for y in range(0, height):
            img3[x][y] = Ip[x][y]





    cv2.imshow('image', img3)
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



def aberrated_wavefront_function(x, y, R, C, A, S):
    # for each x,y
    # S, C are the sphere and cylinder values from prescription in diopters
    # A is cylinder axis in degrees
    # S is related to defocus abberation 
    # R is the radius of the pupil in mm
    # the coeffcients are in μm

    # oblique astigmatism where m, n = -2, 2
    oblique_coefficient = ((R ** 2) * C * np.sin(2*A))/(4 * np.sqrt(6))
    oblique_polynomial = (2 * np.sqrt(6) * x * y)

    # defocus astigmatism where m,n = 0, 2
    defocus_coefficient = -((R ** 2 * (S + C/2))/(4 * np.sqrt(3)))
    defocus_polynomial = np.sqrt(3) * (2 * x ** 2 + 2 * y ** 2 - 1)

    #vertical astigmatism where m, n = 2, 2
    vertical_coefficient = ((R ** 2) * C * np.cos(2*A))/(4 * np.sqrt(6))
    vertical_polynomial = np.sqrt(6) * (x ** 2 - y ** 2)
    
    W = oblique_coefficient*oblique_polynomial + defocus_coefficient*defocus_polynomial + vertical_coefficient*vertical_polynomial   
    return W
    

def S_calc(Sm, d):
    if 1/d < abs(Sm):
        S = Sm + 1/d
    elif (8 - 1/d) < abs(Sm):
        S = Sm - (8 - 1/d)
    else:
        S = 0
    return S

def pupil_function(x, y):
    # binary function that evaluates to 1 inside the projection and 0 outside it.
    # so for now just return 1
    return 1

def generalised_pupil_function(x, y, R, C, A, Sm, d, wavelength):
    # P is the result of the pupil function
    # k is the spherical wavenumber 
    # W is the wavefront aberration function

    P = pupil_function(x, y)
    # wavelength = 550 # used in spherical wavenumber, needs to be performed for each channel
    k = (2*np.pi) / wavelength
    j = np.sqrt(-1+0j)
    W = aberrated_wavefront_function(x, y, R, C, A, S_calc(Sm, d))
    GPF = P * np.exp(j * k * W)
    # print(generalised_pupil_function)
    # point_spread_function = abs(np.fft.fft(list(generalised_pupil_function))) ** 2
    # print(point_spread_function)
    return GPF

# print(PSF(1,2,3,4,5,6, 5))
image()
#look at deconvolution processes probably implement a few 

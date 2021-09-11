# -*- coding: utf-8 -*-

# This script does not introduce the normalization to the average gain function.
# Also, it was calculated for an equal superposition (|+>, |->).
from numpy import sin, cos, exp, pi, sqrt, absolute
from numpy import inf
from numpy.math cimport INFINITY
from scipy.integrate import nquad
import sys
import time
import warnings
warnings.filterwarnings("error")

DEF sigmas_difference = 9000 # [meV] = 9 [eV] = sigma_e - sigma_h
DEF h_bar = 0.6582119569 # [meV⋅ps] = 6.582119569 * 1e−16 [eV⋅s]
DEF N = 1 # [dimensionless constant]
DEF l_z = 1 # [nm]
DEF l_xy = 5 # [nm]
DEF c = 5.1 # [nm/ps]
DEF rho = 33454.4886 # [meV⋅ps/nm] = 5360 [kg/m^3]
DEF k_B = 0.08617333262145 # [meV/K]
DEF E = 1000 # + 0.01046250062870621 # [meV] = 1 eV + sum |g_k|^2 <---- RECALCULATE ADDED VALUE!!!

DEF limit = 10000 # Limit for the number of subdivisions during integration. This value is way too big, but it doesn't increase the runtime.
options={'limit': limit}
# print("##### LIMIT: ", limit, " #####")

# It's called "adjusted", because it is a result of multiplication of the constant comming from f_k
# and constant coming from replacement of a sum with an integral.
# WARNING 1: It also takes into account the power of 2, coming from |g_k|^2 = |f_k / (h_bar * omega_k)|^2!
# WARNING 2: gaussian_squared() also returns the square of the gaussian, so it shuouldn't be squared while calculating the integrals!
cdef inline double adjusted_constant():
    return (N/((2*pi)**2)) * (1/(2*rho*h_bar*(c**3))) * (sigmas_difference**2)

cdef inline double gaussian_squared(double theta, double r):
    return exp(-(1/2) * ((l_z**2)*(r**2)*(cos(theta)**2) + (l_xy**2)*(r**2)*(sin(theta)**2)))

cdef inline double n_k(float T, double r):
    cdef double beta, z
    beta = (1/(k_B*T))
    z = exp(-beta*c*r*h_bar)
    return z/(1-z)

def calc_W_one_part(float t, float T):
    # Calculate the imaginary integral
    cdef double complex integral_imag
    cdef float error_bound_imag
    print("Before")
    integral_imag, error_bound_imag = nquad(lambda theta, r: sin(theta) * r * exp(-(1/2) * ((l_z**2)*(r**2)*(cos(theta)**2) + (l_xy**2)*(r**2)*(sin(theta)**2))) * sin(c*r*t), [[0, pi], [0, inf]], opts=[options,options])
    # integral_imag, error_bound_imag = nquad(lambda theta, r: sin(theta) * r * gaussian_squared(theta, r) * sin(c*r*t), [[0, pi], [0, 100]], opts=[options,options])
    print("integral_imag: ", integral_imag)
    print("After")
    """
    integral_imag *= 1j
    integral_imag *= adjusted_constant()

    # Calculate the real integral
    integral_real, error_bound_real = nquad(lambda theta, r: sin(theta) * r * gaussian_squared(theta, r) * (cos(c*r*t) - 1)*(1 + (2*n_k(T, r))), [[0, pi], [0, inf]], opts=[options,options])
    integral_real *= adjusted_constant()

    # print("### W_one_part ###")
    # print("integral_imag: ", integral_imag)
    # print("integral_real: ", integral_real)

    return exp(integral_imag)*exp(integral_real)
    """
    return 0

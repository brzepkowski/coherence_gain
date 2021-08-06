# This script does not introduce the normalization to the average gain function.
# Also, it was calculated for an equal superposition (|+>, |->).
from numpy import sin, cos, exp, pi, inf, sqrt, absolute
from scipy.integrate import nquad
import sys
import time
import warnings
warnings.filterwarnings("error")

# Newest constants
# sigmas_difference = 9000 # [meV] = 9 [eV] = sigma_e - sigma_h
# h_bar = 0.6582119569 # [meV⋅ps] = 6.582119569 * 1e−16 [eV⋅s]
# # N = 1 # [dimensionless constant]
# l_z = 1 # [nm]
# l_xy = 5 # [nm]
# c = 5.1 # [nm/ps]
# rho = 33454.4886 # [meV⋅ps/nm] = 5360 [kg/m^3]
# k_B = 0.08617333262145 # [meV/K]
# E = 1000 # + 0.01046250062870621 # [meV] = 1 eV + sum |g_k|^2 <---- RECALCULATE ADDED VALUE!!!

# Constants from paper "“Which path” decoherence in quantum dot experiments"
h_bar = 0.6582119569 # [meV⋅ps] = 6.582119569 * 1e−16 [eV⋅s]
sigmas_difference = 9500 # [meV] = 9 [eV] = sigma_e - sigma_h
l_z = 1 # [nm]
l_xy = 4 # [nm]
c = 5.15 # [nm/ps]
k_B = 0.08617333262145 # [meV/K]
rho = 33079.9818 # [meV⋅ps/nm] = 5360 [kg/m^3]
E = 1000 # + 0.01046250062870621 # [meV] = 1 eV + sum |g_k|^2 <---- RECALCULATE ADDED VALUE!!!

limit = 10000 # Limit for the number of subdivisions during integration. This value is way too big, but it doesn't increase the runtime.
options={'limit': limit}

# It's called "adjusted", because it is a result of multiplication of the constant comming from f_k
# and constant coming from replacement of a sum with an integral.
# WARNING 1: It also takes into account the power of 2, coming from |g_k|^2 = |f_k / (h_bar * omega_k)|^2!
# WARNING 2: gaussian_squared() also returns the square of the gaussian, so it shouldn't be squared while calculating the integrals!
def adjusted_constant():
    return (1/((2*pi)**2)) * (1/(2*rho*h_bar*(c**5))) * (sigmas_difference**2)

def gaussian_squared(theta, omega):
    return exp(-(1/(2*(c**2))) * ((l_z**2)*(omega**2)*(cos(theta)**2) + (l_xy**2)*(omega**2)*(sin(theta)**2)))

def n_k(T, omega):
    beta = (1/(k_B*T))
    z = exp(-beta*omega*h_bar)
    return z/(1-z)

def calc_W_one_part(t, T):
    # Calculate the imaginary integral
    integral_imag, error_bound_imag = nquad(lambda theta, omega: sin(theta) * omega * gaussian_squared(theta, omega) * sin(omega*t), [[0, pi], [0, inf]], opts=[options,options])
    integral_imag *= -1j
    integral_imag *= adjusted_constant()

    # Calculate the real integral
    integral_real, error_bound_real = nquad(lambda theta, omega: sin(theta) * omega * gaussian_squared(theta, omega) * (cos(omega*t) - 1)*(1 + (2*n_k(T, omega))), [[0, pi], [0, inf]], opts=[options,options])
    integral_real *= adjusted_constant()

    return exp(integral_imag)*exp(integral_real)

# if __name__ == "__main__":
#     print(calc_W_one_part(3, 100))

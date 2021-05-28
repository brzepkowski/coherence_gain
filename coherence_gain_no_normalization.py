# This script does not introduce the normalization to the average gain function.
# Also, it was calculated for an equal superposition (|+>, |->).
from numpy import sin, cos, exp, pi, inf, sqrt, absolute
from scipy.integrate import nquad
import sys
import time
import warnings
warnings.filterwarnings("error")

sigmas_difference = 9000 # [meV] = 9 [eV] = sigma_e - sigma_h
h_bar = 0.6582119569 # [meV⋅ps] = 6.582119569 * 1e−16 [eV⋅s]
N = 1 # [dimensionless constant]
l_z = 1 # [nm]
l_xy = 5 # [nm]
c = 5.1 # [nm/ps]
rho = 33454.4886 # [meV⋅ps/nm] = 5360 [kg/m^3]
k_B = 0.08617333262145 # [meV/K]
E = 1000 + 0.01046250062870621 # [meV] = 1 eV + sum |g_k|^2

limit = 10000 # Limit for the number of subdivisions during integration. This value is way too big, but it doesn't increase the runtime.
options={'limit': limit}
print("##### LIMIT: ", limit, " #####")

# It's called "adjusted", because it is a result of multiplication of the constant comming from f_k
# and constant coming from replacement of a sum with an integral.
# WARNING 1: It also takes into account the power of 2, coming from |g_k|^2 = |f_k / (h_bar * omega_k)|^2!
# WARNING 2: gaussian_squared() also returns the square of the gaussian, so it shuouldn't be squared while calculating the integrals!
def adjusted_constant():
    return (N/((2*pi)**3)) * (1/(2*rho*h_bar*(c**3))) * (sigmas_difference**2)

def gaussian_squared(theta, r):
    return exp(-(1/4) * ((l_z**2)*(r**2)*(cos(theta)**2) + (l_xy**2)*(r**2)*(sin(theta)**2)))**2

def calc_E():
    integral, error = nquad(lambda theta, r: sin(theta) * r * gaussian_squared(theta, r), [[0, pi], [0, inf]], opts=[options,options])
    integral *= adjusted_constant()
    # print("integral: ", integral)

def calc_W_one_part(t, T):
    # Calculate the imaginary integral
    integral_imag, error_bound_imag = nquad(lambda theta, r: sin(theta) * r * gaussian_squared(theta, r) * sin(c*r*t), [[0, pi], [0, inf]], opts=[options,options])
    integral_imag *= 1j
    integral_imag *= adjusted_constant()

    # Calculate the real integral
    n_k = (1 / (k_B * T))
    integral_real, error_bound_real = nquad(lambda theta, r: sin(theta) * r * gaussian_squared(theta, r) * (cos(c*r*t) - 1)*(1 + (2*n_k)), [[0, pi], [0, inf]], opts=[options,options])
    integral_real *= adjusted_constant()

    # print("### W_one_part ###")
    # print("integral_imag: ", integral_imag)
    # print("integral_real: ", integral_real)

    return exp(integral_imag)*exp(integral_real)

def calc_W_two_parts(t, tau, T):
    # Calculate the imaginary integral
    integral_imag, error_bound_imag = nquad(lambda theta, r: sin(theta) * r * gaussian_squared(theta, r) * ((2*sin(c*r*tau)) + sin(t-(2*tau))), [[0, pi], [0, inf]], opts=[options,options])
    integral_imag *= 1j
    integral_imag *= adjusted_constant()

    # Calculate the real integral
    n_k = (1 / (k_B * T))
    integral_real, error_bound_real = nquad(lambda theta, r: sin(theta) * r * gaussian_squared(theta, r) * ((2*cos(c*r*tau)) + (2*cos(c*r*(t-tau))) - cos(c*r*(t-(2*tau))) - 3)*(1 + (2*n_k)), [[0, pi], [0, inf]], opts=[options,options])
    integral_real *= adjusted_constant()

    # print("### W_two_parts ###")
    # print("integral_imag: ", integral_imag)
    # print("integral_real: ", integral_real)

    return exp(integral_imag)*exp(integral_real)

def calc_W_three_parts(t, tau, T):
    # Calculate the imaginary integral
    integral_imag, error_bound_imag = nquad(lambda theta, r: sin(theta) * r * gaussian_squared(theta, r) * ((2*sin(c*r*tau)) + (2*sin(c*r*(t-(2*tau)))) - sin(c*r*(t-tau))), [[0, pi], [0, inf]], opts=[options,options])
    integral_imag *= 1j
    integral_imag *= adjusted_constant()

    # Calculate the real integral
    n_k = (1 / (k_B * T))
    integral_real, error_bound_real = nquad(lambda theta, r: sin(theta) * r * gaussian_squared(theta, r) * (cos(c*r*(t-tau)) - 1)*(1 + (2*n_k)), [[0, pi], [0, inf]], opts=[options,options])
    integral_real *= adjusted_constant()

    # print("### W_three_parts ###")
    # print("integral_imag: ", integral_imag)
    # print("integral_real: ", integral_real)

    return exp(integral_imag)*exp(integral_real)

def calc_g_av(t, tau, T):
    try:
        W_basic = calc_W_one_part(t, T)

        W_0 = calc_W_one_part(t-tau, T)
        W_1 = calc_W_one_part(t-(2*tau), T)
        W_2 = calc_W_two_parts(t, tau, T)
        W_3 = calc_W_three_parts(t, tau, T)

        first_part = absolute((1/4)*(W_0 + exp(-1j*E*tau/h_bar)*W_1 + exp(1j*E*tau/h_bar)*W_2 + W_3))
        second_part = absolute((1/4)*(-W_0 + exp(-1j*E*tau/h_bar)*W_1 + exp(1j*E*tau/h_bar)*W_2 - W_3))
        third_part = absolute(W_basic)
    except Warning as w:
        print(w)
        print("Warning detected - execution terminated.")
        sys.exit()

    return first_part + second_part - third_part


if __name__ == "__main__":
    # Test calculation of g_av
    print(calc_g_av(10, 1, 1))

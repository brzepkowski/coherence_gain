# This script does not introduce the normalization to the average gain function.
# Also, it was calculated for an equal superposition (|+>, |->).

from numpy import sin, cos, exp, pi, inf, sqrt, absolute
from scipy.integrate import dblquad
import sys

sigma_e = 2
sigma_h = 1
h_bar = 1
rho = 1
N = 1
V = 1
c = 500
k_B = 1.380649
l_e = 1
l_z = 4
E = 2

# It's called "adjusted", because it is a result of multiplication of the constant comming from f_k
# and constant coming from replacement of a sum with an integral.
# WARNING 1: It also takes into account the power of 2, coming from |g_k|^2 = |f_k / (h_bar * omega_k)|^2!
# WARNING 2: gaussian_squared() also returns the square of the gaussian, so it shuouldn't be squared while calculating the integrals!
def adjusted_constant():
    return (N/((2*pi)**3)) * (1/(2*rho*h_bar*(c**3))) * ((sigma_e - sigma_h)**2)

def gaussian_squared(theta, r):
    return exp(-(1/4) * ((l_z**2)*(r**2)*(cos(theta)**2) + (l_e**2)*(r**2)*(sin(theta)**2)))**2

def calc_W_one_part(t, T):
    # Calculate the imaginary integral
    integral_imag, error_bound_imag = dblquad(lambda theta, r: cos(theta) * r * gaussian_squared(theta, r) * sin(c*r*t), -pi/2, pi/2, lambda r: 0, lambda r: inf)
    integral_imag *= 1j
    integral_imag *= adjusted_constant()

    # Calculate the real integral
    n_k = (1 / (k_B * T))
    integral_real, error_bound_real = dblquad(lambda theta, r: cos(theta) * r * gaussian_squared(theta, r) * (cos(c*r*t) - 1)*(1 + (2*n_k)), -pi/2, pi/2, lambda r: 0, lambda r: inf)
    integral_real *= adjusted_constant()

    print("### W_0 ###")
    print("integral_imag: ", integral_imag)
    print("integral_real: ", integral_real)

    return exp(integral_imag)*exp(integral_real)

def calc_W_two_parts(t, tau, T):
    # Calculate the imaginary integral
    integral_imag, error_bound_imag = dblquad(lambda theta, r: cos(theta) * r * gaussian_squared(theta, r) * ((2*sin(c*r*tau)) + sin(t-(2*tau))), -pi/2, pi/2, lambda r: 0, lambda r: inf)
    integral_imag *= 1j
    integral_imag *= adjusted_constant()

    # Calculate the real integral
    n_k = (1 / (k_B * T))
    integral_real, error_bound_real = dblquad(lambda theta, r: cos(theta) * r * gaussian_squared(theta, r) * ((2*cos(c*r*tau)) + (2*cos(c*r*(t-tau))) - cos(c*r*(t-(2*tau))) - 3)*(1 + (2*n_k)), -pi/2, pi/2, lambda r: 0, lambda r: inf)
    integral_real *= adjusted_constant()

    print("### W_0 ###")
    print("integral_imag: ", integral_imag)
    print("integral_real: ", integral_real)

    return exp(integral_imag)*exp(integral_real)

def calc_W_three_parts(t, tau, T):
    # Calculate the imaginary integral
    integral_imag, error_bound_imag = dblquad(lambda theta, r: cos(theta) * r * gaussian_squared(theta, r) * ((2*sin(c*r*tau)) + (2*sin(c*r*(t-(2*tau)))) - sin(c*r*(t-tau))), -pi/2, pi/2, lambda r: 0, lambda r: inf)
    integral_imag *= 1j
    integral_imag *= adjusted_constant()

    # Calculate the real integral
    n_k = (1 / (k_B * T))
    integral_real, error_bound_real = dblquad(lambda theta, r: cos(theta) * r * gaussian_squared(theta, r) * (cos(c*r*(t-tau)) - 1)*(1 + (2*n_k)), -pi/2, pi/2, lambda r: 0, lambda r: inf)
    integral_real *= adjusted_constant()

    print("### W_0 ###")
    print("integral_imag: ", integral_imag)
    print("integral_real: ", integral_real)

    return exp(integral_imag)*exp(integral_real)

def g_av(t, tau, T):
    W_basic = calc_W_one_part(t, T)

    W_0 = calc_W_one_part(t-tau, T)
    W_1 = calc_W_one_part(t-(2*tau), T)
    W_2 = calc_W_two_parts(t, tau, T)
    W_3 = calc_W_three_parts(t, tau, T)

    first_part = absolute((1/4)*(W_0 + exp(-1j*E*tau/h_bar)*W_1 + exp(1j*E*tau/h_bar)*W_2 + W_3))
    second_part = absolute((1/4)*(-W_0 + exp(-1j*E*tau/h_bar)*W_1 + exp(1j*E*tau/h_bar)*W_2 - W_3))
    third_part = absolute(W_basic)

    return first_part + second_part - third_part

print("g_av: ", g_av(2,1,1))

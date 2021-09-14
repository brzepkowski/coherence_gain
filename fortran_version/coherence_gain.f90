module coherence_gain
    implicit none
    real, parameter :: sigmas_difference = 9000 ! [meV] = 9 [eV] = sigma_e - sigma_h
    real, parameter :: h_bar = 0.6582119569 ! [meV⋅ps] = 6.582119569 * 1e−16 [eV⋅s]
    real, parameter :: N = 1 ! [dimensionless constant]
    real, parameter :: l_z = 1 ! [nm]
    real, parameter :: l_xy = 5 ! [nm]
    real, parameter :: c = 5.1 ! [nm/ps]
    real, parameter :: rho = 33454.4886 ! [meV⋅ps/nm] = 5360 [kg/m^3]
    real, parameter :: k_B = 0.08617333262145 ! [meV/K]
    real, parameter :: E = 1000 ! + 0.01046250062870621 # [meV] = 1 eV + sum |g_k|^2 <---- RECALCULATE ADDED VALUE!!!
    real(kind=8), parameter :: pi=4.D0*datan(1.D0)

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!   Subroutines for integration   !!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine integral_finite(f, a, b, result)
    implicit none
    real :: a, abserr, b, epsabs, epsrel, f, res, work, result
    integer :: ier, iwork, key, lenw, limit, neval, last
    parameter(epsabs=0.1, epsrel=0.1)
    parameter(key=2, limit=10000)
    parameter(lenw=4*limit+1)
    dimension :: iwork(limit), work(lenw)
    external f

    call qag(f,a,b,epsabs,epsrel,key,res,abserr,neval,ier,limit,lenw,last,iwork,work)

    result=res
end subroutine integral_finite


subroutine integral_infinite(f, bound, inf, result)
	implicit none
	real :: abserr, epsabs, epsrel, f, res, work, result, bound
	integer :: ier, iwork, key, lenw, limit, neval, inf, last
	! parameter(epsabs=0.000001, epsrel=0.000001)
    parameter(epsabs=0.001, epsrel=0.001)
	parameter(key=2, limit=10000)
	parameter(lenw=4*limit+1)
	dimension :: iwork(limit), work(lenw)
	external f

	call qagi(f,bound,inf,epsabs,epsrel,res,abserr,neval,ier,limit,lenw,last,iwork,work)

	result=res
end subroutine integral_infinite


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!   Basic functions   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! It's called "adjusted", because it is a result of multiplication of the constant comming from f_k
! and constant coming from replacement of a sum with an integral.
! WARNING 1: It also takes into account the power of 2, coming from |g_k|^2 = |f_k / (h_bar * omega_k)|^2!
! WARNING 2: gaussian_squared() also returns the square of the gaussian, so it shuouldn't be squared while calculating the integrals!
real function adjusted_constant()
    implicit none
    adjusted_constant = (N/((2*pi)**2)) * (1/(2*rho*h_bar*(c**3))) * (sigmas_difference**2)
end function adjusted_constant


real function gaussian_squared(theta, r)
    implicit none
    real :: theta, r

    gaussian_squared = exp(-0.5 * ((l_z**2)*(r**2)*(cos(theta)**2) + (l_xy**2)*(r**2)*(sin(theta)**2)))
end function gaussian_squared


real function n_k(T, r)
    implicit none
    real :: T, r
    real :: beta, z

    beta = (1/(k_B*T))
    z = exp(-beta*c*r*h_bar)
    n_k = z/(1-z)
end function n_k


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!   Weyl operators   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! W_one_part
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

double complex function W_one_part(t_time, T_temp)
    real :: t_time, T_temp
    real :: theta, r_val, t_time_val, T_temp_val
    real :: integral_imag_partial, integral_real
    double complex :: integral_imag
    double complex, parameter :: i = cmplx(0, 1)

    integral_imag_partial = W_one_part_imag(t_time)
    integral_imag = integral_imag_partial * i
    integral_imag = integral_imag * adjusted_constant()
    ! print *, integral_imag

    integral_real = W_one_part_real(t_time, T_temp)
    integral_real = integral_real * adjusted_constant()
    ! print *, integral_real

    W_one_part = exp(integral_imag)*exp(integral_real)
end function W_one_part

! Imaginary helper functions
real function W_one_part_imag(t_time)
    implicit none
    real :: t_time, r_val, t_time_val, result
    real, parameter :: r_min = 0.
    integer, parameter :: r_max = 1 ! In fact this corresponds to infinity

    common /one_part_imag_values/r_val
    common /one_part_imag_values/t_time_val

    t_time_val = t_time

    call integral_infinite(W_one_part_imag_2, r_min, r_max, result)
    W_one_part_imag = result
end function W_one_part_imag

real function W_one_part_imag_2(r)
    implicit none
    real :: r, r_val, t_time_val, result
    real, parameter :: theta_min = 0
    real, parameter :: theta_max = pi

    common /one_part_imag_values/r_val
    common /one_part_imag_values/t_time_val

    r_val = r

    call integral_finite(W_one_part_imag_1, theta_min, theta_max, result)
    W_one_part_imag_2 = result
end function W_one_part_imag_2

real function W_one_part_imag_1(theta)
    implicit none
    real :: theta, r_val, t_time_val

    common /one_part_imag_values/r_val
    common /one_part_imag_values/t_time_val

    W_one_part_imag_1 = W_one_part_imag_0(theta, r_val, t_time_val)
end function W_one_part_imag_1

real function W_one_part_imag_0(theta, r, t_time)
    implicit none
    real :: theta, r, t_time

    W_one_part_imag_0 = sin(theta) * r * gaussian_squared(theta, r) * sin(c*r*t_time)
end function W_one_part_imag_0

! Real helper functions
real function W_one_part_real(t_time, T_temp)
    implicit none
    real :: t_time, T_temp, r_val, t_time_val, T_temp_val, result
    real, parameter :: r_min = 0.
    integer, parameter :: r_max = 1 ! In fact this corresponds to infinity

    common /one_part_real_values/r_val ! This line is placed only to avoid the 'Warning: Named COMMON block ‘values’ at (1) shall be of the same size as elsewhere (8 vs 12 bytes)' warning.
    common /one_part_real_values/t_time_val
    common /one_part_real_values/T_temp_val

    t_time_val = t_time
    T_temp_val = T_temp

    call integral_infinite(W_one_part_real_2, r_min, r_max, result)
    W_one_part_real = result
end function W_one_part_real

real function W_one_part_real_2(r)
    implicit none
    real :: r, r_val, t_time_val, T_temp_val, result
    real, parameter :: theta_min = 0
    real, parameter :: theta_max = pi

    common /one_part_real_values/r_val
    common /one_part_real_values/t_time_val
    common /one_part_real_values/T_temp_val

    r_val = r

    call integral_finite(W_one_part_real_1, theta_min, theta_max, result)
    W_one_part_real_2 = result
end function W_one_part_real_2

real function W_one_part_real_1(theta)
    implicit none
    real :: theta, r_val, t_time_val, T_temp_val

    common /one_part_real_values/r_val
    common /one_part_real_values/t_time_val
    common /one_part_real_values/T_temp_val

    W_one_part_real_1 = W_one_part_real_0(theta, r_val, t_time_val, T_temp_val)
end function W_one_part_real_1

real function W_one_part_real_0(theta, r, t_time, T_temp)
    implicit none
    real :: theta, r, t_time, T_temp

    W_one_part_real_0 = sin(theta) * r * gaussian_squared(theta, r) * (cos(c*r*t_time) - 1)*(1 + (2*n_k(T_temp, r)))
end function W_one_part_real_0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! W_two_parts
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

double complex function W_two_parts(t_time, tau_time, T_temp)
    real :: t_time, tau_time, T_temp
    real :: theta, r_val, t_time_val, tau_time_val, T_temp_val
    real :: integral_imag_partial, integral_real
    double complex :: integral_imag
    double complex, parameter :: i = cmplx(0, 1)

    integral_imag_partial = W_two_parts_imag(t_time, tau_time)
    integral_imag = integral_imag_partial * i
    integral_imag = integral_imag * adjusted_constant()
    ! print *, integral_imag_partial

    integral_real = W_two_parts_real(t_time, tau_time, T_temp)
    integral_real = integral_real * adjusted_constant()
    ! print *, integral_real

    W_two_parts = exp(integral_imag)*exp(integral_real)
end function W_two_parts

! Imaginary helper functions
real function W_two_parts_imag(t_time, tau_time)
    implicit none
    real :: t_time, tau_time, r_val, t_time_val, tau_time_val, result
    real, parameter :: r_min = 0.
    integer, parameter :: r_max = 1 ! In fact this corresponds to infinity

    common /two_parts_imag_values/r_val
    common /two_parts_imag_values/t_time_val
    common /two_parts_imag_values/tau_time_val

    t_time_val = t_time
    tau_time_val = tau_time

    call integral_infinite(W_two_parts_imag_2, r_min, r_max, result)
    W_two_parts_imag = result
end function W_two_parts_imag

real function W_two_parts_imag_2(r)
    implicit none
    real :: r, r_val, t_time_val, tau_time_val, result
    real, parameter :: theta_min = 0
    real, parameter :: theta_max = pi

    common /two_parts_imag_values/r_val
    common /two_parts_imag_values/t_time_val
    common /two_parts_imag_values/tau_time_val

    r_val = r

    call integral_finite(W_two_parts_imag_1, theta_min, theta_max, result)
    W_two_parts_imag_2 = result
end function W_two_parts_imag_2

real function W_two_parts_imag_1(theta)
    implicit none
    real :: theta, r_val, t_time_val, tau_time_val

    common /two_parts_imag_values/r_val
    common /two_parts_imag_values/t_time_val
    common /two_parts_imag_values/tau_time_val

    W_two_parts_imag_1 = W_two_parts_imag_0(theta, r_val, t_time_val, tau_time_val)
end function W_two_parts_imag_1

real function W_two_parts_imag_0(theta, r, t_time, tau_time)
    implicit none
    real :: theta, r, t_time, tau_time

    W_two_parts_imag_0 =  sin(theta) * r * gaussian_squared(theta, r) * ((2*sin(c*r*tau_time)) + sin(t_time-(2*tau_time)))
end function W_two_parts_imag_0

! Real helper functions
real function W_two_parts_real(t_time, tau_time, T_temp)
    implicit none
    real :: t_time, tau_time, T_temp, r_val, t_time_val, tau_time_val, T_temp_val, result
    real, parameter :: r_min = 0.
    integer, parameter :: r_max = 1 ! In fact this corresponds to infinity

    common /two_parts_real_values/r_val ! This line is placed only to avoid the 'Warning: Named COMMON block ‘values’ at (1) shall be of the same size as elsewhere (8 vs 12 bytes)' warning.
    common /two_parts_real_values/t_time_val
    common /two_parts_real_values/tau_time_val
    common /two_parts_real_values/T_temp_val

    t_time_val = t_time
    tau_time_val = tau_time
    T_temp_val = T_temp

    call integral_infinite(W_two_parts_real_2, r_min, r_max, result)
    W_two_parts_real = result
end function W_two_parts_real

real function W_two_parts_real_2(r)
    implicit none
    real :: r, r_val, t_time_val, tau_time_val, T_temp_val, result
    real, parameter :: theta_min = 0
    real, parameter :: theta_max = pi

    common /two_parts_real_values/r_val
    common /two_parts_real_values/t_time_val
    common /two_parts_real_values/tau_time_val
    common /two_parts_real_values/T_temp_val

    r_val = r

    call integral_finite(W_two_parts_real_1, theta_min, theta_max, result)
    W_two_parts_real_2 = result
end function W_two_parts_real_2

real function W_two_parts_real_1(theta)
    implicit none
    real :: theta, r_val, t_time_val, tau_time_val, T_temp_val

    common /two_parts_real_values/r_val
    common /two_parts_real_values/t_time_val
    common /two_parts_real_values/tau_time_val
    common /two_parts_real_values/T_temp_val

    W_two_parts_real_1 = W_two_parts_real_0(theta, r_val, t_time_val, tau_time_val, T_temp_val)
end function W_two_parts_real_1

real function W_two_parts_real_0(theta, r, t_time, tau_time, T_temp)
    implicit none
    real :: theta, r, t_time, tau_time, T_temp

    W_two_parts_real_0 = sin(theta) * r * gaussian_squared(theta, r)
    W_two_parts_real_0=W_two_parts_real_0*((2*cos(c*r*tau_time))+(2*cos(c*r*(t_time-tau_time)))-cos(c*r*(t_time-(2*tau_time)))-3)
    W_two_parts_real_0 = W_two_parts_real_0 * (1 + (2*n_k(T_temp, r)))
end function W_two_parts_real_0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! W_three_parts
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

double complex function W_three_parts(t_time, tau_time, T_temp)
    real :: t_time, tau_time, T_temp
    real :: theta, r_val, t_time_val, tau_time_val, T_temp_val
    real :: integral_imag_partial, integral_real
    double complex :: integral_imag
    double complex, parameter :: i = cmplx(0, 1)

    integral_imag_partial = W_three_parts_imag(t_time, tau_time)
    integral_imag = integral_imag_partial * i
    integral_imag = integral_imag * adjusted_constant()
    ! print *, integral_imag_partial

    integral_real = W_three_parts_real(t_time, tau_time, T_temp)
    integral_real = integral_real * adjusted_constant()
    ! print *, integral_real

    W_three_parts = exp(integral_imag)*exp(integral_real)
end function W_three_parts

! Imaginary helper functions
real function W_three_parts_imag(t_time, tau_time)
    implicit none
    real :: t_time, tau_time, r_val, t_time_val, tau_time_val, result
    real, parameter :: r_min = 0.
    integer, parameter :: r_max = 1 ! In fact this corresponds to infinity

    common /three_parts_imag_values/r_val
    common /three_parts_imag_values/t_time_val
    common /three_parts_imag_values/tau_time_val

    t_time_val = t_time
    tau_time_val = tau_time

    call integral_infinite(W_three_parts_imag_2, r_min, r_max, result)
    W_three_parts_imag = result
end function W_three_parts_imag

real function W_three_parts_imag_2(r)
    implicit none
    real :: r, r_val, t_time_val, tau_time_val, result
    real, parameter :: theta_min = 0
    real, parameter :: theta_max = pi

    common /three_parts_imag_values/r_val
    common /three_parts_imag_values/t_time_val
    common /three_parts_imag_values/tau_time_val

    r_val = r

    call integral_finite(W_three_parts_imag_1, theta_min, theta_max, result)
    W_three_parts_imag_2 = result
end function W_three_parts_imag_2

real function W_three_parts_imag_1(theta)
    implicit none
    real :: theta, r_val, t_time_val, tau_time_val

    common /three_parts_imag_values/r_val
    common /three_parts_imag_values/t_time_val
    common /three_parts_imag_values/tau_time_val

    W_three_parts_imag_1 = W_three_parts_imag_0(theta, r_val, t_time_val, tau_time_val)
end function W_three_parts_imag_1

real function W_three_parts_imag_0(theta, r, t_time, tau_time)
    implicit none
    real :: theta, r, t_time, tau_time

    W_three_parts_imag_0 = sin(theta) * r * gaussian_squared(theta, r)
    W_three_parts_imag_0=W_three_parts_imag_0*((2*sin(c*r*tau_time))+(2*sin(c*r*(t_time-(2*tau_time))))-sin(c*r*(t_time-tau_time)))
end function W_three_parts_imag_0

! Real helper functions
real function W_three_parts_real(t_time, tau_time, T_temp)
    implicit none
    real :: t_time, tau_time, T_temp, r_val, t_time_val, tau_time_val, T_temp_val, result
    real, parameter :: r_min = 0.
    integer, parameter :: r_max = 1 ! In fact this corresponds to infinity

    common /three_parts_real_values/r_val ! This line is placed only to avoid the 'Warning: Named COMMON block ‘values’ at (1) shall be of the same size as elsewhere (8 vs 12 bytes)' warning.
    common /three_parts_real_values/t_time_val
    common /three_parts_real_values/tau_time_val
    common /three_parts_real_values/T_temp_val

    t_time_val = t_time
    tau_time_val = tau_time
    T_temp_val = T_temp

    call integral_infinite(W_three_parts_real_2, r_min, r_max, result)
    W_three_parts_real = result
end function W_three_parts_real

real function W_three_parts_real_2(r)
    implicit none
    real :: r, r_val, t_time_val, tau_time_val, T_temp_val, result
    real, parameter :: theta_min = 0
    real, parameter :: theta_max = pi

    common /three_parts_real_values/r_val
    common /three_parts_real_values/t_time_val
    common /three_parts_real_values/tau_time_val
    common /three_parts_real_values/T_temp_val

    r_val = r

    call integral_finite(W_three_parts_real_1, theta_min, theta_max, result)
    W_three_parts_real_2 = result
end function W_three_parts_real_2

real function W_three_parts_real_1(theta)
    implicit none
    real :: theta, r_val, t_time_val, tau_time_val, T_temp_val

    common /three_parts_real_values/r_val
    common /three_parts_real_values/t_time_val
    common /three_parts_real_values/tau_time_val
    common /three_parts_real_values/T_temp_val

    W_three_parts_real_1 = W_three_parts_real_0(theta, r_val, t_time_val, tau_time_val, T_temp_val)
end function W_three_parts_real_1

real function W_three_parts_real_0(theta, r, t_time, tau_time, T_temp)
    implicit none
    real :: theta, r, t_time, tau_time, T_temp

    W_three_parts_real_0 = sin(theta) * r * gaussian_squared(theta, r) * (cos(c*r*(t_time-tau_time)) - 1)*(1 + (2*n_k(T_temp, r)))
end function W_three_parts_real_0


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   g_av    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine g_average(t_time,tau_time,T_temp,W_0,W_1,W_1_phased,W_2,W_2_phased,W_3,phase,D_plus,D_minus,D_t,p_plus,p_minus,g_av)
    implicit none
    real :: t_time, tau_time, T_temp
    real :: p_plus, p_minus, denominator_plus, denominator_minus
    real :: D_plus, D_minus, D_t
    real :: g_plus, g_minus, g_av
    double complex :: W_t, W_tau, W_0, W_1, W_1_phased, W_2, W_2_phased, W_3, phase
    double complex, parameter :: i = cmplx(0, 1)

    W_t = W_one_part(t_time, T_temp)
    W_tau = W_one_part(tau_time, T_temp)

    W_0 = W_one_part(t_time-tau_time, T_temp)
    W_1 = W_one_part(t_time-(2*tau_time), T_temp)
		W_1_phased = exp(-i*E*tau_time/h_bar)*W_1
    W_2 = W_two_parts(t_time, tau_time, T_temp)
		W_2_phased = exp(i*E*tau_time/h_bar)*W_2
    W_3 = W_three_parts(t_time, tau_time, T_temp)
		phase = exp(i*E*tau_time/h_bar)

    p_plus = 0.5*(1 + real(exp(i*E*tau_time/h_bar)*W_tau))
    p_minus = 0.5*(1 - real(exp(i*E*tau_time/h_bar)*W_tau))
    denominator_plus = abs(p_plus)
    denominator_minus = abs(p_minus)

    D_plus = abs(0.25*(W_0 + exp(-i*E*tau_time/h_bar)*W_1 + exp(i*E*tau_time/h_bar)*W_2 + W_3))/denominator_plus
    D_minus = abs(0.25*(-W_0 + exp(-i*E*tau_time/h_bar)*W_1 + exp(i*E*tau_time/h_bar)*W_2 - W_3))/denominator_minus
    D_t = abs(W_t)

    g_plus = (D_plus - D_t)/(1-D_t)
    g_minus = (D_minus - D_t)/(1-D_t)

    g_av = (p_plus*g_plus) + (p_minus*g_minus)
end subroutine g_average


end module coherence_gain

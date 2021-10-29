program main
    use coherence_gain
    implicit none
    real :: tau_time, tau_min_time, tau_max_time ! Taus corresponding to minimal and maximal values of g_av (calculated for t = 20. and starting point of tau = 6.0)
    real :: t_time, t_time_min, t_time_step, T_temp
    real :: p_plus, p_minus
    real :: D_plus, D_minus, D_t
    real :: g_av, period
    double complex :: W_t, W_tau, W_0, W_1, W_1_phased, W_2, W_2_phased, W_3, phase
    integer :: iterator, number_of_time_steps, taus_iterator, num_of_different_taus
    character(40) :: filename
    character(10) :: tau_time_string
    character(4) :: T_temp_string

    num_of_different_taus = 2
    number_of_time_steps = 10000
    period = (2 * pi * h_bar) / E
    print *, "period: ", period

    tau_min_time = 6.00502014
    tau_max_time = 6.00604010
    t_time_step = 0.0001
    T_temp = 70.0
    write(T_temp_string, '(f4.1)' )  T_temp

    ! Generate data for minimas of the g_av
    tau_time = tau_min_time


    t_time = tau_time

    call g_average(t_time,tau_time,T_temp,W_0,W_1,W_1_phased,W_2,W_2_phased, &
                    W_3,phase,D_plus,D_minus,D_t,p_plus,p_minus,g_av)

    print *, "g_av: "
    print *, g_av

    ! do taus_iterator = 0, num_of_different_taus-1
    !     tau_time = tau_time + taus_iterator*500*period
    !     t_time_min = tau_time + t_time_step
    !
    !     write(tau_time_string, '(f10.8)' )  tau_time
    !     filename = "min_g_av_T="//T_temp_string//"_tau="//tau_time_string//".dat"
    !     open(1, file=filename, status='replace')
    !
    !     do iterator = 0, number_of_time_steps
    !         print *, iterator, " / ", number_of_time_steps
    !         t_time = t_time_min + t_time_step*iterator
    !         call g_average(t_time,tau_time,T_temp,W_0,W_1,W_1_phased,W_2,W_2_phased, &
    !                         W_3,phase,D_plus,D_minus,D_t,p_plus,p_minus,g_av)
    !         write(1,*) t_time, W_0, W_1, W_1_phased, W_2, W_2_phased, W_3, phase, D_plus, D_minus, D_t, p_plus, p_minus, g_av
    !     end do
    ! end do
    !
    ! ! Generate data for maximas of the g_av
    ! tau_time = tau_max_time
    !
    ! do taus_iterator = 0, num_of_different_taus-1
    !     tau_time = tau_time + taus_iterator*500*period
    !     t_time_min = tau_time + t_time_step
    !
    !     write(tau_time_string, '(f10.8)' )  tau_time
    !     filename = "max_g_av_T="//T_temp_string//"_tau="//tau_time_string//".dat"
    !     open(1, file=filename, status='replace')
    !
    !     do iterator = 0, number_of_time_steps
    !         print *, iterator, " / ", number_of_time_steps
    !         t_time = t_time_min + t_time_step*iterator
    !         call g_average(t_time,tau_time,T_temp,W_0,W_1,W_1_phased,W_2,W_2_phased, &
    !                         W_3,phase,D_plus,D_minus,D_t,p_plus,p_minus,g_av)
    !         write(1,*) t_time, W_0, W_1, W_1_phased, W_2, W_2_phased, W_3, phase, D_plus, D_minus, D_t, p_plus, p_minus, g_av
    !     end do
    ! end do
end program main

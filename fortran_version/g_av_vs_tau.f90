program main
    use coherence_gain
    implicit none
    real :: t_time, tau_time, tau_time_min, tau_time_step, T_temp
    real :: p_plus, p_minus
    real :: D_plus, D_minus, D_t
    real :: g_av
    double complex :: W_t, W_tau, W_0, W_1, W_1_phased, W_2, W_2_phased, W_3, phase
    integer :: iterator, number_of_iterations
    real :: max_g_av, min_g_av, tau_of_g_av_max, tau_of_g_av_min
    character(30) :: filename
    character(4) :: T_temp_string

    max_g_av = -10000.
    min_g_av = 10000.

    t_time = 20
    tau_time_min = 6.
    tau_time_step = 0.00001
    T_temp = 70

    write(T_temp_string, '(f4.1)' )  T_temp
    filename = "g_av_T="//T_temp_string//".dat"
    open(1, file=filename, status='replace')

    number_of_iterations = 1000
    do iterator = 0, number_of_iterations
        print *, iterator, " / ", number_of_iterations
        tau_time = tau_time_min + tau_time_step*iterator
        call g_average(t_time,tau_time,T_temp,W_0,W_1,W_1_phased,W_2,W_2_phased,W_3,phase,D_plus,D_minus,D_t,p_plus,p_minus,g_av)
        if (g_av > max_g_av) then
            max_g_av = g_av
            tau_of_g_av_max = tau_time
        end if

        if (g_av < min_g_av .and. g_av > 0.) then
            min_g_av = g_av
            tau_of_g_av_min = tau_time
        end if

        write(1,*) tau_time, g_av
    end do
    print *, "tau_of_g_av_max: ", tau_of_g_av_max
    print *, "tau_of_g_av_min: ", tau_of_g_av_min

end program main

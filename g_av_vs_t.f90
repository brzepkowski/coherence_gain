program main
    use coherence_gain
    implicit none
    real :: tau_time
    real :: t_time, t_time_start, t_time_end, t_time_step
    real :: T_temp
    real :: p_plus, p_minus
    real :: D_plus, D_minus, D_t_minus_tau
    real :: g_av, period
    double complex :: W_t, W_tau, W_0, W_1, W_1_phased, W_2, W_2_phased, W_3, phase
    integer :: iterator, number_of_iterations
    character(100) :: filename
    character(12) :: tau_time_string
    character(12) :: t_time_start_string, t_time_end_string
    character(4) :: T_temp_string
    character(100) :: number_of_iterations_string
    character(3) :: min_or_max_type_string ! It can be either 'MIN' or 'MAX'

    period = (2 * pi * h_bar) / E
    t_time_step = 0.00001

    print *, period

    !First, make sure the right number of inputs have been provided
    if (COMMAND_ARGUMENT_COUNT() .ne. 5) then
      write(*,*)'Error, 5 command-line arguments required.'
      write(*,*)'Provide them in the following way:'
      write(*,*)'- MIN/MAX,'
      write(*,*)'- tau_time,'
      write(*,*)'- t_time_start,'
      write(*,*)'- T_temperature,'
      write(*,*)'- number_of_iterations (not smaller than 1000).'
      stop
    endif

    call GET_COMMAND_ARGUMENT(1, min_or_max_type_string)
    call GET_COMMAND_ARGUMENT(2, tau_time_string)
    call GET_COMMAND_ARGUMENT(3, t_time_start_string)
    call GET_COMMAND_ARGUMENT(4, T_temp_string)
    call GET_COMMAND_ARGUMENT(5, number_of_iterations_string)

    ! convert to REALs
    read(tau_time_string,*) tau_time
    read(t_time_start_string,*) t_time_start
    read(T_temp_string,*) T_temp
    read(number_of_iterations_string,*) number_of_iterations

    if (number_of_iterations < 1000) then
      print *, "Too small number_of_iterations! Stopping execution."
      stop
    endif

    write(t_time_end_string, '(f11.9)' )  t_time_start + (number_of_iterations*t_time_step)
    filename = "g_av_vs_t_"//min_or_max_type_string//"_T="//T_temp_string//"_"// &
      tau_time_string//"_"//t_time_start_string//"_"//t_time_end_string//".dat"

    print *, filename
    open(1, file=filename, status='replace')

    do iterator = 0, number_of_iterations
        print *, iterator, " / ", number_of_iterations
        t_time = t_time_start + t_time_step*iterator
        call g_average(t_time,tau_time,T_temp,W_0,W_1,W_1_phased,W_2,W_2_phased, &
                        W_3,phase,D_plus,D_minus,D_t_minus_tau,p_plus,p_minus,g_av)

        pure_phase = exp(i*E*tau_time/h_bar)
        write(1,*) t_time, g_av, W_0, W_1, W_1_phased, W_2, W_2_phased, W_3, D_plus, D_minus, D_t_minus_tau, p_plus, p_minus, pure_phase
    end do


end program main

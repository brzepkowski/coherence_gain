program main
    use coherence_gain
    implicit none
    real :: t_time, tau_time, tau_time_start, tau_time_step, T_temp
    real :: p_plus, p_minus
    real :: D_plus, D_minus, D_t
    real :: g_av
    double complex :: W_0, W_1, W_1_phased, W_2, W_2_phased, W_3, phase
    double complex :: pure_phase
    double complex, parameter :: i = cmplx(0, 1)
    integer :: iterator, number_of_iterations
    real :: max_g_av, min_g_av, tau_of_g_av_max, tau_of_g_av_min
    character(100) :: filename
    character(4) :: tau_time_start_string, tau_time_end_string
    character(4) :: t_time_string
    character(4) :: T_temp_string
    character(100) :: number_of_iterations_string

    max_g_av = -10000.
    min_g_av = 10000.
    tau_time_step = 0.00001

    !First, make sure the right number of inputs have been provided
    if (COMMAND_ARGUMENT_COUNT() .ne. 4) then
      write(*,*)'Error, 4 command-line arguments required.'
      write(*,*)'Provide them in the following way:'
      write(*,*)'- tau_time_start,'
      write(*,*)'- t_time,'
      write(*,*)'- T_temperature,'
      write(*,*)'- number_of_iterations (not smaller than 1000).'
      stop
    endif

    call GET_COMMAND_ARGUMENT(1, tau_time_start_string)   !first, read in the two values
    call GET_COMMAND_ARGUMENT(2, t_time_string)
    call GET_COMMAND_ARGUMENT(3, T_temp_string)
    call GET_COMMAND_ARGUMENT(4, number_of_iterations_string)

    read(tau_time_start_string,*) tau_time_start                    !then, convert them to REALs
    read(t_time_string,*) t_time
    read(T_temp_string,*) T_temp
    read(number_of_iterations_string,*) number_of_iterations

    if (number_of_iterations < 1000) then
      print *, "Too small number_of_iterations! Stopping execution."
      stop
    endif

    write(tau_time_end_string, '(f4.2)' )  tau_time_start + (number_of_iterations*tau_time_step)
    filename = "g_av_vs_tau_T="//T_temp_string//"_"//t_time_string//"_"//tau_time_start_string//"_"//tau_time_end_string//".dat"

    print *, filename
    open(1, file=filename, status='replace')

    do iterator = 0, number_of_iterations
        print *, iterator, " / ", number_of_iterations
        tau_time = tau_time_start + tau_time_step*iterator
        call g_average(t_time,tau_time,T_temp,W_0,W_1,W_1_phased,W_2,W_2_phased,W_3,phase, &
            D_plus,D_minus,D_t,p_plus,p_minus,g_av)

        pure_phase = exp(i*E*tau_time/h_bar)
        write(1,*) tau_time, g_av, D_plus, D_minus, D_t, p_plus, p_minus, pure_phase
    end do

end program main

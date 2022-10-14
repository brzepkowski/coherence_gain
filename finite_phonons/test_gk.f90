program main
    use coherence_gain_finite_phonons
    ! use coherence_gain
    implicit none
    real :: result, sum_finite_phonons, a

    result = sum_gk_squared()

    ! print *, "sum |gk|^2|, k={-inf, inf} (just integration) = ", result
    ! print *, "adjusted constant: ", adjusted_constant()
    a = adjusted_constant() * result
    print *, "sum |gk|^2|, k={-inf, inf} (with adj. const.) = ", a

    sum_finite_phonons = 0
    ! do j = 0, 25000
    !     k = real(j)/10000
    do j = 1, size(ks)
        k = ks(j)
        ! print *, "k = ", k, "|gk|^2 = (just integration)", sum_gk_squared_2(k)
        print *, k, " | ", adjusted_constant() * sum_gk_squared_2(k)
        ! print *, sum_gk_squared_2(k)
        ! sum_finite_phonons = sum_finite_phonons + (adjusted_constant() * sum_gk_squared_2(k) * 0.1282) ! Last value is the width of a rectangle / period od taken ks
        sum_finite_phonons = sum_finite_phonons + (adjusted_constant() * sum_gk_squared_2(k)) ! Last value is the width of a rectangle / period od taken ks
    end do
    print *, "sum finite phonons (with adj. const.) = ", sum_finite_phonons
    !
    alpha = a / sum_finite_phonons
    !
    print *, "alpha = ", alpha
    !
    print *, "calculate_alpha() = ", calculate_alpha()
    !
    print *, sum_finite_phonons * alpha
    !
    ! print *, W_one_part(20.0, 34.0)


end program main

program scratch
    use fasten_module
    implicit none

    real(fasten_rk) :: mu 
    real(fasten_rk), dimension(5) :: x_lag, y_lag
    integer :: ie, i

    do while (1 .eq. 1)
        write (*, fmt="(a)", advance='no') "Enter mass ratio >> "
        read(*, *, iostat=ie) mu 
        if (ie .ne. 0) then 
            exit
        end if
        call lagrange_points(mu, x_lag, y_lag)
        write (*, *) x_lag(1), x_lag(2), x_lag(3)
    end do 



end program scratch

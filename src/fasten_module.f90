module fasten_module
    use, intrinsic :: iso_fortran_env
    implicit none 

    private
    integer, parameter, public :: fasten_rk = real64
    integer, parameter :: wp = fasten_rk

    public :: mass_ratio, primary_distance, pseudopotential, jacobi_constant
    public :: hamiltonian, lagrange_points

    contains 

        !> Compute the three body mass ratio from primary masses/gms
        !! @param gm1 gravitational parameter or mass of the first primary
        !! @param gm2 gravitational parameter or mass of the second primary
        !! @note because this is a nondimensional quantity, units must only be consistent
        pure function mass_ratio(gm1, gm2) result(mu)
            real(wp), intent(in) :: gm1, gm2
            real(wp) :: mu 
            mu = gm2 / (gm1 + gm2)
        end function mass_ratio

        !> Compute the distance from a position to the two three body primaries 
        !! @param mu mass ratio of the system 
        !! @param x position/state from which to compute 
        !! @param r13 the distance to the first primary
        !! @param r23 the distance to the second primary
        pure subroutine primary_distance(mu, x, r13, r23)
            real(wp), intent(in) :: mu 
            real(wp), intent(in) :: x(:)
            real(wp), intent(out) :: r13
            real(wp), intent(out) :: r23

            r13 = x(2)**2 + x(3)**2
            r23 = r13 
            r13 = sqrt(r13 + (x(1) + 0 + mu)**2)
            r23 = sqrt(r23 + (x(1) - 1 + mu)**2)
        end subroutine primary_distance

        pure function pseudopotential(mu, x) result(omega)
            real(wp), intent(in) :: mu 
            real(wp), intent(in) :: x(:)
            real(wp) :: omega, r13, r23 

            call primary_distance(mu, x, r13, r23)
            omega = (1.0_wp - mu) / r13 + mu / r23 + 0.5_wp * (x(1)**2 + x(2)**2)
        end function pseudopotential

        pure function jacobi_constant(mu, x) result(jc)
            real(wp), intent(in) :: mu 
            real(wp), intent(in) :: x(:)
            real(wp) :: jc, v2

            v2 = x(4)**2 + x(5)**2 + x(6)**2
            jc = 2.0_wp * pseudopotential(mu, x) - v2
        end function jacobi_constant

        pure function hamiltonian(mu, x) result(h)
            real(wp), intent(in) :: mu 
            real(wp), intent(in) :: x(:)
            real(wp) :: h
            h = -jacobi_constant(mu, x) / 2.0_wp
        end function hamiltonian

        pure function horner(coeffs, x) result(v)
            real(wp), dimension(0:), intent(in) ::coeffs
            real(wp), intent(in) :: x
            real(wp) :: v 
            integer :: i

            v = coeffs(ubound(coeffs, 1))
            do i = ubound(coeffs, 1)-1,0,-1
                v = v * x + coeffs(i)
            end do 
        end function horner 

        pure subroutine szebehely(mu, x)
            real(wp), intent(in) :: mu 
            real(wp), intent(out) :: x(:)
            real(wp) :: eta, sigma
            real(wp), dimension(6) :: c1, c2, c3 

            c1 = [1.0_wp, -1.0_wp / 3.0_wp, -1.0_wp / 9.0_wp, -23.0_wp / 81.0_wp,  151.0_wp / 243.0_wp, -1.0_wp / 9.0_wp]
            c2 = [1.0_wp,  1.0_wp / 3.0_wp, -1.0_wp / 9.0_wp, -31.0_wp / 81.0_wp, -119.0_wp / 243.0_wp, -1.0_wp / 9.0_wp]
            c3 = [1.0_wp, 0.0_wp, 23.0_wp / 84.0_wp,  761.0_wp / 2352.0_wp, 3163.0_wp / 7056.0_wp, 30703.0_wp / 49392.0_wp]

            eta = (mu / 3.0_wp) ** (1.0_wp / 3.0_wp)
            sigma = 7.0_wp * mu / 12.0_wp

            x(1) =  1.0_wp - mu -   eta * horner(c1, eta)
            x(2) =  1.0_wp - mu +   eta * horner(c2, eta)
            x(3) = -1.0_wp - mu + sigma * horner(c3, sigma)
        end subroutine szebehely

        pure subroutine lp_fdf(fc, dfc, x, f, df)
            real(wp), intent(in) :: fc(:), dfc(:)
            real(wp), intent(in) :: x
            real(wp), intent(out) :: f, df

            f = horner(fc, x)
            df = horner(dfc, x)
        end subroutine

        pure function lp_newton_raphson(fc, dfc, x0, tol) result(x)
            real(wp), intent(in) :: fc(:), dfc(:)
            real(wp), intent(in) :: x0 
            real(wp), intent(in) :: tol
            real(wp) :: x, err, derr

            x = x0
            call lp_fdf(fc, dfc, x, err, derr)

            do while (abs(err) .gt. tol)
                x = x - err / derr
                call lp_fdf(fc, dfc, x, err, derr)
            end do 
        end function lp_newton_raphson

        pure subroutine lp_x_to_g(mu, x)
            real(wp), intent(in) :: mu 
            real(wp), intent(inout) :: x(:)

            x(1) = 1.0_wp - mu - x(1)
            x(2) = x(2) - 1.0_wp + mu
            x(3) = 0.0_wp - mu - x(3)
        end subroutine lp_x_to_g

        pure subroutine lp_g_to_x(mu, g)
            real(wp), intent(in) :: mu 
            real(wp), intent(inout) :: g(:)

            g(1) = 1.0_wp - mu - g(1)
            g(2) = 1.0_wp - mu + g(2)
            g(3) = 0.0_wp - mu - g(3)
        end subroutine 


        pure subroutine lagrange_points(mu, x, y)
            real(wp), intent(in) :: mu 
            real(wp), intent(out) :: x(:), y(:)
            real(wp), dimension(6) :: p1, p2, p3
            real(wp), dimension(5) :: dp1, dp2, dp3
            
            p1  = [-mu, 2.0_wp * mu, -mu, 3.0_wp - 2.0_wp * mu, mu - 3.0_wp, 1.0_wp]
            p2  = [-mu, -2.0_wp * mu, -mu, 3.0_wp - 2.0_wp * mu, 3.0_wp - mu, 1.0_wp] 
            p3  = [mu - 1.0_wp, 2.0_wp * mu -2.0_wp, mu - 1.0_wp, 2.0_wp * mu + 1.0_wp, mu + 2.0_wp, 1.0_wp]
            dp1 = [2.0_wp * mu, -2.0_wp * mu, 9.0_wp - 6.0_wp * mu, 4.0_wp * mu - 12.0_wp, 5.0_wp]
            dp2 = [-2.0_wp * mu, -2.0_wp * mu, 9.0_wp - 6.0_wp * mu, 12.0_wp - 4.0_wp * mu, 5.0_wp]
            dp3 = [ 2.0_wp * mu - 2.0_wp, 2.0_wp * mu - 2.0_wp, 6.0_wp * mu + 3.0_wp,  4.0_wp * mu + 8.0_wp, 5.0_wp]

            call szebehely(mu, x)

            call lp_x_to_g(mu, x)
            x(1) = lp_newton_raphson(p1, dp1, x(1), 1e-12_wp)
            x(2) = lp_newton_raphson(p2, dp2, x(2), 1e-12_wp)
            x(3) = lp_newton_raphson(p3, dp3, x(3), 1e-12_wp)
            call lp_g_to_x(mu, x)

            x(4) = 0.5_wp - mu 
            x(5) = x(4)
            y(1:3) = 0.0_wp
            y(4) = sqrt(3.0_wp) / 2.0_wp
            y(5) = -y(4)

        end subroutine lagrange_points


end module fasten_module

module U1_functions

  use iso_fortran_env, only : dp => real64, i4 => int32
  use parameters, only : Lx, Ly, Lz
  implicit none

contains

  function plaquette(u,x,mu,nu)
    complex(dp) :: plaquette
    complex(dp), dimension(3,0:Lx+1,0:Ly+1,0:Lz+1), intent(in) :: u
    integer(i4), intent(in) :: x(3), mu, nu
    integer(i4), dimension(3) :: x2, x3

    x2 = x
    x3 = x
    
    x2(mu) = x(mu) + 1
    x3(nu) = x(nu) + 1
        
    plaquette = U(mu,x(1),x(2),x(3)) * U(nu,x2(1),x2(2),x2(3)) * &
          conjg(U(mu,x3(1),x3(2),x3(3))) * conjg(U(nu,x(1),x(2),x(3)))
    
  end function plaquette

  function staples(u,x,mu)
    complex(dp) :: staples
    complex(dp), dimension(3,0:Lx+1,0:Ly+1,0:Lz+1), intent(in) :: u
    integer(i4), intent(in) :: x(3), mu
    integer(i4) :: nu
    integer(i4), dimension(3) ::  x2, x3, x4, x5, x6
    
    staples = 0.0_dp
    x3 = x
    x3(mu) = x(mu) + 1
    do nu = 1, 3
       if( nu == mu ) cycle
       x2 = x; x3 = x; x4 = x
       x2(nu) = x(nu) + 1
       x4(nu) = x(nu) - 1
       x6 = x3
       x6(nu) = x3(nu) - 1
       
       staples = staples + u(nu,x(1),x(2),x(3)) * u(mu,x2(1), x2(2), x2(3)) * conjg( u(nu,x3(1), x3(2), x3(3)) ) + &
         conjg( u(nu,x4(1),x4(2),x4(3)) ) * u(mu,x4(1), x4(2), x4(3)) * u(nu,x6(1), x6(2), x6(3))
    end do
    
  end function staples

  function plaquette_value(u)
    real(dp) :: plaquette_value
    complex(dp), dimension(3,0:Lx+1,0:Ly+1,0:Lz+1), intent(in) :: u
    !integer(i4), parameter :: d = 3
    integer(i4) :: x,y,z, mu, nu

    plaquette_value = 0.0_dp
    
    do x = 1, Lx
       do y = 1, Ly
          do z = 1, Lz
             do mu = 1, 2
                do nu = mu+1, 3
                   plaquette_value = plaquette_value + real(plaquette(u,[x,y,z],mu,nu))
                end do
             end do
          end do
       end do
    end do
    
  end function plaquette_value
  
end module U1_functions


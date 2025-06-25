#include "../input/include"
module U1_functions

  use iso_fortran_env, only : dp => real64, i4 => int32
#ifdef PARALLEL
  use parameters, only : Lx, Ly, Lz
#endif
#ifdef SERIAL
  use pbc
#endif
  
  implicit none

contains

  function plaquette(u,x,mu,nu)
    complex(dp) :: plaquette
#ifdef PARALLEL
    complex(dp), dimension(3,0:Lx+1,0:Ly+1,0:Lz+1), intent(in) :: u
#endif
#ifdef SERIAL
    complex(dp), dimension(:,:,:,:), intent(in) :: u
#endif
    integer(i4), intent(in) :: x(3), mu, nu
    integer(i4), dimension(3) :: x2, x3

#ifdef PARALLEL
    x2 = x
    x3 = x
    
    x2(mu) = x(mu) + 1
    x3(nu) = x(nu) + 1
#endif

#ifdef SERIAL
    x2 = ip(x,mu)
    x3 = ip(x,nu)
#endif

    
    plaquette = U(mu,x(1),x(2),x(3)) * U(nu,x2(1),x2(2),x2(3)) * &
          conjg(U(mu,x3(1),x3(2),x3(3))) * conjg(U(nu,x(1),x(2),x(3)))
    
  end function plaquette

  function staples(u,x,mu)
    complex(dp) :: staples
#ifdef PARALLEL
    complex(dp), dimension(3,0:Lx+1,0:Ly+1,0:Lz+1), intent(in) :: u
#endif
#ifdef SERIAL
    complex(dp), dimension(:,:,:,:), intent(in) :: u
#endif
    integer(i4), intent(in) :: x(3), mu
    integer(i4) :: nu
    integer(i4), dimension(3) ::  x2, x3, x4, x5, x6
    
    staples = 0.0_dp
#ifdef PARALLEL
    x3 = x
    x3(mu) = x(mu) + 1
#endif
#ifdef SERIAL
    x3 = ip(x,mu)
#endif
    do nu = 1, 3
       if( nu == mu ) cycle
#ifdef PARALLEL
       x2 = x; x3 = x; x4 = x
       x2(nu) = x(nu) + 1
       x4(nu) = x(nu) - 1
       x6 = x3
       x6(nu) = x3(nu) - 1
#endif
#ifdef SERIAL
       x2 = ip(x,nu)
       x4 = im(x,nu)
       x6 = im(x3,nu)
#endif
       
       staples = staples + u(nu,x(1),x(2),x(3)) * u(mu,x2(1), x2(2), x2(3)) * conjg( u(nu,x3(1), x3(2), x3(3)) ) + &
         conjg( u(nu,x4(1),x4(2),x4(3)) ) * u(mu,x4(1), x4(2), x4(3)) * u(nu,x6(1), x6(2), x6(3))
    end do
    
  end function staples

  function plaquette_value(u)
    real(dp) :: plaquette_value
#ifdef PARALLEL
    complex(dp), dimension(3,0:Lx+1,0:Ly+1,0:Lz+1), intent(in) :: u
#endif
#ifdef SERIAL
    complex(dp), dimension(:,:,:,:), intent(in) :: u
    integer(i4) :: Lx, Ly, Lz
#endif
    
    !integer(i4), parameter :: d = 3
    integer(i4) :: x,y,z, mu, nu
    
    plaquette_value = 0.0_dp

#ifdef SERIAL
    Lx = size(u(1,:,1,1))
    Ly = size(u(1,1,:,1))
    Lz = size(u(1,1,1,:))
#endif
    
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


#include "../input/include"
module U1_functions

  use iso_fortran_env, only : dp => real64, i4 => int32
  use pbc, only : ip,im
  use indices
  
  implicit none

contains

  function plaquette_value(u)
    real(dp) :: plaquette_value
#ifdef SERIAL
    complex(dp), dimension(:,:,:,:), intent(in) :: u
#endif
#ifdef PARALLEL
    complex(dp), dimension(:,:,:,:), intent(in) :: u[*]
#endif
    integer(i4) :: Lx, Ly, Lz
    integer(i4) :: x,y,z, mu, nu
    
    plaquette_value = 0.0_dp

    Lx = size(u(1,:,1,1))
    Ly = size(u(1,1,:,1))
    Lz = size(u(1,1,1,:))
    
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

#ifdef SERIAL
  function plaquette(u,x,mu,nu)
    complex(dp) :: plaquette
    complex(dp), dimension(:,:,:,:), intent(in) :: u
    integer(i4), intent(in) :: x(3), mu, nu
    integer(i4), dimension(3) :: x2, x3 
   
    x2 = ip(x,mu)!x(mu) + 1
    x3 = ip(x,nu)!x(nu) + 1

    plaquette = U(mu,x(1),x(2),x(3)) * U(nu,x2(1),x2(2),x2(3)) * &
          conjg(U(mu,x3(1),x3(2),x3(3))) * conjg(U(nu,x(1),x(2),x(3)))
    
  end function plaquette

  
  function staples(u,x,mu)
    complex(dp) :: staples

    complex(dp), dimension(:,:,:,:), intent(in) :: u
    integer(i4), intent(in) :: x(3), mu
    integer(i4) :: nu
    integer(i4), dimension(3) ::  x2, x3, x4, x5, x6
    
    staples = 0.0_dp
    x3 = ip(x,mu)
    do nu = 1, 3
       if( nu == mu ) cycle
       x2 = ip(x,nu)
       x4 = im(x,nu)
       x6 = im(x3,nu)
       
       staples = staples + u(nu,x(1),x(2),x(3)) * u(mu,x2(1),x2(2),x2(3)) * conjg( u(nu,x3(1),x3(2),x3(3)) ) + &
            conjg( u(nu,x4(1),x4(2),x4(3)) ) * u(mu,x4(1),x4(2),x4(3)) * u(nu,x6(1),x6(2),x6(3))
    end do
    
  end function staples
  
#endif

  
#ifdef PARALLEL
  function plaquette(u,x,mu,nu)
    complex(dp) :: plaquette
    complex(dp), dimension(:,:,:,:), intent(in) :: u[*]
    integer(i4), intent(in) :: x(3), mu, nu
    integer(i4), dimension(0:3) :: x2, x3 


    x2 = ip([this_image(),x],mu)
    x3 = ip([this_image(),x],nu)

    plaquette = U(mu,x(1),x(2),x(3)) * U(nu,x2(1),x2(2),x2(3))[x2(0)] * &
          conjg(U(mu,x3(1),x3(2),x3(3))[x3(0)]) * conjg(U(nu,x(1),x(2),x(3)))
    
  end function plaquette

  function staples(u,x,mu)
    use parameters, only : d
    complex(dp) :: staples
    complex(dp), dimension(:,:,:,:), intent(in) :: u[*]
    integer(i4), intent(in) :: x(3), mu
    integer(i4) :: nu

    integer(i4), dimension(0:3) :: x2, x3, x4, x5, x6
    
    staples = 0.0_dp
   
    x3 = ip([this_image(),x],mu)
    
    do nu = 1, 3
       if( nu == mu ) cycle
       x2 = ip([this_image(),x],nu)
       x4 = im([this_image(),x],nu)
       x6 = im(x3,nu)

       staples = staples + u(nu,x(1),x(2),x(3)) * u(mu,x2(1), x2(2), x2(3))[x2(0)] &
            * conjg( u(nu,x3(1), x3(2), x3(3))[x3(0)] ) + &
            conjg( u(nu,x4(1),x4(2),x4(3))[x4(0)] ) * u(mu,x4(1), x4(2), x4(3))[x4(0)] &
            * u(nu,x6(1), x6(2), x6(3))[x6(0)]
    end do
    
  end function staples

#endif

  
  
end module U1_functions


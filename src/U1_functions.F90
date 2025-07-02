module U1_functions

  use iso_fortran_env, only : dp => real64, i4 => int32
  use pbc, only : ip,im
  use constants, only : twopi, pi2
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
#ifdef SERIAL 
                   plaquette_value = plaquette_value + real(plaquette(u,[x,y,z],mu,nu))
#endif
#ifdef PARALLEL
                   plaquette_value = plaquette_value + real(plaquette(u,[this_image(),x,y,z],mu,nu))
#endif
                end do
             end do
          end do
       end do
    end do
    
  end function plaquette_value

  function plaquette(u,x,mu,nu)
    complex(dp) :: plaquette
#ifdef SERIAL
    complex(dp), dimension(:,:,:,:), intent(in) :: u
    integer(i4), intent(in) :: x(3), mu, nu
    integer(i4), dimension(3) :: x2, x3
#endif
#ifdef PARALLEL
    complex(dp), dimension(:,:,:,:), intent(in) :: u[*]
    integer(i4), intent(in) :: x(0:3), mu, nu
    integer(i4), dimension(0:3) :: x2, x3 
#endif

    x2 = ip(x,mu)
    x3 = ip(x,nu)

#ifdef SERIAL    
    plaquette = U(mu,x(1),x(2),x(3)) * U(nu,x2(1),x2(2),x2(3)) * &
          conjg(U(mu,x3(1),x3(2),x3(3))) * conjg(U(nu,x(1),x(2),x(3)))
#endif

#ifdef PARALLEL
    plaquette = U(mu,x(1),x(2),x(3))[x(0)] * U(nu,x2(1),x2(2),x2(3))[x2(0)] * &
          conjg(U(mu,x3(1),x3(2),x3(3))[x3(0)]) * conjg(U(nu,x(1),x(2),x(3))[x(0)])
#endif
    
  end function plaquette

  
  function staples(u,x,mu)
    complex(dp) :: staples
#ifdef SERIAL
    complex(dp), dimension(:,:,:,:), intent(in) :: u
    integer(i4), intent(in) :: x(3), mu
    integer(i4), dimension(3) ::  x2, x3, x4, x5, x6
#endif
#ifdef PARALLEL
    complex(dp), dimension(:,:,:,:), intent(in) :: u[*]
    integer(i4), intent(in) :: x(0:3), mu
    integer(i4), dimension(0:3) ::  x2, x3, x4, x5, x6
#endif
    integer(i4) :: nu

    staples = 0.0_dp
    x3 = ip(x,mu)
    do nu = 1, 3
       if( nu == mu ) cycle

       x2 = ip(x,nu)
       x4 = im(x,nu)       
       x6 = im(x3,nu)

#ifdef SERIAL
       staples = staples + u(nu,x(1),x(2),x(3)) * u(mu,x2(1),x2(2),x2(3)) * conjg( u(nu,x3(1),x3(2),x3(3)) ) + &
            conjg( u(nu,x4(1),x4(2),x4(3)) ) * u(mu,x4(1),x4(2),x4(3)) * u(nu,x6(1),x6(2),x6(3))
#endif
#ifdef PARALLEL
       staples = staples + u(nu,x(1),x(2),x(3))[x(0)] * u(mu,x2(1),x2(2),x2(3))[x2(0)] &
            * conjg( u(nu,x3(1),x3(2),x3(3))[x3(0)] ) + &
            conjg( u(nu,x4(1),x4(2),x4(3))[x4(0)] ) * u(mu,x4(1),x4(2),x4(3))[x4(0)] &
            * u(nu,x6(1), x6(2), x6(3))[x6(0)]
#endif
    end do
    
  end function staples

  function topological_charge_density(u)
    real(dp) :: topological_charge_density
#ifdef SERIAL
    complex(dp), dimension(:,:,:,:), intent(in) :: u
    integer(i4), dimension(3) :: point, x4, x5, x6
#endif
#ifdef PARALLEL
    complex(dp), dimension(:,:,:,:), intent(in) :: u[*]
    integer(i4), dimension(0:3) :: point, x4, x5, x6
#endif
    integer(i4) :: x, y, z, Lx, Ly, Lz
    complex(dp) :: plq1, plq2, plq3, plq4, plq5, plq6
    real(dp) :: theta1, theta2, theta3, theta4, theta5, theta6, theta_x

    Lx = size(u(1,:,1,1))
    Ly = size(u(1,1,:,1))
    Lz = size(u(1,1,1,:))
    
    
    topological_charge_density = 0.0_dp
    do x = 1, Lx
       do y = 1, Ly
          do z = 1, Lz
#ifdef SERIAL
             point = [x,y,z]
#endif
#ifdef PARALLEL
             point = [this_image(),x,y,z]
#endif
             plq1 = conjg(plaquette(u,point,1,2))
             plq2 = plaquette(u,point,1,3)
             plq3 = conjg(plaquette(u,point,2,3))

             x4 = ip(point,1)
             x5 = ip(point,2)
             x6 = ip(point,3)

             plq4 = plaquette(u,x4,2,3)
             plq5 = conjg(plaquette(u,x5,1,3))
             plq6 = plaquette(u,x6,1,2)

             theta1 = atan2(plq1%im,plq1%re)
             theta2 = atan2(plq2%im,plq2%re)
             theta3 = atan2(plq3%im,plq3%re)
             theta4 = atan2(plq4%im,plq4%re)
             theta5 = atan2(plq5%im,plq5%re)
             theta6 = atan2(plq6%im,plq6%re)

             theta_x = theta1 + theta2 + theta3 + theta4 + theta5 + theta6
             topological_charge_density = topological_charge_density + abs(theta_x)
                       
          end do
       end do
    end do
    
       
    
  end function topological_charge_density
  
  
end module U1_functions


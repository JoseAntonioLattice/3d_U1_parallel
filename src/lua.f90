module lua

  use iso_fortran_env, only : dp => real64, i4 => int32
  use parameters, only : Lx, Ly, Lz
  use U1_functions
  implicit none
  complex(dp), parameter :: i = (0.0_dp,1.0_dp)
  real(dp), parameter :: pi = acos(-1.0_dp)
  
contains

  subroutine metropolis(u,x,mu,beta)
    complex(dp), dimension(3,0:Lx+1,0:Ly+1,0:Lz+1), intent(inout) :: u
    integer(i4), intent(in) :: x(3), mu
    real(dp), intent(in) :: beta
    
    complex(dp) :: u_new
    real(dp) :: phi
    real(dp) :: p, deltaS
    real(dp) :: r
    
    call random_number(phi)
    phi = 2*pi*phi
    u_new = exp(i*phi)

    deltaS = DS(u(mu,x(1),x(2),x(3)),u_new,beta,staples(u,x,mu))

    call random_number(r)
    p = min(1.0_dp,exp(-DeltaS))
    if ( r <= p )then
       u(mu,x(1),x(2),x(3)) = u_new
    end if
    
  end subroutine metropolis

  function DS(uold, unew, beta,stp)
    real(dp) :: DS,beta
    complex(dp) :: uold, unew, stp

    DS = -beta * real( (unew - uold) * conjg(stp) )
    
  end function DS

  
end module lua

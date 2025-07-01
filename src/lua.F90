module lua

  use iso_fortran_env, only : dp => real64, i4 => int32
  use U1_functions
  use random
  use constants, only : pi, twopi, i => ii
  implicit none

  abstract interface
     subroutine lua_function(u,x,mu,beta)
       use iso_fortran_env, only : dp => real64, i4 => int32
       implicit none
#ifdef PARALLEL
       complex(dp), dimension(:,:,:,:), intent(inout) :: u[*]
#endif
#ifdef SERIAL
       complex(dp), dimension(:,:,:,:), intent(inout) :: u
#endif
       integer(i4), intent(in) :: x(3), mu
       real(dp), intent(in) :: beta      
     end subroutine lua_function
  end interface
  
contains
  
  subroutine metropolis(u,x,mu,beta)
#ifdef PARALLEL
    complex(dp), dimension(:,:,:,:), intent(inout) :: u[*]
#endif
#ifdef SERIAL
    complex(dp), dimension(:,:,:,:), intent(inout) :: u
#endif
    integer(i4), intent(in) :: x(3), mu
    real(dp), intent(in) :: beta
    
    complex(dp) :: u_new
    real(dp) :: phi
    real(dp) :: p, deltaS
    real(dp) :: r
    
    call random_number(phi)
    phi = twopi*phi
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

  subroutine heatbath(u,x,mu,beta)
#ifdef PARALLEL
    complex(dp), dimension(:,:,:,:), intent(inout) :: u[*]
#endif
#ifdef SERIAL
    complex(dp), dimension(:,:,:,:), intent(inout) :: u
#endif
    integer(i4), intent(in) :: x(3), mu
    real(dp), intent(in) :: beta

    complex(dp) :: v
    real(dp) :: beta_absv, gamma, phi_p, Z
    logical :: condition
    
    v = conjg(staples(u,x,mu))
    beta_absv = beta*abs(v)
    gamma = atan2(v%im,v%re)
    if(gamma < 0.0_dp) gamma = gamma + 2*pi

    condition = .false.
    do while(.not. condition)
       phi_p = uniform(0.0_dp,twopi)
       Z = uniform(0.0_dp,exp(beta_absv))
       if (Z <= exp(beta_absv*cos(phi_p+gamma))) condition = .true.
    end do

    u(mu,x(1),x(2),x(3)) = exp(i*phi_p)
    
  end subroutine heatbath

subroutine glauber(u,x,mu,beta)
#ifdef PARALLEL
    complex(dp), dimension(:,:,:,:), intent(inout) :: u[*]
#endif
#ifdef SERIAL
    complex(dp), dimension(:,:,:,:), intent(inout) :: u
#endif
    integer(i4), intent(in) :: x(3), mu
    real(dp), intent(in) :: beta
    
    complex(dp) :: u_new
    real(dp) :: phi
    real(dp) :: p, deltaS
    real(dp) :: r
    
    call random_number(phi)
    phi = twopi*phi
    u_new = exp(i*phi)

    deltaS = DS(u(mu,x(1),x(2),x(3)),u_new,beta,staples(u,x,mu))

    call random_number(r)
    p = 1/(1.0_dp + exp(deltaS))
    if ( r <= p )then
       u(mu,x(1),x(2),x(3)) = u_new
    end if
    
  end subroutine glauber
  
  subroutine choose_algorithm(algorithm,u,x,mu,beta)
#ifdef PARALLEL
    complex(dp), dimension(:,:,:,:), intent(inout) :: u[*]
#endif
#ifdef SERIAL
    complex(dp), dimension(:,:,:,:), intent(inout) :: u
#endif
    integer(i4), intent(in) :: x(3), mu
    real(dp), intent(in) :: beta
    character(*) :: algorithm

    select case(algorithm)
    case('glauber')
       call glauber(u,x,mu,beta)
    case('heatbath')
       call heatbath(u,x,mu,beta)
    case('metropolis')
       call metropolis(u,x,mu,beta)
    end select
    
  end subroutine choose_algorithm
  
end module lua

module parameters

  use iso_fortran_env, only : dp => real64, i4 => int32

  implicit none

#ifdef PARALLEL
  integer(i4) :: d(3)
#endif
  logical :: inCluster
  integer(i4) :: L(3)
  integer(i4) :: N_thermalization
  integer(i4) :: N_measurements
  integer(i4) :: N_skip
  logical :: isbeta
  real(dp) :: beta_i, beta_f
  integer(i4) :: n_beta
  character(30) :: algorithm
  integer(i4) :: Nhmc
  real(dp) :: Thmc
  character(10) :: start
  logical :: equilibrium
  integer(i4) :: tau_Q
  integer(i4) :: inunit,outunit
  character(99) :: inputfilename, outputfilename
    
  
  namelist /parametersfile/ inCluster,L,N_thermalization,N_measurements,N_skip, &
       isbeta, beta_i, beta_f, n_beta, algorithm,Nhmc,Thmc, start, equilibrium, tau_Q
contains

  subroutine read_input()

#ifdef PARALLEL
    if(this_image() == 1) then
       read(*,*) d
       print*,"Enter cores array: ", d
#endif
       read(*,'(a)') inputfilename
       print*, 'Enter input parameters file: ', trim(inputfilename)
       open(newunit = inunit,file = trim(inputfilename), status = 'old', action = "read")
       read(inunit, nml = parametersfile)

       if(InCluster) then
          read(*,'(a)') outputfilename
          print*, 'Enter output file: ', trim(outputfilename)
       endif

       write(*,nml = parametersfile)

       if( any(L<=0) ) stop "All elements in L must be > 0"
       if( N_measurements <= 1 ) stop "N_measurements must be > 1"
       if( N_thermalization <= 0 ) stop "N_thermalization must be > 0"
       if( N_skip <= 0 ) stop "N_skip must be > 0"
       if( beta_i < 0.0_dp ) stop "beta_i must be >= 0"
       if( beta_f < 0.0_dp ) stop "beta_f must be >= 0"
       if( N_beta <= 1 ) stop "N_beta must be > 1"
       select case(algorithm)
       case("metropolis","glauber","heatbath","hmc")
       case default
          stop "algorithm must be 'metropolis', 'glauber', 'heatbath' or 'hmc'"
       end select
       if( Nhmc < 3 ) stop "Nhmc mut be > 2"
       if( Thmc <= 0.0_dp ) stop "Thmc mut be > 0"
       if( tau_Q <= 0 ) stop "tau_Q mut be > 0"
       
#ifdef PARALLEL    
    end if

    call co_broadcast(d,source_image = 1)
    call co_broadcast(L,source_image = 1)
    call co_broadcast(N_measurements,source_image = 1)
    call co_broadcast(N_thermalization,source_image = 1)
    call co_broadcast(N_skip,source_image = 1)
    call co_broadcast(isbeta,source_image = 1)
    call co_broadcast(N_beta,source_image = 1)
    call co_broadcast(beta_i, source_image = 1)
    call co_broadcast(beta_f, source_image = 1)
    call co_broadcast(algorithm, source_image = 1)
    call co_broadcast(start, source_image = 1)
    call co_broadcast(equilibrium, source_image = 1)
    call co_broadcast(tau_Q, source_image = 1)
    call co_broadcast(Nhmc, source_image = 1)
    call co_broadcast(Thmc, source_image = 1)

#endif
  end subroutine read_input

  
end module parameters

module parameters

  use iso_fortran_env, only : dp => real64, i4 => int32

  implicit none

#ifdef PARALLEL
  integer(i4) :: d(3)
#endif 
  integer(i4) :: L(3)
  integer(i4) :: N_thermalization
  integer(i4) :: N_measurements
  integer(i4) :: N_skip
  real(dp) :: beta_i, beta_f
  integer(i4) :: n_beta
  character(30) :: algorithm
  character(10) :: start
  logical :: equilibrium
  integer(i4) :: inunit,outunit
  character(99) :: inputfilename, outputfilename
    
  
  namelist /parametersfile/ L,N_thermalization,N_measurements,N_skip, &
       beta_i, beta_f, n_beta, algorithm, start, equilibrium
contains

  subroutine read_input()

#ifdef PARALLEL
    if(this_image() == 1) then
       read(*,*) d
       print*,"Enter cores array: ", d
#endif
       read(*,'(a)') inputfilename
       print*, 'Enter input parameters file: ', trim(inputfilename)
       read(*,'(a)') outputfilename
       print*, 'Enter output file: ', trim(outputfilename)
       open(newunit = inunit,file = trim(inputfilename), status = 'old', action = "read")
       read(inunit, nml = parametersfile)
       write(*,nml = parametersfile)
#ifdef PARALLEL    
    end if

    call co_broadcast(d,source_image = 1)
    call co_broadcast(L,source_image = 1)
    call co_broadcast(N_measurements,source_image = 1)
    call co_broadcast(N_thermalization,source_image = 1)
    call co_broadcast(N_skip,source_image = 1)
    call co_broadcast(N_beta,source_image = 1)
    call co_broadcast(beta_i, source_image = 1)
    call co_broadcast(beta_f, source_image = 1)
    call co_broadcast(algorithm, source_image = 1)
    call co_broadcast(start, source_image = 1)
    call co_broadcast(equilibrium, source_image = 1)
#endif
  end subroutine read_input

  
end module parameters

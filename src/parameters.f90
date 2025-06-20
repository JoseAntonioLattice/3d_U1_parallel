module parameters

  use iso_fortran_env, only : dp => real64, i4 => int32

  implicit none

  integer(i4) :: d(3)
  integer(i4) :: L(3), Lx, Ly, Lz
  integer(i4) :: N_thermalization
  integer(i4) :: N_measurements
  integer(i4) :: N_skip
  real(dp) :: beta_i, beta_f
  integer(i4) :: n_beta

  integer(i4) :: inunit,outunit
  character(99) :: inputfilename, outputfilename
    
  
  namelist /parametersfile/ L,N_thermalization,N_measurements,N_skip, &
       beta_i, beta_f, n_beta
contains

  subroutine read_input()

    if(this_image() == 1) then
       read(*,*) d
       print*,"Enter cores array: ", d
       read(*,'(a)') inputfilename
       print*, 'Enter input parameters file: ', trim(inputfilename)
       read(*,'(a)') outputfilename
       print*, 'Enter output file: ', trim(outputfilename)
       open(newunit = inunit,file = trim(inputfilename), status = 'old', action = "read")
       read(inunit, nml = parametersfile)
       write(*,nml = parametersfile)
    end if

    call co_broadcast(d,source_image = 1)
    call co_broadcast(L,source_image = 1)
    call co_broadcast(N_measurements,source_image = 1)
    call co_broadcast(N_thermalization,source_image = 1)
    call co_broadcast(N_skip,source_image = 1)
    call co_broadcast(N_beta,source_image = 1)
    call co_broadcast(beta_i, source_image = 1)
    call co_broadcast(beta_f, source_image = 1)
    
  end subroutine read_input

  
end module parameters

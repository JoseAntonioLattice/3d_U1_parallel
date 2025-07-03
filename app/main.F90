program U1_3d

  use arrays
  use parameters
  use dynamics
  use files
  use number2string
  implicit none

  integer :: ib
  character(:), allocatable :: filename
  character(100), allocatable :: directories(:)
#ifdef PARALLEL
  if(this_image() == 1 )print*, "3d U(1). Parallel"
#endif

#ifdef SERIAL
  print*, "3d U(1). Serial"
#endif

 
  call read_input()
  call set_memory(u,beta,N_measurements,N_beta,beta_i,beta_f,equilibrium,tau_Q)
  if(equilibrium)then
     directories = [character(100) :: "data","equilibrium",algorithm,"L="//int2str(L(1)) ]
  else
     directories = [character(100) :: "data","out_equilibrium",algorithm,"L="//int2str(L(1)),"tau_Q="//int2str(tau_Q) ]
  end if
#ifdef PARALLEL
  if( this_image() == 1)then
     call create_files(directories,filename)
     print*, "FILENAME: ", filename
  end if
  sync all
#endif
  
#ifdef SERIAL
  call create_files(directories,filename)
  print*, filename
#endif
  

#ifdef PARALLEL
  if(this_image() == 1) then
     open(newunit = outunit, file = filename, status = "unknown", action = "write")
  end if
#endif
#ifdef SERIAL
  open(newunit = outunit, file = filename, status = "unknown", action = "write")
#endif
  if(equilibrium) then
     call eq(start,algorithm,u,beta,N_thermalization,N_skip,N_measurements,outunit)
  else
     call out_eq(start,algorithm,u,beta,tau_Q, N_thermalization,N_measurements,outunit)
  end if
end program U1_3d

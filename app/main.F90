program U1_3d

  use arrays
  use parameters
  use dynamics
  use files
  use number2string
  implicit none

  integer :: ib
  integer(8) :: ti, tf, rate, max
  character(:), allocatable :: filename
  character(100), allocatable :: directories(:)

  call random_init(.false.,.true.)
  
#ifdef PARALLEL
  if(this_image() == 1 ) print*, "3d U(1). Parallel"
#elif SERIAL
  print*, "3d U(1). Serial"
#endif
 
  call read_input()
  call set_memory(u,beta,N_measurements,N_beta,beta_i,beta_f,equilibrium,tau_Q)
  

#ifdef PARALLEL
  if(this_image() == 1) then
#endif
     directories = [character(100) :: "data","equilibrium",algorithm, &
          "Lx="//int2str(L(1)), "Ly="//int2str(L(2)), "Lz="//int2str(L(3))]
     if(.not. equilibrium) then
        directories(2) = "out_equilibrium"
        directories = [character(100) :: directories,"tau_Q="//int2str(tau_Q) ]
     end if
     
     if(inCluster) filename = outputfilename
     call create_files(directories,filename,inCluster)
     print*, "FILENAME: ", filename
     open(newunit = outunit, file = filename, status = "unknown", action = "write")
#ifdef PARALLEL
  end if
  sync all
#endif
  CALL SYSTEM_CLOCK(COUNT_RATE=rate, COUNT_MAX=max)
  CALL SYSTEM_CLOCK(COUNT=ti)
  if(equilibrium) then
     call eq(start,algorithm,u,beta,N_thermalization,N_skip,N_measurements,outunit,isbeta)
  else
     call out_eq(start,algorithm,u,beta,tau_Q, N_thermalization,N_measurements,outunit,isbeta)
  end if
  CALL SYSTEM_CLOCK(COUNT=tf)
  
#ifdef PARALLEL
  if(this_image() == 1) then
#endif
     print*, "elapsed time: ", real(tf - ti)/rate
     write(outunit,"(3/,a,f15.8)") "elapsed time: ", real(tf - ti)/rate
#ifdef PARALLEL
     write(outunit,*) "threads: ", num_images()
     write(outunit,*) "d = ", d
#endif
     close(outunit)
#ifdef PARALLEL
  end if
#endif
end program U1_3d

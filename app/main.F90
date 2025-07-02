program U1_3d

  use arrays
  use parameters
  use dynamics
  use files
  use number2string
  use constants, only: twopi
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
  call set_memory(u,plq,top_den,beta,N_measurements,N_beta,beta_i,beta_f)
  directories = [character(100) :: "data",algorithm,"L="//int2str(L(1)) ]
#ifdef PARALLEL
  if( this_image() == 1)then
     call create_files(directories,filename)
     print*, "FILENAME: ", filename
  end if
#endif
  
#ifdef SERIAL
  call create_files(directories,filename)
  print*, filename
#endif
  
  select case(start)
  case('cold')
     call cold_start(u)
  case("hot")
     call hot_start(u)
  end select
#ifdef PARALLEL
  if(this_image() == 1) then
     open(newunit = outunit, file = filename, status = "unknown", action = "write")
  end if
#endif
#ifdef SERIAL
  open(newunit = outunit, file = filename, status = "unknown", action = "write")
#endif
  do ib = 1, n_beta
     call thermalization(trim(algorithm),u,1/beta(ib),N_thermalization)
     call measurements(trim(algorithm),u,1/beta(ib),N_measurements,N_skip,plq,top_den)
#ifdef PARALLEL
     if(this_image() == 1)then
        print*, beta(ib), (sum(plq)/size(plq))/(3*product(L)), sum(top_den)/(size(top_den)*twopi*product(L))
        write(outunit,*) beta(ib), (sum(plq)/size(plq))/(3*product(L)), sum(top_den)/(size(top_den)*twopi*product(L))
        flush(outunit)
     end if
#endif
#ifdef SERIAL
        print*, beta(ib), (sum(plq)/size(plq))/(3*product(L)), sum(top_den)/(size(top_den)*twopi*product(L))
        write(outunit,*) beta(ib), (sum(plq)/size(plq))/(3*product(L)), sum(top_den)/(size(top_den)*twopi*product(L))
        flush(outunit)
#endif
  end do
end program U1_3d

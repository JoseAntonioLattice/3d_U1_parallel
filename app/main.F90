program U1_3d

  use arrays
  use parameters
  use dynamics

  implicit none

  integer :: ib

#ifdef PARALLEL
  if(this_image() == 1 )print*, "3d U(1). Parallel"
#endif

#ifdef SERIAL
  print*, "3d U(1). Serial"
#endif

  call read_input()
  call set_memory(u,plq,beta,N_measurements,N_beta,beta_i,beta_f)

  select case(start)
  case('cold')
     call cold_start(u)
  case("hot")
     call hot_start(u)
  end select
  
  do ib = 1, n_beta
     call thermalization(trim(algorithm),u,1/beta(ib),N_thermalization)
     call measurements(trim(algorithm),u,1/beta(ib),N_measurements,N_skip,plq)
#ifdef PARALLEL
     if(this_image() == 1) print*, beta(ib), (sum(plq)/size(plq))/(3*product(L))
#endif
#ifdef SERIAL
     print*, beta(ib), (sum(plq)/size(plq))/(3*product(L))
#endif
  end do
end program U1_3d

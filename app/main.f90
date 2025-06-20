program U1_3d_parallel

  use arrays
  use parameters
  use dynamics

  implicit none

  integer :: ib
  
  if(this_image() == 1) print*, "3d U(1) model"

  call read_input()
  call set_memory(u,plq,beta,N_measurements,N_beta,beta_i,beta_f)
  call cold_start(u)
  do ib = 1, n_beta
     call thermalization(u,1/beta(ib),N_thermalization)
     call measurements(u,1/beta(ib),N_measurements,N_skip,plq)
     if(this_image() == 1) print*, beta(ib), (sum(plq)/size(plq))/(3*product(L))
  end do
end program U1_3d_parallel

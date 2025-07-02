module dynamics
  use iso_fortran_env, only : dp => real64, i4 => int32
  use lua
  use constants, only : twopi
  implicit none

contains

  subroutine cold_start(u)
    complex(dp), intent(out), dimension(:,:,:,:) :: u
    u = 1.0_dp
  end subroutine cold_start

  subroutine hot_start(u)
    use constants, only : ii, twopi
    complex(dp), intent(out), dimension(:,:,:,:) :: u
    real(dp), dimension(size(u(:,1,1,1)), size(u(1,:,1,1)), size(u(1,1,:,1)), size(u(1,1,1,:))) :: phi
    u = exp(ii*twopi*phi)
  end subroutine hot_start
  
#ifdef PARALLEL  
  subroutine set_memory(u,plq,top_den,beta,N_measurements,Nbeta,betai,betaf)
    use indices
    use pbc
    use parameters, only : L, d
    complex(dp), intent(inout), allocatable, dimension(:,:,:,:) :: u[:]
    real(dp), intent(inout), allocatable, dimension(:) :: plq[:],top_den[:], beta
    integer(i4), intent(in) :: N_measurements, Nbeta
    real(dp) :: betai, betaf
    
    integer(i4) :: i, Lx, Ly, Lz
    integer(i4), dimension(3) :: a
    integer(i4) :: left, right, up, down, front, back 
   
    integer(i4) :: xc(3,2)

    Lx = L(1)/d(1)
    Ly = L(2)/d(2)
    Lz = L(3)/d(3)

    if(this_image() == 1) print*, Lx, Ly, Lz, num_images()

    allocate(u(3,Lx,Ly,Lz)[*])
    allocate(plq(N_measurements)[*])
    allocate(top_den(N_measurements)[*])
    allocate(beta(Nbeta))

    beta = [(betai + (i-1)*(betaf-betai)/(Nbeta-1), i=1, Nbeta)]
    
   
    a = get_index_array(this_image(),3,d)

    ip1_c = [(mod(i,d(1))+1, i=1,d(1))]; im1_c = [(i-1, i=1,d(1))]; im1_c(1) = d(1)
    ip2_c = [(mod(i,d(2))+1, i=1,d(2))]; im2_c = [(i-1, i=1,d(2))]; im2_c(1) = d(2)
    ip3_c = [(mod(i,d(3))+1, i=1,d(3))]; im3_c = [(i-1, i=1,d(3))]; im3_c(1) = d(3)

    left  = get_index([im1_c(a(1)),a(2),a(3)],3,d)
    right = get_index([ip1_c(a(1)),a(2),a(3)],3,d)

    front = get_index([a(1),im2_c(a(2)),a(3)],3,d)
    back  = get_index([a(1),ip2_c(a(2)),a(3)],3,d)

    up   = get_index([a(1),a(2),im3_c(a(3))],3,d)
    down = get_index([a(1),a(2),ip3_c(a(3))],3,d)

    if(this_image() == 1) print*, "Before initialize_pbc"
    call initialize_pbc([Lx,Ly,Lz],d,this_image(),[left,front,down], [right,back,up])
    if(this_image() == 1) print*, "After initialize_pbc"
    sync all
    
  end subroutine set_memory
#endif

#ifdef SERIAL
  subroutine set_memory(u,plq,top_den,beta,N_measurements,Nbeta,betai,betaf,equilibrium,tau_Q)
    use pbc
    use parameters, only : L

    complex(dp), intent(inout), allocatable, dimension(:,:,:,:) :: u
    real(dp), intent(inout), allocatable, dimension(:) :: plq,top_den, beta
    integer(i4), intent(in) :: N_measurements, Nbeta, tau_Q
    logical :: equilibrium
    real(dp) :: betai, betaf
    
    integer(i4) :: i

    allocate(u(3,L(1),L(2),L(3)))
    allocate(plq(N_measurements))
    allocate(top_den(N_measurements))
    
    if(equilibrium)then
       allocate(beta(Nbeta))
       beta = [(betai + (i-1)*(betaf-betai)/(Nbeta-1), i=1, Nbeta)]
    else
       allocate(beta(-tau_Q:tau_Q))
       beta = 0.5*[(betai +betaf + i*(betaf-betai)/tau_Q, i=-tau_Q, tau_Q)]
    end if
    call initialize_pbc(L)
    
  end subroutine set_memory
#endif

  subroutine thermalization(algorithm,u,beta,N_thermalization)
    character(*), intent(in) :: algorithm
#ifdef PARALLEL
    complex(dp), intent(inout) :: u(:,:,:,:)[*]
#endif
#ifdef SERIAL
    complex(dp), intent(inout) :: u(:,:,:,:)
#endif
    real(dp), intent(in) :: beta
    integer(i4) :: N_thermalization
    integer(i4) :: i_sweeps

    do i_sweeps = 1, N_thermalization
       call sweeps(trim(algorithm),u,beta)
    end do

  end subroutine thermalization

  subroutine measurements(algorithm,u,beta,N_measurements,Nskip,plq,top_den)
    use U1_functions
    character(*), intent(in) :: algorithm
#ifdef PARALLEL
    complex(dp), intent(inout) :: u(:,:,:,:)[*]
    real(dp) :: plq(:)[*], top_den(:)[*]
#endif
#ifdef SERIAL
    complex(dp), intent(inout) :: u(:,:,:,:)
    real(dp) :: plq(:), top_den(:)
#endif
    real(dp), intent(in) :: beta
    integer(i4),intent(in) :: N_measurements, Nskip
    integer(i4) :: i_sweeps, iskip

    do i_sweeps = 1, N_measurements
       do iskip = 1, Nskip
          call sweeps(trim(algorithm),u,beta)
       end do
       plq(i_sweeps) = plaquette_value(u)
       top_den(i_sweeps) = topological_charge_density(u)
       
#ifdef PARALLEL
       sync all
       call co_sum(plq(i_sweeps),result_image = 1)
       call co_sum(top_den(i_sweeps),result_image = 1)
       !if(this_image() == 1) print*, top_den(i_sweeps)/twopi
#endif
    end do
    
  end subroutine measurements
  
  subroutine sweeps(algorithm,u,beta)
    character(*), intent(in) :: algorithm
#ifdef PARALLEL
    complex(dp), intent(inout) :: u(:,:,:,:)[*]
#endif
#ifdef SERIAL
    complex(dp), intent(inout) :: u(:,:,:,:)
#endif
    real(dp), intent(in) :: beta
    integer(i4) :: x,y,z,mu
    integer(i4) :: Lx, Ly, Lz
    
    Lx = size(u(1,:,1,1))
    Ly = size(u(1,1,:,1))
    Lz = size(u(1,1,1,:))


#ifdef PARALLEL
    if(this_image() == 1) print*, Lx, Ly, Lz
#endif
    
    do x = 1, Lx
       do y = 1, Ly
          do z = 1, Lz
             do mu = 1, 3
                call choose_algorithm(trim(algorithm),u,[x,y,z],mu,beta)
             end do
          end do
       end do
    end do   

#ifdef PARALLEL
    sync all
#endif
    
  end subroutine sweeps

  
  subroutine eq(algorithm,u,beta, N_thermalization,Nskip,N_measurements,outunit)
    use parameters, only : L
    use statistics
    character(*), intent(in) :: algorithm
#ifdef PARALLEL
    complex(dp), intent(inout) :: u(:,:,:,:)[*]
#endif
#ifdef SERIAL
    complex(dp), intent(inout) :: u(:,:,:,:)
#endif
    integer(i4), intent(in) :: N_thermalization, Nskip, N_measurements, outunit
    real(dp), intent(in) :: beta(:)
    integer(i4) :: ib, i_sweeps
#ifdef PARALLEL    
    real(dp), dimension(N_measurements), codimension[*] :: plq, top_den
#endif
#ifdef SERIAL    
    real(dp), dimension(N_measurements) :: plq, top_den
#endif
    
    do ib = 1, size(beta)
       call thermalization(trim(algorithm),u,1/beta(ib),N_thermalization)
       call measurements(trim(algorithm),u,1/beta(ib),N_measurements,Nskip,plq,top_den)
#ifdef PARALLEL
       if(this_image() == 1)then
          print*, beta(ib), avr(plq)/(3*product(L)), avr(top_den)/(twopi*product(L))
          write(outunit,*) beta(ib), avr(plq)/(3*product(L)), avr(top_den)/(twopi*product(L))
          flush(outunit)
       end if
#endif
#ifdef SERIAL
       print*, beta(ib), (sum(plq)/size(plq))/(3*product(L)), sum(top_den)/(size(top_den)*twopi*product(L))
       write(outunit,*) beta(ib), (sum(plq)/size(plq))/(3*product(L)), sum(top_den)/(size(top_den)*twopi*product(L))
       flush(outunit)
#endif
    end do
    
  end subroutine eq
  
  subroutine out_eq(algorithm,u,beta, tau_Q, N_thermalization, N_measurements,outunit)
    use parameters, only : L
    use statistics
    character(*), intent(in) :: algorithm
#ifdef PARALLEL
    complex(dp), intent(inout) :: u(:,:,:,:)[*]
#endif
#ifdef SERIAL
    complex(dp), intent(inout) :: u(:,:,:,:)
#endif
    integer(i4), intent(in) :: tau_Q, N_thermalization, N_measurements,outunit
    real(dp), intent(in) :: beta(-tau_Q:tau_Q)
    integer(i4) :: ib, i_sweeps
    real(dp), dimension(-tau_Q:tau_Q,N_measurements) :: plq, top_den
    real(dp) :: err_plq(2), err_top(2)
    do i_sweeps = 1, N_measurements
       call thermalization(algorithm,u,1/beta(-tau_Q),N_thermalization)
       do ib = -tau_Q, tau_Q
          call sweeps(algorithm,u,1/beta(ib))
          plq(ib,i_sweeps) = plaquette_value(u)
          top_den(ib,i_sweeps) = topological_charge_density(u)
#ifdef PARALLEL
          sync all
          call co_sum(plq(ib,i_sweeps),result_image = 1)
          call co_sum(top_den(ib,i_sweeps),result_image = 1)
#endif
       end do
    end do

    do ib = -tau_Q, tau_Q
#ifdef SERIAL
       err_plq = jackknife_max(plq(ib,:))/(3*product(L))
       err_top = jackknife_max(top_den(ib,:))/(twopi*product(L))
       print*, ib, avr(plq(ib,:))/(3*product(L)), err_plq(1), &
            avr(top_den(ib,:))/(twopi*product(L)), err_top(1)
       write(outunit,*) ib, avr(plq(ib,:))/(3*product(L)), err_plq(1), &
            avr(top_den(ib,:))/(twopi*product(L)), err_top(1)
#endif
    end do
  end subroutine out_eq
  
end module dynamics

module dynamics
  use iso_fortran_env, only : dp => real64, i4 => int32
  use lua
  use constants, only : twopi
#if defined(SERIAL) || PARALLEL == 1
  use hybridMC, only : hmc
#endif
  implicit none

  integer(i4) :: left, right, up, down, front, back, &
       left_up, left_down, left_front, left_back, &
       right_up, right_down, right_front, right_back, &
       back_up, back_down, front_up, front_down, &
       left_front_up, left_front_down, &
       left_back_up, left_back_down, &
       right_front_up, right_front_down, &
       right_back_up, right_back_down, &
       Lx, Ly, Lz
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

  subroutine select_start(u,start)
    complex(dp), intent(out), dimension(:,:,:,:) :: u
    character(*), intent(in) :: start
    
    select case(start)
    case('cold')
       call cold_start(u)
    case("hot")
       call hot_start(u)
    end select
    
  end subroutine select_start
  
#ifdef PARALLEL  
  subroutine set_memory(u,beta,N_measurements,Nbeta,betai,betaf,equilibrium,tau_Q)
    use indices
    use pbc
    use parameters, only : L, d
    complex(dp), intent(inout), allocatable, dimension(:,:,:,:) :: u[:]
    real(dp), intent(inout), allocatable, dimension(:) :: beta
    integer(i4), intent(in) :: N_measurements, Nbeta, tau_Q
    logical, intent(in) :: equilibrium
    real(dp) :: betai, betaf
    
    integer(i4) :: i
    integer(i4), dimension(3) :: a
     
   
    integer(i4) :: xc(3,2)

    Lx = L(1)/d(1)
    Ly = L(2)/d(2)
    Lz = L(3)/d(3)

    if(this_image() == 1) print*, Lx, Ly, Lz, num_images()
#if PARALLEL == 1
    allocate(u(3,Lx,Ly,Lz)[*])
#elif PARALLEL == 2
    allocate(u(3,0:Lx+1,0:Ly+1,0:Lz+1)[*])
#endif

    if(equilibrium)then
       allocate(beta(Nbeta))
       beta = [(betai + (i-1)*(betaf-betai)/(Nbeta-1), i=1, Nbeta)]
    else
       allocate(beta(-tau_Q:tau_Q))
       beta = 0.5*[(betai +betaf + i*(betaf-betai)/tau_Q, i=-tau_Q, tau_Q)]
    end if

       
    a = get_index_array(this_image(),d)

    allocate(ip1_c(d(1)),im1_c(d(1)))
    allocate(ip2_c(d(2)),im2_c(d(2)))
    allocate(ip3_c(d(3)),im3_c(d(3)))
    
    ip1_c = [(mod(i,d(1))+1,i=1,d(1))]; im1_c =[(mod(i+d(1)-2,d(1))+1,i=1,d(1))]
    ip2_c = [(mod(i,d(2))+1,i=1,d(2))]; im2_c =[(mod(i+d(2)-2,d(2))+1,i=1,d(2))]
    ip3_c = [(mod(i,d(3))+1,i=1,d(3))]; im3_c =[(mod(i+d(3)-2,d(3))+1,i=1,d(3))]
           
    left  = get_index(im_core(a,1),d)
    right = get_index(ip_core(a,1),d)

    front = get_index(im_core(a,2),d)
    back  = get_index(ip_core(a,2),d)

    up   = get_index(im_core(a,3),d)
    down = get_index(ip_core(a,3),d)
!#if DEBUG == 1
    do i = 1, 3
       print"('image = ',i0,x,', im_',i0,'[',2(i1,','),i1,'] = [',2(i1,','),i1,']' &
            x,', ip_',i0,'[',2(i1,','),i1,'] = [',2(i1,','),i1,']')"&
            ,this_image(), i,a,im_core(a,i),i,a,ip_core(a,i)
    end do
!#endif
    !print*, this_image(),left,right,front,back,up, down
#if PARALLEL == 2
    left_front = get_index(im_core(im_core(a,2),1),d)
    left_down  = get_index(im_core(im_core(a,3),1),d)
    left_back  = get_index(im_core(ip_core(a,2),1),d)
    left_up    = get_index(im_core(ip_core(a,3),1),d)

    right_front = get_index(ip_core(im_core(a,2),1),d)
    right_down  = get_index(ip_core(im_core(a,3),1),d)
    right_back  = get_index(ip_core(ip_core(a,2),1),d)
    right_up    = get_index(ip_core(ip_core(a,3),1),d)

    front_down = get_index(im_core(im_core(a,3),2),d)
    front_up   = get_index(im_core(ip_core(a,3),2),d)
    
    back_down = get_index(ip_core(im_core(a,3),2),d)
    back_up   = get_index(ip_core(ip_core(a,3),2),d)

    left_front_down = get_index(im_core(im_core(im_core(a,3),2),1),d)
    left_front_up   = get_index(im_core(im_core(ip_core(a,3),2),1),d)
    left_back_down  = get_index(im_core(ip_core(im_core(a,3),2),1),d)
    left_back_up    = get_index(im_core(ip_core(ip_core(a,3),2),1),d)

    right_front_down = get_index(ip_core(im_core(im_core(a,3),2),1),d)
    right_front_up   = get_index(ip_core(im_core(ip_core(a,3),2),1),d)
    right_back_down  = get_index(ip_core(ip_core(im_core(a,3),2),1),d)
    right_back_up    = get_index(ip_core(ip_core(ip_core(a,3),2),1),d)

!#if DEBUG == 1
    print'(*(i1,x))', left, right, up, down, front, back, &
         left_up, left_down, left_front, left_back, &
         right_up, right_down, right_front, right_back, &
         back_up, back_down, front_up, front_down, &
         left_front_up, left_front_down, &
         left_back_up, left_back_down, &
         right_front_up, right_front_down, &
         right_back_up, right_back_down
!#endif
#endif

    call initialize_pbc([Lx,Ly,Lz],d,this_image(),[left,front,down], [right,back,up])
    sync all
    
  end subroutine set_memory
#endif

#ifdef SERIAL
  subroutine set_memory(u,beta,N_measurements,Nbeta,betai,betaf,equilibrium,tau_Q)
    use pbc
    use parameters, only : L

    complex(dp), intent(inout), allocatable, dimension(:,:,:,:) :: u
    real(dp), intent(inout), allocatable, dimension(:) :: beta
    integer(i4), intent(in) :: N_measurements, Nbeta, tau_Q
    logical :: equilibrium
    real(dp) :: betai, betaf
    
    integer(i4) :: i

    allocate(u(3,L(1),L(2),L(3)))

    Lx = L(1)
    Ly = L(2)
    Lz = L(3)
    
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

  subroutine thermalization(start,algorithm,u,beta,N_thermalization,isbeta)
    character(*), intent(in) :: start, algorithm
#ifdef SERIAL
    complex(dp), intent(inout) :: u(:,:,:,:)
#elif PARALLEL
    complex(dp), intent(inout) :: u(:,:,:,:)[*]
#endif
    real(dp), intent(in) :: beta
    logical, intent(in) :: isbeta 
    integer(i4) :: N_thermalization
    integer(i4) :: i_sweeps
    

    call select_start(u,start)
    do i_sweeps = 1, N_thermalization
       call sweeps(trim(algorithm),u,beta,isbeta)
    end do

  end subroutine thermalization

  subroutine measurements(algorithm,u,beta,N_measurements,Nskip,plq,top_den,isbeta)
    use U1_functions
    character(*), intent(in) :: algorithm
#ifdef SERIAL
    complex(dp), intent(inout) :: u(:,:,:,:)
    real(dp) :: plq(:), top_den(:)
#elif PARALLEL
    complex(dp), intent(inout) :: u(:,:,:,:)[*]
    real(dp) :: plq(:)[*], top_den(:)[*]
#endif
    real(dp), intent(in) :: beta
    integer(i4),intent(in) :: N_measurements, Nskip
    logical, intent(in) :: isbeta
    integer(i4) :: i_sweeps, iskip

    do i_sweeps = 1, N_measurements
       do iskip = 1, Nskip
          call sweeps(trim(algorithm),u,beta,isbeta)
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
  
  subroutine sweeps(algorithm,u,beta,isbeta)
    use parameters, only : Nhmc, Thmc
    character(*), intent(in) :: algorithm

#ifdef SERIAL
    complex(dp), intent(inout) :: u(:,:,:,:)
#elif PARALLEL
    complex(dp), intent(inout) :: u(:,:,:,:)[*]
#endif
    real(dp), intent(in) :: beta
    real(dp) :: b 
    logical :: isbeta 

    b = beta
    if (.not. isbeta .and. b > epsilon(1.0_dp) ) b = 1/b
    
    select case(algorithm)
    case("metropolis")
       if( .not. isbeta .and. b < epsilon(1.0_dp) )then
          call sweeps_alg(zero_temp_alg,u,b)
       else
          call sweeps_alg(metropolis,u,b)
       end if
    case("glauber")
       if( .not. isbeta .and. 1/b < epsilon(1.0_dp) ) then
          call sweeps_alg(zero_temp_alg,u,b)
       else
          call sweeps_alg(glauber,u,b)
       end if
    case("heatbath")
       if( .not. isbeta .and. b < epsilon(1.0_dp) ) then
          !call sweeps_alg(zero_temp_heatbath,u,b)
       else
          call sweeps_alg(heatbath,u,b)
       end if
#if defined(SERIAL) || PARALLEL == 1
    case("hmc")
       !if( .not. isbeta .and. b < epsilon(1.0_dp) ) error stop 'temperature zero not supported'
       call hmc(u,b,Thmc,Nhmc)
#endif
    end select
#ifdef PARALLEL
    sync all
#endif
    
  end subroutine sweeps

  subroutine sweeps_alg(alg,u,beta)
    procedure(lua_function) :: alg
#ifdef SERIAL
    complex(dp), intent(inout) :: u(:,:,:,:)
    integer(i4) :: point(3)
#elif PARALLEL == 1
    complex(dp), intent(inout) :: u(:,:,:,:)[*]
    integer(i4) :: point(0:3)
    integer :: thisimage
#elif PARALLEL == 2
    complex(dp), intent(inout) :: u(:,0:,0:,0:)[*]
    integer(i4) :: point(3)
#endif
    real(dp), intent(in) :: beta
    integer(i4) :: x,y,z,mu

#if PARALLEL == 1
    thisimage = this_image()
#endif

#if defined(SERIAL) || PARALLEL == 1
    do mu = 1, 3
       do x = 1, Lx
          do y = 1, Ly
             do z = 1, Lz
#if defined(SERIAL) 
                point = [x,y,z]
#elif PARALLEL == 1
                point = [thisimage,x,y,z]
#endif
                call alg(u,point,mu,beta)
             end do
          end do
       end do
    end do
#elif PARALLEL == 2
    !print*, "inside sweeps"
    do mu = 1, 3
       do x = 1, Lx - 1
          do y = 1, Ly - 1
             do z = 1, Lz - 1
                call alg(u,[x,y,z],mu,beta)
             end do
          end do
       end do
    end do
    !print*, "Finished 1st loop"
    !sync all
    u(:,Lx+1,1:Ly-1,1:Lz-1)[left] = u(:,1,1:Ly-1,1:Lz-1)
    u(:,1:Lx-1,Ly+1,1:Lz-1)[front]= u(:,1:Lx-1,1,1:Lz-1)
    u(:,1:Lx-1,1:Ly-1,Lz+1)[down] = u(:,1:Lx-1,1:Ly-1,1)

    u(:,Lx+1,Ly+1,1:Lz-1)[left_front]= u(:,1,1,1:Lz-1)
    u(:,Lx+1,1:Ly-1,Lz+1)[left_down] = u(:,1,1:Ly-1,1)
    u(:,1:Lx-1,Ly+1,Lz+1)[front_down]= u(:,1:Lx-1,1,1)
    
    u(:,Lx+1,Ly+1,Lz+1)[left_front_down] = u(:,1,1,1)
    sync all
    do mu = 1, 3
       do y = 1, Ly - 1
          do z = 1, Lz - 1
             call alg(u,[Lx,y,z],mu,beta)
          end do
       end do
    end do
    u(:,0,1:Ly-1,1:Lz-1)[right] = u(:,Lx,1:Ly-1,1:Lz-1)
    
    u(:,Lx,Ly+1,1:Lz-1)[front]      = u(:,Lx,1,1:Lz-1)
    u(:,0,Ly+1,1:Lz-1)[right_front] = u(:,Lx,1,1:Lz-1)
    
    u(:,Lx,1:Ly-1,Lz+1)[down]      = u(:,Lx,1:Ly-1,1) 
    u(:,0,1:Ly-1,Lz+1)[right_down] = u(:,Lx,1:Ly-1,1)

    u(:,Lx,Ly+1,Lz+1)[front_down]     = u(:,Lx,1,1)
    u(:,0,Ly+1,Lz+1)[right_front_down]= u(:,Lx,1,1)
    sync all

    do mu = 1, 3
       do x = 1, Lx - 1
          do z = 1, Lz - 1
             call alg(u,[x,Ly,z],mu,beta)
          end do
       end do
    end do
    u(:,1:Lx-1,0,1:Lz-1)[back] = u(:,1:Lx-1,Ly,1:Lz-1)

    u(:,Lx+1,Ly,1:Lz-1)[left]     = u(:,Lx,Ly,1:Lz-1)
    u(:,Lx+1,0,1:Lz-1)[left_back] = u(:,Lx,Ly,1:Lz-1)
    
    u(:,1:Lx-1,Ly,Lz+1)[down]     = u(:,1:Lx-1,Ly,1)
    u(:,1:Lx-1,0,Lz+1)[back_down] = u(:,1:Lx-1,Ly,1)

    u(:,Lx+1,Ly,Lz+1)[left_down]    = u(:,1,Ly,1)
    u(:,Lx+1,0,Lz+1)[left_back_down]= u(:,1,Ly,1)
    sync all

    do mu = 1, 3
       do x = 1, Lx -1
          do y = 1, Ly - 1
             call alg(u,[x,y,Lz],mu,beta)
          end do
       end do
    end do
    u(:,1:Lx-1,1:Ly-1,0)[up] = u(:,1:Lx-1,1:Ly-1,Lz)

    u(:,1:Lx-1,Ly+1,Lz)[front]   = u(:,1:Lx-1,1,Lz)
    u(:,1:Lx-1,Ly+1,0)[front_up] = u(:,1:Lx-1,1,Lz)

    u(:,Lx+1,1:Ly-1,Lz)[left] = u(:,1,1:Ly-1,Lz)
    u(:,0,1:Ly-1,Lz+1)[left_up] = u(:,1,1:Ly-1,Lz)

    u(:,Lx+1,Ly+1,Lz)[left_front]  = u(:,1,1,Lz)
    u(:,Lx+1,Ly+1,0)[left_front_up]= u(:,1,1,Lz)
    sync all

    do mu = 1, 3
       do z = 1, Lz - 1
          call alg(u,[Lx,Ly,z],mu,beta)
       end do
    end do
    u(:,0,Ly,1:Lz-1)[right]    = u(:,Lx,Ly,1:Lz-1)
    u(:,Lx,0,1:Lz-1)[back]     = u(:,Lx,Ly,1:Lz-1)
    u(:,0,0,1:Lz-1)[right_back]= u(:,Lx,Ly,1:Lz-1)

    u(:,Lx,Ly,Lz+1)[down]         = u(:,Lx,Ly,1)
    u(:,0,Ly,Lz+1)[right_down]    = u(:,Lx,Ly,1)
    u(:,Lx,0,Lz+1)[back_down]     = u(:,Lx,Ly,1)
    u(:,0,0,Lz+1)[right_back_down]= u(:,Lx,Ly,1)
    sync all

    do mu = 1, 3
       do y = 1, Ly - 1
          call alg(u,[Lx,y,Lz],mu,beta)
       end do
    end do
    u(:,0,1:Ly-1,Lz)[right]  = u(:,Lx,1:Ly-1,Lz)
    u(:,Lx,1:Ly-1,0)[up]     = u(:,Lx,1:Ly-1,Lz)
    u(:,0,1:Ly-1,0)[right_up]= u(:,Lx,1:Ly-1,Lz)

    u(:,Lx,Ly+1,Lz)[front]       = u(:,Lx,1,Lz)
    u(:,0,Ly+1,Lz)[right_front]  = u(:,Lx,1,Lz)
    u(:,Lx,Ly,0)[front_up]       = u(:,Lx,1,Lz)
    u(:,0,Ly+1,0)[right_front_up]= u(:,Lx,1,Lz)
    sync all

    do mu = 1, 3
       do x = 1, Lx - 1
          call alg(u,[x,Ly,Lz],mu,beta)
       end do
    end do
    u(:,1:Lx-1,0,Lz)[back]  = u(:,1:Lx-1,Ly,Lz)
    u(:,1:Lx-1,Ly,0)[up]    = u(:,1:Lx-1,Ly,Lz)
    u(:,1:Lx-1,0,0)[back_up]= u(:,1:Lx-1,Ly,Lz)
    
    u(:,Lx+1,Ly,Lz)[left]      = u(:,1,Ly,Lz)
    u(:,Lx+1,Ly,0)[left_up]    = u(:,1,Ly,Lz)
    u(:,Lx+1,0,Lz)[left_back]  = u(:,1,Ly,Lz)
    u(:,Lx+1,0,0)[left_back_up]= u(:,1,Ly,Lz)
    sync all

    do mu = 1, 3
       call alg(u,[Lx,Ly,Lz],mu,beta)
    end do
    u(:,0,Ly,Lz)[right]= u(:,Lx,Ly,Lz)
    u(:,Lx,0,Lz)[back] = u(:,Lx,Ly,Lz)
    u(:,Lx,Ly,0)[up]   = u(:,Lx,Ly,Lz)

    u(:,0,0,Lz)[right_back]= u(:,Lx,Ly,Lz)
    u(:,0,Ly,0)[right_up]  = u(:,Lx,Ly,Lz)
    u(:,Lx,0,0)[back_up]   = u(:,Lx,Ly,Lz)

    u(:,0,0,0)[right_back_up] = u(:,Lx,Ly,Lz)
    sync all
#endif
    
  end subroutine sweeps_alg
  
  
  subroutine eq(start,algorithm,u,beta, N_thermalization,Nskip,N_measurements,outunit,isbeta)
    use parameters, only : L
    use statistics

    character(*), intent(in) :: start, algorithm

#ifdef SERIAL
    complex(dp), intent(inout) :: u(:,:,:,:)
    real(dp), dimension(N_measurements) :: plq, top_den
#elif PARALLEL
    complex(dp), intent(inout) :: u(:,:,:,:)[*]
    real(dp), dimension(:), codimension[:], allocatable :: plq, top_den
#endif
    integer(i4), intent(in) :: N_thermalization, Nskip, N_measurements, outunit
    real(dp), intent(in) :: beta(:)
    logical, intent(in) :: isbeta
    integer(i4) :: ib, i_sweeps
    real(dp) :: err_plq(2), err_top(2)

#ifdef PARALLEL
    allocate(plq(N_measurements)[*])
    allocate(top_den(N_measurements)[*])
#endif
    
    do ib = 1, size(beta)
       call thermalization(trim(start),trim(algorithm),u,beta(ib),N_thermalization,isbeta)
       call measurements(trim(algorithm),u,beta(ib),N_measurements,Nskip,plq,top_den,isbeta)
#ifdef PARALLEL
       if(this_image() == 1)then
#endif
          err_plq = jackknife_max(plq)/(3*product(L))
          err_top = jackknife_max(top_den)/(twopi*product(L))
          print*, beta(ib), avr(plq)/(3*product(L)), err_plq(1), &
               avr(top_den)/(twopi*product(L)), err_top(1)
          write(outunit,*) beta(ib), avr(plq)/(3*product(L)), err_plq(1), &
               avr(top_den)/(twopi*product(L)), err_top(1)
          flush(outunit)
#ifdef PARALLEL
       end if
#endif
    end do
    
  end subroutine eq
  
  subroutine out_eq(start,algorithm,u,beta, tau_Q, N_thermalization, N_measurements,outunit,isbeta)
    use parameters, only : L
    use statistics
    character(*), intent(in) :: start, algorithm

#ifdef SERIAL
    complex(dp), intent(inout) :: u(:,:,:,:)
    real(dp), dimension(-tau_Q:tau_Q,N_measurements) :: plq, top_den
#elif PARALLEL
    complex(dp), intent(inout) :: u(:,:,:,:)[*]
    real(dp), dimension(:,:), codimension[:], allocatable :: plq, top_den
#endif
    integer(i4), intent(in) :: tau_Q, N_thermalization, N_measurements,outunit
    real(dp), intent(in) :: beta(-tau_Q:tau_Q)
    logical, intent(in) :: isbeta
    integer(i4) :: ib, i_sweeps
    
    real(dp) :: err_plq(2), err_top(2)

#ifdef PARALLEL
    allocate(plq(-tau_Q:tau_Q,N_measurements)[*])
    allocate(top_den(-tau_Q:tau_Q,N_measurements)[*])
#endif
    
    do i_sweeps = 1, N_measurements
       call thermalization(start,algorithm,u,beta(-tau_Q),N_thermalization,isbeta)
       do ib = -tau_Q, tau_Q
          call sweeps(algorithm,u,beta(ib),isbeta)
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
#ifdef PARALLEL
       if(this_image() == 1) then
#endif
          err_plq = jackknife_max(plq(ib,:))/(3*product(L))
          err_top = jackknife_max(top_den(ib,:))/(twopi*product(L))
          print*, beta(ib),ib, avr(plq(ib,:))/(3*product(L)), err_plq(1), &
               avr(top_den(ib,:))/(twopi*product(L)), err_top(1)
          write(outunit,*) ib, avr(plq(ib,:))/(3*product(L)), err_plq(1), &
               avr(top_den(ib,:))/(twopi*product(L)), err_top(1)
#ifdef PARALLEL
       endif
#endif
    end do
  end subroutine out_eq
  
end module dynamics

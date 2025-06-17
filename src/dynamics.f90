module dynamics
  use iso_fortran_env, only : dp => real64, i4 => int32
  use parameters, only : Lx, Ly, Lz, L, d
  use lua
  implicit none

  integer(i4) :: left,right, up, down, back, front, corner(8)
contains

  subroutine cold_start(u)
    complex(dp), intent(out), dimension(:,:,:,:) :: u
    u = 1.0_dp
  end subroutine cold_start
  
  subroutine set_memory(u,plq,beta,N_measurements,Nbeta,betai,betaf)
    use indices
    complex(dp), intent(inout), allocatable, dimension(:,:,:,:) :: u[:]
    real(dp), intent(inout), allocatable, dimension(:) :: plq[:], beta
    integer(i4), intent(in) :: N_measurements, Nbeta
    real(dp) :: betai, betaf
    integer(i4), dimension(:), allocatable :: ip1, im1, ip2, im2, ip3, im3
    integer(i4) :: i, a(3)

    Lx = L(1)/d(1)
    Ly = L(2)/d(2)
    Lz = L(3)/d(3)

    if(this_image() == 1) print*, Lx, Ly, Lz
    
    allocate(ip1(d(1)),im1(d(1)))
    allocate(ip2(d(2)),im2(d(2)))
    allocate(ip3(d(3)),im3(d(3)))
    allocate(u(3,0:Lx+1,0:Ly+1,0:Lz+1)[*])
    allocate(plq(N_measurements)[*])
    allocate(beta(Nbeta))

    beta = [(betai + (i-1)*(betaf-betai)/(Nbeta-1), i=1, Nbeta)]
    
    ip1 = [(mod(i,d(1))+1, i=1,d(1))]; im1 = [(i-1, i=1,d(1))]; im1(1) = d(1)
    ip2 = [(mod(i,d(2))+1, i=1,d(2))]; im2 = [(i-1, i=1,d(2))]; im2(1) = d(2)
    ip3 = [(mod(i,d(3))+1, i=1,d(3))]; im3 = [(i-1, i=1,d(3))]; im3(1) = d(3)

    a = get_index_array(this_image(),3,d)
    left = get_index([im1(a(1)),a(2),a(3)],3,d)
    right= get_index([ip1(a(1)),a(2),a(3)],3,d)
    up   = get_index([a(1),a(2),ip3(a(3))],3,d)
    down = get_index([a(1),a(2),im3(a(3))],3,d)
    front= get_index([a(1),im2(a(2)),a(3)],3,d)
    back = get_index([a(1),ip2(a(2)),a(3)],3,d)

    print*, this_image(), left, right, back, front, up, down
    
    corner(1) = get_index([im1(a(1)),im2(a(2)),ip3(a(3))],3,d)
    corner(2) = get_index([ip1(a(1)),im2(a(2)),ip3(a(3))],3,d)
    corner(3) = get_index([im1(a(1)),im2(a(2)),im3(a(3))],3,d)
    corner(4) = get_index([ip1(a(1)),im2(a(2)),im3(a(3))],3,d)
    corner(5) = get_index([im1(a(1)),ip2(a(2)),ip3(a(3))],3,d)
    corner(6) = get_index([ip1(a(1)),ip2(a(2)),ip3(a(3))],3,d)
    corner(7) = get_index([im1(a(1)),ip2(a(2)),im3(a(3))],3,d)
    corner(8) = get_index([ip1(a(1)),ip2(a(2)),im3(a(3))],3,d)


    print*, this_image(), corner
   
   
    
  end subroutine set_memory

  subroutine thermalization(u,beta,N_thermalization)
    complex(dp), intent(inout) :: u(3,0:Lx+1,0:Ly+1,0:Lz+1)[*]
    real(dp), intent(in) :: beta
    integer(i4) :: N_thermalization
    integer(i4) :: i_sweeps

    do i_sweeps = 1, N_thermalization
       call sweeps(u,beta)
    end do

  end subroutine thermalization

  subroutine measurements(u,beta,N_measurements,Nskip,plq)
    use U1_functions
    complex(dp), intent(inout) :: u(3,0:Lx+1,0:Ly+1,0:Lz+1)[*]
    real(dp), intent(in) :: beta
    integer(i4) :: N_measurements, Nskip
    real(dp) :: plq(:)[*]
    integer(i4) :: i_sweeps, iskip

    do i_sweeps = 1, N_measurements
       do iskip = 1, Nskip
          call sweeps(u,beta)
       end do
       plq(i_sweeps) = plaquette_value(u)
       call co_sum(plq(i_sweeps),result_image = 1)
    end do
    
  end subroutine measurements
  
  subroutine sweeps(u,beta)
    complex(dp), intent(inout) :: u(3,0:Lx+1,0:Ly+1,0:Lz+1)[*]
    real(dp), intent(in) :: beta

    integer(i4) :: x,y,z,mu

    ! Update interior
    do x = 2, Lx-1
       do y = 2, Ly-1
          do z = 2, Lz-1
             do mu = 1, 3
                call metropolis(u,[x,y,z],mu,beta)
             end do
          end do
       end do
    end do   

    ! Update lower band. Done
    do x = 2, Lx-1
       do y = 2, Ly-1
          do mu = 1, 3
             call metropolis(u,[x,y,1],mu,beta)
          end do
       end do
    end do
    u(:,2:Lx-1,2:Ly-1,Lz+1)[down] = u(:,2:Lx-1,2:Ly-1,1)
    sync all
    
    ! Update upper band. Done
    do x = 2, Lx-1
       do y = 2, Ly-1
          do mu = 1, 3
             call metropolis(u,[x,y,Lz],mu,beta)
          end do
       end do
    end do
    u(:,2:Lx-1,2:Ly-1,0)[up] = u(:,2:Lx-1,2:Ly-1,Lz)
    sync all
    
    ! Update left band. Done
    do y = 2, Ly-1
       do z = 2, Lz-1
          do mu = 1, 3
             call metropolis(u,[1,y,z],mu,beta)
          end do
       end do
    end do
    u(:,Lx+1,2:Ly-1,2:Lz-1)[left] = u(:,1,2:Ly-1,2:Lz-1)
    sync all
    
    ! Update right band. Done
    do y = 2, Ly-1
       do z = 2, Lz-1
          do mu = 1, 3
             call metropolis(u,[Lx,y,z],mu,beta)
          end do
       end do
    end do
    u(:,0,2:Ly-1,2:Lz-1)[right] = u(:,Lx,2:Ly-1,2:Lz-1)
    sync all

    ! Update Front band. Done
    do x = 2, Lx-1
       do z = 2, Lz-1
          do mu = 1, 3
             call metropolis(u,[x,1,z],mu,beta)
          end do
       end do
    end do
    u(:,2:Lx-1,Ly+1,2:Lz-1)[front] = u(:,2:Lx-1,1,2:Lz-1)
    sync all
    
    ! Update back band. Done
    do x = 2, Lx-1
       do z = 2, Lz-1
          do mu = 1, 3
             call metropolis(u,[x,Ly,z],mu,beta)
          end do
       end do
    end do
    u(:,2:Lx-1,0,2:Lz-1)[back] = u(:,2:Lx-1,Ly,2:Lz-1)
    sync all

    
    ! Update corners
    ! Up left front. Done
    do mu = 1, 3
       call metropolis(u,[1,1,Lz],mu,beta)
    end do
    u(:,1,1,0)[up] = u(:,1,1,Lz)
    u(:,Lx+1,1,Lz)[left] = u(:,1,1,Lz)
    u(:,1,Ly+1,Lz)[front] = u(:,1,1,Lz)
    u(:,Lx+1,Ly+1,0)[corner(1)] = u(:,1,1,Lz) 
    sync all
    
    ! Up right front. Done
    do mu = 1, 3
       call metropolis(u,[Lx,1,Lz],mu,beta)
    end do
    u(:,Lx,1,0)[up] = u(:,Lx,1,Lz)
    u(:,0,1,Lz)[right] = u(:,Lx,1,Lz)
    u(:,Lx,Ly+1,Lz)[front] = u(:,Lx,1,Lz)
    u(:,0,Ly+1,0)[corner(2)] = u(:,Lx,1,Lz) 
    sync all
    
    ! Down left front. Done
    do mu = 1, 3
       call metropolis(u,[1,1,1],mu,beta)
    end do
    u(:,Lx+1,1,1)[left] = u(:,1,1,1)
    u(:,1,1,Lz+1)[down] = u(:,1,1,1)
    u(:,1,Ly+1,1)[front] = u(:,1,1,1)
    u(:,Lx+1,Ly+1,Lz+1)[corner(3)] = u(:,1,1,1) 
    sync all
    
    ! Down right front. done
    do mu = 1, 3
       call metropolis(u,[Lx,1,1],mu,beta)
    end do
    u(:,Lx,1,Lz+1)[down] = u(:,Lx,1,1)
    u(:,0,1,1)[right] = u(:,Lx,1,1)
    u(:,Lx,Ly+1,1)[front] = u(:,Lx,1,1)
    u(:,0,Ly+1,Lz+1)[corner(4)] = u(:,Lx,1,1) 
    sync all

    ! Up left back. Done
    do mu = 1, 3
       call metropolis(u,[1,Ly,Lz],mu,beta)
    end do
    u(:,1,Ly,0)[up] = u(:,1,Ly,Lz)
    u(:,Lx+1,Ly,Lz)[left] = u(:,1,Ly,Lz)
    u(:,1,0,Lz)[back] = u(:,1,Ly,Lz)
    u(:,Lx+1,0,0)[corner(5)] = u(:,1,Ly,Lz) 
    sync all
    
    ! Up right back. Done
    do mu = 1, 3
       call metropolis(u,[Lx,Ly,Lz],mu,beta)
    end do
    u(:,Lx,Ly,0)[up] = u(:,Lx,Ly,Lz)
    u(:,0,Ly,Lz)[right] = u(:,Lx,Ly,Lz)
    u(:,Lx,0,Lz)[back] = u(:,Lx,Ly,Lz)
    u(:,0,0,0)[corner(6)] = u(:,Lx,Ly,Lz) 
    sync all
    
    ! Down left back. Done
    do mu = 1, 3
       call metropolis(u,[1,Ly,1],mu,beta)
    end do
    u(:,Lx+1,Ly,1)[left] = u(:,1,Ly,1)
    u(:,1,Ly,Lz+1)[down] = u(:,1,Ly,1)
    u(:,1,0,1)[back] = u(:,1,Ly,1)
    u(:,Lx+1,0,Lz+1)[corner(7)] = u(:,1,Ly,1) 
    sync all
    
    ! Down right back. Done
    do mu = 1, 3
       call metropolis(u,[Lx,Ly,1],mu,beta)
    end do
    u(:,Lx,Ly,Lz+1)[down] = u(:,Lx,Ly,1)
    u(:,0,Ly,1)[right] = u(:,Lx,Ly,1)
    u(:,Lx,0,1)[back] = u(:,Lx,Ly,1)
    u(:,0,0,Lz+1)[corner(8)] = u(:,Lx,Ly,1) 
    sync all
    
  end subroutine sweeps




end module dynamics

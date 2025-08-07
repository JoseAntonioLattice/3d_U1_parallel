#define POINT(a,b,c) [a,b,c]

#if defined(SERIAL)

#define SYNCIM 

#elif PARALLEL == 1

#undef POINT
#define POINT(a,b,c) [thisimage,a,b,c]
#define SYNCIM sync all

#elif PARALLEL == 2

#define SYNCIM call sync_sublattice(unew)

#endif

module hybridMC

  use iso_fortran_env, only : dp => real64, i4 => int32
  use random
  use constants
  use U1_functions
  use parallel_utils
  implicit none

  integer(i4) :: Lx, Ly, Lz

contains

  subroutine hmc(u,beta,time,N)
#ifdef SERIAL
    complex(dp), intent(inout), dimension(:,:,:,:) :: u
    complex(dp), dimension(size(u(:,1,1,1)),size(u(1,:,1,1)),size(u(1,1,:,1)),size(u(1,1,1,:))) :: unew
    real(dp) :: DeltaH
    real(dp), dimension(size(u(:,1,1,1)),size(u(1,:,1,1)),size(u(1,1,:,1)),size(u(1,1,1,:))) :: p, F, pnew
#elif PARALLEL == 1
    complex(dp), intent(inout), dimension(:,:,:,:) :: u[*]
    complex(dp), dimension(:,:,:,:), allocatable :: unew[:]
    real(dp), codimension[:], allocatable :: DeltaH
    real(dp), dimension(size(u(:,1,1,1)),size(u(1,:,1,1)),size(u(1,1,:,1)),size(u(1,1,1,:))) :: p, F, pnew
#elif PARALLEL == 2
    complex(dp), intent(inout), dimension(:,0:,0:,0:) :: u[*]
    complex(dp), dimension(:,:,:,:), allocatable :: unew[:]
    real(dp), codimension[:], allocatable :: DeltaH
    real(dp), dimension(3,size(u(1,:,1,1))-2,size(u(1,1,:,1))-2,size(u(1,1,1,:))-2) :: p, F, pnew
#endif
    real(dp), intent(in) :: beta,time
    integer(i4), intent(in) :: N
    integer(i4) :: k
    real(dp) :: dt, r
#if defined(SERIAL)
    logical :: condition
#elif PARALLEL
    logical, allocatable :: condition[:]
#endif

    Lx = size(u(1,:,1,1))
    Ly = size(u(1,1,:,1))
    Lz = size(u(1,1,1,:))

#if PARALLEL == 1
    allocate(unew(3,Lx,Ly,Lz)[*])
    allocate(DeltaH[*],condition[*])
#elif PARALLEL == 2
    Lx = Lx - 2
    Ly = Ly - 2
    Lz = Lz - 2
    allocate(unew(3,0:Lx+1,0:Ly+1,0:Lz+1)[*])
    allocate(DeltaH[*],condition[*])
#endif
    dt = time/N
    p = momentum([Lx,Ly,Lz])

    unew = u
    pnew = p

    unew(:,1:Lx,1:Ly,1:Lz) = unew(:,1:Lx,1:Ly,1:Lz)*exp(0.5*ii*dt*pnew)
    SYNCIM
    F = force(unew,beta)
    
    do k = 1, N - 2
       pnew = pnew + dt*F
       unew(:,1:Lx,1:Ly,1:Lz) = unew(:,1:Lx,1:Ly,1:Lz)*exp(ii*dt*pnew)
       SYNCIM
       F = force(unew,beta)
    end do
    
    pnew = pnew + dt*F
    unew(:,1:Lx,1:Ly,1:Lz) = unew(:,1:Lx,1:Ly,1:Lz)*exp(0.5*ii*dt*pnew)
    SYNCIM
    
    DeltaH = DH(u,unew,p,pnew,beta)

#if defined(PARALLEL)
    sync all
    call co_sum(DeltaH,result_image = 1)
    if(this_image() == 1) then 
#endif
       call random_number(r)
       condition = r <= exp(-DeltaH)
#if defined(PARALLEL)
    end if
    call co_broadcast(condition,source_image=1)
#endif
    if( condition ) u = unew

  end subroutine hmc

#if PARALLEL == 2
  subroutine sync_sublattice(u)
    complex(dp), dimension(:,0:,0:,0:) :: u[*]

    !sync all
    ! Send faces
    u(:, 1:Lx  , 1:Ly  ,   Lz+1)[down]  = u(:, 1:Lx, 1:Ly,  1   )
    u(:, 1:Lx  ,   Ly+1, 1:Lz  )[front] = u(:, 1:Lx, 1   ,  1:Lz)
    u(:,   Lx+1, 1:Ly  , 1:Lz  )[left]  = u(:, 1   , 1:Ly,  1:Lz)
    u(:, 1:Lx  , 1:Ly  , 0     )[up]    = u(:, 1:Lx, 1:Ly,    Lz)
    u(:, 1:Lx  , 0     , 1:Lz  )[back]  = u(:, 1:Lx,   Ly,  1:Lz)
    u(:, 0     , 1:Ly  , 1:Lz  )[right] = u(:,   Lx, 1:Ly,  1:Lz)
    !sync all
    
    ! Send edges
    u(:, 1:Lx, Ly+1, Lz+1)[front_down] = u(:, 1:Lx, 1 , 1 )
    u(:, 1:Lx, Ly+1, 0   )[front_up]   = u(:, 1:Lx, 1 , Lz)
    u(:, 1:Lx, 0   , Lz+1)[back_down]  = u(:, 1:Lx, Ly, 1 )
    u(:, 1:Lx, 0   , 0   )[back_up]    = u(:, 1:Lx, Ly, Lz)

    u(:, Lx+1, 1:Ly, Lz+1)[left_down] = u(:, 1 , 1:Ly ,1 )
    u(:, Lx+1, 1:Ly, 0   )[left_up]   = u(:, 1 , 1:Ly ,Lz)
    u(:, 0   , 1:Ly, Lz+1)[right_down]= u(:, Lx, 1:Ly ,1 )
    u(:, 0   , 1:Ly, 0   )[right_up]  = u(:, Lx, 1:Ly ,Lz)
       
    u(:, Lx+1, Ly+1, 1:Lz)[left_front] = u(:, 1 , 1 , 1:Lz)
    u(:, Lx+1, 0   , 1:Lz)[left_back]  = u(:, 1 , Ly, 1:Lz)
    u(:, 0   , Ly+1, 1:Lz)[right_front]= u(:, Lx, 1 , 1:Lz)
    u(:, 0   , 0   , 1:Lz)[right_back] = u(:, Lx, Ly, 1:Lz)
    
    !sync all
    ! Send corners
    u(:, Lx+1, Ly+1, Lz+1)[left_front_down]  = u(:, 1 , 1, 1)
    u(:, 0   , Ly+1, Lz+1)[right_front_down] = u(:, Lx, 1, 1)
    u(:, Lx+1, Ly+1, 0   )[left_front_up]    = u(:, 1 , 1, Lz)
    u(:, 0   , Ly+1, 0   )[right_front_up]   = u(:, Lx, 1, Lz)
    
    u(:, Lx+1, 0, Lz+1)[left_back_down] = u(:, 1 , Ly, 1)
    u(:, 0   , 0, Lz+1)[right_back_down]= u(:, Lx, Ly, 1)
    u(:, Lx+1, 0, 0   )[left_back_up]   = u(:, 1 , Ly, Lz)
    u(:, 0   , 0, 0   )[right_back_up]  = u(:, Lx, Ly, Lz)
    sync all
  end subroutine sync_sublattice
#endif

  function force(u,beta)
#ifdef SERIAL
    complex(dp), intent(in), dimension(:,:,:,:) :: u
    complex(i4), dimension(size(u(:,1,1,1)),size(u(1,:,1,1)),size(u(1,1,:,1)),size(u(1,1,1,:))) :: force
#elif PARALLEL == 1
    complex(dp), intent(in), dimension(:,:,:,:) :: u[*]
    integer(i4) :: thisimage
    complex(i4), dimension(size(u(:,1,1,1)),size(u(1,:,1,1)),size(u(1,1,:,1)),size(u(1,1,1,:))) :: force
#elif PARALLEL == 2
    complex(dp), intent(in), dimension(:,0:,0:,0:) :: u[*]
    complex(i4), dimension(3,size(u(1,:,1,1))-2,size(u(1,1,:,1))-2,size(u(1,1,1,:))-2) :: force
#endif
    real(dp), intent(in) :: beta
    integer(i4) :: x,y,z,mu
    complex(dp) :: stp

#if PARALLEL == 1
    thisimage = this_image()
#endif
    
    do mu = 1, 3
       do x = 1, Lx
          do y = 1, Ly
             do z = 1, Lz
                stp = u(mu,x,y,z)*conjg(staples(u,POINT(x,y,z),mu))
                force(mu,x,y,z) = -beta*stp%im
             end do
          end do
       end do
    end do
    
  end function force

  function momentum(dim)
    integer(i4), intent(in) :: dim(:)
    real(dp), dimension(size(dim),dim(1),dim(2),dim(3)) :: momentum
    real(dp) :: u1, u2
    integer(i4) :: x, y, z

    do x = 1, dim(1)
       do y = 1, dim(2)
          do z = 1, dim(3)
             momentum(1,x,y,z) = normal()
             momentum(2,x,y,z) = normal()
             momentum(3,x,y,z) = normal()
          end do
       end do
    end do
    
  end function momentum

  function DH(U,Unew,P,Pnew,beta)
    real(dp) :: DH
#ifdef SERIAL
    complex(dp), intent(in), dimension(:,:,:,:) :: u, Unew
#elif PARALLEL == 1
    complex(dp), intent(in), dimension(:,:,:,:) :: u[*], Unew[*]
    integer(i4) :: thisimage
#elif PARALLEL == 2
    complex(dp), intent(in), dimension(:,0:,0:,0:) :: u[*],  Unew[*]
#endif
    real(dp), dimension(:,:,:,:), intent(in) :: P, Pnew
    real(dp), intent(in) :: beta
    integer(i4) :: x,y,z,mu,nu
    real(dp) :: DeltaS

#if PARALLEL == 1
    thisimage = this_image()
#endif
        
    DH = 0.0_dp
    DeltaS = 0.0_dp
    do x = 1, Lx
       do y = 1, Ly
          do z = 1, Lz
             do mu = 1, 2
                DH = DH + (pnew(mu,x,y,z))**2 - (p(mu,x,y,z))**2
                do nu = mu+1, 3
                   DeltaS = DeltaS + real(plaquette(u,POINT(x,y,z),mu,nu) - plaquette(unew,POINT(x,y,z),mu,nu))
                end do
             end do
             DH = DH + (pnew(3,x,y,z))**2 - (p(3,x,y,z))**2
          end do
       end do
    end do
    
    DeltaS = beta*DeltaS
    DH = 0.5*DH + DeltaS
    
  end function DH

end module hybridMC

module hybridMC

  use iso_fortran_env, only : dp => real64, i4 => int32
  use random
  use constants
  use U1_functions
  implicit none

contains

  subroutine hmc(u,beta,time,N)
#ifdef SERIAL
    complex(dp), intent(inout), dimension(:,:,:,:) :: u
    complex(dp), dimension(size(u(:,1,1,1)),size(u(1,:,1,1)),size(u(1,1,:,1)),size(u(1,1,1,:))) :: unew
    real(dp) :: DeltaH, r
#elif PARALLEL == 1
    complex(dp), intent(inout), dimension(:,:,:,:) :: u[*]
    complex(dp), dimension(:,:,:,:), allocatable :: unew[:]
    real(dp), codimension[:], allocatable :: DeltaH, r
#endif
    real(dp), intent(in) :: beta,time
    integer(i4), intent(in) :: N
    integer(i4) :: k, Lx, Ly, Lz
    real(dp) :: dt
    
    real(dp), dimension(size(u(:,1,1,1)),size(u(1,:,1,1)),size(u(1,1,:,1)),size(u(1,1,1,:))) :: p, F, pnew

    Lx = size(u(1,:,1,1))
    Ly = size(u(1,1,:,1))
    Lz = size(u(1,1,1,:))
    
    dt = time/N
    p = momentum([Lx,Ly,Lz])
#if PARALLEL == 1
    allocate(unew(3,Lx,Ly,Lz)[*])
    allocate(DeltaH[*],r[*])
#endif
    unew = u
    pnew = p

    unew = unew*exp(0.5*ii*dt*pnew)
    F = force(unew,beta)
    
    do k = 1, N - 2
       pnew = pnew + dt*F
       unew = unew*exp(ii*dt*pnew)
#if PARALLEL == 1
       sync all
#endif
       F = force(unew,beta)
    end do
    
    pnew = pnew + dt*F
    unew = unew*exp(0.5*ii*dt*pnew)
#if PARALLEL == 1
       sync all
#endif
    DeltaH = DH(u,unew,p,pnew,beta)
    call random_number(r)
#if PARALLEL == 1
    call co_broadcast(DeltaH,source_image=1)
    call co_broadcast(r,source_image=1)
#endif
    if( r <= exp(-DeltaH)) u = unew

  end subroutine hmc

  function force(u,beta)
#ifdef SERIAL
    complex(dp), intent(inout), dimension(:,:,:,:) :: u
    integer(i4) :: point(3)
#elif PARALLEL == 1
    complex(dp), intent(inout), dimension(:,:,:,:) :: u[*]
    integer(i4) :: point(0:3)
    integer(i4) :: thisimage
#endif
    real(dp), intent(in) :: beta
    complex(i4), dimension(size(u(:,1,1,1)),size(u(1,:,1,1)),size(u(1,1,:,1)),size(u(1,1,1,:))) :: force
    integer(i4) :: Lx,Ly,Lz, x, y,z,mu
    complex(dp) :: stp
    
    Lx = size(u(1,:,1,1))
    Ly = size(u(1,1,:,1))
    Lz = size(u(1,1,1,:))

#if PARALLEL == 1
    thisimage = this_image()
#endif
    
    do mu = 1, 3
       do x = 1, Lx
          do y = 1, Ly
             do z = 1, Lz
#ifdef SERIAL
                point = [x,y,z]
#elif PARALLEL == 1
                point = [thisimage,x,y,z]
#endif
                stp = u(mu,x,y,z)*conjg(staples(u,point,mu))
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
    integer(i4) :: Lx,Ly,Lz
    
    Lx = dim(1)
    Ly = dim(2)
    Lz = dim(3)
    
    do x = 1, Lx
       do y = 1, Ly
          do z = 1, Lz
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
    real(dp) :: DHAM
    complex(dp), intent(inout), dimension(:,:,:,:) :: u
    complex(dp), dimension(:,:,:,:), intent(in) :: Unew
    integer(i4) :: point(3) 
#elif PARALLEL == 1
    real(dp), allocatable :: DHAM[:]
    complex(dp), intent(inout), dimension(:,:,:,:) :: u[*]
    complex(dp), dimension(:,:,:,:), intent(in) :: Unew[*]
    integer(i4) :: point(0:3), thisimage
#endif
    
    real(dp), dimension(:,:,:,:), intent(in) :: P, Pnew
    real(dp), intent(in) :: beta
    integer(i4) :: x,y,z,mu,nu, Lx,Ly,Lz
    real(dp) :: DeltaS

    Lx = size(u(1,:,1,1))
    Ly = size(u(1,1,:,1))
    Lz = size(u(1,1,1,:))
#if PARALLEL == 1
    allocate(DHAM[*])
    thisimage = this_image()
#endif
        
    DHAM = 0.0_dp
    DeltaS = 0.0_dp
    do x = 1, Lx
       do y = 1, Ly
          do z = 1, Lz
             do mu = 1, 2
                DHAM = DHAM + (pnew(mu,x,y,z))**2 - (p(mu,x,y,z))**2
                do nu = mu+1, 3
#ifdef SERIAL
                   point = [x,y,z]
#elif PARALLEL == 1
                   point = [thisimage,x,y,z]
#endif
                   DeltaS = DeltaS + real(plaquette(u,point,mu,nu) - plaquette(unew,point,mu,nu))
                end do
             end do
             DHAM = DHAM + (pnew(3,x,y,z))**2 - (p(3,x,y,z))**2
          end do
       end do
    end do
    
    DeltaS = beta*DeltaS
    DHAM = 0.5*DHAM + DeltaS
#if PARALLEL == 1
    sync all
    call co_sum(DHAM,result_image=1)
#endif
    DH = DHAM
  end function DH

end module hybridMC

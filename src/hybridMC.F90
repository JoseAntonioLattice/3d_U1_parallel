module hybridMC

  use iso_fortran_env, only : dp => real64, i4 => int32
  use random
  use constants
  use U1_functions
  implicit none

contains

  subroutine hmc(u,beta,time,N)
    complex(dp), intent(inout), dimension(:,:,:,:) :: u
    real(dp), intent(in) :: beta,time
    integer(i4), intent(in) :: N
    integer(i4) :: k, Lx, Ly, Lz
    real(dp) :: dt, DeltaH, r
    complex(dp), dimension(size(u(:,1,1,1)),size(u(1,:,1,1)),size(u(1,1,:,1)),size(u(1,1,1,:))) :: unew
    real(dp), dimension(size(u(:,1,1,1)),size(u(1,:,1,1)),size(u(1,1,:,1)),size(u(1,1,1,:))) :: p, F, pnew

    Lx = size(u(1,:,1,1))
    Ly = size(u(1,1,:,1))
    Lz = size(u(1,1,1,:))
    
    dt = time/N
    p = momentum([Lx,Ly,Lz])
     
    unew = u
    pnew = p

    unew = unew*exp(0.5*ii*dt*pnew)
    F = force(unew,beta)
    
    do k = 1, N - 2
       pnew = pnew + dt*F
       unew = unew*exp(ii*dt*pnew)
       F = force(unew,beta)
    end do
    
    pnew = pnew + dt*F
    unew = unew*exp(0.5*ii*dt*pnew)
    
    call random_number(r)
    
    DeltaH = DH(u,unew,p,pnew,beta)
    if( r <= exp(-DeltaH)) u = unew

  end subroutine hmc

  function force(u,beta)
    complex(dp), intent(inout), dimension(:,:,:,:) :: u
    real(dp), intent(in) :: beta
    complex(i4), dimension(size(u(:,1,1,1)),size(u(1,:,1,1)),size(u(1,1,:,1)),size(u(1,1,1,:))) :: force
    integer(i4) :: Lx,Ly,Lz, x, y,z,mu
    complex(dp) :: stp
    
    Lx = size(u(1,:,1,1))
    Ly = size(u(1,1,:,1))
    Lz = size(u(1,1,1,:))
    
    do x = 1, Lx
       do y = 1, Ly
          do z = 1, Lz
             do mu = 1, 3
                stp = u(mu,x,y,z)*conjg(staples(u,[x,y,z],mu))
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
    complex(dp), dimension(:,:,:,:), intent(in) :: U, Unew
    real(dp), dimension(:,:,:,:), intent(in) :: P, Pnew
    real(dp), intent(in) :: beta
    integer(i4) :: x,y,z,mu,nu, Lx,Ly,Lz
    real(dp) :: DeltaS

    Lx = size(u(1,:,1,1))
    Ly = size(u(1,1,:,1))
    Lz = size(u(1,1,1,:))
    
    
    DH = 0.0_dp
    DeltaS = 0.0_dp
    do x = 1, Lx
       do y = 1, Ly
          do z = 1, Lz
             do mu = 1, 2
                DH = DH + (pnew(mu,x,y,z))**2 - (p(mu,x,y,z))**2
                do nu = mu+1, 3 
                   DeltaS = DeltaS + real(plaquette(u,[x,y,z],mu,nu) - plaquette(unew,[x,y,z],mu,nu))
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

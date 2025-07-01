module arrays

  use iso_fortran_env, only : dp => real64, i4 => int32
  implicit none
  
#ifdef SERIAL 
  complex(dp), allocatable, dimension(:,:,:,:) :: u
  real(dp), allocatable, dimension(:) :: plq, beta
  
#endif

#ifdef PARALLEL  
  complex(dp), allocatable, dimension(:,:,:,:) :: U[:]
  real(dp), allocatable, dimension(:) :: plq[:], beta
#endif
  
end module arrays

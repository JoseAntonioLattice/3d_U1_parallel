module arrays

  use iso_fortran_env, only : dp => real64, i4 => int32
  implicit none
  
#ifdef SERIAL 
  complex(dp), allocatable, dimension(:,:,:,:) :: u
#endif
#ifdef PARALLEL  
  complex(dp), allocatable, dimension(:,:,:,:) :: U[:]
#endif
  real(dp), allocatable, dimension(:) :: beta
  
end module arrays

module arrays

  use iso_fortran_env, only : dp => real64, i4 => int32
  implicit none

  complex(dp), allocatable, dimension(:,:,:,:) :: U[:]
  real(dp), allocatable, dimension(:) :: plq[:], beta

end module arrays

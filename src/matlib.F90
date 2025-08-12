module matlib

  use iso_fortran_env, only : dp => real64, i4 => int32
  implicit none

contains

  function linspace(xi,xf,n)
    real(dp) :: linspace(n)
    integer(i4), intent(in) :: n
    real(dp), intent(in) :: xi, xf
    integer :: i
    
    linspace = [(xi + (i-1)*(xf-xi)/(n-1), i = 1, n)]
    
  end function linspace
  

end module matlib

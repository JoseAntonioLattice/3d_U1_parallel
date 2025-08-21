subroutine progress_bar(x)
  use iso_fortran_env, only : sp => real32
  implicit none
  real(sp) :: x

  write(*,'(a,"Progress: ",f6.2,"%")', advance = 'no') char(13), x*100
  
end subroutine progress_bar

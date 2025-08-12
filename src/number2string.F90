module number2string

    use iso_fortran_env, only : dp => real64, i4 => int32

    implicit none

    private

    public :: int2str, real2str

contains

    function int2str(k)
      !"Convert an integer to string."
      integer(i4), intent(in) :: k
      character(:), allocatable :: int2str
      character(20) :: var
      write (var, *) k
      int2str = trim(adjustl(var))
    end function int2str

    function real2str(r,int,decimal)
      real(dp), intent(in) :: r
      character(20) :: str
      integer(i4) :: int, decimal
      character(:), allocatable :: real2str

      
      write(str,"(f"//int2str(int+decimal+2)//"."//int2str(decimal)//")") r
      real2str = trim(adjustl(str))
    end function real2str

  end module number2string

module pbc

  use iso_fortran_env, only : i4 => int32

  implicit none

  integer, allocatable, dimension(:) :: ip_array, im_array, ip_t, im_t

contains

  subroutine set_periodic_bounds(L,Lt)

    integer(i4), intent(in) :: L, Lt

    allocate(ip_array(L),im_array(L),ip_t(Lt),im_t(Lt))
    call initialize(L,Lt)

  end subroutine set_periodic_bounds

  subroutine initialize(L,Lt)

    integer(i4), intent(in) :: L,Lt
    integer(i4) :: i

    do i = 1, L
       ip_array(i) = i + 1
       im_array(i) = i - 1
    end do
    ip_array(L) = 1
    im_array(1) = L

    
    do i = 1, Lt
       ip_t(i) = i + 1
       im_t(i) = i - 1
    end do
    ip_t(Lt) = 1
    im_t(1) = Lt
    
  end subroutine initialize

  pure function ip(vector,mu)
    integer, dimension(:), intent(in) :: vector
    integer, intent(in) :: mu

    integer, dimension(size(vector)) :: ip


    ip = vector

    if (mu ==  size(vector)) then
       ip(mu) = ip_t(vector(mu))
    else
       ip(mu) = ip_array(vector(mu))
    end if
  end function ip


  pure function im(vector,mu)

    integer, dimension(:), intent(in) :: vector
    integer, intent(in) :: mu

    integer, dimension(size(vector)) :: im


    im = vector

    
    if (mu == size(vector))then
       im(mu) = im_t(vector(mu))
    else
       im(mu) = im_array(vector(mu))
    end if
  end function im

end module pbc

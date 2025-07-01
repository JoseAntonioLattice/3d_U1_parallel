module pbc

  use iso_fortran_env, only : i4 => int32
  implicit none

#ifdef SERIAL
  integer(i4), allocatable, dimension(:) :: ip1,ip2,ip3, im1,im2,im3
#endif

  
#ifdef PARALLEL
  integer(i4), allocatable, dimension(:,:), codimension[:] :: ip1,ip2,ip3, im1,im2,im3
  integer(i4), dimension(3) :: ip1_c, ip2_c, ip3_c, im1_c, im2_c, im3_c
#endif
  
contains

#ifdef SERIAL

  subroutine initialize_pbc(L)
    integer(i4), intent(in) :: L(3)

    allocate(ip1(L(1)),im1(L(1)))
    allocate(ip2(L(2)),im2(L(2)))
    allocate(ip3(L(3)),im3(L(3)))

    call set_periodic_bounds(ip1,im1,L(1))
    call set_periodic_bounds(ip2,im2,L(2))
    call set_periodic_bounds(ip3,im3,L(3))
    
  end subroutine initialize_pbc

  
subroutine set_periodic_bounds(ip_array,im_array,L)
    integer(i4), intent(in) :: L
    integer, dimension(L), intent(out) :: ip_array, im_array
    integer(i4) :: i

    do i = 1, L
       ip_array(i) = i + 1
       im_array(i) = i - 1
    end do
    ip_array(L) = 1
    im_array(1) = L
    
  end subroutine set_periodic_bounds
  
  function ip(x,mu)
    integer(i4), dimension(3) :: ip
    integer(i4), intent(in) :: x(3), mu

    ip = x
    select case(mu)
    case(1)
       ip(mu) = ip1(x(mu))
    case(2)
       ip(mu) = ip2(x(mu))
    case(3)
       ip(mu) = ip3(x(mu))
    end select
    
  end function ip

  function im(x,mu)
    integer(i4), dimension(3) :: im
    integer(i4), intent(in) :: x(3), mu

    im = x
    select case(mu)
    case(1)
       im(mu) = im1(x(mu))
    case(2)
       im(mu) = im2(x(mu))
    case(3)
       im(mu) = im3(x(mu))
    end select
    
  end function im
#endif
  
#ifdef PARALLEL

    subroutine initialize_pbc(L,d,this_core,left,right)
    integer(i4), intent(in) :: L(3),d(3), right(3), left(3), this_core
    integer :: i
    
    allocate(ip1(L(1),2)[*],im1(L(1),2)[*])
    allocate(ip2(L(2),2)[*],im2(L(2),2)[*])
    allocate(ip3(L(3),2)[*],im3(L(3),2)[*])

    
    call set_periodic_bounds(ip1,im1,L(1),this_core,left(1),right(1))
    call set_periodic_bounds(ip2,im2,L(2),this_core,left(2),right(2))
    call set_periodic_bounds(ip3,im3,L(3),this_core,left(3),right(3))
    
  end subroutine initialize_pbc

  subroutine set_periodic_bounds(ip_array,im_array,L,this_core,left,right)
    integer(i4), intent(in) :: L, right, left, this_core
    integer, dimension(L,2), intent(out), codimension[*] :: ip_array, im_array

    
   integer(i4) :: i

    do i = 1, L
       ip_array(i,1) = i + 1
       im_array(i,1) = i - 1
       ip_array(i,2) = this_core
       im_array(i,2) = this_core
    end do
    ip_array(L,1) = 1
    ip_array(L,2) = right
    im_array(1,1) = L
    im_array(1,2) = left

  end subroutine set_periodic_bounds

  function ip(x,mu)
    integer(i4), dimension(0:3) :: ip
    integer(i4), intent(in) :: x(0:3), mu
    
    ip = x
    
    select case(mu)
    case(1)
       ip(mu) = ip1(x(mu),1)
       ip(0)  = ip1(x(mu),2)[x(0)]
    case(2)
       ip(mu) = ip2(x(mu),1)
       ip(0)  = ip2(x(mu),2)[x(0)]
    case(3)
       ip(mu) = ip3(x(mu),1)
       ip(0)  = ip3(x(mu),2)[x(0)]
    end select
    
  end function ip

  function im(x,mu)
    integer(i4), dimension(0:3) :: im
    integer(i4), intent(in) :: x(0:3), mu
    
    im = x
    
    select case(mu)
    case(1)
       im(mu) = im1(x(mu),1)
       im(0)  = im1(x(mu),2)[x(0)]
    case(2)
       im(mu) = im2(x(mu),1)
       im(0)  = im2(x(mu),2)[x(0)]
    case(3)
       im(mu) = im3(x(mu),1)
       im(0)  = im3(x(mu),2)[x(0)]
    end select
    
  end function im

  
  function im_cores(x,mu)
    integer(i4) :: im_cores(3)
    integer(i4), intent(in) :: x(3), mu
    
    im_cores = x
    
    select case(mu)
    case(1)
       im_cores(mu) = im1_c(x(mu))
    case(2)
       im_cores(mu) = im2_c(x(mu))
    case(3)
       im_cores(mu) = im3_c(x(mu))
    end select
    
  end function im_cores
#endif
end module pbc

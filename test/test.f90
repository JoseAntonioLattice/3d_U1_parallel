program test

  !use indices
  use pbc
  implicit none

  integer, parameter :: n = 10
  integer, dimension(n,2) :: ip, im

  integer :: i
  
  !print*, get_index([2,1,1],3,[2,3,1])
  !print*, get_index_array(5,3,[2,3,1])

  call set_periodic_bounds(ip,im,n,this_image(),this_image()-1,this_image()+1)

  print*, ip(n,:)
  print*, im(1,:)
  
end program test

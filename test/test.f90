program test

  use indices
  implicit none

  print*, get_index([2,1,1],3,[2,3,1])
  print*, get_index_array(5,3,[2,3,1])

end program test

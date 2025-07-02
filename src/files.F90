module files
  use check_files_directories
  implicit none

contains

  subroutine create_files(directories,file)
    character(100), dimension(:), intent(in) :: directories
    character(:), allocatable, intent(out) :: file
    character(:), allocatable :: path
    
    call check_directory(directories,path)
    call numbered_files(path,"measurements",".dat",file)
    
  end subroutine create_files

end module files

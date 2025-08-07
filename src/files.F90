module files
  use check_files_directories
  implicit none

contains

  subroutine create_files(directories,file,inCluster)
    character(100), dimension(:), intent(in) :: directories
    character(:), allocatable, intent(inout) :: file
    logical, intent(in) :: inCluster
    character(:), allocatable :: path, filename
    
    call check_directory(directories,path)

    if(inCluster) then
       file = trim(path)//trim(file)
    else
       call numbered_files(path,"measurements",".dat",file)
    end if
    
  end subroutine create_files

end module files

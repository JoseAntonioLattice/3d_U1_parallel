module check_files_directories

  use number2string
  implicit none

contains
  
  subroutine check_directory(directory,route)
    character(*), intent(in) :: directory(:)
    character(:), allocatable, intent(out), optional :: route
    logical ::  dir_exists
    character(:), allocatable :: path
    integer :: i
    
    path = ''
    do i = 1, size(directory)
       path = path //trim(directory(i))//"/"
       inquire(file=path, exist=dir_exists) ! ask wether the directory exists or not
       if(dir_exists .eqv. .false.) then
          call execute_command_line('mkdir '//path) ! If not it creates
          !print*, "Directory "//path//" created."
       else
          !print*, "Directory "//path//" already exists."
       end if
    end do
    if(present(route)) route = path
  end subroutine check_directory
  
  subroutine check_file(filepath,exists)
    character(*), intent(in) :: filepath
    logical, optional, intent(out) :: exists
    logical :: file_exists
    
    inquire(file = filepath, exist = file_exists)
    if(present(exists)) exists = file_exists
    if(file_exists .eqv. .false.)then
       call execute_command_line('touch '//filepath)
       !print*, "File "//filepath//" created."
    else
       !print*, "File "//filepath//" already exists."
    end if
  end subroutine check_file

  subroutine numbered_files(path,name,ext,file)
    character(*), intent(in) :: path, name, ext
    character(:), allocatable, intent(out), optional :: file
    integer :: i
    logical :: file_exists
    i = 0
    do
       i = i + 1
       call check_file(trim(path)//trim(name)//"_"//int2str(i)//trim(ext), file_exists)
       if(.not. file_exists) exit
    end do
    if(present(file)) file = trim(path)//trim(name)//"_"//int2str(i)//trim(ext)
  end subroutine numbered_files
  
end module check_files_directories

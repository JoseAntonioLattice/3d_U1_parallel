module save_configurations

  use iso_fortran_env, only : dp => real64, i4 => int32
  use check_files_directories
  implicit none

contains

  subroutine save_configuration(u,beta,isbeta,L)
    complex(dp), dimension(:,:,:,:), intent(in) :: u
    real(dp), intent(in) :: beta
    logical, intent(in) :: isbeta
    integer(i4), intent(in) :: L(:)
    character(:), allocatable :: tmp, path, filename
    
    tmp = 'beta'
    if(.not.isbeta) tmp = 'T'
    
    call check_directory(['data','configurations', &
         'Lx='//int2str(L(1)),'Ly='//int2str(L(2)),'Lz='//int2str(L(3)) , &
         tmp,trim(real2str(beta,1,6))],path)
    
    call numberd_files(path,"conf",".bin",filename) 
    open(newunit = outunit, file = filename, form = "unformatted", access = 'sequential')
    write(outunit) u
    close(outunit)
  end subroutine save_configuration

  subroutine read_configuration(u,filename)
    complex(dp), dimension(:,:,:,:), intent(in) :: u
    character(*), intent(in) :: filename
    
    open(newunit = outunit, file = filename, form = "unformatted", access = 'sequential')
    read(outunit) u
    close(outunit)
  end subroutine read_configuration

  
end module save_configurations

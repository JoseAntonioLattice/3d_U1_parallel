module statistics

  use iso_fortran_env, only : dp => real64, i4 => int32
  implicit none
  
contains

  function avr(x)
    real(dp) :: avr
    real(dp), intent(in) :: x(:)

    avr = sum(x)/size(x)

  end function avr

  function var(x)
    real(dp) :: var
    real(dp), intent(in) :: x(:)

    var = sum((x-avr(x))**2)/(size(x)-1)

  end function var


  function std_err(x)
    real(dp) :: std_Err
    real(dp), intent(in) :: x(:)

    std_Err = sqrt(var(x)/size(x))

  end function std_err

  function jackknife(x,bins)
    real(dp) :: jackknife
    real(dp), intent(in) :: x(:)
    integer(i4), intent(in) :: bins

    integer(i4) :: MM, NN, i
    real(dp) :: xbar, sum_x
    real(dp) :: x_m(bins)
    
    
    NN = size(x)
    MM = NN/bins
    
    xbar = avr(x)
    sum_x = sum(x)
    x_m = 1.0_dp/(NN-MM) * [(sum_x - sum(x(MM*(i-1)+1:MM*i)),i=1,bins)]

    jackknife = sqrt( real(bins - 1,dp)/bins * sum( (x_m - xbar)**2) )
    
  end function jackknife

  function jackknife_max(x)
    real(dp) :: jackknife_max(2)
    real(dp), intent(in) :: x(:)
    
    integer(i4) :: i
    real(dp) :: err(size(x))

    err = 0
    do i = 2, size(x)
       if( mod(size(x),i) /= 0) cycle
       err(i) = jackknife(x,i) 
    end do
    jackknife_max(1) = maxval(err)
    jackknife_max(2) = maxloc(err,dim=1)
  end function jackknife_max
  
end module statistics

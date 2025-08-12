module random

  use iso_fortran_env, only : dp => real64, i4 => int32
  use constants, only : pi

  implicit none

  private
  public :: uniform, normal, exponential
  public :: random_choice
  
contains

  function random_choice(n)
    integer(i4) :: random_choice
    integer(i4), intent(in) :: n
    real(dp) :: r

    call random_number(r)
    random_choice = floor(n*r) + 1
    
  end function random_choice
  
  function uniform(a,b)
    real(dp), optional :: a, b
    real(dp) :: uniform, r, a_, b_

    
    if(.not. present(a) .and. .not. present(b))then
       a_ = 0.0_dp
       b_ = 1.0_dp
    else if(present(a) .and. present(b)) then
       a_ = a
       b_ = b
    end if

    call random_number(r)
    uniform = (b_-a_) * r  + a_

  end function uniform

  function normal(mu, sigma)
    real(dp), optional :: mu, sigma
    real(dp) :: normal, u1, u2, mu_, sigma_

    call random_number(u1)
    call random_number(u2)

    if(.not. present(mu) .and. .not. present(sigma))then
       mu_ = 0.0_dp
       sigma_ = 1.0_dp
    else if(present(mu) .and. present(sigma))then
       mu_ = mu
       sigma_ = sigma
    end if
    normal = sigma_ * sqrt(-2*log(1.0_dp-u1))*cos(2*pi*u2) + mu_
    
  end function normal

  function exponential(lambda)
    real(dp), optional :: lambda
    real(dp) :: exponential, u, lambda_

    if(.not. present(lambda))then
       lambda_ = 1.0_dp
    else
       lambda_ = lambda
    end if

    call random_number(u)
    exponential = -log(1.0_dp-u)/lambda_
    
  end function exponential

end module random

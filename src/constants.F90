module constants

  use iso_fortran_env, only : i4 => int32, dp => real64
  implicit none

  !=========================================
  !             math constants
  !=========================================
  real(dp), parameter :: pi = acos(-1.0_dp)
  real(dp), parameter :: twopi = 2*pi
  real(dp), parameter :: pi2 = pi**2
  real(dp), parameter :: e = exp(1.0_dp)
  real(dp), parameter :: euler_mascheroni = 0.57721566490153286060_dp

  complex(dp), parameter :: ii = (0.0_dp,1.0_dp)
  
  !=========================================
  !                 Physics                  
  !=========================================

  ! speeed of light
  real(dp), parameter :: c = 299792458.0_dp ! m/s

  ! Gravitational constant
  real(dp), parameter :: G_newton = 6.67430e-11_dp ! N m^2 kg^-2

  ! Planck's constant h. Reduced Planck's constant hbar
  real(dp), parameter :: h = 6.62607015e-34_dp          ! J s
  real(dp), parameter :: h_MeVs = 4.135667696e-21_dp    ! MeV s
  real(dp), parameter :: hbar = 1.054571817e-34_dp      ! J s
  real(dp), parameter :: hbar_MeVs = 6.582119569e-22_dp ! MeV s

  ! electron-volts
  real(dp), parameter :: eV = 1.602176463 ! J

  ! Natural units
  ! =============
  ! hbar * c
  real(dp), parameter :: hbarC_MeVfm = 197.32698044_dp  ! MeV fm
  
end module constants

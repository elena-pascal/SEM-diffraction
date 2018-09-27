!******************************************************************************
!*                       Code generated with sympy 1.0                        *
!*                                                                            *
!*              See http://www.sympy.org/ for more information.               *
!*                                                                            *
!*                This file is part of 'Screwbeta_generation'                 *
!******************************************************************************

subroutine beta_ECCI(x, y, z, betaH)
implicit none
REAL*8, intent(in) :: x
REAL*8, intent(in) :: y
REAL*8, intent(in) :: z
REAL*8, intent(out) :: betaH

betaH = -0.318309886183791d0*(-5.77659207740314d0*(1.22464679914735d-16* &
      y - 1.0d0*z)*((1.22464679914735d-16*y - 1.0d0*z)**2 + ( &
      -0.880063298291132d0*x - 0.474856389870595d0*y - &
      5.81531357909691d-17*z)**2 + (0.474856389870595d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)**2)**(-0.5d0) + &
      5.77659207740314d0)*(0.474856389870595d0*x - 0.880063298291132d0* &
      y - 1.0777667012993d-16*z)/(-1.22464679914735d-16*y + 1.0d0*z + &
      sqrt((1.22464679914735d-16*y - 1.0d0*z)**2 + ( &
      -0.880063298291132d0*x - 0.474856389870595d0*y - &
      5.81531357909691d-17*z)**2 + (0.474856389870595d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)**2))**2 + &
      0.318309886183791d0*(-2.63298709806897d0*(1.22464679914735d-16*y &
      - 1.0d0*z)*((1.22464679914735d-16*y - 1.0d0*z)**2 + ( &
      -0.880063298291132d0*x - 0.474856389870595d0*y - &
      5.81531357909691d-17*z)**2 + (0.474856389870595d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)**2)**(-0.5d0) + &
      2.63298709806897d0)*(-0.880063298291132d0*x - 0.474856389870595d0 &
      *y - 5.81531357909691d-17*z)/(-1.22464679914735d-16*y + 1.0d0*z + &
      sqrt((1.22464679914735d-16*y - 1.0d0*z)**2 + ( &
      -0.880063298291132d0*x - 0.474856389870595d0*y - &
      5.81531357909691d-17*z)**2 + (0.474856389870595d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)**2))**2 + &
      0.467624807244861d0*(-0.880063298291132d0*x - 0.474856389870595d0 &
      *y - 5.81531357909691d-17*z)**2*((1.22464679914735d-16*y - 1.0d0* &
      z)**2 + (-0.880063298291132d0*x - 0.474856389870595d0*y - &
      5.81531357909691d-17*z)**2 + (0.474856389870595d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)**2)**(-0.5d0)/( &
      -1.22464679914735d-16*y + 1.0d0*z + sqrt((1.22464679914735d-16*y &
      - 1.0d0*z)**2 + (-0.880063298291132d0*x - 0.474856389870595d0*y - &
      5.81531357909691d-17*z)**2 + (0.474856389870595d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)**2))**2 - &
      0.159275747491687d0*(-0.880063298291132d0*x - 0.474856389870595d0 &
      *y - 5.81531357909691d-17*z)*(0.474856389870595d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)*(( &
      1.22464679914735d-16*y - 1.0d0*z)**2 + (-0.880063298291132d0*x - &
      0.474856389870595d0*y - 5.81531357909691d-17*z)**2 + ( &
      0.474856389870595d0*x - 0.880063298291132d0*y - &
      1.0777667012993d-16*z)**2)**(-0.5d0)/(-1.22464679914735d-16*y + &
      1.0d0*z + sqrt((1.22464679914735d-16*y - 1.0d0*z)**2 + ( &
      -0.880063298291132d0*x - 0.474856389870595d0*y - &
      5.81531357909691d-17*z)**2 + (0.474856389870595d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)**2))**2 - &
      1.90139406858438d0*(0.474856389870595d0*x - 0.880063298291132d0*y &
      - 1.0777667012993d-16*z)**2*((1.22464679914735d-16*y - 1.0d0*z)** &
      2 + (-0.880063298291132d0*x - 0.474856389870595d0*y - &
      5.81531357909691d-17*z)**2 + (0.474856389870595d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)**2)**(-0.5d0)/( &
      -1.22464679914735d-16*y + 1.0d0*z + sqrt((1.22464679914735d-16*y &
      - 1.0d0*z)**2 + (-0.880063298291132d0*x - 0.474856389870595d0*y - &
      5.81531357909691d-17*z)**2 + (0.474856389870595d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)**2))**2 + &
      4.50431897836711d0/(-3.84734138744357d-16*y + 3.14159265358979d0* &
      z + 3.14159265358979d0*sqrt((1.22464679914735d-16*y - 1.0d0*z)**2 &
      + (-0.880063298291132d0*x - 0.474856389870595d0*y - &
      5.81531357909691d-17*z)**2 + (0.474856389870595d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)**2)) - ( &
      -0.884730161365563d0*x - 0.477374492552354d0*y - &
      5.84615144298831d-17*z)/((3.14159265358979d0*( &
      -0.880063298291132d0*x - 0.474856389870595d0*y - &
      5.81531357909691d-17*z)**2/(0.474856389870595d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)**2 + &
      3.14159265358979d0)*(0.474856389870595d0*x - 0.880063298291132d0* &
      y - 1.0777667012993d-16*z)**2) + 0.542431997197586d0/(( &
      3.14159265358979d0*(-0.880063298291132d0*x - 0.474856389870595d0* &
      y - 5.81531357909691d-17*z)**2/(0.474856389870595d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)**2 + &
      3.14159265358979d0)*(0.474856389870595d0*x - 0.880063298291132d0* &
      y - 1.0777667012993d-16*z))

end subroutine

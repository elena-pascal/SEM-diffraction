!******************************************************************************
!*                       Code generated with sympy 1.0                        *
!*                                                                            *
!*              See http://www.sympy.org/ for more information.               *
!*                                                                            *
!*                This file is part of 'mixedbeta_generation'                 *
!******************************************************************************

subroutine beta_ECCI(x, y, z, betaH)
implicit none
REAL*8, intent(in) :: x
REAL*8, intent(in) :: y
REAL*8, intent(in) :: z
REAL*8, intent(out) :: betaH

betaH = -0.309454457911282d0*(1.22464679914735d-16*y - 1.0d0*z)**2*( &
      -0.474856389870595d0*x + 0.880063298291132d0*y + &
      1.0777667012993d-16*z)*(0.880063298291132d0*x + &
      0.474856389870595d0*y + 5.81531357909692d-17*z)*(( &
      1.22464679914735d-16*y - 1.0d0*z)**2 + (-0.474856389870595d0*x + &
      0.880063298291132d0*y + 1.0777667012993d-16*z)**2 + ( &
      0.880063298291132d0*x + 0.474856389870595d0*y + &
      5.81531357909692d-17*z)**2)**(-1.5d0)/(-1.22464679914735d-16*y + &
      1.0d0*z + sqrt((1.22464679914735d-16*y - 1.0d0*z)**2 + ( &
      -0.474856389870595d0*x + 0.880063298291132d0*y + &
      1.0777667012993d-16*z)**2 + (0.880063298291132d0*x + &
      0.474856389870595d0*y + 5.81531357909692d-17*z)**2))**2 - &
      0.21694867484351d0*(1.22464679914735d-16*y - 1.0d0*z)*(-2.4d0*( &
      1.22464679914735d-16*y - 1.0d0*z)*((1.22464679914735d-16*y - &
      1.0d0*z)**2 + (-0.474856389870595d0*x + 0.880063298291132d0*y + &
      1.0777667012993d-16*z)**2 + (0.880063298291132d0*x + &
      0.474856389870595d0*y + 5.81531357909692d-17*z)**2)**(-0.5d0) + &
      2.4d0)/(-1.22464679914735d-16*y + 1.0d0*z + sqrt(( &
      1.22464679914735d-16*y - 1.0d0*z)**2 + (-0.474856389870595d0*x + &
      0.880063298291132d0*y + 1.0777667012993d-16*z)**2 + ( &
      0.880063298291132d0*x + 0.474856389870595d0*y + &
      5.81531357909692d-17*z)**2))**2 + 0.309454457911282d0*( &
      1.22464679914735d-16*y - 1.0d0*z)*(-2.0d0*(1.22464679914735d-16*y &
      - 1.0d0*z)*((1.22464679914735d-16*y - 1.0d0*z)**2 + ( &
      -0.474856389870595d0*x + 0.880063298291132d0*y + &
      1.0777667012993d-16*z)**2 + (0.880063298291132d0*x + &
      0.474856389870595d0*y + 5.81531357909692d-17*z)**2)**(-0.5d0) + 2 &
      )*(-0.474856389870595d0*x + 0.880063298291132d0*y + &
      1.0777667012993d-16*z)*(0.880063298291132d0*x + &
      0.474856389870595d0*y + 5.81531357909692d-17*z)*(( &
      1.22464679914735d-16*y - 1.0d0*z)**2 + (-0.474856389870595d0*x + &
      0.880063298291132d0*y + 1.0777667012993d-16*z)**2 + ( &
      0.880063298291132d0*x + 0.474856389870595d0*y + &
      5.81531357909692d-17*z)**2)**(-0.5d0)/(-1.22464679914735d-16*y + &
      1.0d0*z + sqrt((1.22464679914735d-16*y - 1.0d0*z)**2 + ( &
      -0.474856389870595d0*x + 0.880063298291132d0*y + &
      1.0777667012993d-16*z)**2 + (0.880063298291132d0*x + &
      0.474856389870595d0*y + 5.81531357909692d-17*z)**2))**3 - &
      0.319997842785255d0*(1.22464679914735d-16*y - 1.0d0*z)*( &
      -0.474856389870595d0*x + 0.880063298291132d0*y + &
      1.0777667012993d-16*z)**2*(0.880063298291132d0*x + &
      0.474856389870595d0*y + 5.81531357909692d-17*z)*(( &
      1.22464679914735d-16*y - 1.0d0*z)**2 + (-0.474856389870595d0*x + &
      0.880063298291132d0*y + 1.0777667012993d-16*z)**2 + ( &
      0.880063298291132d0*x + 0.474856389870595d0*y + &
      5.81531357909692d-17*z)**2)**(-1.5d0)/(-1.22464679914735d-16*y + &
      1.0d0*z + sqrt((1.22464679914735d-16*y - 1.0d0*z)**2 + ( &
      -0.474856389870595d0*x + 0.880063298291132d0*y + &
      1.0777667012993d-16*z)**2 + (0.880063298291132d0*x + &
      0.474856389870595d0*y + 5.81531357909692d-17*z)**2))**2 - &
      0.639995685570511d0*(1.22464679914735d-16*y - 1.0d0*z)*( &
      -0.474856389870595d0*x + 0.880063298291132d0*y + &
      1.0777667012993d-16*z)**2*(0.880063298291132d0*x + &
      0.474856389870595d0*y + 5.81531357909692d-17*z)*1.0/(( &
      1.22464679914735d-16*y - 1.0d0*z)**2 + (-0.474856389870595d0*x + &
      0.880063298291132d0*y + 1.0777667012993d-16*z)**2 + ( &
      0.880063298291132d0*x + 0.474856389870595d0*y + &
      5.81531357909692d-17*z)**2)/(-1.22464679914735d-16*y + 1.0d0*z + &
      sqrt((1.22464679914735d-16*y - 1.0d0*z)**2 + ( &
      -0.474856389870595d0*x + 0.880063298291132d0*y + &
      1.0777667012993d-16*z)**2 + (0.880063298291132d0*x + &
      0.474856389870595d0*y + 5.81531357909692d-17*z)**2))**3 - &
      0.17266146729041d0*(1.22464679914735d-16*y - 1.0d0*z)*( &
      -0.474856389870595d0*x + 0.880063298291132d0*y + &
      1.0777667012993d-16*z)*(0.880063298291132d0*x + &
      0.474856389870595d0*y + 5.81531357909692d-17*z)**2*(( &
      1.22464679914735d-16*y - 1.0d0*z)**2 + (-0.474856389870595d0*x + &
      0.880063298291132d0*y + 1.0777667012993d-16*z)**2 + ( &
      0.880063298291132d0*x + 0.474856389870595d0*y + &
      5.81531357909692d-17*z)**2)**(-1.5d0)/(-1.22464679914735d-16*y + &
      1.0d0*z + sqrt((1.22464679914735d-16*y - 1.0d0*z)**2 + ( &
      -0.474856389870595d0*x + 0.880063298291132d0*y + &
      1.0777667012993d-16*z)**2 + (0.880063298291132d0*x + &
      0.474856389870595d0*y + 5.81531357909692d-17*z)**2))**2 - &
      0.34532293458082d0*(1.22464679914735d-16*y - 1.0d0*z)*( &
      -0.474856389870595d0*x + 0.880063298291132d0*y + &
      1.0777667012993d-16*z)*(0.880063298291132d0*x + &
      0.474856389870595d0*y + 5.81531357909692d-17*z)**2*1.0/(( &
      1.22464679914735d-16*y - 1.0d0*z)**2 + (-0.474856389870595d0*x + &
      0.880063298291132d0*y + 1.0777667012993d-16*z)**2 + ( &
      0.880063298291132d0*x + 0.474856389870595d0*y + &
      5.81531357909692d-17*z)**2)/(-1.22464679914735d-16*y + 1.0d0*z + &
      sqrt((1.22464679914735d-16*y - 1.0d0*z)**2 + ( &
      -0.474856389870595d0*x + 0.880063298291132d0*y + &
      1.0777667012993d-16*z)**2 + (0.880063298291132d0*x + &
      0.474856389870595d0*y + 5.81531357909692d-17*z)**2))**3 + &
      0.711078202963581d0*(1.22464679914735d-16*y - 1.0d0*z)*( &
      -0.474856389870595d0*x + 0.880063298291132d0*y + &
      1.0777667012993d-16*z)*((1.22464679914735d-16*y - 1.0d0*z)**2 + ( &
      -0.474856389870595d0*x + 0.880063298291132d0*y + &
      1.0777667012993d-16*z)**2 + (0.880063298291132d0*x + &
      0.474856389870595d0*y + 5.81531357909692d-17*z)**2)**(-0.5d0)/( &
      -1.22464679914735d-16*y + 1.0d0*z + sqrt((1.22464679914735d-16*y &
      - 1.0d0*z)**2 + (-0.474856389870595d0*x + 0.880063298291132d0*y + &
      1.0777667012993d-16*z)**2 + (0.880063298291132d0*x + &
      0.474856389870595d0*y + 5.81531357909692d-17*z)**2))**2 + &
      0.433897349687021d0*(1.22464679914735d-16*y - 1.0d0*z)*( &
      0.880063298291132d0*x + 0.474856389870595d0*y + &
      5.81531357909692d-17*z)**2*((1.22464679914735d-16*y - 1.0d0*z)**2 &
      + (-0.474856389870595d0*x + 0.880063298291132d0*y + &
      1.0777667012993d-16*z)**2 + (0.880063298291132d0*x + &
      0.474856389870595d0*y + 5.81531357909692d-17*z)**2)**(-1.5d0)/( &
      -1.22464679914735d-16*y + 1.0d0*z + sqrt((1.22464679914735d-16*y &
      - 1.0d0*z)**2 + (-0.474856389870595d0*x + 0.880063298291132d0*y + &
      1.0777667012993d-16*z)**2 + (0.880063298291132d0*x + &
      0.474856389870595d0*y + 5.81531357909692d-17*z)**2)) + &
      0.610511749960027d0*(1.22464679914735d-16*y - 1.0d0*z)*( &
      0.880063298291132d0*x + 0.474856389870595d0*y + &
      5.81531357909692d-17*z)*((1.22464679914735d-16*y - 1.0d0*z)**2 + &
      (-0.474856389870595d0*x + 0.880063298291132d0*y + &
      1.0777667012993d-16*z)**2 + (0.880063298291132d0*x + &
      0.474856389870595d0*y + 5.81531357909692d-17*z)**2)**(-0.5d0)/( &
      -1.22464679914735d-16*y + 1.0d0*z + sqrt((1.22464679914735d-16*y &
      - 1.0d0*z)**2 + (-0.474856389870595d0*x + 0.880063298291132d0*y + &
      1.0777667012993d-16*z)**2 + (0.880063298291132d0*x + &
      0.474856389870595d0*y + 5.81531357909692d-17*z)**2))**2 + &
      0.318309886183791d0*(-5.16852975346597d0*(1.22464679914735d-16*y &
      - 1.0d0*z)*((1.22464679914735d-16*y - 1.0d0*z)**2 + ( &
      -0.880063298291132d0*x - 0.474856389870595d0*y - &
      5.81531357909691d-17*z)**2 + (0.474856389870595d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)**2)**(-0.5d0) + &
      5.16852975346597d0)*(0.474856389870595d0*x - 0.880063298291132d0* &
      y - 1.0777667012993d-16*z)/(-1.22464679914735d-16*y + 1.0d0*z + &
      sqrt((1.22464679914735d-16*y - 1.0d0*z)**2 + ( &
      -0.880063298291132d0*x - 0.474856389870595d0*y - &
      5.81531357909691d-17*z)**2 + (0.474856389870595d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)**2))**2 - &
      0.318309886183791d0*(-3.68618193729656d0*(1.22464679914735d-16*y &
      - 1.0d0*z)*((1.22464679914735d-16*y - 1.0d0*z)**2 + ( &
      -0.880063298291132d0*x - 0.474856389870595d0*y - &
      5.81531357909691d-17*z)**2 + (0.474856389870595d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)**2)**(-0.5d0) + &
      3.68618193729656d0)*(-0.880063298291132d0*x - 0.474856389870595d0 &
      *y - 5.81531357909691d-17*z)/(-1.22464679914735d-16*y + 1.0d0*z + &
      sqrt((1.22464679914735d-16*y - 1.0d0*z)**2 + ( &
      -0.880063298291132d0*x - 0.474856389870595d0*y - &
      5.81531357909691d-17*z)**2 + (0.474856389870595d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)**2))**2 + &
      0.318309886183791d0*(-0.256399081739261d0*(1.22464679914735d-16*y &
      - 1.0d0*z)*((1.22464679914735d-16*y - 1.0d0*z)**2 + ( &
      -0.474856389870595d0*x + 0.880063298291132d0*y + &
      1.0777667012993d-16*z)**2 + (0.880063298291132d0*x + &
      0.474856389870595d0*y + 5.81531357909692d-17*z)**2)**(-1.5d0) + &
      0.256399081739261d0*(-0.4d0*(1.22464679914735d-16*y - 1.0d0*z)*(( &
      1.22464679914735d-16*y - 1.0d0*z)**2 + (-0.474856389870595d0*x + &
      0.880063298291132d0*y + 1.0777667012993d-16*z)**2 + ( &
      0.880063298291132d0*x + 0.474856389870595d0*y + &
      5.81531357909692d-17*z)**2)**(-0.5d0) + 0.4d0)/( &
      -1.22464679914735d-16*y + 1.0d0*z + sqrt((1.22464679914735d-16*y &
      - 1.0d0*z)**2 + (-0.474856389870595d0*x + 0.880063298291132d0*y + &
      1.0777667012993d-16*z)**2 + (0.880063298291132d0*x + &
      0.474856389870595d0*y + 5.81531357909692d-17*z)**2))**2)*( &
      0.880063298291132d0*x + 0.474856389870595d0*y + &
      5.81531357909692d-17*z) + 0.520676819624425d0*(-2.0d0*( &
      1.22464679914735d-16*y - 1.0d0*z)*((1.22464679914735d-16*y - &
      1.0d0*z)**2 + (-0.474856389870595d0*x + 0.880063298291132d0*y + &
      1.0777667012993d-16*z)**2 + (0.880063298291132d0*x + &
      0.474856389870595d0*y + 5.81531357909692d-17*z)**2)**(-0.5d0) + 2 &
      )*(0.880063298291132d0*x + 0.474856389870595d0*y + &
      5.81531357909692d-17*z)**2/(-1.22464679914735d-16*y + 1.0d0*z + &
      sqrt((1.22464679914735d-16*y - 1.0d0*z)**2 + ( &
      -0.474856389870595d0*x + 0.880063298291132d0*y + &
      1.0777667012993d-16*z)**2 + (0.880063298291132d0*x + &
      0.474856389870595d0*y + 5.81531357909692d-17*z)**2))**3 - &
      0.433897349687021d0*(-1.0d0*(1.22464679914735d-16*y - 1.0d0*z)*(( &
      1.22464679914735d-16*y - 1.0d0*z)**2 + (-0.474856389870595d0*x + &
      0.880063298291132d0*y + 1.0777667012993d-16*z)**2 + ( &
      0.880063298291132d0*x + 0.474856389870595d0*y + &
      5.81531357909692d-17*z)**2)**(-0.5d0) + 1)*(0.880063298291132d0*x &
      + 0.474856389870595d0*y + 5.81531357909692d-17*z)**2*(( &
      1.22464679914735d-16*y - 1.0d0*z)**2 + (-0.474856389870595d0*x + &
      0.880063298291132d0*y + 1.0777667012993d-16*z)**2 + ( &
      0.880063298291132d0*x + 0.474856389870595d0*y + &
      5.81531357909692d-17*z)**2)**(-0.5d0)/(-1.22464679914735d-16*y + &
      1.0d0*z + sqrt((1.22464679914735d-16*y - 1.0d0*z)**2 + ( &
      -0.474856389870595d0*x + 0.880063298291132d0*y + &
      1.0777667012993d-16*z)**2 + (0.880063298291132d0*x + &
      0.474856389870595d0*y + 5.81531357909692d-17*z)**2))**2 + &
      0.154727228955641d0*(-0.8d0*(1.22464679914735d-16*y - 1.0d0*z)*(( &
      1.22464679914735d-16*y - 1.0d0*z)**2 + (-0.474856389870595d0*x + &
      0.880063298291132d0*y + 1.0777667012993d-16*z)**2 + ( &
      0.880063298291132d0*x + 0.474856389870595d0*y + &
      5.81531357909692d-17*z)**2)**(-0.5d0) + 0.8d0)*( &
      -0.474856389870595d0*x + 0.880063298291132d0*y + &
      1.0777667012993d-16*z)*(0.880063298291132d0*x + &
      0.474856389870595d0*y + 5.81531357909692d-17*z)/( &
      -1.22464679914735d-16*y + 1.0d0*z + sqrt((1.22464679914735d-16*y &
      - 1.0d0*z)**2 + (-0.474856389870595d0*x + 0.880063298291132d0*y + &
      1.0777667012993d-16*z)**2 + (0.880063298291132d0*x + &
      0.474856389870595d0*y + 5.81531357909692d-17*z)**2))**3 + &
      0.21694867484351d0*(0.4d0*(1.22464679914735d-16*y - 1.0d0*z)*(( &
      1.22464679914735d-16*y - 1.0d0*z)**2 + (-0.474856389870595d0*x + &
      0.880063298291132d0*y + 1.0777667012993d-16*z)**2 + ( &
      0.880063298291132d0*x + 0.474856389870595d0*y + &
      5.81531357909692d-17*z)**2)**(-0.5d0) - 0.4d0)/( &
      -1.22464679914735d-16*y + 1.0d0*z + sqrt((1.22464679914735d-16*y &
      - 1.0d0*z)**2 + (-0.474856389870595d0*x + 0.880063298291132d0*y + &
      1.0777667012993d-16*z)**2 + (0.880063298291132d0*x + &
      0.474856389870595d0*y + 5.81531357909692d-17*z)**2)) + &
      0.318309886183791d0*(-0.265134823400107d0*(-0.474856389870595d0*x &
      + 0.880063298291132d0*y + 1.0777667012993d-16*z)*(( &
      1.22464679914735d-16*y - 1.0d0*z)**2 + (-0.474856389870595d0*x + &
      0.880063298291132d0*y + 1.0777667012993d-16*z)**2 + ( &
      0.880063298291132d0*x + 0.474856389870595d0*y + &
      5.81531357909692d-17*z)**2)**(-1.5d0) - 0.106053929360043d0*( &
      -0.474856389870595d0*x + 0.880063298291132d0*y + &
      1.0777667012993d-16*z)*((1.22464679914735d-16*y - 1.0d0*z)**2 + ( &
      -0.474856389870595d0*x + 0.880063298291132d0*y + &
      1.0777667012993d-16*z)**2 + (0.880063298291132d0*x + &
      0.474856389870595d0*y + 5.81531357909692d-17*z)**2)**(-0.5d0)/( &
      -1.22464679914735d-16*y + 1.0d0*z + sqrt((1.22464679914735d-16*y &
      - 1.0d0*z)**2 + (-0.474856389870595d0*x + 0.880063298291132d0*y + &
      1.0777667012993d-16*z)**2 + (0.880063298291132d0*x + &
      0.474856389870595d0*y + 5.81531357909692d-17*z)**2))**2)*( &
      0.880063298291132d0*x + 0.474856389870595d0*y + &
      5.81531357909692d-17*z) + 0.318309886183791d0*( &
      -0.143058988271891d0*(0.880063298291132d0*x + 0.474856389870595d0 &
      *y + 5.81531357909692d-17*z)*((1.22464679914735d-16*y - 1.0d0*z) &
      **2 + (-0.474856389870595d0*x + 0.880063298291132d0*y + &
      1.0777667012993d-16*z)**2 + (0.880063298291132d0*x + &
      0.474856389870595d0*y + 5.81531357909692d-17*z)**2)**(-1.5d0) - &
      0.0572235953087564d0*(0.880063298291132d0*x + 0.474856389870595d0 &
      *y + 5.81531357909692d-17*z)*((1.22464679914735d-16*y - 1.0d0*z) &
      **2 + (-0.474856389870595d0*x + 0.880063298291132d0*y + &
      1.0777667012993d-16*z)**2 + (0.880063298291132d0*x + &
      0.474856389870595d0*y + 5.81531357909692d-17*z)**2)**(-0.5d0)/( &
      -1.22464679914735d-16*y + 1.0d0*z + sqrt((1.22464679914735d-16*y &
      - 1.0d0*z)**2 + (-0.474856389870595d0*x + 0.880063298291132d0*y + &
      1.0777667012993d-16*z)**2 + (0.880063298291132d0*x + &
      0.474856389870595d0*y + 5.81531357909692d-17*z)**2))**2)*( &
      0.880063298291132d0*x + 0.474856389870595d0*y + &
      5.81531357909692d-17*z) - 0.654674730142803d0*( &
      -0.880063298291132d0*x - 0.474856389870595d0*y - &
      5.81531357909691d-17*z)**2*((1.22464679914735d-16*y - 1.0d0*z)**2 &
      + (-0.880063298291132d0*x - 0.474856389870595d0*y - &
      5.81531357909691d-17*z)**2 + (0.474856389870595d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)**2)**(-0.5d0)/( &
      -1.22464679914735d-16*y + 1.0d0*z + sqrt((1.22464679914735d-16*y &
      - 1.0d0*z)**2 + (-0.880063298291132d0*x - 0.474856389870595d0*y - &
      5.81531357909691d-17*z)**2 + (0.474856389870595d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)**2))**2 - &
      0.295381905529366d0*(-0.880063298291132d0*x - 0.474856389870595d0 &
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
      0.127999137114102d0*(-0.474856389870595d0*x + 0.880063298291132d0 &
      *y + 1.0777667012993d-16*z)**2*(0.880063298291132d0*x + &
      0.474856389870595d0*y + 5.81531357909692d-17*z)*(( &
      1.22464679914735d-16*y - 1.0d0*z)**2 + (-0.474856389870595d0*x + &
      0.880063298291132d0*y + 1.0777667012993d-16*z)**2 + ( &
      0.880063298291132d0*x + 0.474856389870595d0*y + &
      5.81531357909692d-17*z)**2)**(-0.5d0)/(-1.22464679914735d-16*y + &
      1.0d0*z + sqrt((1.22464679914735d-16*y - 1.0d0*z)**2 + ( &
      -0.474856389870595d0*x + 0.880063298291132d0*y + &
      1.0777667012993d-16*z)**2 + (0.880063298291132d0*x + &
      0.474856389870595d0*y + 5.81531357909692d-17*z)**2))**3 + &
      0.448680613060976d0*(-0.474856389870595d0*x + 0.880063298291132d0 &
      *y + 1.0777667012993d-16*z)*(0.880063298291132d0*x + &
      0.474856389870595d0*y + 5.81531357909692d-17*z)**2*(( &
      1.22464679914735d-16*y - 1.0d0*z)**2 + (-0.474856389870595d0*x + &
      0.880063298291132d0*y + 1.0777667012993d-16*z)**2 + ( &
      0.880063298291132d0*x + 0.474856389870595d0*y + &
      5.81531357909692d-17*z)**2)**(-1.5d0)/(-1.22464679914735d-16*y + &
      1.0d0*z + sqrt((1.22464679914735d-16*y - 1.0d0*z)**2 + ( &
      -0.474856389870595d0*x + 0.880063298291132d0*y + &
      1.0777667012993d-16*z)**2 + (0.880063298291132d0*x + &
      0.474856389870595d0*y + 5.81531357909692d-17*z)**2)) + &
      0.448680613060976d0*(-0.474856389870595d0*x + 0.880063298291132d0 &
      *y + 1.0777667012993d-16*z)*(0.880063298291132d0*x + &
      0.474856389870595d0*y + 5.81531357909692d-17*z)**2*1.0/(( &
      1.22464679914735d-16*y - 1.0d0*z)**2 + (-0.474856389870595d0*x + &
      0.880063298291132d0*y + 1.0777667012993d-16*z)**2 + ( &
      0.880063298291132d0*x + 0.474856389870595d0*y + &
      5.81531357909692d-17*z)**2)/(-1.22464679914735d-16*y + 1.0d0*z + &
      sqrt((1.22464679914735d-16*y - 1.0d0*z)**2 + ( &
      -0.474856389870595d0*x + 0.880063298291132d0*y + &
      1.0777667012993d-16*z)**2 + (0.880063298291132d0*x + &
      0.474856389870595d0*y + 5.81531357909692d-17*z)**2))**2 - &
      1.14589805826251d0*(-0.474856389870595d0*x + 0.880063298291132d0* &
      y + 1.0777667012993d-16*z)*(0.880063298291132d0*x + &
      0.474856389870595d0*y + 5.81531357909692d-17*z)**2*(( &
      1.22464679914735d-16*y - 1.0d0*z)**2 + (-0.474856389870595d0*x + &
      0.880063298291132d0*y + 1.0777667012993d-16*z)**2 + ( &
      0.880063298291132d0*x + 0.474856389870595d0*y + &
      5.81531357909692d-17*z)**2)**(-0.5d0)/(-1.22464679914735d-16*y + &
      1.0d0*z + sqrt((1.22464679914735d-16*y - 1.0d0*z)**2 + ( &
      -0.474856389870595d0*x + 0.880063298291132d0*y + &
      1.0777667012993d-16*z)**2 + (0.880063298291132d0*x + &
      0.474856389870595d0*y + 5.81531357909692d-17*z)**2))**3 + &
      0.309454457911282d0*(-0.474856389870595d0*x + 0.880063298291132d0 &
      *y + 1.0777667012993d-16*z)*(0.880063298291132d0*x + &
      0.474856389870595d0*y + 5.81531357909692d-17*z)*(( &
      1.22464679914735d-16*y - 1.0d0*z)**2 + (-0.474856389870595d0*x + &
      0.880063298291132d0*y + 1.0777667012993d-16*z)**2 + ( &
      0.880063298291132d0*x + 0.474856389870595d0*y + &
      5.81531357909692d-17*z)**2)**(-0.5d0)/(-1.22464679914735d-16*y + &
      1.0d0*z + sqrt((1.22464679914735d-16*y - 1.0d0*z)**2 + ( &
      -0.474856389870595d0*x + 0.880063298291132d0*y + &
      1.0777667012993d-16*z)**2 + (0.880063298291132d0*x + &
      0.474856389870595d0*y + 5.81531357909692d-17*z)**2))**2 + &
      0.0897361226121953d0*(-0.474856389870595d0*x + &
      0.880063298291132d0*y + 1.0777667012993d-16*z)*(( &
      1.22464679914735d-16*y - 1.0d0*z)**2 + (-0.474856389870595d0*x + &
      0.880063298291132d0*y + 1.0777667012993d-16*z)**2 + ( &
      0.880063298291132d0*x + 0.474856389870595d0*y + &
      5.81531357909692d-17*z)**2)**(-0.5d0)/(-1.22464679914735d-16*y + &
      1.0d0*z + sqrt((1.22464679914735d-16*y - 1.0d0*z)**2 + ( &
      -0.474856389870595d0*x + 0.880063298291132d0*y + &
      1.0777667012993d-16*z)**2 + (0.880063298291132d0*x + &
      0.474856389870595d0*y + 5.81531357909692d-17*z)**2)) + &
      0.0863307336452049d0*(-0.189942555948238d0*x + &
      0.352025319316453d0*y + 4.3110668051972d-17*z)/( &
      -1.22464679914735d-16*y + 1.0d0*z + sqrt((1.22464679914735d-16*y &
      - 1.0d0*z)**2 + (-0.474856389870595d0*x + 0.880063298291132d0*y + &
      1.0777667012993d-16*z)**2 + (0.880063298291132d0*x + &
      0.474856389870595d0*y + 5.81531357909692d-17*z)**2))**2 + &
      0.159998921392628d0*(0.352025319316453d0*x + 0.189942555948238d0* &
      y + 2.32612543163877d-17*z)/(-1.22464679914735d-16*y + 1.0d0*z + &
      sqrt((1.22464679914735d-16*y - 1.0d0*z)**2 + ( &
      -0.474856389870595d0*x + 0.880063298291132d0*y + &
      1.0777667012993d-16*z)**2 + (0.880063298291132d0*x + &
      0.474856389870595d0*y + 5.81531357909692d-17*z)**2))**2 + &
      1.70124732452287d0*(0.474856389870595d0*x - 0.880063298291132d0*y &
      - 1.0777667012993d-16*z)**2*((1.22464679914735d-16*y - 1.0d0*z)** &
      2 + (-0.880063298291132d0*x - 0.474856389870595d0*y - &
      5.81531357909691d-17*z)**2 + (0.474856389870595d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)**2)**(-0.5d0)/( &
      -1.22464679914735d-16*y + 1.0d0*z + sqrt((1.22464679914735d-16*y &
      - 1.0d0*z)**2 + (-0.880063298291132d0*x - 0.474856389870595d0*y - &
      5.81531357909691d-17*z)**2 + (0.474856389870595d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)**2))**2 + &
      0.242094922645643d0*(0.880063298291132d0*x + 0.474856389870595d0* &
      y + 5.81531357909692d-17*z)**3*((1.22464679914735d-16*y - 1.0d0*z &
      )**2 + (-0.474856389870595d0*x + 0.880063298291132d0*y + &
      1.0777667012993d-16*z)**2 + (0.880063298291132d0*x + &
      0.474856389870595d0*y + 5.81531357909692d-17*z)**2)**(-1.5d0)/( &
      -1.22464679914735d-16*y + 1.0d0*z + sqrt((1.22464679914735d-16*y &
      - 1.0d0*z)**2 + (-0.474856389870595d0*x + 0.880063298291132d0*y + &
      1.0777667012993d-16*z)**2 + (0.880063298291132d0*x + &
      0.474856389870595d0*y + 5.81531357909692d-17*z)**2)) + &
      0.242094922645643d0*(0.880063298291132d0*x + 0.474856389870595d0* &
      y + 5.81531357909692d-17*z)**3*1.0/((1.22464679914735d-16*y - &
      1.0d0*z)**2 + (-0.474856389870595d0*x + 0.880063298291132d0*y + &
      1.0777667012993d-16*z)**2 + (0.880063298291132d0*x + &
      0.474856389870595d0*y + 5.81531357909692d-17*z)**2)/( &
      -1.22464679914735d-16*y + 1.0d0*z + sqrt((1.22464679914735d-16*y &
      - 1.0d0*z)**2 + (-0.474856389870595d0*x + 0.880063298291132d0*y + &
      1.0777667012993d-16*z)**2 + (0.880063298291132d0*x + &
      0.474856389870595d0*y + 5.81531357909692d-17*z)**2))**2 - &
      0.581027814349542d0*(0.880063298291132d0*x + 0.474856389870595d0* &
      y + 5.81531357909692d-17*z)**3*((1.22464679914735d-16*y - 1.0d0*z &
      )**2 + (-0.474856389870595d0*x + 0.880063298291132d0*y + &
      1.0777667012993d-16*z)**2 + (0.880063298291132d0*x + &
      0.474856389870595d0*y + 5.81531357909692d-17*z)**2)**(-0.5d0)/( &
      -1.22464679914735d-16*y + 1.0d0*z + sqrt((1.22464679914735d-16*y &
      - 1.0d0*z)**2 + (-0.474856389870595d0*x + 0.880063298291132d0*y + &
      1.0777667012993d-16*z)**2 + (0.880063298291132d0*x + &
      0.474856389870595d0*y + 5.81531357909692d-17*z)**2))**3 - &
      0.435770860762157d0*(0.880063298291132d0*x + 0.474856389870595d0* &
      y + 5.81531357909692d-17*z)*((1.22464679914735d-16*y - 1.0d0*z)** &
      2 + (-0.474856389870595d0*x + 0.880063298291132d0*y + &
      1.0777667012993d-16*z)**2 + (0.880063298291132d0*x + &
      0.474856389870595d0*y + 5.81531357909692d-17*z)**2)**(-0.5d0)/( &
      -1.22464679914735d-16*y + 1.0d0*z + sqrt((1.22464679914735d-16*y &
      - 1.0d0*z)**2 + (-0.474856389870595d0*x + 0.880063298291132d0*y + &
      1.0777667012993d-16*z)**2 + (0.880063298291132d0*x + &
      0.474856389870595d0*y + 5.81531357909692d-17*z)**2)) + &
      0.121047461322821d0*(4.22430383179743d0*x + 2.27931067137886d0*y &
      + 2.79135051796652d-16*z)/(-1.22464679914735d-16*y + 1.0d0*z + &
      sqrt((1.22464679914735d-16*y - 1.0d0*z)**2 + ( &
      -0.474856389870595d0*x + 0.880063298291132d0*y + &
      1.0777667012993d-16*z)**2 + (0.880063298291132d0*x + &
      0.474856389870595d0*y + 5.81531357909692d-17*z)**2))**2 + &
      0.0455370902743939d0*((1.22464679914735d-16*y - 1.0d0*z)**2 + ( &
      -0.474856389870595d0*x + 0.880063298291132d0*y + &
      1.0777667012993d-16*z)**2 + (0.880063298291132d0*x + &
      0.474856389870595d0*y + 5.81531357909692d-17*z)**2)**(-0.5d0) - &
      0.502461983514667d0/(-1.22464679914735d-16*y + 1.0d0*z + sqrt(( &
      1.22464679914735d-16*y - 1.0d0*z)**2 + (-0.474856389870595d0*x + &
      0.880063298291132d0*y + 1.0777667012993d-16*z)**2 + ( &
      0.880063298291132d0*x + 0.474856389870595d0*y + &
      5.81531357909692d-17*z)**2)) - 3.28790477395282d0/( &
      -3.84734138744357d-16*y + 3.14159265358979d0*z + &
      3.14159265358979d0*sqrt((1.22464679914735d-16*y - 1.0d0*z)**2 + ( &
      -0.880063298291132d0*x - 0.474856389870595d0*y - &
      5.81531357909691d-17*z)**2 + (0.474856389870595d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)**2)) - &
      0.564888152839834d0*(-1.76012659658226d0*x - 0.94971277974119d0*y &
      - 1.16306271581938d-16*z)/(2.8d0*(-0.474856389870595d0*x + &
      0.880063298291132d0*y + 1.0777667012993d-16*z)**2 + 2.8d0*( &
      0.880063298291132d0*x + 0.474856389870595d0*y + &
      5.81531357909692d-17*z)**2) - 1.04692143047561d0*( &
      -0.94971277974119d0*x + 1.76012659658226d0*y + &
      2.1555334025986d-16*z)/(2.8d0*(-0.474856389870595d0*x + &
      0.880063298291132d0*y + 1.0777667012993d-16*z)**2 + 2.8d0*( &
      0.880063298291132d0*x + 0.474856389870595d0*y + &
      5.81531357909692d-17*z)**2) + 1.04692143047561d0*(5.6d0*( &
      -0.474856389870595d0*x + 0.880063298291132d0*y + &
      1.0777667012993d-16*z)**2 - 5.6d0*(0.880063298291132d0*x + &
      0.474856389870595d0*y + 5.81531357909692d-17*z)**2)*( &
      -0.474856389870595d0*x + 0.880063298291132d0*y + &
      1.0777667012993d-16*z)/(2.8d0*(-0.474856389870595d0*x + &
      0.880063298291132d0*y + 1.0777667012993d-16*z)**2 + 2.8d0*( &
      0.880063298291132d0*x + 0.474856389870595d0*y + &
      5.81531357909692d-17*z)**2)**2 + 0.564888152839834d0*(5.6d0*( &
      -0.474856389870595d0*x + 0.880063298291132d0*y + &
      1.0777667012993d-16*z)**2 - 5.6d0*(0.880063298291132d0*x + &
      0.474856389870595d0*y + 5.81531357909692d-17*z)**2)*( &
      0.880063298291132d0*x + 0.474856389870595d0*y + &
      5.81531357909692d-17*z)/(2.8d0*(-0.474856389870595d0*x + &
      0.880063298291132d0*y + 1.0777667012993d-16*z)**2 + 2.8d0*( &
      0.880063298291132d0*x + 0.474856389870595d0*y + &
      5.81531357909692d-17*z)**2)**2 + 0.402876757010956d0*( &
      -0.474856389870595d0*x + 0.880063298291132d0*y + &
      1.0777667012993d-16*z)/(1.4d0*(-0.474856389870595d0*x + &
      0.880063298291132d0*y + 1.0777667012993d-16*z)**2 + 1.4d0*( &
      0.880063298291132d0*x + 0.474856389870595d0*y + &
      5.81531357909692d-17*z)**2) + 0.746661633165595d0*( &
      0.880063298291132d0*x + 0.474856389870595d0*y + &
      5.81531357909692d-17*z)/(1.4d0*(-0.474856389870595d0*x + &
      0.880063298291132d0*y + 1.0777667012993d-16*z)**2 + 1.4d0*( &
      0.880063298291132d0*x + 0.474856389870595d0*y + &
      5.81531357909692d-17*z)**2) - 2.09065257286367d0*( &
      -0.474856389870595d0*x + 0.880063298291132d0*y + &
      1.0777667012993d-16*z)**2*(0.880063298291132d0*x + &
      0.474856389870595d0*y + 5.81531357909692d-17*z)/(1.4d0*( &
      -0.474856389870595d0*x + 0.880063298291132d0*y + &
      1.0777667012993d-16*z)**2 + 1.4d0*(0.880063298291132d0*x + &
      0.474856389870595d0*y + 5.81531357909692d-17*z)**2)**2 - &
      1.12805491963068d0*(-0.474856389870595d0*x + 0.880063298291132d0* &
      y + 1.0777667012993d-16*z)*(0.880063298291132d0*x + &
      0.474856389870595d0*y + 5.81531357909692d-17*z)**2/(1.4d0*( &
      -0.474856389870595d0*x + 0.880063298291132d0*y + &
      1.0777667012993d-16*z)**2 + 1.4d0*(0.880063298291132d0*x + &
      0.474856389870595d0*y + 5.81531357909692d-17*z)**2)**2 - &
      1.04692143047561d0*(-0.135673254248742d0*x + 0.251446656654609d0* &
      y + 3.07933343228372d-17*z)/((-0.474856389870595d0*x + &
      0.880063298291132d0*y + 1.0777667012993d-16*z)**2 + ( &
      0.880063298291132d0*x + 0.474856389870595d0*y + &
      5.81531357909692d-17*z)**2) - 0.564888152839834d0*( &
      0.251446656654609d0*x + 0.135673254248742d0*y + &
      1.66151816545626d-17*z)/((-0.474856389870595d0*x + &
      0.880063298291132d0*y + 1.0777667012993d-16*z)**2 + ( &
      0.880063298291132d0*x + 0.474856389870595d0*y + &
      5.81531357909692d-17*z)**2) - (-0.884730161365563d0*x - &
      0.477374492552354d0*y - 5.84615144298831d-17*z)/(( &
      3.14159265358979d0*(-0.880063298291132d0*x - 0.474856389870595d0* &
      y - 5.81531357909691d-17*z)**2/(0.474856389870595d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)**2 + &
      3.14159265358979d0)*(0.474856389870595d0*x - 0.880063298291132d0* &
      y - 1.0777667012993d-16*z)**2) + 0.542431997197586d0/(( &
      3.14159265358979d0*(-0.880063298291132d0*x - 0.474856389870595d0* &
      y - 5.81531357909691d-17*z)**2/(0.474856389870595d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)**2 + &
      3.14159265358979d0)*(0.474856389870595d0*x - 0.880063298291132d0* &
      y - 1.0777667012993d-16*z)) + 0.402876757010956d0/((1 + ( &
      0.880063298291132d0*x + 0.474856389870595d0*y + &
      5.81531357909692d-17*z)**2/(-0.474856389870595d0*x + &
      0.880063298291132d0*y + 1.0777667012993d-16*z)**2)*( &
      -0.474856389870595d0*x + 0.880063298291132d0*y + &
      1.0777667012993d-16*z)) - 0.746661633165595d0*( &
      0.880063298291132d0*x + 0.474856389870595d0*y + &
      5.81531357909692d-17*z)/((1 + (0.880063298291132d0*x + &
      0.474856389870595d0*y + 5.81531357909692d-17*z)**2/( &
      -0.474856389870595d0*x + 0.880063298291132d0*y + &
      1.0777667012993d-16*z)**2)*(-0.474856389870595d0*x + &
      0.880063298291132d0*y + 1.0777667012993d-16*z)**2)

end subroutine
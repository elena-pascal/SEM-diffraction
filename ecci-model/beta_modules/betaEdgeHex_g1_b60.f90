!******************************************************************************
!*                       Code generated with sympy 1.0                        *
!*                                                                            *
!*              See http://www.sympy.org/ for more information.               *
!*                                                                            *
!*                 This file is part of 'Edgebeta_generation'                 *
!******************************************************************************

subroutine beta_ECCI(x, y, z, betaH)
implicit none
REAL*8, intent(in) :: x
REAL*8, intent(in) :: y
REAL*8, intent(in) :: z
REAL*8, intent(out) :: betaH

betaH = 0.530493356419341d0*(1.22464679914735d-16*y - 1.0d0*z)**2*( &
      -0.880063298291132d0*x + 0.474856389870594d0*y + &
      5.81531357909691d-17*z)*(-0.474856389870594d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)*(( &
      1.22464679914735d-16*y - 1.0d0*z)**2 + (-0.880063298291132d0*x + &
      0.474856389870594d0*y + 5.81531357909691d-17*z)**2 + ( &
      -0.474856389870594d0*x - 0.880063298291132d0*y - &
      1.0777667012993d-16*z)**2)**(-1.5d0)/(-1.22464679914735d-16*y + &
      1.0d0*z + sqrt((1.22464679914735d-16*y - 1.0d0*z)**2 + ( &
      -0.880063298291132d0*x + 0.474856389870594d0*y + &
      5.81531357909691d-17*z)**2 + (-0.474856389870594d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)**2))**2 - &
      0.0255233735110012d0*(1.22464679914735d-16*y - 1.0d0*z)*(-2.4d0*( &
      1.22464679914735d-16*y - 1.0d0*z)*((1.22464679914735d-16*y - &
      1.0d0*z)**2 + (-0.880063298291132d0*x + 0.474856389870594d0*y + &
      5.81531357909691d-17*z)**2 + (-0.474856389870594d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)**2)**(-0.5d0) + &
      2.4d0)/(-1.22464679914735d-16*y + 1.0d0*z + sqrt(( &
      1.22464679914735d-16*y - 1.0d0*z)**2 + (-0.880063298291132d0*x + &
      0.474856389870594d0*y + 5.81531357909691d-17*z)**2 + ( &
      -0.474856389870594d0*x - 0.880063298291132d0*y - &
      1.0777667012993d-16*z)**2))**2 - 0.530493356419341d0*( &
      1.22464679914735d-16*y - 1.0d0*z)*(-2.0d0*(1.22464679914735d-16*y &
      - 1.0d0*z)*((1.22464679914735d-16*y - 1.0d0*z)**2 + ( &
      -0.880063298291132d0*x + 0.474856389870594d0*y + &
      5.81531357909691d-17*z)**2 + (-0.474856389870594d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)**2)**(-0.5d0) + 2) &
      *(-0.880063298291132d0*x + 0.474856389870594d0*y + &
      5.81531357909691d-17*z)*(-0.474856389870594d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)*(( &
      1.22464679914735d-16*y - 1.0d0*z)**2 + (-0.880063298291132d0*x + &
      0.474856389870594d0*y + 5.81531357909691d-17*z)**2 + ( &
      -0.474856389870594d0*x - 0.880063298291132d0*y - &
      1.0777667012993d-16*z)**2)**(-0.5d0)/(-1.22464679914735d-16*y + &
      1.0d0*z + sqrt((1.22464679914735d-16*y - 1.0d0*z)**2 + ( &
      -0.880063298291132d0*x + 0.474856389870594d0*y + &
      5.81531357909691d-17*z)**2 + (-0.474856389870594d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)**2))**3 + &
      0.29599108678356d0*(1.22464679914735d-16*y - 1.0d0*z)*( &
      -0.880063298291132d0*x + 0.474856389870594d0*y + &
      5.81531357909691d-17*z)**2*(-0.474856389870594d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)*(( &
      1.22464679914735d-16*y - 1.0d0*z)**2 + (-0.880063298291132d0*x + &
      0.474856389870594d0*y + 5.81531357909691d-17*z)**2 + ( &
      -0.474856389870594d0*x - 0.880063298291132d0*y - &
      1.0777667012993d-16*z)**2)**(-1.5d0)/(-1.22464679914735d-16*y + &
      1.0d0*z + sqrt((1.22464679914735d-16*y - 1.0d0*z)**2 + ( &
      -0.880063298291132d0*x + 0.474856389870594d0*y + &
      5.81531357909691d-17*z)**2 + (-0.474856389870594d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)**2))**2 + &
      0.591982173567119d0*(1.22464679914735d-16*y - 1.0d0*z)*( &
      -0.880063298291132d0*x + 0.474856389870594d0*y + &
      5.81531357909691d-17*z)**2*(-0.474856389870594d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)*1.0/(( &
      1.22464679914735d-16*y - 1.0d0*z)**2 + (-0.880063298291132d0*x + &
      0.474856389870594d0*y + 5.81531357909691d-17*z)**2 + ( &
      -0.474856389870594d0*x - 0.880063298291132d0*y - &
      1.0777667012993d-16*z)**2)/(-1.22464679914735d-16*y + 1.0d0*z + &
      sqrt((1.22464679914735d-16*y - 1.0d0*z)**2 + ( &
      -0.880063298291132d0*x + 0.474856389870594d0*y + &
      5.81531357909691d-17*z)**2 + (-0.474856389870594d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)**2))**3 + &
      0.0510467470220023d0*(1.22464679914735d-16*y - 1.0d0*z)*( &
      -0.880063298291132d0*x + 0.474856389870594d0*y + &
      5.81531357909691d-17*z)**2*((1.22464679914735d-16*y - 1.0d0*z)**2 &
      + (-0.880063298291132d0*x + 0.474856389870594d0*y + &
      5.81531357909691d-17*z)**2 + (-0.474856389870594d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)**2)**(-1.5d0)/( &
      -1.22464679914735d-16*y + 1.0d0*z + sqrt((1.22464679914735d-16*y &
      - 1.0d0*z)**2 + (-0.880063298291132d0*x + 0.474856389870594d0*y + &
      5.81531357909691d-17*z)**2 + (-0.474856389870594d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)**2)) - &
      0.54856773048901d0*(1.22464679914735d-16*y - 1.0d0*z)*( &
      -0.880063298291132d0*x + 0.474856389870594d0*y + &
      5.81531357909691d-17*z)*(-0.474856389870594d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)**2*(( &
      1.22464679914735d-16*y - 1.0d0*z)**2 + (-0.880063298291132d0*x + &
      0.474856389870594d0*y + 5.81531357909691d-17*z)**2 + ( &
      -0.474856389870594d0*x - 0.880063298291132d0*y - &
      1.0777667012993d-16*z)**2)**(-1.5d0)/(-1.22464679914735d-16*y + &
      1.0d0*z + sqrt((1.22464679914735d-16*y - 1.0d0*z)**2 + ( &
      -0.880063298291132d0*x + 0.474856389870594d0*y + &
      5.81531357909691d-17*z)**2 + (-0.474856389870594d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)**2))**2 - &
      1.09713546097802d0*(1.22464679914735d-16*y - 1.0d0*z)*( &
      -0.880063298291132d0*x + 0.474856389870594d0*y + &
      5.81531357909691d-17*z)*(-0.474856389870594d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)**2*1.0/(( &
      1.22464679914735d-16*y - 1.0d0*z)**2 + (-0.880063298291132d0*x + &
      0.474856389870594d0*y + 5.81531357909691d-17*z)**2 + ( &
      -0.474856389870594d0*x - 0.880063298291132d0*y - &
      1.0777667012993d-16*z)**2)/(-1.22464679914735d-16*y + 1.0d0*z + &
      sqrt((1.22464679914735d-16*y - 1.0d0*z)**2 + ( &
      -0.880063298291132d0*x + 0.474856389870594d0*y + &
      5.81531357909691d-17*z)**2 + (-0.474856389870594d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)**2))**3 + &
      0.582745837215453d0*(1.22464679914735d-16*y - 1.0d0*z)*( &
      -0.880063298291132d0*x + 0.474856389870594d0*y + &
      5.81531357909691d-17*z)*((1.22464679914735d-16*y - 1.0d0*z)**2 + &
      (-0.880063298291132d0*x + 0.474856389870594d0*y + &
      5.81531357909691d-17*z)**2 + (-0.474856389870594d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)**2)**(-0.5d0)/( &
      -1.22464679914735d-16*y + 1.0d0*z + sqrt((1.22464679914735d-16*y &
      - 1.0d0*z)**2 + (-0.880063298291132d0*x + 0.474856389870594d0*y + &
      5.81531357909691d-17*z)**2 + (-0.474856389870594d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)**2))**2 - &
      0.359334232156874d0*(1.22464679914735d-16*y - 1.0d0*z)*( &
      -0.474856389870594d0*x - 0.880063298291132d0*y - &
      1.0777667012993d-16*z)*((1.22464679914735d-16*y - 1.0d0*z)**2 + ( &
      -0.880063298291132d0*x + 0.474856389870594d0*y + &
      5.81531357909691d-17*z)**2 + (-0.474856389870594d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)**2)**(-0.5d0)/( &
      -1.22464679914735d-16*y + 1.0d0*z + sqrt((1.22464679914735d-16*y &
      - 1.0d0*z)**2 + (-0.880063298291132d0*x + 0.474856389870594d0*y + &
      5.81531357909691d-17*z)**2 + (-0.474856389870594d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)**2))**2 + &
      0.318309886183791d0*(-0.256399081739261d0*(1.22464679914735d-16*y &
      - 1.0d0*z)*((1.22464679914735d-16*y - 1.0d0*z)**2 + ( &
      -0.880063298291132d0*x + 0.474856389870594d0*y + &
      5.81531357909691d-17*z)**2 + (-0.474856389870594d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)**2)**(-1.5d0) + &
      0.256399081739261d0*(-0.4d0*(1.22464679914735d-16*y - 1.0d0*z)*(( &
      1.22464679914735d-16*y - 1.0d0*z)**2 + (-0.880063298291132d0*x + &
      0.474856389870594d0*y + 5.81531357909691d-17*z)**2 + ( &
      -0.474856389870594d0*x - 0.880063298291132d0*y - &
      1.0777667012993d-16*z)**2)**(-0.5d0) + 0.4d0)/( &
      -1.22464679914735d-16*y + 1.0d0*z + sqrt((1.22464679914735d-16*y &
      - 1.0d0*z)**2 + (-0.880063298291132d0*x + 0.474856389870594d0*y + &
      5.81531357909691d-17*z)**2 + (-0.474856389870594d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)**2))**2)*( &
      -0.880063298291132d0*x + 0.474856389870594d0*y + &
      5.81531357909691d-17*z) + 0.0612560964264028d0*(-2.0d0*( &
      1.22464679914735d-16*y - 1.0d0*z)*((1.22464679914735d-16*y - &
      1.0d0*z)**2 + (-0.880063298291132d0*x + 0.474856389870594d0*y + &
      5.81531357909691d-17*z)**2 + (-0.474856389870594d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)**2)**(-0.5d0) + 2) &
      *(-0.880063298291132d0*x + 0.474856389870594d0*y + &
      5.81531357909691d-17*z)**2/(-1.22464679914735d-16*y + 1.0d0*z + &
      sqrt((1.22464679914735d-16*y - 1.0d0*z)**2 + ( &
      -0.880063298291132d0*x + 0.474856389870594d0*y + &
      5.81531357909691d-17*z)**2 + (-0.474856389870594d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)**2))**3 - &
      0.0510467470220023d0*(-1.0d0*(1.22464679914735d-16*y - 1.0d0*z)* &
      ((1.22464679914735d-16*y - 1.0d0*z)**2 + (-0.880063298291132d0*x &
      + 0.474856389870594d0*y + 5.81531357909691d-17*z)**2 + ( &
      -0.474856389870594d0*x - 0.880063298291132d0*y - &
      1.0777667012993d-16*z)**2)**(-0.5d0) + 1)*(-0.880063298291132d0*x &
      + 0.474856389870594d0*y + 5.81531357909691d-17*z)**2*(( &
      1.22464679914735d-16*y - 1.0d0*z)**2 + (-0.880063298291132d0*x + &
      0.474856389870594d0*y + 5.81531357909691d-17*z)**2 + ( &
      -0.474856389870594d0*x - 0.880063298291132d0*y - &
      1.0777667012993d-16*z)**2)**(-0.5d0)/(-1.22464679914735d-16*y + &
      1.0d0*z + sqrt((1.22464679914735d-16*y - 1.0d0*z)**2 + ( &
      -0.880063298291132d0*x + 0.474856389870594d0*y + &
      5.81531357909691d-17*z)**2 + (-0.474856389870594d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)**2))**2 - &
      0.265246678209671d0*(-0.8d0*(1.22464679914735d-16*y - 1.0d0*z)*(( &
      1.22464679914735d-16*y - 1.0d0*z)**2 + (-0.880063298291132d0*x + &
      0.474856389870594d0*y + 5.81531357909691d-17*z)**2 + ( &
      -0.474856389870594d0*x - 0.880063298291132d0*y - &
      1.0777667012993d-16*z)**2)**(-0.5d0) + 0.8d0)*( &
      -0.880063298291132d0*x + 0.474856389870594d0*y + &
      5.81531357909691d-17*z)*(-0.474856389870594d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)/( &
      -1.22464679914735d-16*y + 1.0d0*z + sqrt((1.22464679914735d-16*y &
      - 1.0d0*z)**2 + (-0.880063298291132d0*x + 0.474856389870594d0*y + &
      5.81531357909691d-17*z)**2 + (-0.474856389870594d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)**2))**3 + &
      0.0255233735110012d0*(0.4d0*(1.22464679914735d-16*y - 1.0d0*z)*(( &
      1.22464679914735d-16*y - 1.0d0*z)**2 + (-0.880063298291132d0*x + &
      0.474856389870594d0*y + 5.81531357909691d-17*z)**2 + ( &
      -0.474856389870594d0*x - 0.880063298291132d0*y - &
      1.0777667012993d-16*z)**2)**(-0.5d0) - 0.4d0)/( &
      -1.22464679914735d-16*y + 1.0d0*z + sqrt((1.22464679914735d-16*y &
      - 1.0d0*z)**2 + (-0.880063298291132d0*x + 0.474856389870594d0*y + &
      5.81531357909691d-17*z)**2 + (-0.474856389870594d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)**2)) + &
      0.318309886183791d0*(-0.143058988271891d0*(-0.880063298291132d0*x &
      + 0.474856389870594d0*y + 5.81531357909691d-17*z)*(( &
      1.22464679914735d-16*y - 1.0d0*z)**2 + (-0.880063298291132d0*x + &
      0.474856389870594d0*y + 5.81531357909691d-17*z)**2 + ( &
      -0.474856389870594d0*x - 0.880063298291132d0*y - &
      1.0777667012993d-16*z)**2)**(-1.5d0) - 0.0572235953087564d0*( &
      -0.880063298291132d0*x + 0.474856389870594d0*y + &
      5.81531357909691d-17*z)*((1.22464679914735d-16*y - 1.0d0*z)**2 + &
      (-0.880063298291132d0*x + 0.474856389870594d0*y + &
      5.81531357909691d-17*z)**2 + (-0.474856389870594d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)**2)**(-0.5d0)/( &
      -1.22464679914735d-16*y + 1.0d0*z + sqrt((1.22464679914735d-16*y &
      - 1.0d0*z)**2 + (-0.880063298291132d0*x + 0.474856389870594d0*y + &
      5.81531357909691d-17*z)**2 + (-0.474856389870594d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)**2))**2)*( &
      -0.880063298291132d0*x + 0.474856389870594d0*y + &
      5.81531357909691d-17*z) - 0.318309886183791d0*( &
      -0.265134823400107d0*(-0.474856389870594d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)*(( &
      1.22464679914735d-16*y - 1.0d0*z)**2 + (-0.880063298291132d0*x + &
      0.474856389870594d0*y + 5.81531357909691d-17*z)**2 + ( &
      -0.474856389870594d0*x - 0.880063298291132d0*y - &
      1.0777667012993d-16*z)**2)**(-1.5d0) - 0.106053929360043d0*( &
      -0.474856389870594d0*x - 0.880063298291132d0*y - &
      1.0777667012993d-16*z)*((1.22464679914735d-16*y - 1.0d0*z)**2 + ( &
      -0.880063298291132d0*x + 0.474856389870594d0*y + &
      5.81531357909691d-17*z)**2 + (-0.474856389870594d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)**2)**(-0.5d0)/( &
      -1.22464679914735d-16*y + 1.0d0*z + sqrt((1.22464679914735d-16*y &
      - 1.0d0*z)**2 + (-0.880063298291132d0*x + 0.474856389870594d0*y + &
      5.81531357909691d-17*z)**2 + (-0.474856389870594d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)**2))**2)*( &
      -0.880063298291132d0*x + 0.474856389870594d0*y + &
      5.81531357909691d-17*z) + 0.0664574297458624d0*(-5.6d0*( &
      -0.880063298291132d0*x + 0.474856389870594d0*y + &
      5.81531357909691d-17*z)**2 + 5.6d0*(-0.474856389870594d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)**2)*( &
      -0.880063298291132d0*x + 0.474856389870594d0*y + &
      5.81531357909691d-17*z)/(2.8d0*(-0.880063298291132d0*x + &
      0.474856389870594d0*y + 5.81531357909691d-17*z)**2 + 2.8d0*( &
      -0.474856389870594d0*x - 0.880063298291132d0*y - &
      1.0777667012993d-16*z)**2)**2 - 0.123167227114778d0*(-5.6d0*( &
      -0.880063298291132d0*x + 0.474856389870594d0*y + &
      5.81531357909691d-17*z)**2 + 5.6d0*(-0.474856389870594d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)**2)*( &
      -0.474856389870594d0*x - 0.880063298291132d0*y - &
      1.0777667012993d-16*z)/(2.8d0*(-0.880063298291132d0*x + &
      0.474856389870594d0*y + 5.81531357909691d-17*z)**2 + 2.8d0*( &
      -0.474856389870594d0*x - 0.880063298291132d0*y - &
      1.0777667012993d-16*z)**2)**2 + 0.0142408778026848d0*( &
      -4.22430383179743d0*x + 2.27931067137885d0*y + &
      2.79135051796652d-16*z)/(-1.22464679914735d-16*y + 1.0d0*z + sqrt &
      ((1.22464679914735d-16*y - 1.0d0*z)**2 + (-0.880063298291132d0*x &
      + 0.474856389870594d0*y + 5.81531357909691d-17*z)**2 + ( &
      -0.474856389870594d0*x - 0.880063298291132d0*y - &
      1.0777667012993d-16*z)**2))**2 + 0.0284817556053696d0*( &
      -0.880063298291132d0*x + 0.474856389870594d0*y + &
      5.81531357909691d-17*z)**3*((1.22464679914735d-16*y - 1.0d0*z)**2 &
      + (-0.880063298291132d0*x + 0.474856389870594d0*y + &
      5.81531357909691d-17*z)**2 + (-0.474856389870594d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)**2)**(-1.5d0)/( &
      -1.22464679914735d-16*y + 1.0d0*z + sqrt((1.22464679914735d-16*y &
      - 1.0d0*z)**2 + (-0.880063298291132d0*x + 0.474856389870594d0*y + &
      5.81531357909691d-17*z)**2 + (-0.474856389870594d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)**2)) + &
      0.0284817556053696d0*(-0.880063298291132d0*x + &
      0.474856389870594d0*y + 5.81531357909691d-17*z)**3*1.0/(( &
      1.22464679914735d-16*y - 1.0d0*z)**2 + (-0.880063298291132d0*x + &
      0.474856389870594d0*y + 5.81531357909691d-17*z)**2 + ( &
      -0.474856389870594d0*x - 0.880063298291132d0*y - &
      1.0777667012993d-16*z)**2)/(-1.22464679914735d-16*y + 1.0d0*z + &
      sqrt((1.22464679914735d-16*y - 1.0d0*z)**2 + ( &
      -0.880063298291132d0*x + 0.474856389870594d0*y + &
      5.81531357909691d-17*z)**2 + (-0.474856389870594d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)**2))**2 - &
      0.0683562134528872d0*(-0.880063298291132d0*x + &
      0.474856389870594d0*y + 5.81531357909691d-17*z)**3*(( &
      1.22464679914735d-16*y - 1.0d0*z)**2 + (-0.880063298291132d0*x + &
      0.474856389870594d0*y + 5.81531357909691d-17*z)**2 + ( &
      -0.474856389870594d0*x - 0.880063298291132d0*y - &
      1.0777667012993d-16*z)**2)**(-0.5d0)/(-1.22464679914735d-16*y + &
      1.0d0*z + sqrt((1.22464679914735d-16*y - 1.0d0*z)**2 + ( &
      -0.880063298291132d0*x + 0.474856389870594d0*y + &
      5.81531357909691d-17*z)**2 + (-0.474856389870594d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)**2))**3 - &
      0.0527859544777619d0*(-0.880063298291132d0*x + &
      0.474856389870594d0*y + 5.81531357909691d-17*z)**2*( &
      -0.474856389870594d0*x - 0.880063298291132d0*y - &
      1.0777667012993d-16*z)*((1.22464679914735d-16*y - 1.0d0*z)**2 + ( &
      -0.880063298291132d0*x + 0.474856389870594d0*y + &
      5.81531357909691d-17*z)**2 + (-0.474856389870594d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)**2)**(-1.5d0)/( &
      -1.22464679914735d-16*y + 1.0d0*z + sqrt((1.22464679914735d-16*y &
      - 1.0d0*z)**2 + (-0.880063298291132d0*x + 0.474856389870594d0*y + &
      5.81531357909691d-17*z)**2 + (-0.474856389870594d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)**2)) - &
      0.0527859544777619d0*(-0.880063298291132d0*x + &
      0.474856389870594d0*y + 5.81531357909691d-17*z)**2*( &
      -0.474856389870594d0*x - 0.880063298291132d0*y - &
      1.0777667012993d-16*z)*1.0/((1.22464679914735d-16*y - 1.0d0*z)**2 &
      + (-0.880063298291132d0*x + 0.474856389870594d0*y + &
      5.81531357909691d-17*z)**2 + (-0.474856389870594d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)**2)/( &
      -1.22464679914735d-16*y + 1.0d0*z + sqrt((1.22464679914735d-16*y &
      - 1.0d0*z)**2 + (-0.880063298291132d0*x + 0.474856389870594d0*y + &
      5.81531357909691d-17*z)**2 + (-0.474856389870594d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)**2))**2 + &
      0.245082725460052d0*(-0.880063298291132d0*x + 0.474856389870594d0 &
      *y + 5.81531357909691d-17*z)**2*(-0.474856389870594d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)*(( &
      1.22464679914735d-16*y - 1.0d0*z)**2 + (-0.880063298291132d0*x + &
      0.474856389870594d0*y + 5.81531357909691d-17*z)**2 + ( &
      -0.474856389870594d0*x - 0.880063298291132d0*y - &
      1.0777667012993d-16*z)**2)**(-0.5d0)/(-1.22464679914735d-16*y + &
      1.0d0*z + sqrt((1.22464679914735d-16*y - 1.0d0*z)**2 + ( &
      -0.880063298291132d0*x + 0.474856389870594d0*y + &
      5.81531357909691d-17*z)**2 + (-0.474856389870594d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)**2))**3 - &
      0.219427092195604d0*(-0.880063298291132d0*x + 0.474856389870594d0 &
      *y + 5.81531357909691d-17*z)*(-0.474856389870594d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)**2*(( &
      1.22464679914735d-16*y - 1.0d0*z)**2 + (-0.880063298291132d0*x + &
      0.474856389870594d0*y + 5.81531357909691d-17*z)**2 + ( &
      -0.474856389870594d0*x - 0.880063298291132d0*y - &
      1.0777667012993d-16*z)**2)**(-0.5d0)/(-1.22464679914735d-16*y + &
      1.0d0*z + sqrt((1.22464679914735d-16*y - 1.0d0*z)**2 + ( &
      -0.880063298291132d0*x + 0.474856389870594d0*y + &
      5.81531357909691d-17*z)**2 + (-0.474856389870594d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)**2))**3 - &
      0.530493356419341d0*(-0.880063298291132d0*x + 0.474856389870594d0 &
      *y + 5.81531357909691d-17*z)*(-0.474856389870594d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)*(( &
      1.22464679914735d-16*y - 1.0d0*z)**2 + (-0.880063298291132d0*x + &
      0.474856389870594d0*y + 5.81531357909691d-17*z)**2 + ( &
      -0.474856389870594d0*x - 0.880063298291132d0*y - &
      1.0777667012993d-16*z)**2)**(-0.5d0)/(-1.22464679914735d-16*y + &
      1.0d0*z + sqrt((1.22464679914735d-16*y - 1.0d0*z)**2 + ( &
      -0.880063298291132d0*x + 0.474856389870594d0*y + &
      5.81531357909691d-17*z)**2 + (-0.474856389870594d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)**2))**2 - &
      0.0512671600896654d0*(-0.880063298291132d0*x + &
      0.474856389870594d0*y + 5.81531357909691d-17*z)*(( &
      1.22464679914735d-16*y - 1.0d0*z)**2 + (-0.880063298291132d0*x + &
      0.474856389870594d0*y + 5.81531357909691d-17*z)**2 + ( &
      -0.474856389870594d0*x - 0.880063298291132d0*y - &
      1.0777667012993d-16*z)**2)**(-0.5d0)/(-1.22464679914735d-16*y + &
      1.0d0*z + sqrt((1.22464679914735d-16*y - 1.0d0*z)**2 + ( &
      -0.880063298291132d0*x + 0.474856389870594d0*y + &
      5.81531357909691d-17*z)**2 + (-0.474856389870594d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)**2)) - &
      0.0105571908955524d0*(-0.474856389870594d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)*(( &
      1.22464679914735d-16*y - 1.0d0*z)**2 + (-0.880063298291132d0*x + &
      0.474856389870594d0*y + 5.81531357909691d-17*z)**2 + ( &
      -0.474856389870594d0*x - 0.880063298291132d0*y - &
      1.0777667012993d-16*z)**2)**(-0.5d0)/(-1.22464679914735d-16*y + &
      1.0d0*z + sqrt((1.22464679914735d-16*y - 1.0d0*z)**2 + ( &
      -0.880063298291132d0*x + 0.474856389870594d0*y + &
      5.81531357909691d-17*z)**2 + (-0.474856389870594d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)**2)) + &
      0.274283865244505d0*(-0.352025319316453d0*x + 0.189942555948238d0 &
      *y + 2.32612543163876d-17*z)/(-1.22464679914735d-16*y + 1.0d0*z + &
      sqrt((1.22464679914735d-16*y - 1.0d0*z)**2 + ( &
      -0.880063298291132d0*x + 0.474856389870594d0*y + &
      5.81531357909691d-17*z)**2 + (-0.474856389870594d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)**2))**2 - &
      0.14799554339178d0*(-0.189942555948238d0*x - 0.352025319316453d0* &
      y - 4.3110668051972d-17*z)/(-1.22464679914735d-16*y + 1.0d0*z + &
      sqrt((1.22464679914735d-16*y - 1.0d0*z)**2 + ( &
      -0.880063298291132d0*x + 0.474856389870594d0*y + &
      5.81531357909691d-17*z)**2 + (-0.474856389870594d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)**2))**2 + &
      0.0455370902743939d0*((1.22464679914735d-16*y - 1.0d0*z)**2 + ( &
      -0.880063298291132d0*x + 0.474856389870594d0*y + &
      5.81531357909691d-17*z)**2 + (-0.474856389870594d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)**2)**(-0.5d0) - &
      0.0430412603166453d0/(-1.22464679914735d-16*y + 1.0d0*z + sqrt(( &
      1.22464679914735d-16*y - 1.0d0*z)**2 + (-0.880063298291132d0*x + &
      0.474856389870594d0*y + 5.81531357909691d-17*z)**2 + ( &
      -0.474856389870594d0*x - 0.880063298291132d0*y - &
      1.0777667012993d-16*z)**2)) + 0.123167227114778d0*( &
      -0.949712779741188d0*x - 1.76012659658226d0*y - &
      2.1555334025986d-16*z)/(2.8d0*(-0.880063298291132d0*x + &
      0.474856389870594d0*y + 5.81531357909691d-17*z)**2 + 2.8d0*( &
      -0.474856389870594d0*x - 0.880063298291132d0*y - &
      1.0777667012993d-16*z)**2) - 0.0664574297458624d0*( &
      1.76012659658226d0*x - 0.949712779741188d0*y - &
      1.16306271581938d-16*z)/(2.8d0*(-0.880063298291132d0*x + &
      0.474856389870594d0*y + 5.81531357909691d-17*z)**2 + 2.8d0*( &
      -0.474856389870594d0*x - 0.880063298291132d0*y - &
      1.0777667012993d-16*z)**2) + 1.27999137114102d0*( &
      -0.880063298291132d0*x + 0.474856389870594d0*y + &
      5.81531357909691d-17*z)/(1.4d0*(-0.880063298291132d0*x + &
      0.474856389870594d0*y + 5.81531357909691d-17*z)**2 + 1.4d0*( &
      -0.474856389870594d0*x - 0.880063298291132d0*y - &
      1.0777667012993d-16*z)**2) - 0.690645869161638d0*( &
      -0.474856389870594d0*x - 0.880063298291132d0*y - &
      1.0777667012993d-16*z)/(1.4d0*(-0.880063298291132d0*x + &
      0.474856389870594d0*y + 5.81531357909691d-17*z)**2 + 1.4d0*( &
      -0.474856389870594d0*x - 0.880063298291132d0*y - &
      1.0777667012993d-16*z)**2) + 1.93380843365259d0*( &
      -0.880063298291132d0*x + 0.474856389870594d0*y + &
      5.81531357909691d-17*z)**2*(-0.474856389870594d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)/(1.4d0*( &
      -0.880063298291132d0*x + 0.474856389870594d0*y + &
      5.81531357909691d-17*z)**2 + 1.4d0*(-0.474856389870594d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)**2)**2 - &
      3.58397583919486d0*(-0.880063298291132d0*x + 0.474856389870594d0* &
      y + 5.81531357909691d-17*z)*(-0.474856389870594d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)**2/(1.4d0*( &
      -0.880063298291132d0*x + 0.474856389870594d0*y + &
      5.81531357909691d-17*z)**2 + 1.4d0*(-0.474856389870594d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)**2)**2 - &
      0.0664574297458624d0*(-0.251446656654609d0*x + &
      0.135673254248741d0*y + 1.66151816545626d-17*z)/(( &
      -0.880063298291132d0*x + 0.474856389870594d0*y + &
      5.81531357909691d-17*z)**2 + (-0.474856389870594d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)**2) + &
      0.123167227114778d0*(-0.135673254248741d0*x - 0.251446656654609d0 &
      *y - 3.07933343228372d-17*z)/((-0.880063298291132d0*x + &
      0.474856389870594d0*y + 5.81531357909691d-17*z)**2 + ( &
      -0.474856389870594d0*x - 0.880063298291132d0*y - &
      1.0777667012993d-16*z)**2) - 1.27999137114102d0*( &
      -0.880063298291132d0*x + 0.474856389870594d0*y + &
      5.81531357909691d-17*z)/(((-0.880063298291132d0*x + &
      0.474856389870594d0*y + 5.81531357909691d-17*z)**2/( &
      -0.474856389870594d0*x - 0.880063298291132d0*y - &
      1.0777667012993d-16*z)**2 + 1)*(-0.474856389870594d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z)**2) - &
      0.690645869161638d0/(((-0.880063298291132d0*x + &
      0.474856389870594d0*y + 5.81531357909691d-17*z)**2/( &
      -0.474856389870594d0*x - 0.880063298291132d0*y - &
      1.0777667012993d-16*z)**2 + 1)*(-0.474856389870594d0*x - &
      0.880063298291132d0*y - 1.0777667012993d-16*z))

end subroutine

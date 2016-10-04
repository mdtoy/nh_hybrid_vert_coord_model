module physical_parameters

!-----------------------------------------------------------------
!	this modules specifices physical parameters and the units
!	  of those parameters (MKS is standard)
!-----------------------------------------------------------------

use kinds

implicit none
save

real (kind=dbl_kind), parameter ::                                &
       pi           = 3.14159265358979323846_dbl_kind,   &
       grav         = 9.8100_dbl_kind,                   &
       gas_const_R  = 287.000_dbl_kind,                  &
       spec_heat_cp = 1005.000_dbl_kind,                 &
       spec_heat_cv = 718.000_dbl_kind,                  &
       gas_const_Rv = 461.52_dbl_kind,                   &
       L_vap        = 2.50E+06_dbl_kind,                 &
       rhowater     = 1000._dbl_kind,                    &
       epsilon      = gas_const_R/gas_const_Rv,          &
       kappa        = gas_const_R/spec_heat_cp,          &
       gamma        = 1.0_dbl_kind/(1.0_dbl_kind-kappa), &
       gamma_mns1   = gamma - 1.0_dbl_kind,              &
       inv_kappa_mns1 = 1.0_dbl_kind/kappa-1.0_dbl_kind, &
       p0_ref       = 1.000E+05_dbl_kind,                &
       f_cor        = 0.000E-04_dbl_kind

real (kind=dbl_kind), parameter ::                                &
                c0        = 0.00000_dbl_kind,                     &
                c1        = 1.00000_dbl_kind,                     &
                c2        = 2.00000_dbl_kind,                     &
                c3        = 3.00000_dbl_kind,                     &
                c4        = 4.00000_dbl_kind,                     &
                p5        = 0.50000_dbl_kind,                     &
                p25       = 0.25000_dbl_kind,                     &
                inv6      = c1/6.0_dbl_kind,                      &
                inv8      = c1/8.0_dbl_kind,                      &
                inv12     = c1/12.0_dbl_kind,                     &
                invgrav    = c1/grav,                             &
                invepsilon = c1/epsilon,                          &
                inv_gas_const_R = c1/gas_const_R,                 &
                inv_spec_heat_cp = c1/spec_heat_cp,               &
                inv_L_vap = c1/L_vap,                             &
                inv_p0_ref = c1/p0_ref


!-------------------------------------------------------------------
!   VARIABLE DEFINITION
!-------------------------------------------------------------------
! pi = pi
! grav = gravitational acceleration (m/s^2)
! gas_const_R = gas constant for dry air (J/ (kg K))
! gas_const_Rv = gas constant for water vapor (J/ (kg K))
! L_vap = latent heat of vaporization (J/kg)
! epsilon = ratio of Rd to Rv
! spec_heat_cp = specific heat at constant pressure (J/ (kg K))
! kappa = gas_const_R/spec_heat_cp
! p0_ref = reference pressure (Pa)
! f_cor = Coriolis parameter
! rhowater = density of liquid water (kg/m^3)
!-------------------------------------------------------------------

end module physical_parameters


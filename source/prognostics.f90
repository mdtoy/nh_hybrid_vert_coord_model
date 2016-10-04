module prognostics

!-----------------------------------------------------------------
!
!   Contains arrays related to the prognostic variables, and
!       those diagnostic variables derived from the prognostic
!       variables.  Variables are initialized to zero here.
!       Surface geopotential also declared here.
!
!-----------------------------------------------------------------


use kinds
use physical_parameters
use model_parameters
use eta_coordinate


implicit none
save


!
! Declare layer-center prognostic variables
!
real (kind = dbl_kind), dimension(im,jm,nlm,ntprog) ::       &
                     u, v,   &   ! horizontal velocity components (m/s)
                     m           ! pseudo-density (kg/m^2/eta)

!
! Declare layer-edge prognostic variables
!                     
real (kind = dbl_kind), dimension(im,jm,1:nlm-1,ntprog) ::   &
                     w_l2,   &     ! vertical velocity (m/s)
                     phi_l2        ! geopotential (J/kg)
                                   ! Note:  w_l2 and phi_l2 at
                                   ! k = 0 (at surface)
                                   ! and k = nlm (model top)  are not
                                   ! prognosed but are determined by
                                   ! the upper and lower boundary
                                   ! conditions

real (kind = dbl_kind), dimension(im,jm,0:nlm,ntprog) ::     &
                     th_l2,  &     ! potential temperature (K)
                     qvpc_l2_m_l2, & ! ( water vapor + cloud mixing ratios )
                                     !   * m_l2 (kg/m^2/eta)
                     qr_l2_m_l2      ! rain water mixing ratio *
                                     ! m_l2 (kg/m^2/eta)


!
! Declare AB3 tendencies of prognostic variables
!
real (kind = dbl_kind), dimension(im,jm,nlm,nttend) ::       &
                     u_f,     &   ! d/dt (u)   (m/s^2)
                     v_f,     &   ! d/dt (v)   (m/s^2)
                     m_f          ! d/dt (m)   (kg/m^2/eta/s)

real (kind = dbl_kind), dimension(im,jm,1:nlm-1,nttend) ::   &
                     w_l2_f       ! d/dt (w_l2)   (m/s^2)

real (kind=dbl_kind), dimension(im,jm,1:nlm-1,nttend) ::     &
                     phi_l2_f, th_l2_f


real (kind = dbl_kind), dimension(im,jm,nttend) ::           &
                     th_l2_f_0, th_l2_f_nlm   ! d/dt (th_l2) at bottom
                                              ! and top layers (K/s)

real (kind=dbl_kind), dimension(im,jm,0:nlm,nttend) ::     &
                     qvpc_l2_m_l2_f, & ! d/dt (m_l2*qvpc_l2) (kg/m^2/eta/s)
                     qr_l2_m_l2_f      ! d/dt (m_l2*qr_l2) (kg/m^2/eta/s)

real (kind=dbl_kind), dimension(im,jm,0:nlm,nttend) ::       &
                     qc_horiz_adv_ab3,  &  ! horiz. cloud water advection
                     qc_vert_adv_ab3


!
! Declare Euler-forward tendencies of prognostic variables
!
real (kind = dbl_kind), dimension(im,jm,nlm) ::              &
                     m_f_ef

real (kind = dbl_kind), dimension(im,jm,0:nlm) ::            &
                     qvpc_l2_m_l2_f_ef,   &
                     qc_vert_adv_ef   ! vert. cloud water advection

real (kind=dbl_kind), dimension(im,jm,1:nlm-1) ::            &
                     th_l2_f_ef


!
! Declare layer-center diagnostic variables
!
real (kind = dbl_kind), dimension(im,jm,nlm) ::              &
                     ke_horiz,  &   ! contribution to kinetic energy
                                    ! from horizontal velocity (J/kg)
                     pv,        &   ! vertical component of potential 
                                    ! vorticity (m^2 eta/kg/s)
                     pl,        &   ! pressure (Pa)
                     T_temp         ! temperature (K)



!
! Declare layer-edge diagnostic variables
!
real (kind = dbl_kind), dimension(im,jm,0:nlm) ::            &
                     m_l2,   &      ! pseudo-density (kg/m^2/eta)
                     qvpc_l2,  &  ! water vapor + cloud mixing ratio (kg/kg)
                     qc_l2,    &  ! cloud-water mix. ratio (kg/kg)
                     qr_l2,    &  ! rain-water mix. ratio (kg/kg)
                     qT_l2        ! total water mix. ratio (kg/kg)


real (kind = dbl_kind), dimension(im,jm,1:nlm-1) ::          &
                     eta_dot_l2, &  ! generalized vert. velocity (eta/s)
                                    ! at layer edges.
                     eta_dot_l2_sigma, eta_dot_l2_theta,             &
                     eta_dot_l2_Q_diab, eta_dot_l2_tau,              &
                     eta_dot_l2_regrid, eta_dot_l2_prime
                                    ! generalized vert. velocity associated
                                    ! with sigma-, theta-, Q_diab-forcings,
                                    ! nudging of F toward eta_l2, 
                                    ! vertical regridding, and setting
                                    ! of F to F_set.
                                    ! Note:  eta_dots equal zero at top
                                    ! and bottom surfaces so are not
                                    ! declared at those levels.

real (kind = dbl_kind), dimension(im,jm,1:nlm-1) ::          &
                     F_set,   &  ! set point for F_th_sgma
                     dF_deta     ! vertical derivative of F_th_sgma

logical (kind = log_kind), dimension(im,jm,1:nlm-1) ::               &
           smooth_point,  &  ! true if grid point smoothed at current time step
           in_cloud          ! true when Q_diabatic (latent htg) >~ 0

!
! Declare the diabatic heating rate (at layer edges)
!
real (kind = dbl_kind), dimension(im,jm,0:nlm-1) ::                  &
           Q_diabatic,   &    ! total diabatic heating rate (W/kg)
           Q_diab_rain_evap   ! diabatic heating rate due to rain
                              ! evaporation (W/kg)



!
! Declare surface value of Exner function.
!
real (kind = dbl_kind), dimension(im,jm) ::                          &
           P_exnr_l2_0   ! Surface Exner function (J/kg/K)



!
! Declare diagnostic vertical velocity at lowest layer edge 
! (i.e. at the lower boundary)  (Note: This may be non-zero if there
! is surface topography).  The upper boundary is assumed flat so the
! vertical velocity there is always zero.
!
real (kind = dbl_kind), dimension(im,jm) ::                  &
                     w_l2_0         ! vertical velocity at surface (m/s)



!
! Declare surface geopotential and model "thickness" (z_top - z_surface)
!
real (kind = dbl_kind), dimension (im,jm) ::                 &
                     phis,  &        ! surface geopotential (J/kg)
                     h_1, inv_h_1    ! z_top - z_surface (m)


!
! Declare surface-precipitation budgets
!
real (kind = dbl_kind), dimension (im,jm) ::                 &
                     RAIN_RATE,  &  ! surface rain rate (mm/day)
                     RAIN_SUM       ! Accumulated rain (mm)

! Rain terminal velocity (m/s)
real (kind = dbl_kind), dimension (im,jm,0:nlm) ::           &
                     Vt




contains


!======================================================================
! BEGINNING OF INIT_PROGNOSTICS
!======================================================================

subroutine init_prognostics

implicit none

! initialize prognostic arrays
u = c0
v = c0
m = c0
w_l2 = c0
phi_l2 = c0
th_l2 = c0
qvpc_l2_m_l2 = c0
qr_l2_m_l2 = c0

u_f = c0
v_f = c0
m_f = c0
w_l2_f = c0
phi_l2_f = c0
th_l2_f = c0
th_l2_f_0 = c0; th_l2_f_nlm = c0
qvpc_l2_m_l2_f = c0
qr_l2_m_l2_f = c0
m_f_ef = c0
qvpc_l2_m_l2_f_ef = c0
th_l2_f_ef = c0
qc_horiz_adv_ab3 = c0
qc_vert_adv_ab3 = c0; qc_vert_adv_ef = c0

! initialize eta_dot_l2 (for diabatic heating calculation)
eta_dot_l2 = c0


end subroutine init_prognostics

!======================================================================
! END OF INIT_PROGNOSTICS
!======================================================================


end module prognostics
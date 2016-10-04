module step

!-----------------------------------------------------------------
!   This module updates the diagnostic variables and then time-
!   steps the prognostic variables.
!   Euler forward is used for the first time step, followed by
!   Adams-Bashforth 2nd order for the second time step, and
!   subsequently by Adams-Bashforth 3rd order time stepping.
!   The "traditional" HPGF horizontal discretization is used.
!-----------------------------------------------------------------

use kinds
use physical_parameters
use model_parameters
use eta_coordinate
use prognostics
use mass_tendency
use eta_dot_diagnosis
use momentum_tendency
use physics
use tracer_transport
use semi_implicit_solvers

implicit none
save


real (kind = dbl_kind), dimension(im,jm,0:nlm) ::                    &
           qc_l2_n4   ! Storage of current qc_l2 for output module




contains


!======================================================================
! BEGINNING OF UPDATE_DIAGNOSTICS
!======================================================================

subroutine update_diagnostics ( tau, w1, w2, w3 )

!----------------------------------------------------------------------
! PURPOSE:
!   Updates diagnostic variables and calls subroutines to calculate
!   tendencies of prognostic variables in preparation for
!   time-stepping the prognostic variables.
!----------------------------------------------------------------------

implicit none

!-------------------------------------------------------------------
! INTENT IN 
!-------------------------------------------------------------------
real (kind = dbl_kind), intent(in) :: tau         ! time in hours
real (kind = dbl_kind), intent(in) :: w1, w2, w3  ! Time stepping
                                                  ! weights

!-------------------------------------------------------------------
! LOCAL
!-------------------------------------------------------------------
integer :: i,j,k

real (kind=dbl_kind), dimension(im,jm) ::                            &
           f1, f2, f3, f4     ! working variables

real (kind=dbl_kind), dimension(im,jm) ::                            &
           d_i_phis, d_j_phis ! differences of phis in i and j directions

real (kind=dbl_kind), dimension(im,jm) ::                            &
           inv_aa11, aa21, aa01, aa12, aa10, bb  ! coefficients used in
                                                 ! surface pressure calc

!
! Declare layer-center diagnostic variables
!
real (kind = dbl_kind), dimension(im,jm,nlm) ::                      &
           th,        &   ! potential temperature (K)
           phi,       &   ! geopotential (J/kg)
           m_eta_dot, &   ! vertical mass flux used for calculation of
                          ! vertical advection of vertical momentum
                          ! (kg/m^2/s)
           rho,       &   ! dry-air density (kg/m^3)
           P_exnr         ! Exner function (J/kg/K)


! Layer-edge diagnostic variable
real (kind=dbl_kind), dimension (im,jm,0:nlm) ::                     &
          m_PI_l2,   &      ! Exner function weighted pseudo-density at layer
                            ! edges (J/m^2/K)
          rho_l2            ! layer-edge dry-air density (kg/m^3)


!
! Declare moisture-related variables
!

real (kind = dbl_kind), dimension(0:nlm) ::                                &
           T_temp_l2_iter,  &  ! layer-edge temperature (K)
           rho_l2_k,        &  ! layer-edge dry-air density (kg/m^3)
           qsat_iter,       &  ! layer-edge saturation mix. ratio (kg/kg)
           esat_l2_iter,    &  ! layer-edge saturation vapor-pressure
           qc_l2_iter,      &  ! layer-edge cloud water mixing ratio (kg/kg)
           th_l2_iter, th_l2_int,  &
           th_l2_n, m_l2_n, qc_l2_n, m_PI_l2_n,  &
           m_l2_np1, qvpc_l2_np1, phi_l2_np1

real (kind = dbl_kind), dimension(im,jm,0:nlm) ::                          &
           esat,    &          ! layer-edge saturation vapor-pressure
           qsat,    &          ! layer-edge saturation mix. ratio (kg/kg)
           pl_l2,   &          ! layer-edge pressure (Pa)
           T_temp_l2           ! layer-edge temperature (K)

real (kind = dbl_kind), dimension(0:nlm-1) ::                              &
           f1_fac1, f1_fac2, f1_fac3

real (kind = dbl_kind), dimension(0:nlm-1) ::                              &
           F1vec, a_elem, step_change_th

real (kind = dbl_kind) ::                                                  &
           F2vec, dF2_dth, dF2_dT, step_change_T

real (kind = dbl_kind) ::                                                  &
           P_exnr_l2,   &   ! Exner function at top model layer
           rho_l2_pt,   &   ! Temporary holder for rho_l2
           fpt1, fpt2, fpt3 ! working variables

real (kind = dbl_kind) ::                                                  &
           dqsatdth
           

real (kind = dbl_kind) ::                                                  &
           max_change_theta, max_change_T

real (kind = dbl_kind), parameter ::                                       &
           converge_cond_theta_T = 1.0E-10,    &
           r_fac_param1 = c1/(r_fac+c1)



! Temporary line
integer :: max_number_iter


!
! Declare mass flux variables to be used in the "3rd-order" Takacs
! advection schemes for continuity and theta advection
!
real (kind = dbl_kind), dimension(im,jm,nlm) ::                      &
           F_x,          &  ! "3rd-order" mass flux in x-direction
                            ! (colocated with u-points) (kg/s/m/eta)
           F_y              ! "3rd-order" mass flux in y-direction
                            ! (colocated with v-points) (kg/s/m/eta)

real (kind = dbl_kind), dimension(im,jm,1:nlm-1) ::                  &
           F_eta_l2_AB3, &  ! vertical mass flux (kg/s/m^2)
           F_eta_l2_EF      ! (colocated with eta_dot_l2 points)
                            ! AB3 and Euler forward fluxes respectively

real (kind=dbl_kind), dimension(im,jm,1:nlm-1) ::                    &
           exp_horiz_phi_advec, &  ! explicitly weighted phi advection
           eta_dot_sigma_coeff     ! coefficient for calculating
                                   ! eta_dot_l2_sigma

real (kind = dbl_kind), dimension(im,jm,nlm) ::                      &
           m_f_temp        ! multi-use mass tendency

real (kind = dbl_kind), dimension(im,jm,1:nlm-1) ::                  &
           w_l2_f_temp     ! multi-use w_l2 tendency

integer :: km0, km1     ! indices for k level (km0 = k-0) and the k-1
                        ! level (km1 = k-1)  (they toggle between 0 and 1)

real (kind=dbl_kind), dimension(im,jm,0:1) ::                        &
          f1_bi_lev     ! working variable



!
! Declare the vertical pressure gradient force to be used in the
! vertical momentum equation and the horizontal momentum equation
! (as part of the horizontal pressure gradient force)
! (lives at w_l2 points)
!
real (kind = dbl_kind), dimension(im,jm,0:nlm-1) ::                  &
           vpgf_l2          ! (m/s^2)



!
! Declare vertical velocity advection for top and bottom half-layers
!
real (kind = dbl_kind), dimension(im,jm) ::                          &
           w_advect_0       ! advection of vertical velocity
                            ! associated with the bottom half-layer

!
! Declare time tendency of vertical velocity at surface
!
real (kind = dbl_kind), dimension(im,jm) ::                          &
           w_l2_0_f         ! (m/s^2)


!
! Declare advection of horizontal velocity plus Coriolis acceleration
! in the bottom layer
!
real (kind = dbl_kind), dimension(im,jm) ::                          &
           u_advect_coriol_1,   &  ! advection of the x-component of
                                   ! velocity plus Coriolis accel.
                                   ! in the bottom layer
           v_advect_coriol_1       ! advection of the y-component of
                                   ! velocity plus Coriolis accel.
                                   ! in the bottom layer


!
! Declare variables dealing with the iterative solution for the
! Exner function at the surface (and moisture -- iter and max_iter).
!
integer :: iter                        ! iteration counter

integer, parameter :: max_iter = 100   ! max. number of iterations allowed

real (kind = dbl_kind) :: max_error    ! maximum error between
                                       ! subsequent iterations

real (kind = dbl_kind), parameter ::                                 &
           converge = 1.0E-14_dbl_kind  ! convergence criterion

real (kind = dbl_kind), dimension(im,jm) ::                          &
           vpgf_l2_0_old    ! Holder for old values for comparison
                            ! with new values in Gauss-Seidel iteration


!
! "New moisture" variables
!
real (kind = dbl_kind), dimension(im,jm,nlm) ::             &
        moist_corr1, &     ! moisture correction for equation of
                           ! state at layer centers (-)
        moist_corr1_hpgf   ! moisture correction for horizontal
                           ! pressure gradient force at layer centers (-)

real (kind = dbl_kind), dimension(im,jm,0:nlm) ::           &
        moist_corr2_l2     ! moisture correction for vertical
                           ! pressure gradient force at layer edges (-)

!
! Kessler microphysics variables
!
real (kind = dbl_kind), dimension(0:nlm) ::                 &
        qr_l2_k,     &     ! Rain water mixing ratio (kg/kg)
        qc_l2_k,     &     ! Cloud water mixing ration (kg/kg)
        sedim_fac,   &     ! Factor used in sedimentation calculation
        dqr_sedim          ! Change in qr_l2 due to sedimentation

real (kind = dbl_kind) ::                                   &
        ppt,              &  ! Surface rain rate (mm/s)
        qrprod_collect,   &  ! Rain production from collection
        qrprod_auto_conv, &  ! Rain production from autoconversion
        qv_l2,            &  ! Water vapor mixing ratio (kg/kg)
        prod,             &  ! Factor related to supersaturation
        evap,             &  ! time-step rain evaporation (kg/kg)
        rcgs, qrr,        &  ! misc. factors
        rho_l2_k0            ! alternative density at bottom level
                             ! used for sedimentation

real (kind = dbl_kind), parameter ::                        &
        qc_crit = 0.001_dbl_kind  ! "Critical" cloud-water mixing ratio
                                  ! for rain formation from autoconversion

! Hole filling variables for rain
real (kind = dbl_kind) :: neg_qr, total_qr, PHI_qr





! Interpolate geopotential to layer centers
phi(:,:,1) = p5 * ( phi_l2(:,:,1,n4) + phis(:,:) )
do k = 2, nlm-1
   phi(:,:,k) = p5 * ( phi_l2(:,:,k,n4) + phi_l2(:,:,k-1,n4) )
end do
phi(:,:,nlm) = p5 * ( grav*z_top + phi_l2(:,:,nlm-1,n4) )


! Interpolate theta to layer centers, diagnose layer-center density
! from pseudo-density, and diagnose layer-center temperature, pressure,
! and Exner function.
! Note:  Density rho is dry-air density
th(:,:,1) = p5 * ( th_l2(:,:,1,n4) + th_l2(:,:,0,n4) )
rho(:,:,1) = m(:,:,1,n4) * grav * d_eta(1) /                            &
             ( phi_l2(:,:,1,n4) - phis(:,:) )
moist_corr1(:,:,1) = c1 + invepsilon*p5*( max(c0,qvpc_l2(:,:,1)) +      &
             max(c0,qvpc_l2(:,:,0)) - qc_l2(:,:,1) - qc_l2(:,:,0) )
T_temp(:,:,1) = th(:,:,1)**gamma * (rho(:,:,1)*gas_const_R*             &
                   moist_corr1(:,:,1)*inv_p0_ref)**gamma_mns1
pl(:,:,1) = rho(:,:,1)*gas_const_R*T_temp(:,:,1)*moist_corr1(:,:,1)
P_exnr(:,:,1) = spec_heat_cp * (pl(:,:,1)*inv_p0_ref) ** kappa
do k = 2,nlm-1
   th(:,:,k) = p5 * ( th_l2(:,:,k,n4) + th_l2(:,:,k-1,n4) )
   rho(:,:,k) = m(:,:,k,n4) * grav * d_eta(k) /                         &
                ( phi_l2(:,:,k,n4) - phi_l2(:,:,k-1,n4) )
   moist_corr1(:,:,k) =   c1 + invepsilon*p5*( max(c0,qvpc_l2(:,:,k)) + &
          max(c0,qvpc_l2(:,:,k-1)) - qc_l2(:,:,k) - qc_l2(:,:,k-1) )
   T_temp(:,:,k) = th(:,:,k)**gamma * (rho(:,:,k)*gas_const_R*          &
                      moist_corr1(:,:,k)*inv_p0_ref)**gamma_mns1
   pl(:,:,k) = rho(:,:,k)*gas_const_R*T_temp(:,:,k)*moist_corr1(:,:,k)
   P_exnr(:,:,k) = spec_heat_cp * (pl(:,:,k)*inv_p0_ref) ** kappa
end do
th(:,:,nlm) = p5 * ( th_l2(:,:,nlm,n4) + th_l2(:,:,nlm-1,n4) )
rho(:,:,nlm) = m(:,:,nlm,n4) * grav * d_eta(nlm) /                      &
             ( grav*z_top - phi_l2(:,:,nlm-1,n4) )
moist_corr1(:,:,nlm) = c1 + invepsilon*p5*( max(c0,qvpc_l2(:,:,nlm)) +  &
       max(c0,qvpc_l2(:,:,nlm-1)) - qc_l2(:,:,nlm) - qc_l2(:,:,nlm-1) )
T_temp(:,:,nlm) = th(:,:,nlm)**gamma * (rho(:,:,nlm)*gas_const_R*       &
                     moist_corr1(:,:,nlm)*inv_p0_ref)**gamma_mns1
pl(:,:,nlm) = rho(:,:,nlm)*gas_const_R*T_temp(:,:,nlm)*                 &
                 moist_corr1(:,:,nlm)
P_exnr(:,:,nlm) = spec_heat_cp * (pl(:,:,nlm)*inv_p0_ref) ** kappa


! Store current value of cloud field
qc_l2_n4 = qc_l2

! Calculate current value of qT_l2
qT_l2(:,:,:) = qvpc_l2(:,:,:) + qr_l2(:,:,:)



! Calculate effect of subgrid scale turbulence on momentum and
! potential temperature tendencies
call calc_subgrid_turbulent_flx_div(th,phi,rho,n4)



! Calculate horizontal advection of cloud water
! Now, F_x and F_y are calculated in subroutine get_horiz_div_mvqc
! located in module tracer_transport
! Note:  mDqc_Dt at this point will equal horiz_div(m*v*qc)
call get_horiz_div_mvqc ( F_x, F_y, qc_horiz_adv_ab3(:,:,:,n3_f) )




! Diagnose eta_dot_l2 and time-step potential temperature and geopotential
! Note:  These are intermediate values in preparation for semi-implicit
!        time stepping.  Does not include smoothing processes.
call get_eta_dot_l2_1 ( w1,w2,w3,P_exnr,th,phi,F_x,F_y,m_PI_l2,      &
                        exp_horiz_phi_advec,eta_dot_sigma_coeff )


! Calculate intermediate values of th_l2 at top and bottom levels
th_l2(:,:,0,n3) = th_l2(:,:,0,n4) + dt * ( w3*th_l2_f_0(:,:,n3_f) +  &
                  w2*th_l2_f_0(:,:,n2_f) + w1*th_l2_f_0(:,:,n1_f) +  &
                     F_turb_th_l2(:,:,0) )

th_l2(:,:,nlm,n3) = th_l2(:,:,nlm,n4) + dt *                         &
           ( w3*th_l2_f_nlm(:,:,n3_f) + w2*th_l2_f_nlm(:,:,n2_f) +   &
                                        w1*th_l2_f_nlm(:,:,n1_f) )



! Obtain intermediate values of mass based on horizontal advection
! and known values of generalized vertical velocity
call get_mf_1 ( eta_dot_l2_sigma,eta_dot_l2_theta+eta_dot_l2_tau,    &
                eta_dot_l2_Q_diab,m_l2,F_x,F_y,F_eta_l2_AB3,         &
                F_eta_l2_EF,m_f(:,:,:,n3_f),m_f_ef,m_f_temp )


! Time-step m based on current tendencies
m(:,:,:,n3) = m(:,:,:,n4) + dt * ( w3*m_f(:,:,:,n3_f) +              &
           w2*m_f(:,:,:,n2_f) + w1*m_f(:,:,:,n1_f) + m_f_ef(:,:,:) + &
           m_f_temp(:,:,:) )




! Diagnose m_eta_dot (vertical mass flux interpolated to layer centers)
! based on known values of generalized vertical velocity
! Note:  This is only used for vertical advection of vertical momentum
!        w_l2.  For now, advection of momentum due to diabatic effects
!         (i.e., eta_dot_l2_Q_diab) are treated with AB3 time stepping
km0 = mod(1,2)
km1 = mod(0,2)
f1_bi_lev(:,:,km0) = p5 * m_l2(:,:,1) *                              &
     ( eta_dot_l2_sigma(:,:,1) + eta_dot_l2_theta(:,:,1) +           &
       eta_dot_l2_Q_diab(:,:,1) + eta_dot_l2_tau(:,:,1) )
m_eta_dot(:,:,1) = f1_bi_lev(:,:,km0)
do k = 2, nlm-1
   km0 = mod(k,2)
   km1 = mod(k-1,2)
   f1_bi_lev(:,:,km0) = p5 * m_l2(:,:,k) *                           &
        ( eta_dot_l2_sigma(:,:,k) + eta_dot_l2_theta(:,:,k) +        &
          eta_dot_l2_Q_diab(:,:,k) + eta_dot_l2_tau(:,:,k) )
   m_eta_dot(:,:,k) = f1_bi_lev(:,:,km0) + f1_bi_lev(:,:,km1)
end do
km0 = mod(nlm,2)
km1 = mod(nlm-1,2)
m_eta_dot(:,:,nlm) = f1_bi_lev(:,:,km1)



! Diagnose vertical velocity at lowest layer edge (i.e. at the lower boundary)
w_l2_0(:,:) = p5*invdx *  ( u(ip1(:),:,1,n4) *                       &
                           (phis(ip1(:),:)-phis(:,:))*invgrav   +    &
                           u(:,:,1,n4) *                             &
                           (phis(:,:)-phis(im1(:),:))*invgrav ) +    &
              p5*invdy *  ( v(:,jp1(:),1,n4) *                       &
                           (phis(:,jp1(:))-phis(:,:))*invgrav   +    &
                           v(:,:,1,n4) *                             &
                           (phis(:,:)-phis(:,jm1(:)))*invgrav )



! Calculate moisture correction for pressure gradient force
! at layer centers and at layer edges
do k = 1, nlm
   moist_corr1_hpgf(:,:,k) = moist_corr1(:,:,k) /                       &
      ( c1 + p5*(max(c0,qT_l2(:,:,k))+max(c0,qT_l2(:,:,k-1))) )
end do
moist_corr2_l2(:,:,:) = (c1+invepsilon*                                 &
     (max(c0,qvpc_l2(:,:,:))-qc_l2(:,:,:))) / (c1+max(c0,qT_l2(:,:,:)))


! Calculate the vertical pressure gradient force for levels 1 thru nlm-1
! (First guess at surface (k=0) to be calculated below)
do k = 1, nlm-1
   vpgf_l2(:,:,k) = grav * moist_corr2_l2(:,:,k) * th_l2(:,:,k,n4) * &
                    ( P_exnr(:,:,k+1) - P_exnr(:,:,k) ) /            &
                    ( phi(:,:,k+1) - phi(:,:,k) )
end do



! Obtain intermediate values of vertical velocity (w) based on horizontal
! advection and known values of generalized vertical velocity
call get_wl2f_1 ( u(:,:,:,n4),v(:,:,:,n4),m_eta_dot,w_l2(:,:,:,n4),     &
           w_l2_0,m_l2,m(:,:,:,n4),th_l2(:,:,:,n4),P_exnr,phi,          &
           vpgf_l2,phi_l2(:,:,:,n4),phis,w_advect_0,w_l2_f(:,:,:,n3_f), &
           w_l2_f_temp )


! Time-step w_l2 based on current tendencies
w_l2(:,:,:,n3) = w_l2(:,:,:,n4) + dt * ( w3*w_l2_f(:,:,:,n3_f) +        &
                 w2*w_l2_f(:,:,:,n2_f) + w1*w_l2_f(:,:,:,n1_f) +        &
                                              w_l2_f_temp(:,:,:) )


! Perform trapezoidal (semi-implicit) time-stepping of m, w_l2
! and phi_l2 for handling vertically propagating acoustic waves
call semi_implicit_solve_m_w_phi ( phi, exp_horiz_phi_advec,            &
                  eta_dot_sigma_coeff, moist_corr1, moist_corr2_l2 )


! Continue diagnosing generalized vertical velocity and finalizing
! geopotential field and adiabatic effects on potential temperature
call get_eta_dot_l2_2 ( w3,P_exnr,th,phi,m_PI_l2 )


! Update m and w_l2 with advection due to newly calculated
! eta_dot_l2_regrid and eta_dot_l2_prime

call get_mf_2 ( eta_dot_l2_prime,eta_dot_l2_regrid,m_l2,w3,          &
                F_eta_l2_AB3,F_eta_l2_EF,m_f(:,:,:,n3_f),            &
                m_f_temp,m_f_ef )

! Update m based on newly computed tendencies
m(:,:,:,n3) = m(:,:,:,n3) + dt * ( w3*m_f_temp(:,:,:) + m_f_ef(:,:,:) )



! Rediagnose m_eta_dot (vertical mass flux interpolated to layer centers)
km0 = mod(1,2)
km1 = mod(0,2)
f1_bi_lev(:,:,km0) = p5 * m_l2(:,:,1) *                              &
     ( eta_dot_l2_regrid(:,:,1) + eta_dot_l2_prime(:,:,1) )
m_eta_dot(:,:,1) = f1_bi_lev(:,:,km0)
do k = 2, nlm-1
   km0 = mod(k,2)
   km1 = mod(k-1,2)
   f1_bi_lev(:,:,km0) = p5 * m_l2(:,:,k) *                           &
        ( eta_dot_l2_regrid(:,:,k) + eta_dot_l2_prime(:,:,k) )
   m_eta_dot(:,:,k) = f1_bi_lev(:,:,km0) + f1_bi_lev(:,:,km1)
end do
km0 = mod(nlm,2)
km1 = mod(nlm-1,2)
m_eta_dot(:,:,nlm) = f1_bi_lev(:,:,km1)


! Use these values of m_eta_dot to update vertical advection of w_l2
call get_wl2f_2 ( m_eta_dot,w_l2(:,:,:,n4),w_l2_0,m_l2,w_advect_0,   &
                  w_l2_f(:,:,:,n3_f),w_l2_f_temp )

! Update w_l2 based on newly updated tendency
w_l2(:,:,:,n3) = w_l2(:,:,:,n3) + dt * w3 * w_l2_f_temp(:,:,:)




!
! Get tendencies of prognostic variables u and v
!
call get_uf_vf ( u(:,:,:,n4),v(:,:,:,n4),eta_dot_l2,m(:,:,:,n4),     &
                 m_l2,th,P_exnr,phi,vpgf_l2,phi_l2(:,:,:,n4),        &
                 phis,moist_corr1_hpgf,pv,ke_horiz,                  &
                 u_advect_coriol_1,v_advect_coriol_1,                &
                 u_f(:,:,:,n3_f),v_f(:,:,:,n3_f) )



!
! Calculate "first guess" of the Exner function at the surface
! assuming the time tendency of the vertical velocity at the
! surface is zero.
! (In the case of no orography, this will be the final answer).
!
P_exnr_l2_0(:,:) = P_exnr(:,:,1) + ( phi(:,:,1) - phis(:,:) ) *      &
                   ( w_advect_0(:,:) + grav ) /                      &
                    (grav * moist_corr2_l2(:,:,0) * th_l2(:,:,0,n4))


! Calculate "first guess" of vertical pressure gradient force at surface
vpgf_l2(:,:,0) = grav * moist_corr2_l2(:,:,0) * th_l2(:,:,0,n4) *    &
                    ( P_exnr(:,:,1) - P_exnr_l2_0(:,:) ) /           &
                    ( phi(:,:,1) - phis(:,:) )


!
! Use Gauss-Seidel relaxation method to determine Exner function
! at the surface.
!

! Calculate differences in surface geopotential
do j = 1, jm
   do i = 1, im
      d_i_phis(i,j) = phis(i,j) - phis(im1(i),j)
      d_j_phis(i,j) = phis(i,j) - phis(i,jm1(j))
   end do
end do

! Calculate coefficients (aa's) and right-hand side of eqn (bb)
do j = 1, jm
   do i = 1, im

      inv_aa11(i,j) = c1 / ( c1 + inv8*invgrav**2 *                  &
         ( invdx**2 * (d_i_phis(ip1(i),j)**2+d_i_phis(i,j)**2) +     &
           invdy**2 * (d_j_phis(i,jp1(j))**2+d_j_phis(i,j)**2) ) )
      aa21(i,j) = inv8*(invgrav*invdx*d_i_phis(ip1(i),j))**2
      aa01(i,j) = inv8*(invgrav*invdx*d_i_phis(i,j))**2
      aa12(i,j) = inv8*(invgrav*invdy*d_j_phis(i,jp1(j)))**2
      aa10(i,j) = inv8*(invgrav*invdy*d_j_phis(i,j))**2

      f1(i,j) = -invdx*p25*(moist_corr1_hpgf(ip1(i),j,1)+                  &
         moist_corr1_hpgf(i,j,1))*(th(ip1(i),j,1)+th(i,j,1))*              &
         (P_exnr(ip1(i),j,1)-P_exnr(i,j,1))+invdx*invgrav*p25*             &
         (vpgf_l2(ip1(i),j,1)+vpgf_l2(i,j,1))*(phi_l2(ip1(i),j,1,n4)-      &
         phi_l2(i,j,1,n4))-u_advect_coriol_1(ip1(i),j)
      f2(i,j) = -invdx*p25*(moist_corr1_hpgf(i,j,1)+                       &
         moist_corr1_hpgf(im1(i),j,1))*(th(i,j,1)+th(im1(i),j,1))*         &
         (P_exnr(i,j,1)-P_exnr(im1(i),j,1))+invdx*invgrav*p25*             &
         (vpgf_l2(i,j,1)+vpgf_l2(im1(i),j,1))*(phi_l2(i,j,1,n4)-           &
         phi_l2(im1(i),j,1,n4))-u_advect_coriol_1(i,j)
      f3(i,j) = -invdy*p25*(moist_corr1_hpgf(i,jp1(j),1)+                  &
         moist_corr1_hpgf(i,j,1))*(th(i,jp1(j),1)+th(i,j,1))*              &
         (P_exnr(i,jp1(j),1)-P_exnr(i,j,1))+invdy*invgrav*p25*             &
         (vpgf_l2(i,jp1(j),1)+vpgf_l2(i,j,1))*(phi_l2(i,jp1(j),1,n4)-      &
         phi_l2(i,j,1,n4))-v_advect_coriol_1(i,jp1(j))
      f4(i,j) = -invdy*p25*(moist_corr1_hpgf(i,j,1)+                       &
         moist_corr1_hpgf(i,jm1(j),1))*(th(i,j,1)+th(i,jm1(j),1))*         &
         (P_exnr(i,j,1)-P_exnr(i,jm1(j),1))+invdy*invgrav*p25*             &
         (vpgf_l2(i,j,1)+vpgf_l2(i,jm1(j),1))*(phi_l2(i,j,1,n4)-           &
         phi_l2(i,jm1(j),1,n4))-v_advect_coriol_1(i,j)
      bb(i,j) = - grav - w_advect_0(i,j) - p5*invgrav*                     &
             ( invdx*(f1(i,j)*d_i_phis(ip1(i),j)+f2(i,j)*d_i_phis(i,j)) +  &
               invdy*(f3(i,j)*d_j_phis(i,jp1(j))+f4(i,j)*d_j_phis(i,j)) )

   end do
end do

!
! Perform Gauss-Seidel iteration
!
iter = 0
max_error = 1.E+10_dbl_kind  ! Note:  Set max_error equal to zero to use
                             !        the assumption that the time-tendency
                             !        of w_l2 at the surface is zero

do while (max_error .ge. converge)

   vpgf_l2_0_old(:,:) = vpgf_l2(:,:,0)
   do j = 1, jm
      do i = 1, im
         vpgf_l2(i,j,0) = inv_aa11(i,j) * ( bb(i,j) -                      &
           aa21(i,j)*vpgf_l2(ip1(i),j,0) - aa01(i,j)*vpgf_l2(im1(i),j,0) - &
           aa12(i,j)*vpgf_l2(i,jp1(j),0) - aa10(i,j)*vpgf_l2(i,jm1(j),0) )
      end do
   end do

   ! Calculate max_error
   max_error = c0
   do j = 1,jm
      do i = 1,im
         max_error = max(max_error,abs(vpgf_l2(i,j,0)-vpgf_l2_0_old(i,j)))
      end do
   end do

   iter = iter + 1
   if (iter .ge. max_iter) then
      print *
      print *,   &
         "Max. iteration count reached in subroutine update diagnostics."
      print *, "Max. error = ", max_error
      print *
      print *, "PROGRAM STOPPING "
      print *
      stop
   end if

end do  ! while (max_error .ge. converge)


! Temporary line
print "(A39,I5)", "Surface pressure solver iter =", iter


! Now that we have vpgf_l2 at the surface, update u and v tendencies
! at k = 1 to include HPGF and calculate surface Exner function
f1(:,:) = - p25 * ( moist_corr1_hpgf(:,:,1) +                        &
   moist_corr1_hpgf(im1(:),:,1) ) * ( th(:,:,1) + th(im1(:),:,1) ) * &
      ( P_exnr(:,:,1) - P_exnr(im1(:),:,1) )

f2(:,:) = invgrav * p5 *                                             &
             ( p5*(vpgf_l2(:,:,1)+vpgf_l2(im1(:),:,1)) *             &
               (phi_l2(:,:,1,n4)-phi_l2(im1(:),:,1,n4))   +          &
               p5*(vpgf_l2(:,:,0)+vpgf_l2(im1(:),:,0)) *             &
               (phis(:,:)-phis(im1(:),:))  )

u_f(:,:,1,n3_f) = - u_advect_coriol_1(:,:) +                         &
                   invdx * ( f1(:,:) + f2(:,:) )


f1(:,:) = - p25 * ( moist_corr1_hpgf(:,:,1) +                        &
   moist_corr1_hpgf(:,jm1(:),1) ) * ( th(:,:,1) + th(:,jm1(:),1) ) * &
      ( P_exnr(:,:,1) - P_exnr(:,jm1(:),1) )

f2(:,:) = invgrav * p5 *                                             &
             ( p5*(vpgf_l2(:,:,1)+vpgf_l2(:,jm1(:),1)) *             &
               (phi_l2(:,:,1,n4)-phi_l2(:,jm1(:),1,n4))   +          &
               p5*(vpgf_l2(:,:,0)+vpgf_l2(:,jm1(:),0)) *             &
               (phis(:,:)-phis(:,jm1(:)))  )

v_f(:,:,1,n3_f) = - v_advect_coriol_1(:,:) +                         &
                   invdy * ( f1(:,:) + f2(:,:) )

! Diagnose the time tendency of the vertical velocity at the surface
w_l2_0_f(:,:) = p5*invdx * ( u_f(ip1(:),:,1,n3_f)*                   &
                            (phis(ip1(:),:)-phis(:,:))*invgrav +     &
            u_f(:,:,1,n3_f)*(phis(:,:)-phis(im1(:),:))*invgrav ) +   &
                p5*invdy * ( v_f(:,jp1(:),1,n3_f)*                   &
                            (phis(:,jp1(:))-phis(:,:))*invgrav +     &
            v_f(:,:,1,n3_f)*(phis(:,:)-phis(:,jm1(:)))*invgrav )

! Calculate Exner function at the surface based on w_l2_0_f
P_exnr_l2_0(:,:) = P_exnr(:,:,1) + ( phi(:,:,1) - phis(:,:) ) *      &
                   ( w_l2_0_f(:,:) + w_advect_0(:,:) + grav ) /      &
                   (grav * moist_corr2_l2(:,:,0) * th_l2(:,:,0,n4))





call get_tracer_tendencies ( F_x, F_y, F_eta_l2_AB3, F_eta_l2_EF,    &
                                eta_dot_l2 )


!
! Time-step qvpc_l2, qr_l2 and th_l2 at model top and bottom based on current
! (i.e., final) values of eta_dot_l2 and Q_diabatic
! Note reversed order of n3 and n4 compared to subroutine step_dynamics
! (n3 is now n+1 time step value)
!
qvpc_l2_m_l2(:,:,:,n3) = qvpc_l2_m_l2(:,:,:,n4) + dt *               &
   ( w3*qvpc_l2_m_l2_f(:,:,:,n3_f) + w2*qvpc_l2_m_l2_f(:,:,:,n2_f) + &
     w1*qvpc_l2_m_l2_f(:,:,:,n1_f) + qvpc_l2_m_l2_f_ef(:,:,:) +      &
        F_turb_qvpc_l2(:,:,:)*m_l2(:,:,:) )
qr_l2_m_l2(:,:,:,n3) = qr_l2_m_l2(:,:,:,n4) + dt *                   &
   ( w3*qr_l2_m_l2_f(:,:,:,n3_f) + w2*qr_l2_m_l2_f(:,:,:,n2_f) +     &
     w1*qr_l2_m_l2_f(:,:,:,n1_f) + F_turb_qr_l2(:,:,:)*m_l2(:,:,:) )


! Calculate current value of Q_diabatic at surface and update surface th_l2
! Now add effects of subgrid-scale cloud diffusion
Q_diabatic(:,:,0) = L_vap * ( ( w3*qc_horiz_adv_ab3(:,:,0,n3_f) +        &
   w2*qc_horiz_adv_ab3(:,:,0,n2_f) + w1*qc_horiz_adv_ab3(:,:,0,n1_f) ) / &
      m_l2(:,:,0) + F_turb_qc_l2(:,:,0) )
th_l2(:,:,0,n3) = th_l2(:,:,0,n3) + dt * Q_diabatic(:,:,0)/P_exnr(:,:,1)


! Calculate diabatic heating associated with horizontal cloud-
! water advection and update theta accordingly
do k = 1, nlm-1
   Q_diabatic(:,:,k) = L_vap * ( ( w3*qc_horiz_adv_ab3(:,:,k,n3_f) +        &
      w2*qc_horiz_adv_ab3(:,:,k,n2_f) + w1*qc_horiz_adv_ab3(:,:,k,n1_f) ) / &
         m_l2(:,:,k) + F_turb_qc_l2(:,:,k) )
   f1(:,:) = m_l2(:,:,k)*Q_diabatic(:,:,k)*d_eta_l2(k) /                    &
                                 m_PI_l2(:,:,k)  ! Q_diabatic_forcing
   th_l2(:,:,k,n3) = th_l2(:,:,k,n3) + dt*f1(:,:)
end do



! Temporary line
max_number_iter = 0

!
! Begin iterative procedure to implicitely calculate effects of
! latent heating
!
do j = 1, jm
   do i = 1, im

      ! Assign *_n (current time-step n) values of prognostic variables
      ! (as well as m_l2, qvpc_l2 and qc_l2)
      m_l2_n(0) = m_l2(i,j,0)
      th_l2_n(0) = th_l2(i,j,0,n4)
      qc_l2_n(0) = qc_l2(i,j,0)
      m_PI_l2_n(0) = m_PI_l2(i,j,0)
      do k = 1, nlm-1
         m_l2_n(k) = m_l2(i,j,k)
         th_l2_n(k) = th_l2(i,j,k,n4)
         qc_l2_n(k) = qc_l2(i,j,k)
         m_PI_l2_n(k) = m_PI_l2(i,j,k)
      end do
      m_l2_n(nlm) = m_l2(i,j,nlm)
      th_l2_n(nlm) = th_l2(i,j,nlm,n4)
      qc_l2_n(nlm) = qc_l2(i,j,nlm)
      m_PI_l2_n(nlm) = m_PI_l2(i,j,nlm)

      ! Assign *_int (intermediate/starred) values of prognostic variables
      th_l2_int(0) = th_l2(i,j,0,n3)
      do k = 1, nlm-1
         th_l2_int(k) = th_l2(i,j,k,n3)
      end do
      th_l2_int(nlm) = th_l2(i,j,nlm,n3)

      ! Assign known *_np1 (n+1) values of prognostic variables
      m_l2_np1(0) = m(i,j,1,n3)
      qvpc_l2_np1(0) = qvpc_l2_m_l2(i,j,0,n3)/m_l2_np1(0)
      phi_l2_np1(0) = phis(i,j)
      do k = 1, nlm-1
         m_l2_np1(k) = p5 * ( m(i,j,k,n3)*d_eta(k) +                       &
                        m(i,j,k+1,n3)*d_eta(k+1) ) * inv_d_eta_l2(k)
         qvpc_l2_np1(k) = qvpc_l2_m_l2(i,j,k,n3)/m_l2_np1(k)
         phi_l2_np1(k) = phi_l2(i,j,k,n3)
      end do
      m_l2_np1(nlm) = m(i,j,nlm,n3)
      qvpc_l2_np1(nlm) = qvpc_l2_m_l2(i,j,nlm,n3)/m_l2_np1(nlm)
      phi_l2_np1(nlm) = grav*z_top



      ! Update Q_diabatic and th_l2_int to include latent heating associated
      ! with vertical cloud water advection by eta_dot_l2
      do k = 0, nlm-1
         fpt1 = L_vap * ( w3*qc_vert_adv_ab3(i,j,k,n3_f)+                  &
           w2*qc_vert_adv_ab3(i,j,k,n2_f)+w1*qc_vert_adv_ab3(i,j,k,n1_f) + &
             qc_vert_adv_ef(i,j,k) ) / m_l2_n(k)
         Q_diabatic(i,j,k) = Q_diabatic(i,j,k) + fpt1
         th_l2_int(k) = th_l2_int(k) + dt*fpt1*m_l2_n(k)*d_eta_l2(k) /     &
                           m_PI_l2_n(k)
      end do
      fpt1 = L_vap * ( w3*qc_vert_adv_ab3(i,j,nlm,n3_f)+                    &
        w2*qc_vert_adv_ab3(i,j,nlm,n2_f)+w1*qc_vert_adv_ab3(i,j,nlm,n1_f) + &
          qc_vert_adv_ef(i,j,nlm) ) / m_l2_n(nlm)
      th_l2_int(nlm) = th_l2_int(nlm) + dt*fpt1*m_l2_n(nlm)*d_eta_l2(nlm) / &
                          m_PI_l2_n(nlm)



      ! Assign initial iteration values (i) values of prognostic variables
      th_l2_iter(:) = th_l2_int(:)


      ! Find initial cloud amount (qc_l2_iter) throughout the column
      k = 0
      fpt1 = c1 + invepsilon*(max(c0,qvpc_l2_np1(k))-qc_l2_n(k))
                                ! moisture correction
                                ! for state eqn. (qc_l2 based on n4 step)
      rho_l2_pt = p0_ref*inv_gas_const_R*(P_exnr_l2_0(i,j)*                &
                     inv_spec_heat_cp)**inv_kappa_mns1/(fpt1*th_l2_n(k))
      rho_l2(i,j,k) = rho_l2_pt   ! Saved for later use
      rho_l2_k(k) = rho_l2_pt     ! Saved for iterative loop
      T_temp_l2_iter(k) = th_l2_iter(k)**gamma*                            &
                     (rho_l2_pt*gas_const_R*fpt1*inv_p0_ref)**gamma_mns1
      esat_l2_iter(k) = 611.2_dbl_kind*exp(17.67_dbl_kind*                 &
                (T_temp_l2_iter(k)-273.15_dbl_kind)/(T_temp_l2_iter(k)-    &
                    29.65_dbl_kind))
      qsat_iter(k) = esat_l2_iter(k) /                                     &
                        (rho_l2_pt*gas_const_Rv*T_temp_l2_iter(k))
      qc_l2_iter(k) = max(c0,qvpc_l2_np1(k)-qsat_iter(k))
      do k = 1, nlm-1
         fpt1 = c1 + invepsilon*(max(c0,qvpc_l2_np1(k))-qc_l2_n(k))
                                  ! moisture correction
                                  ! for state eqn. (qc_l2 based on n4 step)
         rho_l2_pt = m_l2_np1(k)*grav*d_eta_l2(k)*c2/                      &
                                      (phi_l2_np1(k+1)-phi_l2_np1(k-1))
         rho_l2(i,j,k) = rho_l2_pt   ! Saved for later use
         rho_l2_k(k) = rho_l2_pt     ! Saved for iterative loop
         T_temp_l2_iter(k) = th_l2_iter(k)**gamma*                         &
                        (rho_l2_pt*gas_const_R*fpt1*inv_p0_ref)**gamma_mns1
         esat_l2_iter(k) = 611.2_dbl_kind*exp(17.67_dbl_kind*              &
                   (T_temp_l2_iter(k)-273.15_dbl_kind)/(T_temp_l2_iter(k)- &
                       29.65_dbl_kind))
         qsat_iter(k) = esat_l2_iter(k) /                                  &
                           (rho_l2_pt*gas_const_Rv*T_temp_l2_iter(k))
         qc_l2_iter(k) = max(c0,qvpc_l2_np1(k)-qsat_iter(k))
      end do
      k = nlm
      ! Estimate P_exnr_l2 at the model top from hydrostatics
      fpt1 = c1 + invepsilon*(max(c0,qvpc_l2_np1(k))-qc_l2_n(k))
                                ! moisture correction
                                ! for state eqn. (qc_l2 based on n4 step)
      fpt2 = fpt1/(c1+max(c0,qT_l2(i,j,k)))   ! moisture correction for vpgf
      P_exnr_l2 = P_exnr(i,j,k) - (grav*z_top-phi(i,j,k))/                 &
                                       (fpt2*th_l2_n(k))
      rho_l2_pt = p0_ref*inv_gas_const_R*(P_exnr_l2*                       &
                     inv_spec_heat_cp)**inv_kappa_mns1/(fpt1*th_l2_n(k))
      rho_l2(i,j,k) = rho_l2_pt   ! Saved for later use
      rho_l2_k(k) = rho_l2_pt     ! Saved for iterative loop
      T_temp_l2_iter(k) = th_l2_iter(k)**gamma*                            &
                     (rho_l2_pt*gas_const_R*fpt1*inv_p0_ref)**gamma_mns1
      esat_l2_iter(k) = 611.2_dbl_kind*exp(17.67_dbl_kind*                 &
                (T_temp_l2_iter(k)-273.15_dbl_kind)/(T_temp_l2_iter(k)-    &
                    29.65_dbl_kind))
      qsat_iter(k) = esat_l2_iter(k) /                                     &
                        (rho_l2_pt*gas_const_Rv*T_temp_l2_iter(k))
      ! NOTE:  This will be the final value of qc_l2_iter(nlm)
      qc_l2_iter(k) = max(c0,qvpc_l2_np1(k)-qsat_iter(k))
      ! Variables used in Kessler microphysics scheme
      pl_l2(i,j,k) = rho_l2(i,j,k)*gas_const_R*T_temp_l2_iter(k)*fpt1
      esat(i,j,k) = esat_l2_iter(k)
      qsat(i,j,k) = qsat_iter(k)
      T_temp_l2(i,j,k) = T_temp_l2_iter(k)


      ! Initialize f1_fac
      do k = 0,  nlm-1
         f1_fac1(k) = m_PI_l2_n(k)*inv_L_vap*inv_d_eta_l2(k)
         f1_fac2(k) = (m_l2_n(k)*qc_l2_n(k) - f1_fac1(k)*th_l2_int(k))/  &
                         m_l2_np1(k)
         f1_fac3(k) = f1_fac1(k)/m_l2_np1(k)
      end do



      ! Note:  This algorithm is based on Newton's method.  In each
      ! iteration, the matrix equation Ay=F is solved for y, which
      ! is the incremental change in the dependent variables (th_l2).
      ! A is the jacobian, and F is the vector-
      ! valued function of nonlinear residuals which is zero when the system
      ! is solved.
      ! NOTE:  In this simplified "quasi-Lagrangian" method, A is a diagonal
      !        matrix.

      !
      ! Begin iteration
      !
      iter = 0
      max_change_theta = 1.E+10_dbl_kind
      max_change_T = 1.E+10_dbl_kind
      do while ( (max_change_theta.ge.converge_cond_theta_T) .and.   &
                 (max_change_T.ge.converge_cond_theta_T) )
        

         ! Initialize max_change
         max_change_theta = c0
         max_change_T = c0


         ! Initialize F vector, i.e., vector-valued function of
         ! nonlinear residuals
         
         ! F1(theta)
         do k = 0, nlm-1
            F1vec(k) = f1_fac2(k) + f1_fac3(k)*th_l2_iter(k) - qc_l2_iter(k)
         end do
         
         
         ! Form Jacobian matrix

         ! D(F1)

         ! form "a" elements
         do k = 0, nlm-1
            a_elem(k) = f1_fac3(k)
            if ( qvpc_l2_np1(k)-qsat_iter(k).gt.c0 ) then
               fpt1 = 4302.645_dbl_kind*T_temp_l2_iter(k)/                 &
                         (T_temp_l2_iter(k)-29.65_dbl_kind)**2-c1
               fpt2 = c1 - gamma_mns1*fpt1/(c1+epsilon/qsat_iter(k))
               dqsatdth = gamma*qsat_iter(k)*fpt1 / (th_l2_iter(k)*fpt2)
               a_elem(k) = a_elem(k) + dqsatdth
            end if
         end do

         ! Solve for the incremental change in the state vector F1
         ! Note:  f_elem is substituted by -F1vec
         ! Note:  Don't increment th_l2_iter yet
         do k = 0, nlm-1
            step_change_th(k) = -F1vec(k) / a_elem(k)
         end do


         ! We're not done yet -- need to iteratively find T_temp_l2_iter
         
         ! F2(theta,T_temp)
         do k = 0, nlm-1
            fpt1 = rho_l2_k(k)*gas_const_R*(c1+invepsilon*                 &
                                  (max(c0,qvpc_l2_np1(k))-qc_l2_iter(k)))
            F2vec = T_temp_l2_iter(k) - th_l2_iter(k)**gamma *             &
                          (fpt1*inv_p0_ref)**gamma_mns1
            dF2_dth = - gamma*T_temp_l2_iter(k)/th_l2_iter(k)
            dF2_dT = c1
            if ( qvpc_l2_np1(k)-qsat_iter(k).gt.c0 ) then
               fpt2 = 4302.645_dbl_kind*esat_l2_iter(k)/                   &
                               (T_temp_l2_iter(k)-29.65_dbl_kind)**2
               dF2_dT = dF2_dT - gamma_mns1*(fpt2-esat_l2_iter(k)/         &
                              T_temp_l2_iter(k))/fpt1
            end if
            step_change_T = - (F2vec + dF2_dth*step_change_th(k)) / dF2_dT
            T_temp_l2_iter(k) = T_temp_l2_iter(k) + step_change_T
            ! Determine max change in delta_T's
            max_change_T = max(max_change_T,abs(step_change_T))
         end do



         do k = 0, nlm-1

            ! Now update th_l2_iter
            th_l2_iter(k) = th_l2_iter(k) + step_change_th(k)

            ! Determine max changes in delta_th's
            max_change_theta = max(max_change_theta,abs(step_change_th(k)))

            ! Update diagnostic variables based on new th_l2_iters
            ! NOTE: With QL-method, rho_l2's do not change
            esat_l2_iter(k) = 611.2_dbl_kind*exp(17.67_dbl_kind*           &
                   (T_temp_l2_iter(k)-273.15_dbl_kind)/(T_temp_l2_iter(k)- &
                       29.65_dbl_kind))
            qsat_iter(k) = esat_l2_iter(k) /                               &
                              (rho_l2_k(k)*gas_const_Rv*T_temp_l2_iter(k))
            qc_l2_iter(k) = max(c0,qvpc_l2_np1(k)-qsat_iter(k))

         end do
         



         iter = iter + 1
         if (iter .ge. max_iter) then
            print *
            print *, "Max. iteration count reached in subroutine &
                     &update diagnostics -- Latent heating iteration."
            print *, "max_change_theta = ", max_change_theta
            print *, "max_change_T = ", max_change_T
            print *, "i =", i, "  j =", j
            print *
            print *, "PROGRAM STOPPING "
            print *
            stop
         end if

      end do   ! while ( (max_change_theta.ge.converge_cond_theta_T) .and.
               !         (max_change_T.ge.converge_cond_theta_T) )




      ! Finalize prognostic and diagnostic variables
      ! Note: Final calculation of m_l2 delayed until end of subroutine
      th_l2(i,j,0,n3) = th_l2_iter(0)
      qc_l2(i,j,0) = qc_l2_iter(0)
      do k = 1, nlm-1
         th_l2(i,j,k,n3) = th_l2_iter(k)
         qc_l2(i,j,k) = qc_l2_iter(k)
      end do
      th_l2(i,j,nlm,n3) = th_l2_iter(nlm)
      qc_l2(i,j,nlm) = qc_l2_iter(nlm)


      ! The following is needed for the Kessler microphysics scheme
      do k = 0, nlm-1
         qsat(i,j,k) = qsat_iter(k)
         esat(i,j,k) = esat_l2_iter(k)
         fpt1 = rho_l2_k(k)*gas_const_R*(c1+invepsilon*              &
                   (max(c0,qvpc_l2_np1(k))-qc_l2_iter(k)))
         pl_l2(i,j,k) = fpt1*T_temp_l2_iter(k)
         T_temp_l2(i,j,k) = T_temp_l2_iter(k)
      end do



      do k = 0, nlm-1
         Q_diabatic(i,j,k) = Q_diabatic(i,j,k) + L_vap * invdt*            &
            (m_l2_np1(k)*qc_l2_iter(k)-m_l2_n(k)*qc_l2_n(k)) / m_l2_n(k)
      end do


      ! Set in_cloud, i.e., if there is diabatic (latent) heating
      in_cloud(i,j,:) = .false.
      do k = 1, nlm-1
         if ( abs(Q_diabatic(i,j,k)).gt.1.E-12_dbl_kind )            &
            in_cloud(i,j,k) = .true.
      end do



      ! Temporary lines
      max_number_iter = max(max_number_iter, iter)


   end do  ! i = 1, im
end do  ! j = 1, jm



! Temporary line
print "(A52,I5)", "In this time step, moisture solver &
                          &max_number_iter =", max_number_iter






! Calculate new (n+1) values of m_l2
m_l2(:,:,0) = m(:,:,1,n3)
do k = 1, nlm-1
   m_l2(:,:,k) = p5 * ( m(:,:,k,n3)*d_eta(k) +                       &
                        m(:,:,k+1,n3)*d_eta(k+1) ) * inv_d_eta_l2(k)
end do
m_l2(:,:,nlm) = m(:,:,nlm,n3)

! Calculate new (n+1) values of qvpc_l2 and qr_l2
qvpc_l2(:,:,:) = qvpc_l2_m_l2(:,:,:,n3) / m_l2(:,:,:)
qr_l2(:,:,:) = qr_l2_m_l2(:,:,:,n3) / m_l2(:,:,:)


!
! Perform hole filling of rain field
!
neg_qr = c0
total_qr = c0
do k = 0, nlm
   do j = 1, jm
      do i = 1, im
         total_qr = total_qr + m_l2(i,j,k)*qr_l2(i,j,k)*d_eta_l2(k)
         if ( qr_l2(i,j,k) .lt. c0 ) then
            neg_qr = neg_qr + m_l2(i,j,k)*qr_l2(i,j,k)*d_eta_l2(k)
            ! Replace negative value of qr_l2 by zero
            qr_l2(i,j,k) = c0
         end if
      end do
   end do
end do
if ( total_qr .eq. c0 ) then
   PHI_qr = c1
   print *, "PHI_qr = ", PHI_qr
else
   PHI_qr = ( total_qr + neg_qr ) / total_qr
   print *, "PHI_qr = ", PHI_qr
   ! Adjust qr_l2
   qr_l2(:,:,:) = PHI_qr * qr_l2(:,:,:)
end if



!
! Kessler microphysics -- Method copied from WRF model
!

do j = 1, jm
   do i = 1, im

      ! Calculate terminal velocities
      ! Vt(:) = 5._dbl_kind
      do k = 0, nlm
         rho_l2_k(k) = rho_l2(i,j,k)
         qr_l2_k(k) = qr_l2(i,j,k)
         qc_l2_k(k) = qc_l2(i,j,k)
         qrr = max(c0,qr_l2_k(k)*0.001_dbl_kind*rho_l2_k(k))
         Vt(i,j,k) = 36.34_dbl_kind*(qrr**0.1364_dbl_kind)*          &
                                  sqrt(rho_l2_k(0)/rho_l2_k(k))
      end do

      ! Sedimentation
      ! First, calculate alternative rho_l2_k0
      rho_l2_k0 = m_l2(i,j,0)*grav*d_eta_l2(0)/(phi(i,j,1)-phis(i,j))
      sedim_fac(0) = 2*grav*dt / ( rho_l2_k0*                        &
                        (phi_l2(i,j,1,n3)-phis(i,j)) )
      sedim_fac(1) = 2*grav*dt / ( rho_l2_k(1)*                      &
                        (phi_l2(i,j,2,n3)-phis(i,j)) )
      do k = 2, nlm-2
         sedim_fac(k) = 2*grav*dt / ( rho_l2_k(k)*                   &
                           (phi_l2(i,j,k+1,n3)-phi_l2(i,j,k-1,n3)) )
      end do
      sedim_fac(nlm-1) = 2*grav*dt / ( rho_l2_k(nlm-1)*              &
                            (grav*z_top-phi_l2(i,j,nlm-2,n3)) )

      ! Use upstream fluxes
      fpt1 = p5*( rho_l2_k(1)*Vt(i,j,1)*qr_l2_k(1) +                  &
                  rho_l2_k0*Vt(i,j,0)*qr_l2_k(0) )    ! rain flux out of
                                                      ! bottom of domain
      dqr_sedim(0) = sedim_fac(0) * ( rho_l2_k(1)*Vt(i,j,1)*         &
             qr_l2_k(1) - fpt1 )
      do k = 1, nlm-1
         dqr_sedim(k) = sedim_fac(k) * ( rho_l2_k(k+1)*Vt(i,j,k+1)*  &
                   qr_l2_k(k+1) - rho_l2_k(k)*Vt(i,j,k)*qr_l2_k(k) )
      end do
      dqr_sedim(nlm) = c0

      ! Calculate surface precipitation
      ppt = 1000*fpt1/rhowater ! surf. rain rate (mm/s)
      RAIN_RATE(i,j) = ppt*86400   ! surface rain rate (mm/day)
      RAIN_SUM(i,j) = RAIN_SUM(i,j) + ppt*dt  ! Accumulated rain (mm)


      ! Calculate rain production from collection and autoconversion
      ! Do this using pre-sedimentation qr_l2's
      ! Then calculate evaporation
      ! NOTE:  For now, no microphysics performed at top layer edge
      do k = 0, nlm-1

         fpt1 = c1 / ( c1 +                                          &
                2.2_dbl_kind*dt*max(c0,qr_l2_k(k))**0.875_dbl_kind )
         qrprod_collect = qc_l2_k(k)*(c1-fpt1)
         qrprod_auto_conv = fpt1*0.001_dbl_kind*dt*                  &
                               max(c0,qc_l2_k(k)-qc_crit)
         ! Modify cloud budget by collection and autoconversion
         qc_l2_k(k) = max(c0,                                        &
                      qc_l2_k(k) - qrprod_collect - qrprod_auto_conv)

         ! Now modify rain-water mixing ratio due to sedimentation
         qr_l2_k(k) = qr_l2_k(k) + dqr_sedim(k)
         ! Modify rain budget by collection and autoconversion
         qr_l2_k(k) = max(c0,                                        &
                      qr_l2_k(k) + qrprod_collect + qrprod_auto_conv)

         ! Calculate rain-drop evaporation
         ! NOTE:  In this model, the WRF variable "product" is not
         !        calculated because it would automatically be zero
         !        since at this point qv=qsat if there is a cloud
         qv_l2 = max(c0,qvpc_l2(i,j,k) - qc_l2(i,j,k))
         fpt1 = 4302.645_dbl_kind*L_vap*inv_spec_heat_cp*            &
                   qsat(i,j,k)/(T_temp_l2(i,j,k)-29.65_dbl_kind)**2
         prod = (qv_l2-qsat(i,j,k)) / (c1+pl_l2(i,j,k)/              &
                  (pl_l2(i,j,k)-esat(i,j,k))*fpt1)
         rcgs = 0.001_dbl_kind*rho_l2_k(k)
         fpt1 = ( 1.6_dbl_kind + 124.9_dbl_kind*                     &
                   (rcgs*qr_l2_k(k))**0.2046_dbl_kind ) *            &
                   (rcgs*qr_l2_k(k))**0.525_dbl_kind
         fpt2 = 2.55E+08_dbl_kind/(pl_l2(i,j,k)*qsat(i,j,k)) +       &
                   5.4E+05_dbl_kind
         fpt3 = max(c0,qsat(i,j,k)-qv_l2)/(rcgs*qsat(i,j,k))
         evap = min( dt*(fpt1/fpt2)*fpt3, max(c0,-prod-qc_l2_k(k)),  &
                      qr_l2_k(k) )

         ! Update all variables
         Q_diab_rain_evap(i,j,k) = - L_vap*evap*invdt
         fpt1 = spec_heat_cp*(pl_l2(i,j,k)*inv_p0_ref)**kappa ! Exner funct.
         ! Latent cooling effect of evaporation on potential temperature
         th_l2(i,j,k,n3) = th_l2(i,j,k,n3) + dt *                    &
                                   Q_diab_rain_evap(i,j,k)/fpt1
         Q_diabatic(i,j,k) = Q_diabatic(i,j,k) + Q_diab_rain_evap(i,j,k)

         ! When updating qvpc_l2, need to account for decrease in qc_l2
         ! by collection and auto conversion of rain drops calculated
         ! above.  Note:  Change in qc_l2 is qc_l2_k(k) - qc_l2(i,j,k)
         qvpc_l2(i,j,k) = qvpc_l2(i,j,k) + evap + qc_l2_k(k) - qc_l2(i,j,k)
         qr_l2_k(k) = qr_l2_k(k) - evap

         qc_l2(i,j,k) = qc_l2_k(k)

         qr_l2(i,j,k) = qr_l2_k(k)
         qr_l2_m_l2(i,j,k,n3) = qr_l2_k(k) * m_l2(i,j,k)

         qvpc_l2_m_l2(i,j,k,n3) = qvpc_l2(i,j,k) * m_l2(i,j,k)

      end do


   end do
end do




end subroutine update_diagnostics

!======================================================================
! END OF UPDATE_DIAGNOSTICS
!======================================================================




!======================================================================
! BEGINNING OF STEP_DYNAMICS
!======================================================================

subroutine step_dynamics ( step_count, w1, w2, w3 )

!----------------------------------------------------------------------
! PURPOSE:
!   Performs the dynamics time stepping using the Adams-Bashworth 3rd.
!   order scheme.
!----------------------------------------------------------------------

implicit none

!-------------------------------------------------------------------
! INTENT IN 
!-------------------------------------------------------------------
integer (kind = int_kind), intent(in) :: step_count

real (kind = dbl_kind), intent(in) :: w1, w2, w3  ! Time stepping
                                                  ! weights




! Advance prognostic time step indices
n4 = mod(step_count+1,2) + 1
n3 = mod(step_count,2) + 1


! Step prognostic variables
u(:,:,:,n4) = u(:,:,:,n3) + dt * ( w3*u_f(:,:,:,n3_f) +              &
                                   w2*u_f(:,:,:,n2_f) +              &
                                   w1*u_f(:,:,:,n1_f) )

v(:,:,:,n4) = v(:,:,:,n3) + dt * ( w3*v_f(:,:,:,n3_f) +              &
                                   w2*v_f(:,:,:,n2_f) +              &
                                   w1*v_f(:,:,:,n1_f) )

! m has already been time-stepped in subroutine update_diagnostics


! w_l2 has already been time-stepped in subroutine update_diagnostics


! qvpc_l2_m_l2 has already been time-stepped in subroutine update_diagnostics


! Note:  th_l2(:,:,:,n4) for k = 1 to nlm-1 was already
! "time-stepped" in subroutine get_eta_dot_l2 and at model
! top and bottom in subroutine update_diagnostics


! Note:  phi_l2(:,:,:,n4) was already "time-stepped" in 
! subroutine get_eta_dot_l2



! Advance tendency time step indices
n3_f = mod(step_count+2,3) + 1
n2_f = mod(step_count+1,3) + 1
n1_f = mod(step_count,3) + 1


end subroutine step_dynamics

!======================================================================
! END OF STEP_DYNAMICS
!======================================================================




end module step

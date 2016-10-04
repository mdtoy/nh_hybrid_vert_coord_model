module momentum_tendency

!-----------------------------------------------------------------------
! PURPOSE: Calculates the tendencies of u and v.
!-----------------------------------------------------------------------

use kinds
use model_parameters
use physical_parameters
use eta_coordinate
use physics
use semi_implicit_solvers


implicit none
save


! Declare variables related to 2D divergence damping
! Note:  For now damping is only calculated on x-z plane
real (kind = dbl_kind), dimension(im,jm,nlm) ::                      &
          divrgnc    ! velocity divergence (du/dx+dw/dz) (s^-1)

real (kind = dbl_kind), parameter ::                                 &
           alpha_d = 2700._dbl_kind  ! divergence damping (bldr_wnd_strm)
!            alpha_d = 540._dbl_kind  ! divergence damping
                                    ! coefficient (m^2/s)



contains


!======================================================================
! BEGINNING OF GET_WL2F_1
!======================================================================

subroutine get_wl2f_1 ( u, v, m_eta_dot, w_l2, w_l2_0, m_l2, m,      &
                      th_l2, P_exnr, phi, vpgf_l2, phi_l2,           &
                      phis, w_advect_0, w_l2_f, w_l2_f_exp_trap )


!---------------------------------------------------------------------------
! PURPOSE:
!   Computes the time tendency of the vertical component of velocity
!   and also the advection of vertical velocity associated with the 
!   top and bottom half-layers for use in the diagnosis of the Exner
!   function (and therefore pressure) at the top and bottom boundaries.
!---------------------------------------------------------------------------

implicit none

!---------------------------------------------------------------------------
! INTENT IN
!---------------------------------------------------------------------------
real (kind=dbl_kind), dimension(im,jm,nlm), intent(in) ::            &
          u, v,  &       ! horiz. velocity components (m/s)
          m_eta_dot      ! vertical mass flux at layer centers (kg/m^2/s)

real (kind=dbl_kind), dimension(im,jm,1:nlm-1), intent(in) ::        &
          w_l2           ! vertical velocity (m/s)

real (kind=dbl_kind), dimension(im,jm), intent(in) ::                &
          w_l2_0         ! vertical velocity at surface (m/s)

real (kind=dbl_kind), dimension(im,jm,0:nlm), intent(in) ::          &
          m_l2           ! pseudo-density at layer edges (kg/m^2/eta)

real (kind=dbl_kind), dimension(im,jm,nlm), intent(in) ::            &
          m              ! pseudo-density at layer centers (kg/m^2/eta)

real (kind=dbl_kind), dimension(im,jm,0:nlm), intent(in) ::          &
          th_l2          ! potential temperature at layer edges (K)

real (kind=dbl_kind), dimension(im,jm,nlm), intent(in) ::            &
          P_exnr,   &    ! Exner function at layer centers (J/kg/K)
          phi            ! layer center geopotential (J/kg)

real (kind = dbl_kind), dimension(im,jm,0:nlm-1), intent(in) ::      &
          vpgf_l2        ! Vertical pressure gradient force (m/s^2)

real (kind=dbl_kind), dimension(im,jm,1:nlm-1), intent(in) ::        &
          phi_l2         ! geopotential at layer edges (J/kg)

real (kind=dbl_kind), dimension(im,jm), intent(in) ::                &
          phis           ! surface geopotential (J/kg)

!---------------------------------------------------------------------------
! INTENT OUT
!---------------------------------------------------------------------------
real (kind=dbl_kind), dimension(im,jm), intent(out) ::               &
          w_advect_0       ! advection of vertical velocity
                           ! associated with the bottom half-layer

real (kind=dbl_kind), dimension(im,jm,1:nlm-1), intent(out) ::       &
          w_l2_f           ! time tendency of w_l2 (m/s^2)

real (kind=dbl_kind), dimension(im,jm,1:nlm-1), intent(out) ::       &
          w_l2_f_exp_trap  ! time tendency of w_l2 (m/s^2) due to
                           ! explicit part of trapezoidal scheme

!---------------------------------------------------------------------------
! LOCAL
!---------------------------------------------------------------------------
integer :: k

real (kind=dbl_kind), dimension(im,jm) ::                            &
       mu_l2,    &     ! x-direction mass flux
                       ! at half-integer levels  (kg/m/eta/s)
       mv_l2           ! y-direction mass flux
                       ! at half-integer levels  (kg/m/eta/s)

real (kind=dbl_kind), parameter ::                                   &
       kappa_w = 6.85E+08_dbl_kind   ! del-4 diffusion coefficient for
                                     ! w_l2 (m^4/s)

real (kind=dbl_kind), dimension(im,jm) :: f1, f2, f3   ! working variables



! Start out by calculating velocity divergence
! Note:  This subroutine (get_wldf) MUST be called before get_uf_vf because
!        velocity divergence is needed in that both subroutines
! This is the corrected divergence which takes into account sloping
! coordinate surfaces (April 5, 2010)
! Corrected on February 8, 2011

f1(:,:) = c1/(phi(:,:,2)-phi(:,:,1))  ! inverse delta-phi for 
                                      ! du/dz and dv/dz calculations

f2(:,:) = invdx*(u(ip1(:),:,1)-u(:,:,1)) -                               &
     invdx*p25*(phi(ip1(:),:,1)-phi(im1(:),:,1))*                        &
           (u(ip1(:),:,2)+u(:,:,2)-u(ip1(:),:,1)-u(:,:,1))*f1(:,:) +     &
       grav*(w_l2(:,:,1)-w_l2_0(:,:))/(phi_l2(:,:,1)-phis(:,:))

divrgnc(:,:,1) = f2(:,:) + invdy*(v(:,jp1(:),1)-v(:,:,1)) -              &
  invdy*p25*(phi(:,jp1(:),1)-phi(:,jm1(:),1))*(v(:,jp1(:),2)+v(:,:,2)-   &
  v(:,jp1(:),1)-v(:,:,1))*f1(:,:)

do k = 2, nlm-1

   f1(:,:) = c1/(phi(:,:,k+1)-phi(:,:,k-1))  ! inverse delta-phi for 
                                             ! du/dz and dv/dz calculations
   f2(:,:) = invdx*(u(ip1(:),:,k)-u(:,:,k)) -                            &
      invdx*p25*(phi(ip1(:),:,k)-phi(im1(:),:,k))*                       &
     (u(ip1(:),:,k+1)+u(:,:,k+1)-u(ip1(:),:,k-1)-u(:,:,k-1))*f1(:,:) +   &
        grav*(w_l2(:,:,k)-w_l2(:,:,k-1))/(phi_l2(:,:,k)-phi_l2(:,:,k-1))

   divrgnc(:,:,k) = f2(:,:) + invdy*(v(:,jp1(:),k)-v(:,:,k)) -           &
     invdy*p25*(phi(:,jp1(:),k)-phi(:,jm1(:),k))*(v(:,jp1(:),k+1)+       &
     v(:,:,k+1)-v(:,jp1(:),k-1)-v(:,:,k-1))*f1(:,:)

end do

f1(:,:) = c1/(phi(:,:,nlm)-phi(:,:,nlm-1))  ! inverse delta-phi for 
                                            ! du/dz and dv/dz calculations

f2(:,:) = invdx*(u(ip1(:),:,nlm)-u(:,:,nlm)) -                           &
   invdx*p25*(phi(ip1(:),:,nlm)-phi(im1(:),:,nlm))*                      &
  (u(ip1(:),:,nlm)+u(:,:,nlm)-u(ip1(:),:,nlm-1)-u(:,:,nlm-1))*f1(:,:) -  &
     grav*w_l2(:,:,nlm-1)/(grav*z_top-phi_l2(:,:,nlm-1))

divrgnc(:,:,nlm) = f2(:,:) + invdy*(v(:,jp1(:),nlm)-v(:,:,nlm)) -     &
  invdy*p25*(phi(:,jp1(:),nlm)-phi(:,jm1(:),nlm))*(v(:,jp1(:),nlm)+   &
  v(:,:,nlm)-v(:,jp1(:),nlm-1)-v(:,:,nlm-1))*f1(:,:)




! Start at layer k = 3/2

! Calculate advection term

mu_l2(:,:) = p5 * ( d_eta(1)*p5*(m(:,:,1)+m(im1(:),:,1))*u(:,:,1) +      &
               d_eta(2)*p5*(m(:,:,2)+m(im1(:),:,2))*u(:,:,2) ) /         &
                                     d_eta_l2(1)

mv_l2(:,:) = p5 * ( d_eta(1)*p5*(m(:,:,1)+m(:,jm1(:),1))*v(:,:,1) +      &
               d_eta(2)*p5*(m(:,:,2)+m(:,jm1(:),2))*v(:,:,2) ) /         &
                                     d_eta_l2(1)

f1(:,:) = (c1/dx) *                                                      &
        (  mu_l2(:,:) * p5 * ( w_l2(:,:,1) - w_l2(im1(:),:,1) ) +        &
           mu_l2(ip1(:),:) * p5 * ( w_l2(ip1(:),:,1) - w_l2(:,:,1) )  )

f2(:,:) = (c1/dy) *                                                      &
        (  mv_l2(:,:) * p5 * ( w_l2(:,:,1) - w_l2(:,jm1(:),1) ) +        &
           mv_l2(:,jp1(:)) * p5 * ( w_l2(:,jp1(:),1) - w_l2(:,:,1) )  )

f3(:,:) = (c1/d_eta_l2(1)) *                                             &
        (  m_eta_dot(:,:,1) * p5 * ( w_l2(:,:,1) - w_l2_0(:,:) ) +       &
           m_eta_dot(:,:,2) * p5 * ( w_l2(:,:,2) - w_l2(:,:,1) )  )

w_l2_f(:,:,1) = - (c1/m_l2(:,:,1)) * ( f1(:,:) + f2(:,:) + f3(:,:) )


! Add contributions from gravity,
! subgrid-scale turbulent momentum flux and divergence damping

w_l2_f(:,:,1) = w_l2_f(:,:,1) - grav + F_turb_w_l2(:,:,1) +          &
   alpha_d*grav*(divrgnc(:,:,2)-divrgnc(:,:,1))/                     &
      (phi(:,:,2)-phi(:,:,1))


! Add del4 diffusion term to w_l2_f
! First, calculate Laplacian of w_l2
f1(:,:) = invdx2 * (w_l2(ip1(:),:,1)-2*w_l2(:,:,1)+                  &
                                           w_l2(im1(:),:,1)) +       &
          invdy2 * (w_l2(:,jp1(:),1)-2*w_l2(:,:,1)+w_l2(:,jm1(:),1))


! Now, calculate Laplacian of Laplacian of w_l2, multiply by
! diffusion coefficient, and subtract from w_l2_f
w_l2_f(:,:,1) = w_l2_f(:,:,1)   -   kappa_w *    (                   &
       invdx2 * (f1(ip1(:),:)-2*f1(:,:)+f1(im1(:),:)) +              &
       invdy2 * (f1(:,jp1(:))-2*f1(:,:)+f1(:,jm1(:)))    )




! Move on up

do k = 2, nlm-2

   ! Calculate advection term

   mu_l2(:,:) = p5 * ( d_eta(k)*p5*(m(:,:,k)+m(im1(:),:,k))*u(:,:,k) +      &
                  d_eta(k+1)*p5*(m(:,:,k+1)+m(im1(:),:,k+1))*u(:,:,k+1) ) / &
                                        d_eta_l2(k)

   mv_l2(:,:) = p5 * ( d_eta(k)*p5*(m(:,:,k)+m(:,jm1(:),k))*v(:,:,k) +      &
                  d_eta(k+1)*p5*(m(:,:,k+1)+m(:,jm1(:),k+1))*v(:,:,k+1) ) / &
                                        d_eta_l2(k)

   f1(:,:) = (c1/dx) *                                                      &
           (  mu_l2(:,:) * p5 * ( w_l2(:,:,k) - w_l2(im1(:),:,k) ) +        &
              mu_l2(ip1(:),:) * p5 * ( w_l2(ip1(:),:,k) - w_l2(:,:,k) )  )

   f2(:,:) = (c1/dy) *                                                      &
           (  mv_l2(:,:) * p5 * ( w_l2(:,:,k) - w_l2(:,jm1(:),k) ) +        &
              mv_l2(:,jp1(:)) * p5 * ( w_l2(:,jp1(:),k) - w_l2(:,:,k) )  )

   f3(:,:) = (c1/d_eta_l2(k)) *                                             &
           (  m_eta_dot(:,:,k) * p5 * ( w_l2(:,:,k) - w_l2(:,:,k-1) ) +     &
              m_eta_dot(:,:,k+1) * p5 * ( w_l2(:,:,k+1) - w_l2(:,:,k) )  )

   w_l2_f(:,:,k) = - (c1/m_l2(:,:,k)) * ( f1(:,:) + f2(:,:) + f3(:,:) )


   ! Add contributions from gravity,
   ! subgrid-scale turbulent momentum flux and divergence damping

   w_l2_f(:,:,k) = w_l2_f(:,:,k) - grav + F_turb_w_l2(:,:,k) +       &
      alpha_d*grav*(divrgnc(:,:,k+1)-divrgnc(:,:,k))/                &
         (phi(:,:,k+1)-phi(:,:,k))


   ! Add del4 diffusion term to w_l2_f
   ! First, calculate Laplacian of w_l2
   f1(:,:) = invdx2 * (w_l2(ip1(:),:,k)-2*w_l2(:,:,k)+                  &
                                              w_l2(im1(:),:,k)) +       &
             invdy2 * (w_l2(:,jp1(:),k)-2*w_l2(:,:,k)+w_l2(:,jm1(:),k))


   ! Now, calculate Laplacian of Laplacian of w_l2, multiply by
   ! diffusion coefficient, and subtract from w_l2_f
   w_l2_f(:,:,k) = w_l2_f(:,:,k)   -   kappa_w *    (                &
          invdx2 * (f1(ip1(:),:)-2*f1(:,:)+f1(im1(:),:)) +           &
          invdy2 * (f1(:,jp1(:))-2*f1(:,:)+f1(:,jm1(:)))    )

end do


! Layer k = nlm - 1/2

! Calculate advection term

mu_l2(:,:) = p5 * ( d_eta(nlm-1)*p5*(m(:,:,nlm-1)+m(im1(:),:,nlm-1))*       &
                                                         u(:,:,nlm-1) +     &
               d_eta(nlm)*p5*(m(:,:,nlm)+m(im1(:),:,nlm))*u(:,:,nlm) ) /    &
                                     d_eta_l2(nlm-1)

mv_l2(:,:) = p5 * ( d_eta(nlm-1)*p5*(m(:,:,nlm-1)+m(:,jm1(:),nlm-1))*       &
                                                        v(:,:,nlm-1) +      &
               d_eta(nlm)*p5*(m(:,:,nlm)+m(:,jm1(:),nlm))*v(:,:,nlm) ) /    &
                                     d_eta_l2(nlm-1)

f1(:,:) = (c1/dx) *                                                         &
    (  mu_l2(:,:) * p5 * ( w_l2(:,:,nlm-1) - w_l2(im1(:),:,nlm-1) ) +       &
       mu_l2(ip1(:),:) * p5 * ( w_l2(ip1(:),:,nlm-1) - w_l2(:,:,nlm-1) )  )

f2(:,:) = (c1/dy) *                                                         &
    (  mv_l2(:,:) * p5 * ( w_l2(:,:,nlm-1) - w_l2(:,jm1(:),nlm-1) ) +       &
       mv_l2(:,jp1(:)) * p5 * ( w_l2(:,jp1(:),nlm-1) - w_l2(:,:,nlm-1) )  )

f3(:,:) = (c1/d_eta_l2(nlm-1)) *                                            &
    (  m_eta_dot(:,:,nlm-1) * p5 * ( w_l2(:,:,nlm-1) - w_l2(:,:,nlm-2) ) +  &
       m_eta_dot(:,:,nlm) * p5 * ( - w_l2(:,:,nlm-1) )  )

w_l2_f(:,:,nlm-1) = - (c1/m_l2(:,:,nlm-1)) * ( f1(:,:) + f2(:,:) + f3(:,:) )


! Add contributions from gravity,
! subgrid-scale turbulent momentum flux and divergence damping

w_l2_f(:,:,nlm-1) = w_l2_f(:,:,nlm-1) - grav + F_turb_w_l2(:,:,nlm-1) +     &
   alpha_d*grav*                                                            &
      (divrgnc(:,:,nlm)-divrgnc(:,:,nlm-1))/(phi(:,:,nlm)-phi(:,:,nlm-1))


! Add del4 diffusion term to w_l2_f
! First, calculate Laplacian of w_l2
f1(:,:) = invdx2 * (w_l2(ip1(:),:,nlm-1)-2*w_l2(:,:,nlm-1)+          &
                                         w_l2(im1(:),:,nlm-1)) +     &
          invdy2 * (w_l2(:,jp1(:),nlm-1)-2*w_l2(:,:,nlm-1)+          &
                                            w_l2(:,jm1(:),nlm-1))


! Now, calculate Laplacian of Laplacian of w_l2, multiply by
! diffusion coefficient, and subtract from w_l2_f
w_l2_f(:,:,nlm-1) = w_l2_f(:,:,nlm-1)   -   kappa_w *    (           &
       invdx2 * (f1(ip1(:),:)-2*f1(:,:)+f1(im1(:),:)) +              &
       invdy2 * (f1(:,jp1(:))-2*f1(:,:)+f1(:,jm1(:)))            )




! Calculate w-tendency due to explicit part of trapezoidal scheme for
! vertical pressure gradient force
do k = 1, nlm-1
   w_l2_f_exp_trap(:,:,k) = - alpha_si*vpgf_l2(:,:,k)
end do



! Calculate advection terms at bottom layer to be used to surface
! Exner function

mu_l2(:,:) = p5*(m(:,:,1)+m(im1(:),:,1))*u(:,:,1)

mv_l2(:,:) = p5*(m(:,:,1)+m(:,jm1(:),1))*v(:,:,1)

f1(:,:) = (c1/dx) *                                                      &
        (  mu_l2(:,:) * p5 * ( w_l2_0(:,:) - w_l2_0(im1(:),:) ) +        &
           mu_l2(ip1(:),:) * p5 * ( w_l2_0(ip1(:),:) - w_l2_0(:,:) )  )

f2(:,:) = (c1/dy) *                                                      &
        (  mv_l2(:,:) * p5 * ( w_l2_0(:,:) - w_l2_0(:,jm1(:)) ) +        &
           mv_l2(:,jp1(:)) * p5 * ( w_l2_0(:,jp1(:)) - w_l2_0(:,:) )  )

f3(:,:) = (c1/d_eta_l2(0)) *                                             &
        (  m_eta_dot(:,:,1) * p5 * ( w_l2(:,:,1) - w_l2_0(:,:) )  )

w_advect_0(:,:) = (c1/m_l2(:,:,0)) * ( f1(:,:) + f2(:,:) + f3(:,:) )



end subroutine get_wl2f_1

!======================================================================
! END OF GET_WL2F_1
!======================================================================




!======================================================================
! BEGINNING OF GET_WL2F_2
!======================================================================

subroutine get_wl2f_2 ( m_eta_dot, w_l2, w_l2_0, m_l2, w_advect_0,   &
                        w_l2_f, w_l2_f_temp )


!---------------------------------------------------------------------------
! PURPOSE:
!   Updates the time tendency of the vertical component of velocity
!   and also the advection of vertical velocity associated with the 
!   top and bottom half-layers for use in the diagnosis of the Exner
!   function (and therefore pressure) at the top and bottom boundaries.
!---------------------------------------------------------------------------

implicit none

!---------------------------------------------------------------------------
! INTENT IN
!---------------------------------------------------------------------------
real (kind=dbl_kind), dimension(im,jm,nlm), intent(in) ::            &
          m_eta_dot      ! vertical mass flux at layer centers (kg/m^2/s)

real (kind=dbl_kind), dimension(im,jm,1:nlm-1), intent(in) ::        &
          w_l2           ! vertical velocity (m/s)

real (kind=dbl_kind), dimension(im,jm), intent(in) ::                &
          w_l2_0         ! vertical velocity at surface (m/s)

real (kind=dbl_kind), dimension(im,jm,0:nlm), intent(in) ::          &
          m_l2           ! pseudo-density at layer edges (kg/m^2/eta)

!---------------------------------------------------------------------------
! INTENT INOUT
!---------------------------------------------------------------------------
real (kind=dbl_kind), dimension(im,jm), intent(inout) ::             &
          w_advect_0       ! advection of vertical velocity
                           ! associated with the bottom half-layer

real (kind=dbl_kind), dimension(im,jm,1:nlm-1), intent(inout) ::     &
          w_l2_f           ! time tendency of w_l2 (m/s^2)

!---------------------------------------------------------------------------
! INTENT OUT
!---------------------------------------------------------------------------
real (kind=dbl_kind), dimension(im,jm,1:nlm-1), intent(out) ::       &
          w_l2_f_temp      ! updated portion of time tendency of w_l2 (m/s^2)

!---------------------------------------------------------------------------
! LOCAL
!---------------------------------------------------------------------------
integer :: k

real (kind=dbl_kind), dimension(im,jm) :: f1    ! working variable



! Start at layer k = 3/2

! Calculate advection term

f1(:,:) = inv_d_eta_l2(1) *                                              &
        (  m_eta_dot(:,:,1) * p5 * ( w_l2(:,:,1) - w_l2_0(:,:) ) +       &
           m_eta_dot(:,:,2) * p5 * ( w_l2(:,:,2) - w_l2(:,:,1) )  )

w_l2_f_temp(:,:,1) = - f1(:,:)/m_l2(:,:,1)
w_l2_f(:,:,1) = w_l2_f(:,:,1) + w_l2_f_temp(:,:,1)


! Move on up

do k = 2, nlm-2

   ! Calculate advection term

   f1(:,:) = inv_d_eta_l2(k) *                                             &
           (  m_eta_dot(:,:,k) * p5 * ( w_l2(:,:,k) - w_l2(:,:,k-1) ) +     &
              m_eta_dot(:,:,k+1) * p5 * ( w_l2(:,:,k+1) - w_l2(:,:,k) )  )

   w_l2_f_temp(:,:,k) = - f1(:,:)/m_l2(:,:,k)
   w_l2_f(:,:,k) = w_l2_f(:,:,k) + w_l2_f_temp(:,:,k)

end do


! Layer k = nlm - 1/2

! Calculate advection term

f1(:,:) = inv_d_eta_l2(nlm-1) *                                            &
    (  m_eta_dot(:,:,nlm-1) * p5 * ( w_l2(:,:,nlm-1) - w_l2(:,:,nlm-2) ) +  &
       m_eta_dot(:,:,nlm) * p5 * ( - w_l2(:,:,nlm-1) )  )

w_l2_f_temp(:,:,nlm-1) = - f1(:,:)/m_l2(:,:,nlm-1)
w_l2_f(:,:,nlm-1) = w_l2_f(:,:,nlm-1) + w_l2_f_temp(:,:,nlm-1)




! Update advection terms at bottom layer to be used to calculate surface
! Exner function

f1(:,:) = inv_d_eta_l2(0) *                                          &
        (  m_eta_dot(:,:,1) * p5 * ( w_l2(:,:,1) - w_l2_0(:,:) )  )

w_advect_0(:,:) = w_advect_0(:,:) + f1(:,:)/m_l2(:,:,0)
               



end subroutine get_wl2f_2

!======================================================================
! END OF GET_WL2F_2
!======================================================================





!======================================================================
! BEGINNING OF GET_UF_VF
!======================================================================

subroutine get_uf_vf ( u, v, eta_dot_l2, m, m_l2, th, P_exnr, phi,   &
                       vpgf_l2, phi_l2, phis, moist_corr1_hpgf,      &
                       pv, ke_horiz, u_advect_coriol_1,              &
                       v_advect_coriol_1, u_f, v_f )

!---------------------------------------------------------------------------
! PURPOSE:
!   Computes the time tendency of the x and y components of velocity,
!   i.e. u and v respectively.  Also returns the advection
!   of horizontal velocity and Coriolis acceleration at the bottom
!   layer for use in subroutine update_diagnostics to compute the
!   Exner function at the surface.
!---------------------------------------------------------------------------

implicit none

!---------------------------------------------------------------------------
! INTENT IN
!---------------------------------------------------------------------------
real (kind=dbl_kind), dimension(im,jm,nlm), intent(in) ::            &
          u, v           ! horiz. velocity components (m/s)

real (kind=dbl_kind), dimension(im,jm,1:nlm-1), intent(in) ::        &
          eta_dot_l2     ! generalized vertical velocity at layer
                         ! edges (eta/s)

real (kind=dbl_kind), dimension(im,jm,nlm), intent(in) ::            &
          m              ! pseudo-density at layer centers (kg/m^2/eta)

real (kind=dbl_kind), dimension(im,jm,0:nlm), intent(in) ::          &
          m_l2           ! pseudo-density at layer edges (kg/m^2/eta)

real (kind=dbl_kind), dimension(im,jm,nlm), intent(in) ::            &
          th,      &     ! potential temperature at layer centers (K)
          P_exnr,  &     ! Exner function at layer centers (J/kg/K)
          phi            ! layer center geopotential (J/kg)

real (kind = dbl_kind), dimension(im,jm,0:nlm-1), intent(in) ::      &
          vpgf_l2        ! Vertical pressure gradient force (m/s^2)

real (kind=dbl_kind), dimension(im,jm,1:nlm-1), intent(in) ::        &
          phi_l2         ! geopotential at layer edges (J/kg)

real (kind=dbl_kind), dimension(im,jm), intent(in) ::                &
          phis      ! surface geopotential (J/kg)

real (kind=dbl_kind), dimension(im,jm,nlm), intent(in) ::            &
          moist_corr1_hpgf  ! moisture correction for pressure gradient
                            ! force at layer centers (-)

!---------------------------------------------------------------------------
! INTENT OUT
!---------------------------------------------------------------------------
real (kind=dbl_kind), dimension(im,jm,nlm), intent(out) ::           &
          pv,      &     ! vert. component of pot. vorticity (m^2 eta/kg/s)
          ke_horiz       ! contribution to kinetic energy
                         ! from horizontal velocity (J/kg)

real (kind=dbl_kind), dimension(im,jm), intent(out) ::               &
          u_advect_coriol_1,   &  ! advection of the x-component of
                                  ! velocity plus Coriolis accel.
                                  ! in the bottom layer
          v_advect_coriol_1       ! advection of the y-component of
                                  ! velocity plus Coriolis accel.
                                  ! in the bottom layer

real (kind=dbl_kind), dimension(im,jm,nlm), intent(out) ::           &
          u_f, v_f       ! tendency of u and v (m/s^2)

!---------------------------------------------------------------------------
! LOCAL
!---------------------------------------------------------------------------
integer :: k

real (kind=dbl_kind), dimension(im,jm) :: f1, f2   ! working variables
real (kind=dbl_kind), dimension(im,jm) :: del2v  ! special working variable

real (kind=dbl_kind), parameter ::                                   &
          inv24 = c1/24.00000_dbl_kind

real (kind=dbl_kind), dimension(im,jm) ::                            &
          zeta,            &                  ! relative vorticity (s^-1)
          alfa, beta, gamm, delt, epsln, fi   ! linear combinations of pv

real (kind=dbl_kind), dimension(im,jm) ::                            &
          m_u,          &         ! mass interpolated to u points
          m_v                     ! mass interpolated to v points

real (kind=dbl_kind), parameter ::                                   &
       kappa_uv = 6.85E+08_dbl_kind     ! del-4 diffusion coefficient for
                                      ! u and v (m^4/s)




!---------------------------------------------------------------------------
! Calculate the (PV k) cross (m times velocity) term.       (term 1 of 5)
! Potential enstrophy and energy conserving scheme of Arakawa and 
! Lamb (1981) is used.
!---------------------------------------------------------------------------

do k = 1, nlm

   ! calculate mass interpolated to u and v points
   m_u(:,:)  =  p5 * (m(:,:,k) + m(im1(:),:,k))
   m_v(:,:)  =  p5 * (m(:,:,k) + m(:,jm1(:),k))


   ! calculate potential vorticity
   zeta(:,:) = (c1/dy) * ( u(:,jm1(:),k) - u(:,:,k) )  +             &
               (c1/dx) * ( v(:,:,k) - v(im1(:),:,k) )
   pv(:,:,k) = ( f_cor + zeta(:,:) ) /                               &
                  ( p25 * (m(:,:,k) + m(im1(:),:,k) +                &
                           m(im1(:),jm1(:),k) + m(:,jm1(:),k)) )



   ! calculate linear combinations of pv ( eqn. (3.34) of AL (1981) )
   alfa(:,:)  = inv24 * ( c2*pv(ip1(:),jp1(:),k) + pv(:,jp1(:),k) +  &
                            c2*pv(:,:,k) + pv(ip1(:),:,k) )
   beta(:,:)  = inv24 * ( pv(:,jp1(:),k) + c2*pv(im1(:),jp1(:),k) +  &
                            pv(im1(:),:,k) + c2*pv(:,:,k) )
   gamm(:,:)  = inv24 * ( c2*pv(:,jp1(:),k) + pv(im1(:),jp1(:),k) +  &
                            c2*pv(im1(:),:,k) + pv(:,:,k) )
   delt(:,:)  = inv24 * ( pv(ip1(:),jp1(:),k) + c2*pv(:,jp1(:),k) +  &
                            pv(:,:,k) + c2*pv(ip1(:),:,k) )
   epsln(:,:) = inv24 * ( pv(ip1(:),jp1(:),k) + pv(:,jp1(:),k) -     &
                            pv(:,:,k) - pv(ip1(:),:,k) )
   fi(:,:)    = inv24 * (-pv(ip1(:),jp1(:),k) + pv(:,jp1(:),k) +     &
                            pv(:,:,k) - pv(ip1(:),:,k) ) 


   !
   ! calculate u and v tendencies -- term 1 of 5 ( see eqns. (3.5) and (3.6)
   !                                 of AL (1981) )

   u_f(:,:,k) = alfa(:,:)  * m_v(:,jp1(:)) * v(:,jp1(:),k)           +   &
                beta(:,:)  * m_v(im1(:),jp1(:)) * v(im1(:),jp1(:),k) +   &
                gamm(:,:)  * m_v(im1(:),:) * v(im1(:),:,k)           +   &
                delt(:,:)  * m_v(:,:) * v(:,:,k)                     -   &
                epsln(:,:) * m_u(ip1(:),:) * u(ip1(:),:,k)           +   &
                epsln(im1(:),:) * m_u(im1(:),:) * u(im1(:),:,k)
   v_f(:,:,k) = -gamm(ip1(:),:) * m_u(ip1(:),:) * u(ip1(:),:,k)      -   &
                 delt(:,:) * m_u(:,:) * u(:,:,k)                     -   &
                 alfa(:,jm1(:)) * m_u(:,jm1(:)) * u(:,jm1(:),k)      -   &
                 beta(ip1(:),jm1(:)) *                                   &
                               m_u(ip1(:),jm1(:)) * u(ip1(:),jm1(:),k) - &
                 fi(:,:) * m_v(:,jp1(:)) * v(:,jp1(:),k)             +   &
                 fi(:,jm1(:)) * m_v(:,jm1(:)) * v(:,jm1(:),k)

end do



!---------------------------------------------------------------------------
! Add contribution of the horiz. gradient of kinetic energy.   (term 2 of 5)
! And also divergence damping!!!
!---------------------------------------------------------------------------


! Calculate contribution to kinetic energy
! from the horizontal velocity  ( see eqn (3.41) of AL (1981) )
! Note:  expect SICK to result from use of this K.E.
ke_horiz(:,:,:) = p5 * (p5*u(:,:,:)**2 + p5*u(ip1(:),:,:)**2) +          &
                  p5 * (p5*v(:,:,:)**2 + p5*v(:,jp1(:),:)**2)

! Add contribution of horiz. gradient of K.E.
! This is the corrected divergence of divergence which takes into account
! sloping coordinate surfaces (April 5, 2010)
! Corrected on February 8, 2011

! k = 1
f1(:,:) = c1/(phi(:,:,2)-phi(:,:,1))  ! inverse delta-phi for 
                                      ! dD/dz calculations
f2(:,:) = invdx*(divrgnc(:,:,1)-divrgnc(im1(:),:,1)) - invdx*            &
   (phi(:,:,1)-phi(im1(:),:,1))*p5*(divrgnc(:,:,2)+                      &
    divrgnc(im1(:),:,2)-divrgnc(:,:,1)-divrgnc(im1(:),:,1))*f1(:,:)
u_f(:,:,1) = u_f(:,:,1) -                                                &
                invdx*(ke_horiz(:,:,1)-ke_horiz(im1(:),:,1)) +           &
                   alpha_d*f2(:,:)
f2(:,:) = invdy*(divrgnc(:,:,1)-divrgnc(:,jm1(:),1)) - invdy*            &
      (phi(:,:,1)-phi(:,jm1(:),1))*p5*(divrgnc(:,:,2)+                   &
       divrgnc(:,jm1(:),2)-divrgnc(:,:,1)-divrgnc(:,jm1(:),1))*f1(:,:)
v_f(:,:,1) = v_f(:,:,1) - invdy *                                        &
      (ke_horiz(:,:,1)-ke_horiz(:,jm1(:),1)) + alpha_d*f2(:,:)

do k = 2, nlm-1

   f1(:,:) = c1/(phi(:,:,k+1)-phi(:,:,k-1))  ! inverse delta-phi for 
                                             ! dD/dz calculations
   f2(:,:) = invdx*(divrgnc(:,:,k)-divrgnc(im1(:),:,k)) - invdx*         &
      (phi(:,:,k)-phi(im1(:),:,k))*p5*(divrgnc(:,:,k+1)+                 &
       divrgnc(im1(:),:,k+1)-divrgnc(:,:,k-1)-divrgnc(im1(:),:,k-1))*    &
          f1(:,:)
   u_f(:,:,k) = u_f(:,:,k) -                                             &
                   invdx*(ke_horiz(:,:,k)-ke_horiz(im1(:),:,k)) +        &
                      alpha_d*f2(:,:)
   f2(:,:) = invdy*(divrgnc(:,:,k)-divrgnc(:,jm1(:),k)) - invdy*         &
      (phi(:,:,k)-phi(:,jm1(:),k))*p5*(divrgnc(:,:,k+1)+                 &
       divrgnc(:,jm1(:),k+1)-divrgnc(:,:,k-1)-divrgnc(:,jm1(:),k-1))*    &
          f1(:,:)
   v_f(:,:,k) = v_f(:,:,k) - invdy *                                     &
      (ke_horiz(:,:,k)-ke_horiz(:,jm1(:),k)) + alpha_d*f2(:,:)

end do

! k = nlm
f1(:,:) = c1/(phi(:,:,nlm)-phi(:,:,nlm-1))  ! inverse delta-phi for 
                                            ! dD/dz calculations
f2(:,:) = invdx*(divrgnc(:,:,nlm)-divrgnc(im1(:),:,nlm)) - invdx*        &
   (phi(:,:,nlm)-phi(im1(:),:,nlm))*p5*(divrgnc(:,:,nlm)+                &
    divrgnc(im1(:),:,nlm)-divrgnc(:,:,nlm-1)-divrgnc(im1(:),:,nlm-1))*   &
       f1(:,:)
u_f(:,:,nlm) = u_f(:,:,nlm) -                                            &
                 invdx*(ke_horiz(:,:,nlm)-ke_horiz(im1(:),:,nlm)) +      &
                    alpha_d*f2(:,:)
f2(:,:) = invdy*(divrgnc(:,:,nlm)-divrgnc(:,jm1(:),nlm)) - invdy*        &
   (phi(:,:,nlm)-phi(:,jm1(:),nlm))*p5*(divrgnc(:,:,nlm)+                &
    divrgnc(:,jm1(:),nlm)-divrgnc(:,:,nlm-1)-divrgnc(:,jm1(:),nlm-1))*   &
       f1(:,:)
v_f(:,:,nlm) = v_f(:,:,nlm) - invdy *                                    &
   (ke_horiz(:,:,nlm)-ke_horiz(:,jm1(:),nlm)) + alpha_d*f2(:,:)



!---------------------------------------------------------------------------
! Add contributions due to vert. advection of horiz. momentum and
! subgrid-scale turbulent momentum flux.            (terms 3 & 4 of 5)
!---------------------------------------------------------------------------


! Bottom level

! calculate mass interpolated to u and v points
m_u(:,:)  =  p5 * (m(:,:,1) + m(im1(:),:,1))
m_v(:,:)  =  p5 * (m(:,:,1) + m(:,jm1(:),1))

u_f(:,:,1) = u_f(:,:,1) - ( p5*(m_l2(:,:,1)*eta_dot_l2(:,:,1) +      &  
                      m_l2(im1(:),:,1)*eta_dot_l2(im1(:),:,1)) *     &
                  p5*(u(:,:,2) - u(:,:,1)) ) / (d_eta(1)*m_u(:,:)) + &
                  F_turb_u(:,:,1)
v_f(:,:,1) = v_f(:,:,1) - ( p5*(m_l2(:,:,1)*eta_dot_l2(:,:,1) +      &
                      m_l2(:,jm1(:),1)*eta_dot_l2(:,jm1(:),1)) *     &
                  p5*(v(:,:,2) - v(:,:,1)) ) / (d_eta(1)*m_v(:,:)) + &
                  F_turb_v(:,:,1)


! Continue up to level nlm-1

do k = 2, nlm-1
   ! calculate mass interpolated to u and v points
   m_u(:,:)  =  p5 * (m(:,:,k) + m(im1(:),:,k))
   m_v(:,:)  =  p5 * (m(:,:,k) + m(:,jm1(:),k))

   u_f(:,:,k) = u_f(:,:,k) - ( p5*(m_l2(:,:,k-1)*eta_dot_l2(:,:,k-1) + &
                      m_l2(im1(:),:,k-1)*eta_dot_l2(im1(:),:,k-1)) *   &
                  p5*(u(:,:,k) - u(:,:,k-1))   +                       &
                  p5*(m_l2(:,:,k)*eta_dot_l2(:,:,k) +                  &
                      m_l2(im1(:),:,k)*eta_dot_l2(im1(:),:,k)) *       &
                  p5*(u(:,:,k+1) - u(:,:,k)) ) / (d_eta(k)*m_u(:,:)) + &
                  F_turb_u(:,:,k)
   v_f(:,:,k) = v_f(:,:,k) - ( p5*(m_l2(:,:,k-1)*eta_dot_l2(:,:,k-1) + &
                      m_l2(:,jm1(:),k-1)*eta_dot_l2(:,jm1(:),k-1)) *   &
                  p5*(v(:,:,k) - v(:,:,k-1))   +                       &
                  p5*(m_l2(:,:,k)*eta_dot_l2(:,:,k) +                  &
                      m_l2(:,jm1(:),k)*eta_dot_l2(:,jm1(:),k)) *       &
                  p5*(v(:,:,k+1) - v(:,:,k)) ) / (d_eta(k)*m_v(:,:)) + &
                  F_turb_v(:,:,k)
end do


! Top level

! calculate mass interpolated to u and v points
m_u(:,:)  =  p5 * (m(:,:,nlm) + m(im1(:),:,nlm))
m_v(:,:)  =  p5 * (m(:,:,nlm) + m(:,jm1(:),nlm))

u_f(:,:,nlm) = u_f(:,:,nlm) - ( p5*(m_l2(:,:,nlm-1)*eta_dot_l2(:,:,nlm-1) + &
                  m_l2(:,:,nlm-1)*eta_dot_l2(im1(:),:,nlm-1)) *             &
              p5*(u(:,:,nlm) - u(:,:,nlm-1)) ) / (d_eta(nlm)*m_u(:,:)) +    &
              F_turb_u(:,:,nlm)
v_f(:,:,nlm) = v_f(:,:,nlm) - ( p5*(m_l2(:,:,nlm-1)*eta_dot_l2(:,:,nlm-1) + &
                  m_l2(:,jm1(:),nlm-1)*eta_dot_l2(:,jm1(:),nlm-1)) *        &
              p5*(v(:,:,nlm) - v(:,:,nlm-1)) ) / (d_eta(nlm)*m_v(:,:)) +    &
              F_turb_v(:,:,nlm)




! Extract the layer 1 advection and Coriolis terms
u_advect_coriol_1(:,:) = - u_f(:,:,1)
v_advect_coriol_1(:,:) = - v_f(:,:,1)



!---------------------------------------------------------------------------
! Add contribution of the horizontal pressure gradient force.  (term 5 of 5)
!---------------------------------------------------------------------------

! Skip bottom layer which is calculated using Gauss-Seidel method back
! in subroutine update_diagnostics


! Move upward

do k = 2, nlm-1

   f1(:,:) = - p25 * ( moist_corr1_hpgf(:,:,k) +                     &
                      moist_corr1_hpgf(im1(:),:,k) ) *               &
                  ( th(:,:,k) + th(im1(:),:,k) ) *                   &
                  ( P_exnr(:,:,k) - P_exnr(im1(:),:,k) )

   f2(:,:) = (c1/grav) * p5 *                                        &
                ( p5*(vpgf_l2(:,:,k)+vpgf_l2(im1(:),:,k)) *          &
                  (phi_l2(:,:,k)-phi_l2(im1(:),:,k))         +       &
                  p5*(vpgf_l2(:,:,k-1)+vpgf_l2(im1(:),:,k-1)) *      &
                  (phi_l2(:,:,k-1)-phi_l2(im1(:),:,k-1))  )
   
   u_f(:,:,k) = u_f(:,:,k) +                                         &
                      (c1/dx) * ( f1(:,:) + f2(:,:) )


   f1(:,:) = - p25 * ( moist_corr1_hpgf(:,:,k) +                     &
                      moist_corr1_hpgf(:,jm1(:),k) ) *               &
                  ( th(:,:,k) + th(:,jm1(:),k) ) *                   &
                  ( P_exnr(:,:,k) - P_exnr(:,jm1(:),k) )

   f2(:,:) = (c1/grav) * p5 *                                        &
                ( p5*(vpgf_l2(:,:,k)+vpgf_l2(:,jm1(:),k)) *          &
                  (phi_l2(:,:,k)-phi_l2(:,jm1(:),k))         +       &
                  p5*(vpgf_l2(:,:,k-1)+vpgf_l2(:,jm1(:),k-1)) *      &
                  (phi_l2(:,:,k-1)-phi_l2(:,jm1(:),k-1))  )

   v_f(:,:,k) = v_f(:,:,k) +                                         &
                      (c1/dy) * ( f1(:,:) + f2(:,:) )

end do

! Top layer treated differently due to z_top = constant

f1(:,:) = - p25 * ( moist_corr1_hpgf(:,:,nlm) +                      &
                   moist_corr1_hpgf(im1(:),:,nlm) ) *                &
               ( th(:,:,nlm) + th(im1(:),:,nlm) ) *                  &
               ( P_exnr(:,:,nlm) - P_exnr(im1(:),:,nlm) )

f2(:,:) = (c1/grav) * p5 *                                           &
             ( p5*(vpgf_l2(:,:,nlm-1)+vpgf_l2(im1(:),:,nlm-1)) *     &
               (phi_l2(:,:,nlm-1)-phi_l2(im1(:),:,nlm-1))  )
   
u_f(:,:,nlm) = u_f(:,:,nlm) +                                        &
                   (c1/dx) * ( f1(:,:) + f2(:,:) )


f1(:,:) = - p25 * ( moist_corr1_hpgf(:,:,nlm) +                      &
                   moist_corr1_hpgf(:,jm1(:),nlm) ) *                &
               ( th(:,:,nlm) + th(:,jm1(:),nlm) ) *                  &
               ( P_exnr(:,:,nlm) - P_exnr(:,jm1(:),nlm) )

f2(:,:) = (c1/grav) * p5 *                                           &
             ( p5*(vpgf_l2(:,:,nlm-1)+vpgf_l2(:,jm1(:),nlm-1)) *     &
               (phi_l2(:,:,nlm-1)-phi_l2(:,jm1(:),nlm-1))  )

v_f(:,:,nlm) = v_f(:,:,nlm) +                                        &
                   (c1/dy) * ( f1(:,:) + f2(:,:) )



!
! Add del4 diffusion to terms u_f and v_f
!
do k = 1, nlm

   ! First, calculate Laplacians of u and v

   f1(:,:) = invdx2 * (u(ip1(:),:,k)-2*u(:,:,k)+u(im1(:),:,k)) +     &
             invdy2 * (u(:,jp1(:),k)-2*u(:,:,k)+u(:,jm1(:),k))
   f2(:,:) = invdx2 * (v(ip1(:),:,k)-2*v(:,:,k)+v(im1(:),:,k)) +     &
             invdy2 * (v(:,jp1(:),k)-2*v(:,:,k)+v(:,jm1(:),k))


   ! Now, calculate Laplacians of Laplacians of u and v, multiply by
   ! diffusion coefficient, and subtract from u_f and v_f

   u_f(:,:,k) = u_f(:,:,k)   -   kappa_uv *    (                     &
          invdx2 * (f1(ip1(:),:)-2*f1(:,:)+f1(im1(:),:)) +           &
          invdy2 * (f1(:,jp1(:))-2*f1(:,:)+f1(:,jm1(:)))   )

   v_f(:,:,k) = v_f(:,:,k)   -   kappa_uv *    (                     &
          invdx2 * (f2(ip1(:),:)-2*f2(:,:)+f2(im1(:),:)) +           &
          invdy2 * (f2(:,jp1(:))-2*f2(:,:)+f2(:,jm1(:)))   )

end do



end subroutine get_uf_vf

!======================================================================
! END OF GET_UF_VF
!======================================================================



end module momentum_tendency

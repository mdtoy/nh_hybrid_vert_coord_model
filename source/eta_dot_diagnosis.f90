module eta_dot_diagnosis

!-----------------------------------------------------------------------
! PURPOSE: Diagnoses vertical velocity eta_dot and in the process,
!          potential temperature and geopotential are predicted.
! 
!-----------------------------------------------------------------------

use kinds
use model_parameters
use physical_parameters
use eta_coordinate
use prognostics
use physics
use semi_implicit_solvers


implicit none
save




real (kind=dbl_kind), dimension (im,jm) ::                   &
     th_l2_f_sigma_0, th_l2_f_sigma_nlm,  &  ! theta tendencies at bottom and
                                             ! top boundaries due to Takacs
                                             ! vertical (sigma) advection
     th_l2_f_theta_0, th_l2_f_theta_nlm      ! theta tendencies at bottom and
                                             ! top boundaries due horizontal
                                             ! advection

real (kind=dbl_kind), dimension(im,jm,1:nlm-1) ::            &
          sigma_l2_star,   &  ! iterated value of sigma_l2
          th_l2_star          ! iterated value of th_l2

real (kind=dbl_kind), parameter ::                                   &
          r_fac_param1 = c1/(r_fac+c1)   ! miscellaneous parameter

real (kind=dbl_kind), parameter ::                           &
!           tau_eta_l2 = 1.e+15_dbl_kind 
          tau_eta_l2 = 1800._dbl_kind    ! relaxation time constant to return
                                         ! F(sgma,theta) to target value
                                         ! i.e., eta_l2  (s)




contains


!======================================================================
! BEGINNING OF GET_ETA_DOT_L2_1
!======================================================================

subroutine get_eta_dot_l2_1 (w1, w2, w3, P_exnr, th, phi, F_x, F_y,  &
                   m_PI_l2, exp_horiz_phi_advec, eta_dot_sigma_coeff)

!---------------------------------------------------------------------------
! PURPOSE:
!   Diagnoses the vertical velocity eta_dot_l2 in a generalized vertical
!   coordinate.  Also potential temperature and geopotential are time
!   stepped.
!   Subroutine 1 of 2.
!---------------------------------------------------------------------------

implicit none

!---------------------------------------------------------------------------
! INTENT IN
!---------------------------------------------------------------------------
real (kind=dbl_kind), intent(in) :: w1, w2, w3    ! Time stepping
                                                  ! weights

real (kind=dbl_kind), dimension(im,jm,nlm), intent(in) ::            &
          P_exnr,  &     ! Exner function at layer centers (J/kg/K)
          th,      &     ! potential temperature at layer centers (K)
          phi            ! geopotential at layer centers (J/kg)

real (kind=dbl_kind), dimension(im,jm,nlm), intent(in) ::            &
           F_x,          &  ! "3rd-order" mass flux in x-direction
                            ! (kg/s/m/eta) (colocated with u-points)
           F_y              ! "3rd-order" mass flux in y-direction
                            ! (kg/s/m/eta) (colocated with v-points)

!---------------------------------------------------------------------------
! INTENT OUT
!---------------------------------------------------------------------------
real (kind=dbl_kind), dimension (im,jm,0:nlm), intent(out) ::        &
          m_PI_l2           ! Exner function weighted pseudo-density at layer
                            ! edges (J/m^2/K)

real (kind=dbl_kind), dimension(im,jm,1:nlm-1), intent(out) ::       &
          exp_horiz_phi_advec, &  ! explicitly weighted phi advection
          eta_dot_sigma_coeff     ! coefficient for calculating
                                  ! eta_dot_l2_sigma

!---------------------------------------------------------------------------
! LOCAL
!---------------------------------------------------------------------------
integer :: i,j,k

real (kind=dbl_kind), dimension (im,jm,0:nlm) ::                     &
          sigma_l2      ! normalized height (-)

real (kind=dbl_kind), dimension (im,jm,1:nlm-1) ::                   &
          dF_dsigma_l2,   & ! partial derivative of F (functional definition
                            ! of vertical coordinate eta) w.r.t. sigma holding
                            ! theta constant (-)
          dF_dtheta_l2      ! partial derivative of F w.r.t. theta holding
                            ! sigma constant (K^-1)


! The "three forcings"
real (kind=dbl_kind), dimension (im,jm,1:nlm-1) ::                   &
          theta_forcing,  &  ! direct effect of horizontal theta
                             ! advection on theta tendency
          Q_diab_forcing     ! direct effect of diabatic heating on theta
                             ! tendency


! Variables and parameters used in Takacs "third-order" advection scheme
real (kind=dbl_kind), dimension(im,jm) ::                            &
          F_thphi_horiz_l2      ! x- & y- "3rd.-order" mass-weighted
                                ! theta/phi-fluxcomponents. (K kg/s/m/eta)/
                                ! (W/m/eta) (located at theta/phi-levels 
                                ! with horizontal positions corresponding 
                                ! to u and v respectively)
real (kind=dbl_kind), dimension(im,jm) ::                            &
          G_flx_contrib

real (kind=dbl_kind), dimension(im,jm,0:1) ::                        &
          F_PI_x_l2, F_PI_y_l2  ! horizontal Exner function weighted mass
                                ! flux interpolated to layer edges

real (kind=dbl_kind), dimension(im,jm) ::                            &
          F_horiz_l2    ! horizontal mass flux interpolated to layer edges

integer :: kp0, kp1     ! indices for k level (kp0 = k+0) and the k+1
                        ! level (kp1 = k+1)  (they toggle between 0 and 1)


! "New generation" variables
real (kind=dbl_kind), dimension(im,jm,0:nlm) ::                      &
          F_th_sgma       ! diagnosed F(theta,sigma)
real (kind=dbl_kind), parameter ::                                   &
          inv_dF_deta_cutoff = c1/dF_deta_cutoff, &  ! inverse dF_deta_cutoff
          inv_dF_deta_cutoff_sqr = inv_dF_deta_cutoff**2   ! square of inverse


real (kind=dbl_kind) :: fpt1                     ! working variable
real (kind=dbl_kind), dimension(im,jm) :: f1,f2  ! more working variables
real (kind=dbl_kind), dimension(im,jm) ::  &
          ftend1,ftend2                          ! even more working variables


real (kind=dbl_kind), dimension(im,jm,1:nlm-1,nttend) ::     &
          horiz_phi_advec       ! horizontal phi advection at multiple
                                ! time levels




!---------------------------------------------------------------------------
! Calculate initial sigma_l2's, dF_dsigma_l2's and dF_dtheta_l2's
!---------------------------------------------------------------------------
sigma_l2(:,:,0) = c0
do k = 1, nlm-1
   sigma_l2(:,:,k) = invgrav*(phi_l2(:,:,k,n4)-phis(:,:))/h_1(:,:)
   f1(:,:) = (c1-sigma_l2(:,:,k))**r_fac
   dF_dsigma_l2(:,:,k) = (th_l2(:,:,k,n4)-theta_min) * r_fac *       &
                            f1(:,:)/(c1-sigma_l2(:,:,k)) -           &
                            dth_dsigma_min * (c1-f1(:,:))
   dF_dtheta_l2(:,:,k) = c1-f1(:,:)
end do
sigma_l2(:,:,nlm) = c1



!---------------------------------------------------------------------------
! Calculate initial F_th_sgma and dF_deta
!---------------------------------------------------------------------------
F_th_sgma(:,:,0) = eta_l2(0)      ! bottom layer edge (invariant)
do k = 1, nlm                     ! move up
   do j = 1, jm
      do i = 1, im
         fpt1 = (c1-sigma_l2(i,j,k))**r_fac
         F_th_sgma(i,j,k) = theta_min * fpt1 + dth_dsigma_min *      &
            (c1-sigma_l2(i,j,k)) * (c1-r_fac_param1*fpt1) +          &
               th_l2(i,j,k,n4) * (c1-fpt1)
      end do
   end do
end do

do k = 1, nlm-1
   do j = 1, jm
      do i = 1, im
         dF_deta(i,j,k) = p5*(F_th_sgma(i,j,k+1)-F_th_sgma(i,j,k-1))*  &
                             inv_d_eta_l2(k)
      end do
   end do
end do





!---------------------------------------------------------------------------
! Calculate Exner function weighted pseudo-density at layer edges
!---------------------------------------------------------------------------
m_PI_l2(:,:,0) = p5 * P_exnr(:,:,1)*m(:,:,1,n4)*d_eta(1)
do k = 1, nlm-1
   m_PI_l2(:,:,k) = p5 * ( P_exnr(:,:,k+1)*m(:,:,k+1,n4)*d_eta(k+1) +  &
                           P_exnr(:,:,k)*m(:,:,k,n4)*d_eta(k) )
end do
m_PI_l2(:,:,nlm) = p5 * P_exnr(:,:,nlm)*m(:,:,nlm,n4)*d_eta(nlm)




!---------------------------------------------------------------------------
! Calculate geopotential and potential temperature tendencies due to
! horizontal advection.  Takacs advection scheme is used.
!
! Also calculate theta-, sigma- and Q_diab-forcings and the initial
! estimates of eta_dot_l2_sigma, eta_dot_l2_theta and eta_dot_l2_Q.
!---------------------------------------------------------------------------


!---------------------------------------------------------------------------
! start at bottom layer edge (surface)
!---------------------------------------------------------------------------

kp0 = mod(0,2)
kp1 = mod(1,2)

! x-contribution (theta)
F_horiz_l2(:,:) = F_x(:,:,1)

G_flx_contrib(:,:) = - inv12 * abs(F_horiz_l2(:,:)) *                   &
   ( 3*(th_l2(:,:,0,n4)-th_l2(im1(:),:,0,n4)) - th_l2(ip1(:),:,0,n4) +  &
     th_l2(im2(:),:,0,n4) )

f1(:,:) = inv12*F_x(:,:,1)*                                             &
   ( 7*(P_exnr(:,:,1)*th_l2(:,:,0,n4)+P_exnr(im1(:),:,1)*               &
    th_l2(im1(:),:,0,n4)) - P_exnr(ip1(:),:,1)*th_l2(ip1(:),:,0,n4) -   &
    P_exnr(im2(:),:,1)*th_l2(im2(:),:,0,n4) ) * d_eta(1)
F_thphi_horiz_l2(:,:) = p5 * f1(:,:)

F_PI_x_l2(:,:,kp1) = inv12 * F_x(:,:,1) *                               &
   ( 7*(P_exnr(:,:,1)+P_exnr(im1(:),:,1)) - P_exnr(ip1(:),:,1) -        &
        P_exnr(im2(:),:,1) ) * d_eta(1)

f1(:,:) = (c1/m_PI_l2(:,:,0)) *                                         &
            ( F_thphi_horiz_l2(ip1(:),:) - F_thphi_horiz_l2(:,:) -      &
   th_l2(:,:,0,n4)*p5*( F_PI_x_l2(ip1(:),:,kp1) - F_PI_x_l2(:,:,kp1) ) )
ftend1(:,:) = f1(:,:) + (c1/m_l2(:,:,0))*                               &
                 ( G_flx_contrib(ip1(:),:) - G_flx_contrib(:,:) )


! y-contribution (theta)
F_horiz_l2(:,:) = F_y(:,:,1)

G_flx_contrib(:,:) = - inv12 * abs(F_horiz_l2(:,:)) *                   &
   ( 3*(th_l2(:,:,0,n4)-th_l2(:,jm1(:),0,n4)) - th_l2(:,jp1(:),0,n4) +  &
     th_l2(:,jm2(:),0,n4) )

f1(:,:) = inv12*F_y(:,:,1)*                                             &
   ( 7*(P_exnr(:,:,1)*th_l2(:,:,0,n4)+P_exnr(:,jm1(:),1)*               &
    th_l2(:,jm1(:),0,n4)) - P_exnr(:,jp1(:),1)*th_l2(:,jp1(:),0,n4) -   &
    P_exnr(:,jm2(:),1)*th_l2(:,jm2(:),0,n4) ) * d_eta(1)
F_thphi_horiz_l2(:,:) = p5 * f1(:,:)

F_PI_y_l2(:,:,kp1) = inv12 * F_y(:,:,1) *                               &
   ( 7*(P_exnr(:,:,1)+P_exnr(:,jm1(:),1)) - P_exnr(:,jp1(:),1) -        &
        P_exnr(:,jm2(:),1) ) * d_eta(1)

f1(:,:) = (c1/m_PI_l2(:,:,0)) *                                         &
      ( F_thphi_horiz_l2(:,jp1(:)) - F_thphi_horiz_l2(:,:) -            &
   th_l2(:,:,0,n4)*p5*( F_PI_y_l2(:,jp1(:),kp1) - F_PI_y_l2(:,:,kp1) ) )
ftend2(:,:) = f1(:,:) + (c1/m_l2(:,:,0))*                               &
                     ( G_flx_contrib(:,jp1(:)) - G_flx_contrib(:,:) )


! Calculate theta tendency
th_l2_f_theta_0(:,:) = - ( invdx*ftend1(:,:) + invdy*ftend2(:,:) )



!---------------------------------------------------------------------------
! move upward
!
! Note:  In this loop we calculate forcings and initial values of
!        eta_dot_l2_x's.
!---------------------------------------------------------------------------

do k = 1, nlm-1

   kp0 = mod(k,2)
   kp1 = mod(k+1,2)

   ! x-contribution (theta)
   F_horiz_l2(:,:) = p5*inv_d_eta_l2(k) *                                  &
                    ( F_x(:,:,k)*d_eta(k) + F_x(:,:,k+1)*d_eta(k+1) )

   G_flx_contrib(:,:) = - inv12 * abs(F_horiz_l2(:,:)) *                   &
      ( 3*(th_l2(:,:,k,n4)-th_l2(im1(:),:,k,n4)) - th_l2(ip1(:),:,k,n4) +  &
        th_l2(im2(:),:,k,n4) )

   f1(:,:) = inv12*F_x(:,:,k+1)*                                           &
      ( 7*(P_exnr(:,:,k+1)*th_l2(:,:,k,n4)+P_exnr(im1(:),:,k+1)*           &
       th_l2(im1(:),:,k,n4)) - P_exnr(ip1(:),:,k+1)*th_l2(ip1(:),:,k,n4) - &
       P_exnr(im2(:),:,k+1)*th_l2(im2(:),:,k,n4) ) * d_eta(k+1)
   f2(:,:) = inv12*F_x(:,:,k)*                                             &
      ( 7*(P_exnr(:,:,k)*th_l2(:,:,k,n4)+P_exnr(im1(:),:,k)*               &
       th_l2(im1(:),:,k,n4)) - P_exnr(ip1(:),:,k)*th_l2(ip1(:),:,k,n4) -   &
       P_exnr(im2(:),:,k)*th_l2(im2(:),:,k,n4) ) * d_eta(k)
   F_thphi_horiz_l2(:,:) = p5 * ( f1(:,:) + f2(:,:) )

   F_PI_x_l2(:,:,kp1) = inv12 * F_x(:,:,k+1) *                            &
      ( 7*(P_exnr(:,:,k+1)+P_exnr(im1(:),:,k+1)) - P_exnr(ip1(:),:,k+1) - &
           P_exnr(im2(:),:,k+1) ) * d_eta(k+1)

   f1(:,:) = (c1/m_PI_l2(:,:,k)) *                                        &
               (  F_thphi_horiz_l2(ip1(:),:) - F_thphi_horiz_l2(:,:) -    &
                 th_l2(:,:,k,n4) * p5 *                                   & 
                    ( (F_PI_x_l2(ip1(:),:,kp1)+F_PI_x_l2(ip1(:),:,kp0)) - &
                      (F_PI_x_l2(:,:,kp1)+F_PI_x_l2(:,:,kp0)) )  )
   theta_forcing(:,:,k) = f1(:,:) + (c1/m_l2(:,:,k))*                     &
                     ( G_flx_contrib(ip1(:),:) - G_flx_contrib(:,:) )


   ! x-contribution (phi)
   G_flx_contrib(:,:) = - inv12 * abs(F_horiz_l2(:,:)) *             &
      ( 3*(phi_l2(:,:,k,n4)-phi_l2(im1(:),:,k,n4)) -                 &
        phi_l2(ip1(:),:,k,n4) + phi_l2(im2(:),:,k,n4) )
   F_thphi_horiz_l2(:,:) = inv12 * F_horiz_l2(:,:) *                 &
      ( 7*(phi_l2(:,:,k,n4)+phi_l2(im1(:),:,k,n4)) -                 &
        phi_l2(ip1(:),:,k,n4) - phi_l2(im2(:),:,k,n4) ) +            &
           G_flx_contrib(:,:)
   horiz_phi_advec(:,:,k,n3_f) = F_thphi_horiz_l2(ip1(:),:) -        &
                            F_thphi_horiz_l2(:,:) -                  &
                         phi_l2(:,:,k,n4) *                          & 
                         ( F_horiz_l2(ip1(:),:) - F_horiz_l2(:,:) )



   ! y-contribution (theta)
   F_horiz_l2(:,:) = (p5*inv_d_eta_l2(k)) *                                &
                     ( F_y(:,:,k)*d_eta(k) + F_y(:,:,k+1)*d_eta(k+1) )

   G_flx_contrib(:,:) = - inv12 * abs(F_horiz_l2(:,:)) *                   &
      ( 3*(th_l2(:,:,k,n4)-th_l2(:,jm1(:),k,n4)) - th_l2(:,jp1(:),k,n4) +  &
        th_l2(:,jm2(:),k,n4) )

   f1(:,:) = inv12*F_y(:,:,k+1)*                                           &
      ( 7*(P_exnr(:,:,k+1)*th_l2(:,:,k,n4)+P_exnr(:,jm1(:),k+1)*           &
       th_l2(:,jm1(:),k,n4)) - P_exnr(:,jp1(:),k+1)*th_l2(:,jp1(:),k,n4) - &
       P_exnr(:,jm2(:),k+1)*th_l2(:,jm2(:),k,n4) ) * d_eta(k+1)
   f2(:,:) = inv12*F_y(:,:,k)*                                             &
      ( 7*(P_exnr(:,:,k)*th_l2(:,:,k,n4)+P_exnr(:,jm1(:),k)*               &
       th_l2(:,jm1(:),k,n4)) - P_exnr(:,jp1(:),k)*th_l2(:,jp1(:),k,n4) -   &
       P_exnr(:,jm2(:),k)*th_l2(:,jm2(:),k,n4) ) * d_eta(k)
   F_thphi_horiz_l2(:,:) = p5 * ( f1(:,:) + f2(:,:) )

   F_PI_y_l2(:,:,kp1) = inv12 * F_y(:,:,k+1) *                            &
      ( 7*(P_exnr(:,:,k+1)+P_exnr(:,jm1(:),k+1)) - P_exnr(:,jp1(:),k+1) - &
           P_exnr(:,jm2(:),k+1) ) * d_eta(k+1)

   f1(:,:) = (c1/m_PI_l2(:,:,k)) *                                        &
               ( F_thphi_horiz_l2(:,jp1(:)) - F_thphi_horiz_l2(:,:) -     &
                 th_l2(:,:,k,n4) * p5 *                                   & 
                    ( (F_PI_y_l2(:,jp1(:),kp1)+F_PI_y_l2(:,jp1(:),kp0)) - &
                      (F_PI_y_l2(:,:,kp1)+F_PI_y_l2(:,:,kp0)) )   )
   ftend1(:,:) = f1(:,:) + (c1/m_l2(:,:,k))*                         &
                     ( G_flx_contrib(:,jp1(:)) - G_flx_contrib(:,:) )


   ! y-contribution (phi)
   G_flx_contrib(:,:) = - inv12 * abs(F_horiz_l2(:,:)) *             &
      ( 3*(phi_l2(:,:,k,n4)-phi_l2(:,jm1(:),k,n4)) -                 &
        phi_l2(:,jp1(:),k,n4) + phi_l2(:,jm2(:),k,n4) )
   F_thphi_horiz_l2(:,:) = inv12 * F_horiz_l2(:,:) *                 &
      ( 7*(phi_l2(:,:,k,n4)+phi_l2(:,jm1(:),k,n4)) -                 &
        phi_l2(:,jp1(:),k,n4) - phi_l2(:,jm2(:),k,n4) ) +            &
           G_flx_contrib(:,:)
   ftend2(:,:) = F_thphi_horiz_l2(:,jp1(:)) - F_thphi_horiz_l2(:,:) -   &
                    phi_l2(:,:,k,n4) *                                  & 
                       ( F_horiz_l2(:,jp1(:)) - F_horiz_l2(:,:) )




   !
   ! Calculate forcings
   !

   ! Update sigma forcing
   horiz_phi_advec(:,:,k,n3_f) = - (c1/m_l2(:,:,k)) *                   &
               (invdx*horiz_phi_advec(:,:,k,n3_f)+invdy*ftend2(:,:))
   exp_horiz_phi_advec(:,:,k) = w3*horiz_phi_advec(:,:,k,n3_f) +        &
        w2*horiz_phi_advec(:,:,k,n2_f) + w1*horiz_phi_advec(:,:,k,n1_f)


   ! Update theta forcing
   theta_forcing(:,:,k) = - (invdx*theta_forcing(:,:,k)+invdy*ftend1(:,:))

   ! Calculate Q_diab forcing which for now only includes contribution
   ! from subgrid-scale turbulent heat flux
   Q_diab_forcing(:,:,k) = F_turb_th_l2(:,:,k)




   !
   ! Calculate first estimate of eta_dot_l2_x's
   !

   do j = 1, jm
      do i = 1, im

         if ( dF_deta(i,j,k) .ge. dF_deta_cutoff ) then

            eta_dot_sigma_coeff(i,j,k) = dF_dsigma_l2(i,j,k)*              &
                                inv_h_1(i,j)*invgrav/dF_deta(i,j,k)
            eta_dot_l2_sigma(i,j,k) = eta_dot_sigma_coeff(i,j,k)*          &
                    ( w_l2(i,j,k,n4)*grav + exp_horiz_phi_advec(i,j,k) )

            eta_dot_l2_theta(i,j,k) = dF_dtheta_l2(i,j,k) *                &
                               theta_forcing(i,j,k) / dF_deta(i,j,k)

            eta_dot_l2_Q_diab(i,j,k) = dF_dtheta_l2(i,j,k) *               &
                               Q_diab_forcing(i,j,k) / dF_deta(i,j,k)

         else if ( dF_deta(i,j,k) .gt. 0 ) then

            ! dF_deta is greater than 0 and less than dF_deta_cutoff
            eta_dot_sigma_coeff(i,j,k) = (c1-dF_deta(i,j,k)*               &
               inv_dF_deta_cutoff)*d_eta_l2(k)/(phi(i,j,k+1)-phi(i,j,k)) + &
               inv_h_1(i,j)*invgrav*dF_deta(i,j,k)*inv_dF_deta_cutoff_sqr* &
               dF_dsigma_l2(i,j,k)
            eta_dot_l2_sigma(i,j,k) = eta_dot_sigma_coeff(i,j,k)*          &
                    ( w_l2(i,j,k,n4)*grav + exp_horiz_phi_advec(i,j,k) )

            eta_dot_l2_theta(i,j,k) = dF_deta(i,j,k)*inv_dF_deta_cutoff_sqr*  &
                               dF_dtheta_l2(i,j,k)*theta_forcing(i,j,k)

            eta_dot_l2_Q_diab(i,j,k) = dF_deta(i,j,k)*inv_dF_deta_cutoff_sqr* &
                               dF_dtheta_l2(i,j,k)*Q_diab_forcing(i,j,k)

         else

            ! dF_deta is less than or equal to 0
            eta_dot_sigma_coeff(i,j,k) = d_eta_l2(k) /                     &
                                            (phi(i,j,k+1)-phi(i,j,k))
            eta_dot_l2_sigma(i,j,k) = eta_dot_sigma_coeff(i,j,k)*          &
                    ( w_l2(i,j,k,n4)*grav + exp_horiz_phi_advec(i,j,k) )

            eta_dot_l2_theta(i,j,k) = c0

            eta_dot_l2_Q_diab(i,j,k) = c0

         end if

      end do
   end do


end do  !  do k = 1, nlm-1



!---------------------------------------------------------------------------
! top layer edge (model top)
!---------------------------------------------------------------------------

kp0 = mod(nlm,2)
kp1 = mod(nlm+1,2)

! x-contribution (theta)
F_horiz_l2(:,:) = F_x(:,:,nlm)

G_flx_contrib(:,:) = - inv12 * abs(F_horiz_l2(:,:)) *                      &
   ( 3*(th_l2(:,:,nlm,n4)-th_l2(im1(:),:,nlm,n4)) -                        &
     th_l2(ip1(:),:,nlm,n4) + th_l2(im2(:),:,nlm,n4) )

f1(:,:) = inv12*F_x(:,:,nlm)*                                              &
   ( 7*(P_exnr(:,:,nlm)*th_l2(:,:,nlm,n4)+P_exnr(im1(:),:,nlm)*            &
   th_l2(im1(:),:,nlm,n4)) - P_exnr(ip1(:),:,nlm)*th_l2(ip1(:),:,nlm,n4) - &
    P_exnr(im2(:),:,nlm)*th_l2(im2(:),:,nlm,n4) ) * d_eta(nlm)
F_thphi_horiz_l2(:,:) = p5 * f1(:,:)

f1(:,:) = (c1/m_PI_l2(:,:,nlm)) *                                          &
     ( F_thphi_horiz_l2(ip1(:),:) - F_thphi_horiz_l2(:,:) -                &
   th_l2(:,:,nlm,n4)*p5*( F_PI_x_l2(ip1(:),:,kp0) - F_PI_x_l2(:,:,kp0) )  )
ftend1(:,:) = f1(:,:) + (c1/m_l2(:,:,nlm))*                                &
                       ( G_flx_contrib(ip1(:),:) - G_flx_contrib(:,:) )


! y-contribution (theta)
F_horiz_l2(:,:) = F_y(:,:,nlm)

G_flx_contrib(:,:) = - inv12 * abs(F_horiz_l2(:,:)) *                      &
   ( 3*(th_l2(:,:,nlm,n4)-th_l2(:,jm1(:),nlm,n4)) -                        &
     th_l2(:,jp1(:),nlm,n4) + th_l2(:,jm2(:),nlm,n4) )

f1(:,:) = inv12*F_y(:,:,nlm)*                                              &
   ( 7*(P_exnr(:,:,nlm)*th_l2(:,:,nlm,n4)+P_exnr(:,jm1(:),nlm)*            &
   th_l2(:,jm1(:),nlm,n4)) - P_exnr(:,jp1(:),nlm)*th_l2(:,jp1(:),nlm,n4) - &
    P_exnr(:,jm2(:),nlm)*th_l2(:,jm2(:),nlm,n4) ) * d_eta(nlm)
F_thphi_horiz_l2(:,:) = p5 * f1(:,:)

f1(:,:) = (c1/m_PI_l2(:,:,nlm)) *                                          &
     ( F_thphi_horiz_l2(:,jp1(:)) - F_thphi_horiz_l2(:,:) -                &
   th_l2(:,:,nlm,n4)*p5*( F_PI_y_l2(:,jp1(:),kp0) - F_PI_y_l2(:,:,kp0) )  )
ftend2(:,:) = f1(:,:) + (c1/m_l2(:,:,nlm))*                                &
                       ( G_flx_contrib(:,jp1(:)) - G_flx_contrib(:,:) )



! Calculate theta tendency
th_l2_f_theta_nlm(:,:) = - ( invdx*ftend1(:,:) + invdy*ftend2(:,:) )




!===========================================================================
!---------------------------------------------------------------------------
!  Step 1:  Calculate eta_dot_l2_tau (for nudging F toward eta_l2)
!---------------------------------------------------------------------------
!===========================================================================

! Note that eta_dot_l2_tau will equal zero everywhere except where
! 0 < dF_deta < dF_deta_cutoff

do k = 1, nlm-1
   do j = 1, jm
      do i = 1, im

         if ( dF_deta(i,j,k) .lt. dF_deta_cutoff ) then

            if ( dF_deta(i,j,k) .gt. c0 ) then
               ! 0 < dF_deta < dF_deta_cutoff
               eta_dot_l2_tau(i,j,k) = dF_deta(i,j,k) *              &
                                        inv_dF_deta_cutoff_sqr *     &
                            (F_th_sgma(i,j,k)-eta_l2(k)) / tau_eta_l2
            else
               ! dF_deta <= 0
               eta_dot_l2_tau(i,j,k) = c0
            end if

         else

            ! dF_deta >= dF_deta_cutoff
            eta_dot_l2_tau(i,j,k) = c0

         end if

      end do
   end do
end do




!===========================================================================
!---------------------------------------------------------------------------
!  Step 2a:  Update phi_l2
!---------------------------------------------------------------------------
!===========================================================================

! Calculate phi tendency due to sigma forcing and vertical difference
! terms involving eta_dot_l2_sigma and eta_dot_l2_tau
! Note:  With semi-implicit time-diffencing, this is only AB3 effects
!        of vertical advection of phi by non-eta_dot_sigma velocity
do k = 1, nlm-1
   phi_l2_f(:,:,k,n3_f) = - (eta_dot_l2_theta(:,:,k)+                &
               eta_dot_l2_Q_diab(:,:,k)+eta_dot_l2_tau(:,:,k)) *     &
                   (phi(:,:,k+1)-phi(:,:,k))*inv_d_eta_l2(k)
end do

! Find resulting intermediate (starred) values of phi_l2.
!
! Note:  All modifications to th_l2 and phi_l2 will be made on the n3 time
!        level.  Values at the n4 time level are preserved.  Final "time-
!        stepped" values of th_l2 and phi_l2 at the end of this subroutine
!        will be at the n3 time level.
!        Also, the following is a temporary assignment of phi_l2, which
!        will be changed after smoothing and iterative eta_dot_l2_prime
!        diagnosis procedures have occurred.
! 
! With semi-implicit time-differencing of acoustic terms, this is a very
! preliminary value of geopotential
do k = 1, nlm-1
   do j = 1, jm
      do i = 1, im
         fpt1 = exp_horiz_phi_advec(i,j,k) + alpha_si*               &
            ( w_l2(i,j,k,n4)*grav - eta_dot_l2_sigma(i,j,k)*         &
              (phi(i,j,k+1)-phi(i,j,k))*inv_d_eta_l2(k) )
         phi_l2(i,j,k,n3) = phi_l2(i,j,k,n4) + dt *                  &
              ( w3*phi_l2_f(i,j,k,n3_f) + w2*phi_l2_f(i,j,k,n2_f) +  &
                w1*phi_l2_f(i,j,k,n1_f) + fpt1 )
      end do
   end do
end do



! Calculate theta tendency due to theta and Q_diab forcings and vertical
! difference terms involving eta_dot_l2_theta, eta_dot_l2_Q_diab and
! eta_dot_l2_tau
do k = 1, nlm-1
   ! AB3
   th_l2_f(:,:,k,n3_f) = theta_forcing(:,:,k) -                          &
                       (eta_dot_l2_theta(:,:,k)+eta_dot_l2_tau(:,:,k)) * &
                                 (th(:,:,k+1)-th(:,:,k))*inv_d_eta_l2(k)
   ! Euler forward
   th_l2_f_ef(:,:,k) = Q_diab_forcing(:,:,k) -                            &
                 eta_dot_l2_Q_diab(:,:,k) *                               &
                                 (th(:,:,k+1)-th(:,:,k))*inv_d_eta_l2(k)
end do


! Update theta tendency due to vertical Takacs advection by eta_dot_l2_sigma
! (AB3)
call vertical_Takacs_theta_advec_1 ( eta_dot_l2_sigma,m_PI_l2,          &
                                         P_exnr,th,th_l2_f(:,:,:,n3_f) )


! Calculate intermediate values of theta based on adiabatic and pre-
! smoothing processes
do k = 1, nlm-1
   th_l2(:,:,k,n3) = th_l2(:,:,k,n4) + dt *                          &
                ( w3*th_l2_f(:,:,k,n3_f) + w2*th_l2_f(:,:,k,n2_f) +  &
                  w1*th_l2_f(:,:,k,n1_f) + th_l2_f_ef(:,:,k) )
end do

!--------------------------------------------------------------------
! Calculate final values of theta tendency at bottom and top
! layers.  For now turbulent effect of turbulent heat flux
! at top layer is neglected
!--------------------------------------------------------------------
th_l2_f_0(:,:,n3_f) = th_l2_f_sigma_0(:,:) + th_l2_f_theta_0(:,:)
th_l2_f_nlm(:,:,n3_f) = th_l2_f_sigma_nlm(:,:) + th_l2_f_theta_nlm(:,:)



! Now continue with the semi-implicit process.
! Smoothing process (Step 2b) will be taken up in subroutine
! get_eta_dot_l2_2.


end subroutine get_eta_dot_l2_1

!======================================================================
! END OF GET_ETA_DOT_L2_1
!======================================================================





!======================================================================
! BEGINNING OF GET_ETA_DOT_L2_2
!======================================================================

subroutine get_eta_dot_l2_2 (w3, P_exnr, th, phi, m_PI_l2)

!---------------------------------------------------------------------------
! PURPOSE:
!   Diagnoses the vertical velocity eta_dot_l2 in a generalized vertical
!   coordinate.  Also potential temperature and geopotential are time
!   stepped.
!   Subroutine 2 of 2.
!---------------------------------------------------------------------------

implicit none

!---------------------------------------------------------------------------
! INTENT IN
!---------------------------------------------------------------------------
real (kind=dbl_kind), intent(in) :: w3    ! Time stepping weight

real (kind=dbl_kind), dimension(im,jm,nlm), intent(in) ::            &
          P_exnr,  &     ! Exner function at layer centers (J/kg/K)
          th,      &     ! potential temperature at layer centers (K)
          phi            ! geopotential at layer centers (J/kg)

real (kind=dbl_kind), dimension (im,jm,0:nlm), intent(in) ::         &
          m_PI_l2           ! Exner function weighted pseudo-density at layer
                            ! edges (J/m^2/K)

!---------------------------------------------------------------------------
! LOCAL
!---------------------------------------------------------------------------
integer :: i,j,k

! Variables and parameters used in vertical regridding procedure
real (kind=dbl_kind), parameter ::                                   &
          del2_phi_max = 1.E-04_dbl_kind*grav, &
          ! del2_phi_max = 2.E-04_dbl_kind*grav, &  ! maximum absolute value
                                  ! of the second spatial deriv of phi_l2 w.r.t.
                                  ! x,y for smoothing purposes
             ! Note:  set to 3.4E-07 (first try)
          del4_phi_max = 3.4E-11_dbl_kind*grav, &  ! maximum absolute value
                                  ! of the fourth spatial deriv of phi_l2 w.r.t.
                                  ! x,y for smoothing purposes
             ! Note:  set to 3.4E-11 for Boulder wind storm
          d2z_dz_max = 0.4_dbl_kind     ! maximum layer thickness ratio
real (kind=dbl_kind), parameter ::                                   &
          epsln_d4z = 0.70_dbl_kind,  &
          epsln_d2z = 0.52_dbl_kind     ! new method smoothing limit
                       ! (0 for "ultimate" smoothing, 1 for no smoothing)
real (kind=dbl_kind), parameter ::                                   &
          spatial_fac = 6*invdx4 + 8*invdx2dy2 + 6*invdy4,  &
          inv_spatial_fac = c1 / spatial_fac,               &
          spatial_fac2 = invdx2+invdy2,                     &
          inv_spatial_fac2 = c1 / spatial_fac2
real (kind=dbl_kind) :: del4_phi,  &   ! 4th. derivative of phi_l2 w.r.t. x,y
                        del2_phi,  &   ! 2nd. derivative of phi_l2 w.r.t. x,y
                        d2z_dz          ! layer thickness ratio
real (kind=dbl_kind), dimension(im,jm,1:nlm-1) ::            &
                        phi_l2_smooth        ! smoothed geopotential field
real (kind=dbl_kind), dimension(1:nlm-1) ::                  &
                        phi_l2_temp_column   ! temporary holder
real (kind=dbl_kind), dimension(im,jm) ::                    &
                        phi_l2_temp_level    ! temporary holder


logical (kind = log_kind) :: dF_deta_test

real (kind=dbl_kind) :: fpt1, fpt2, fpt3, fpt4   ! working variables



!===========================================================================
!---------------------------------------------------------------------------
!  Step 2b:  Smooth geopotential field (phi_l2) where necessary
!---------------------------------------------------------------------------
!===========================================================================

!---------------------------------------------------------------------------
!  Now we conditionally smooth the geopotential height of coordinate surfaces
!  and determine any new "smoothing points".  The conditional smoothing
!  is such that the horizontal curvature of phi_l2 (d2phi_dx2 -- the
!  second derivative in x and d2phi_dx2 -- the fourth derivative in x)
!  are kept below prescribed maximum values.
!---------------------------------------------------------------------------

!--------------------------------------------------------------------
! Initialize regridding contribution to vertical velocity
!--------------------------------------------------------------------
eta_dot_l2_regrid(:,:,:) = c0

!--------------------------------------------------------------------
! Perform smoothing:
!--------------------------------------------------------------------

! Check where abs(d2phi_dx2) exceeds limit and smooth if necessary
! Smoothing only along x for now (limited to 2D for now)
! Also, smooth in vertical


! Initialize "smoothed" geopotential field and logical
! "smooth_point" variable
phi_l2_smooth(:,:,:) = phi_l2(:,:,:,n3)
smooth_point(:,:,:) = .false.


!
! Perform vertical smoothing
!
do j = 1, jm
   do i = 1, im

      phi_l2_temp_column(:) = phi_l2_smooth(i,j,:)

      fpt1 = phi_l2_temp_column(2) + phis(i,j)
      fpt2 = p5 * ( phi_l2_temp_column(2) - phis(i,j) )
      d2z_dz = ( fpt1 - 2*phi_l2_temp_column(1) ) / fpt2
      if ( abs(d2z_dz) .gt. d2z_dz_max ) then
         smooth_point(i,j,1) = .true.   ! set to smooth point
                                        ! (if not already)
         ! perform vertical smoothing (new method)
         phi_l2_smooth(i,j,1) = p5 * ( fpt1 -                           &
                         d2z_dz_max*fpt2*abs(d2z_dz)/d2z_dz )
      end if
 
      do k = 2, nlm-2
         ! Calculate layer thickness ratio
         fpt1 = phi_l2_temp_column(k+1) + phi_l2_temp_column(k-1)
         fpt2 = p5 * ( phi_l2_temp_column(k+1) -                        &
                          phi_l2_temp_column(k-1) )
         d2z_dz = ( fpt1 - 2*phi_l2_temp_column(k) ) / fpt2
         if ( abs(d2z_dz) .gt. d2z_dz_max ) then
            smooth_point(i,j,k) = .true.   ! set to smooth point
                                           ! (if not already)
            ! perform vertical smoothing (new method)
            phi_l2_smooth(i,j,k) = p5 * ( fpt1 -                        &
                            d2z_dz_max*fpt2*abs(d2z_dz)/d2z_dz )
         end if
      end do

      fpt1 = grav*z_top + phi_l2_temp_column(nlm-2)
      fpt2 = p5 * (grav*z_top - phi_l2_temp_column(nlm-2))
      d2z_dz = ( fpt1 - 2*phi_l2_temp_column(nlm-1) ) / fpt2
      if ( abs(d2z_dz) .gt. d2z_dz_max ) then
         smooth_point(i,j,nlm-1) = .true.   ! set to smooth point
                                            ! (if not already)
         ! perform vertical smoothing (new method)
         phi_l2_smooth(i,j,nlm-1) = p5 * ( fpt1 -                       &
                         d2z_dz_max*fpt2*abs(d2z_dz)/d2z_dz )
      end if

   end do   ! i = 1, im
end do   ! j = 1, im


!
! Perform del-4 horizontal smoothing
!
do k = 1, nlm-1
   phi_l2_temp_level(:,:) = phi_l2_smooth(:,:,k)
   do j = 1, jm
      do i = 1, im
         ! Calculate fourth derivative of phi_l2 w.r.t. x
         ! (based on results from d2z_dz smoothing)
         fpt1 = phi_l2_temp_level(ip2(i),j) + phi_l2_temp_level(im2(i),j) - &
            4*( phi_l2_temp_level(ip1(i),j) + phi_l2_temp_level(im1(i),j) )
         fpt2 = phi_l2_temp_level(i,jp2(j)) + phi_l2_temp_level(i,jm2(j)) - &
            4*( phi_l2_temp_level(i,jp1(j)) + phi_l2_temp_level(i,jm1(j)) )
         fpt3 = phi_l2_temp_level(ip1(i),jp1(j))+                           &
           phi_l2_temp_level(ip1(i),jm1(j))-2*phi_l2_temp_level(ip1(i),j) - &
           2*(phi_l2_temp_level(i,jp1(j))+phi_l2_temp_level(i,jm1(j))) +    &
           phi_l2_temp_level(im1(i),jp1(j))-2*phi_l2_temp_level(im1(i),j)+  &
           phi_l2_temp_level(im1(i),jm1(j))
         del4_phi = spatial_fac * phi_l2_temp_level(i,j) + invdx4*fpt1 +    &
                       2*invdx2dy2*fpt3 + invdy4*fpt2
         if ( abs(del4_phi) .gt. del4_phi_max ) then
            smooth_point(i,j,k) = .true.   ! set to smooth point
                                           ! (if not already)
            ! perform horizontal del-4 smoothing (new method)
            fpt4 = del4_phi_max + epsln_d4z*(abs(del4_phi)-del4_phi_max)
                        ! fpt4 is modified smoothing "target"
            phi_l2_smooth(i,j,k) = inv_spatial_fac *                        &
               ( fpt4*abs(del4_phi)/del4_phi - invdx4*fpt1 - invdy4*fpt2 -  &
                 2*invdx2dy2*fpt3 )
         end if
      end do
   end do
end do


!
! Perform del-2 horizontal smoothing
!
do k = 1, nlm-1
   phi_l2_temp_level(:,:) = phi_l2_smooth(:,:,k)
   do j = 1, jm
      do i = 1, im
         ! Calculate second derivative of phi_l2 w.r.t. x
         ! (based on results from d2z_dz and del-4 smoothing)
         fpt1 = phi_l2_temp_level(ip1(i),j) + phi_l2_temp_level(im1(i),j)
         fpt2 = phi_l2_temp_level(i,jp1(j)) + phi_l2_temp_level(i,jm1(j))
         del2_phi = invdx2*fpt1 + invdy2*fpt2 - 2*spatial_fac2*             &
                       phi_l2_temp_level(i,j)
         if ( abs(del2_phi) .gt. del2_phi_max ) then
            smooth_point(i,j,k) = .true.   ! set to smooth point
                                           ! (if not already)
            ! perform horizontal del-2 smoothing (new method)
            fpt4 = del2_phi_max + epsln_d2z*(abs(del2_phi)-del2_phi_max)
                        ! fpt4 is modified smoothing "target"
            phi_l2_smooth(i,j,k) = p5*inv_spatial_fac2*                     &
               ( invdx2*fpt1 + invdy2*fpt2 - fpt4*abs(del2_phi)/del2_phi )
         end if
      end do
   end do
end do


! Now that smoothing field has been found, set actual field
! and calculate eta_dot_l2s for smoothed points
do k = 1, nlm-1
   do j = 1, jm
      do i = 1, im
         if ( smooth_point(i,j,k) ) then
            ! Note, this is an Euler-forward-like process
            eta_dot_l2_regrid(i,j,k) = -invdt*                          &
               (phi_l2_smooth(i,j,k)-phi_l2(i,j,k,n3)) * d_eta_l2(k) /  &
               (phi(i,j,k+1)-phi(i,j,k))
            ! Set final value of phi_l2
            phi_l2(i,j,k,n3) = phi_l2_smooth(i,j,k)
         end if
      end do
   end do
end do





!===========================================================================
!---------------------------------------------------------------------------
!  Step 3:  Update th_l2 based on forcings and now-known vertical velocities
!           and prepare for final steps (i.e., iterative eta_dot_l2_prime
!           diagnosis procedure)
!---------------------------------------------------------------------------
!===========================================================================


! Calculate theta tendency due to vertical Takacs advection by
! eta_dot_l2_regrid (Euler Forward)
call vertical_Takacs_theta_advec_2 ( eta_dot_l2_regrid,m_PI_l2,      &
                                         P_exnr,th,th_l2_f_ef )


! In preparation for final steps, find intermediate (starred) values of
! sigma_l2 and theta_l2 resulting from processes so far.
!
! Note:  All modifications to th_l2 and phi_l2 will be made on the n3 time
!        level.  Values at the n4 time level are preserved.  Final "time-
!        stepped" values of th_l2 and phi_l2 at the end of this subroutine
!        will be at the n3 time level.

do k = 1, nlm-1

   ! Calculate sigma_l2_star from updated phi_l2
   sigma_l2_star(:,:,k) = invgrav*(phi_l2(:,:,k,n3)-phis(:,:))/h_1(:,:)

   th_l2_star(:,:,k) = th_l2(:,:,k,n3) + dt * th_l2_f_ef(:,:,k)

end do



!===========================================================================
!---------------------------------------------------------------------------
!  Step 4:  Iterative eta_dot_l2_prime diagnosis and final determination
!           of phi_l2 and th_l2
!---------------------------------------------------------------------------
!===========================================================================

! Calculate F_set, the value to set F to in current time step
do k = 1, nlm-1
   do j = 1, jm
      do i = 1, im

         dF_deta_test = dF_deta(i,j,k).lt.dF_deta_cutoff
         if ( smooth_point(i,j,k) .or. in_cloud(i,j,k) .or.          &
              dF_deta_test ) then
            ! set F_set to current value of F
            ! ( this causes eta_dot_l2_prime to equal zero, unless
            ! dF_deta >= dF_deta_cutoff
            ! (see lines below) )
            fpt1 = (c1-sigma_l2_star(i,j,k))**r_fac
            F_set(i,j,k) = theta_min * fpt1 + dth_dsigma_min *       &
               (c1-sigma_l2_star(i,j,k)) * (c1-r_fac_param1*fpt1) +  &
                  th_l2_star(i,j,k) * (c1-fpt1)
         end if
         if ( .not.dF_deta_test ) then
            ! dF_deta >= dF_deta_cutoff:  Nudge F_set toward eta_l2
            ! (Note:  For 0 < dF_deta < dF_deta_cutoff, this is done
            !         with eta_dot_l2_tau)
            F_set(i,j,k) = F_set(i,j,k) + (eta_l2(k)-F_set(i,j,k))*   &
                                                          dt/tau_eta_l2
         end if

      end do
   end do
end do




! Now find iterative solution for eta_dot_l2_prime
call vert_mass_flx_iterative_solver ( phi, th, w3, eta_dot_l2_prime )



!-----------------------------------------------------------------
! Now that we have eta_dot_l2_prime, perform final update
! of tendencies of phi_l2 and th_l2
!-----------------------------------------------------------------
do k = 1, nlm-1

   phi_l2_f(:,:,k,n3_f) = phi_l2_f(:,:,k,n3_f)  -                      &
                           eta_dot_l2_prime(:,:,k) *                   &
                            (phi(:,:,k+1)-phi(:,:,k)) * inv_d_eta_l2(k)

   th_l2_f(:,:,k,n3_f) = th_l2_f(:,:,k,n3_f)  -                        &
                          eta_dot_l2_prime(:,:,k) *                    &
                           (th(:,:,k+1)-th(:,:,k)) * inv_d_eta_l2(k)


end do  !  k = 1, nlm-1




!===========================================================================
!---------------------------------------------------------------------------
!  Wrap-up steps 1 - 4.
!---------------------------------------------------------------------------
!===========================================================================

!--------------------------------------------------------------------
! Calculate final value of phi_l2 from final sigma_l2_star and
! update th_l2
!--------------------------------------------------------------------
do k = 1, nlm-1
   phi_l2(:,:,k,n3) = phis(:,:) + grav*h_1(:,:)*sigma_l2_star(:,:,k)
   th_l2(:,:,k,n3) = th_l2_star(:,:,k)
end do

!--------------------------------------------------------------------
! Finally, calculate eta_dot_l2 and vertical velocity used for
! "natural" Konor and Arakawa  calculations of mass and
! velocity tendencies and provisional values for smoothing points
!--------------------------------------------------------------------
do k = 1, nlm-1
   eta_dot_l2(:,:,k) =                                               &
          eta_dot_l2_sigma(:,:,k)  + eta_dot_l2_theta(:,:,k) +       &
          eta_dot_l2_Q_diab(:,:,k) + eta_dot_l2_tau(:,:,k)   +       &
          eta_dot_l2_regrid(:,:,k) + eta_dot_l2_prime(:,:,k)
end do




end subroutine get_eta_dot_l2_2

!======================================================================
! END OF GET_ETA_DOT_L2_2
!======================================================================








!======================================================================
! SUBROUTINES USED BY SUBROUTINE GET_ETA_DOT_L2
!======================================================================



!***************************************************************************
subroutine vertical_Takacs_theta_advec_1 (eta_dot_l2, m_PI_l2,       &
                              P_exnr, th, th_l2_f )

! Calculates tendency of theta due to vertical advection using
! Takacs scheme.
! This subroutine also calculates th_l2_f_sigma_0 and th_l2_f_sigma_nlm

implicit none

!---------------------------------------------------------------------------
! INTENT IN 
!---------------------------------------------------------------------------
real (kind=dbl_kind), dimension(im,jm,1:nlm-1), intent(in) ::      &
          eta_dot_l2     ! (AB3) generalized vert. velocity (eta/s)
real (kind=dbl_kind), dimension (im,jm,0:nlm), intent(in) ::       &
          m_PI_l2        ! Exner function weighted pseudo-density at layer
                         ! edges (J/m^2/K)
real (kind=dbl_kind), dimension(im,jm,nlm), intent(in) ::          &
          P_exnr,  &     ! Exner function at layer centers (J/kg/K)
          th             ! potential temperature at layer centers (K)

!---------------------------------------------------------------------------
! INTENT INOUT
!---------------------------------------------------------------------------
real (kind=dbl_kind), dimension(im,jm,1:nlm-1), intent(inout) ::     &
          th_l2_f        ! (AB3) Pot. temp. tendency due to vert. advection

!---------------------------------------------------------------------------
! LOCAL
!---------------------------------------------------------------------------
integer :: i, j, k
real (kind=dbl_kind), dimension(nlm) ::                            &
          F_PI_th_eta,  &  ! vertical "3rd.-order" mass- and PI-weighted 
                           ! theta-flux  (W/m^2)
          F_eta            ! vertical mass flux interpolated to layer centers
real (kind=dbl_kind), dimension(nlm) :: G_flx_contrib_vert
real (kind=dbl_kind) :: fpt1    ! working variable




do j = 1, jm
   do i = 1, im


      !
      ! Adams-Bashforth 3rd-order (AB3) Tendencies
      !

      !---------------------------------------------------------------------
      ! start at bottom layer edge (surface)
      !---------------------------------------------------------------------

      ! Calculate F_eta at all levels
      ! (do this only once)
      F_eta(1) = p5*m_l2(i,j,1)*eta_dot_l2(i,j,1)
      do k = 2,nlm-1
         F_eta(k) = p5*( m_l2(i,j,k)*eta_dot_l2(i,j,k) +             &
                         m_l2(i,j,k-1)*eta_dot_l2(i,j,k-1) )
      end do
      F_eta(nlm) = p5*m_l2(i,j,nlm-1)*eta_dot_l2(i,j,nlm-1)


      G_flx_contrib_vert(1) = c0

      F_PI_th_eta(1) = p5*F_eta(1)*P_exnr(i,j,1)*                    &
            (th_l2(i,j,1,n4)+th_l2(i,j,0,n4))

      fpt1 = (c1/m_PI_l2(i,j,0)) *                                   &
            ( F_PI_th_eta(1) -                                       &
                         th_l2(i,j,0,n4)*P_exnr(i,j,1)*F_eta(1)   )

      ! Calculate theta tendency
      th_l2_f_sigma_0(i,j) = - fpt1 - (c1/m_l2(i,j,0)) *             &
                        G_flx_contrib_vert(1) * inv_d_eta_l2(0)


      !---------------------------------------------------------------------
      ! move upward
      !---------------------------------------------------------------------
      
      do k = 1, nlm-2

         G_flx_contrib_vert(k+1) = - inv12 * abs(F_eta(k+1)) *       &
            ( 3*(th_l2(i,j,k+1,n4)-th_l2(i,j,k,n4)) -                &
              th_l2(i,j,k+2,n4) + th_l2(i,j,k-1,n4) )

         F_PI_th_eta(k+1) = inv12*F_eta(k+1)*P_exnr(i,j,k+1)*        &
            ( 7*(th_l2(i,j,k+1,n4)+th_l2(i,j,k,n4)) -                &
              th_l2(i,j,k+2,n4) - th_l2(i,j,k-1,n4) )

         fpt1 = (c1/m_PI_l2(i,j,k)) *                                &
                     ( F_PI_th_eta(k+1) - F_PI_th_eta(k) -           &
                  th_l2(i,j,k,n4) * ( P_exnr(i,j,k+1)*F_eta(k+1) -   &
                                      P_exnr(i,j,k)*F_eta(k) )   )

         ! Calculate theta tendency
         th_l2_f(i,j,k) = th_l2_f(i,j,k) - fpt1 - (c1/m_l2(i,j,k))*         &
             (G_flx_contrib_vert(k+1)-G_flx_contrib_vert(k))*inv_d_eta_l2(k)

      end do


      !---------------------------------------------------------------------
      ! second-from-the-top layer edge
      !---------------------------------------------------------------------

      G_flx_contrib_vert(nlm) = c0

      F_PI_th_eta(nlm) = p5*F_eta(nlm)*P_exnr(i,j,nlm)*              &
         (th_l2(i,j,nlm,n4)+th_l2(i,j,nlm-1,n4))

      fpt1 = (c1/m_PI_l2(i,j,nlm-1)) *                               &
           ( F_PI_th_eta(nlm) - F_PI_th_eta(nlm-1) -                 &
             th_l2(i,j,nlm-1,n4) * ( P_exnr(i,j,nlm)*F_eta(nlm) -    &
                P_exnr(i,j,nlm-1)*F_eta(nlm-1) )   )

      ! Calculate theta tendency
      th_l2_f(i,j,nlm-1) = th_l2_f(i,j,nlm-1) - fpt1 -               &
         (c1/m_l2(i,j,nlm-1)) *                                      &
           (G_flx_contrib_vert(nlm)-G_flx_contrib_vert(nlm-1)) *     &
                                                 inv_d_eta_l2(nlm-1)



      !---------------------------------------------------------------------
      ! top layer edge (model top)
      !---------------------------------------------------------------------

      fpt1 = (c1/m_PI_l2(i,j,nlm)) *                                 &
                 ( - F_PI_th_eta(nlm) +                              &
                    th_l2(i,j,nlm,n4)*P_exnr(i,j,nlm)*F_eta(nlm)   )

      ! Calculate theta tendency
      th_l2_f_sigma_nlm(i,j) = - fpt1 - (c1/m_l2(i,j,nlm)) *         &
                     ( - G_flx_contrib_vert(nlm) ) * inv_d_eta_l2(nlm)


   end do  !  i = 1, im
end do     !  j = 1, jm



end subroutine vertical_Takacs_theta_advec_1
!***************************************************************************





!***************************************************************************
subroutine vertical_Takacs_theta_advec_2 (eta_dot_l2, m_PI_l2,       &
                              P_exnr, th, th_l2_f)

! Calculates tendency of theta due to vertical advection using
! Takacs scheme.

implicit none

!---------------------------------------------------------------------------
! INTENT IN 
!---------------------------------------------------------------------------
real (kind=dbl_kind), dimension(im,jm,1:nlm-1), intent(in) ::      &
          eta_dot_l2     !  (EF) generalized vert. velocity (eta/s)
real (kind=dbl_kind), dimension (im,jm,0:nlm), intent(in) ::       &
          m_PI_l2        ! Exner function weighted pseudo-density at layer
                         ! edges (J/m^2/K)
real (kind=dbl_kind), dimension(im,jm,nlm), intent(in) ::          &
          P_exnr,  &     ! Exner function at layer centers (J/kg/K)
          th             ! potential temperature at layer centers (K)

!---------------------------------------------------------------------------
! INTENT OUT
!---------------------------------------------------------------------------
real (kind=dbl_kind), dimension(im,jm,1:nlm-1), intent(out) ::     &
          th_l2_f      !  (EF) Pot. temp. tendency due to vert. advection

!---------------------------------------------------------------------------
! LOCAL
!---------------------------------------------------------------------------
integer :: i, j, k
real (kind=dbl_kind), dimension(nlm) ::                            &
          F_PI_th_eta,  &  ! vertical "3rd.-order" mass- and PI-weighted 
                           ! theta-flux  (W/m^2)
          F_eta            ! vertical mass flux interpolated to layer centers
real (kind=dbl_kind), dimension(nlm) :: G_flx_contrib_vert
real (kind=dbl_kind) :: fpt1    ! working variable




do j = 1, jm
   do i = 1, im


      !
      ! Euler-forward (EF) Tendencies
      !

      !---------------------------------------------------------------------
      ! start at bottom layer edge (surface)
      !---------------------------------------------------------------------

      ! Calculate F_eta at all levels
      ! (do this only once)
      F_eta(1) = p5*m_l2(i,j,1)*eta_dot_l2(i,j,1)
      do k = 2,nlm-1
         F_eta(k) = p5*( m_l2(i,j,k)*eta_dot_l2(i,j,k) +             &
                         m_l2(i,j,k-1)*eta_dot_l2(i,j,k-1) )
      end do
      F_eta(nlm) = p5*m_l2(i,j,nlm-1)*eta_dot_l2(i,j,nlm-1)


      G_flx_contrib_vert(1) = c0

      F_PI_th_eta(1) = p5*F_eta(1)*P_exnr(i,j,1)*                    &
            (th_l2(i,j,1,n4)+th_l2(i,j,0,n4))


      !---------------------------------------------------------------------
      ! move upward
      !---------------------------------------------------------------------
      
      do k = 1, nlm-2

         G_flx_contrib_vert(k+1) = - inv12 * abs(F_eta(k+1)) *       &
            ( 3*(th_l2(i,j,k+1,n4)-th_l2(i,j,k,n4)) -                &
              th_l2(i,j,k+2,n4) + th_l2(i,j,k-1,n4) )

         F_PI_th_eta(k+1) = inv12*F_eta(k+1)*P_exnr(i,j,k+1)*        &
            ( 7*(th_l2(i,j,k+1,n4)+th_l2(i,j,k,n4)) -                &
              th_l2(i,j,k+2,n4) - th_l2(i,j,k-1,n4) )

         fpt1 = (c1/m_PI_l2(i,j,k)) *                                &
                     ( F_PI_th_eta(k+1) - F_PI_th_eta(k) -           &
                  th_l2(i,j,k,n4) * ( P_exnr(i,j,k+1)*F_eta(k+1) -   &
                                      P_exnr(i,j,k)*F_eta(k) )   )

         ! Calculate theta tendency
         th_l2_f(i,j,k) = - fpt1 - (c1/m_l2(i,j,k)) *                &
          (G_flx_contrib_vert(k+1)-G_flx_contrib_vert(k))*inv_d_eta_l2(k)

      end do


      !---------------------------------------------------------------------
      ! second-from-the-top layer edge
      !---------------------------------------------------------------------

      G_flx_contrib_vert(nlm) = c0

      F_PI_th_eta(nlm) = p5*F_eta(nlm)*P_exnr(i,j,nlm)*              &
         (th_l2(i,j,nlm,n4)+th_l2(i,j,nlm-1,n4))

      fpt1 = (c1/m_PI_l2(i,j,nlm-1)) *                               &
           ( F_PI_th_eta(nlm) - F_PI_th_eta(nlm-1) -                 &
             th_l2(i,j,nlm-1,n4) * ( P_exnr(i,j,nlm)*F_eta(nlm) -    &
                P_exnr(i,j,nlm-1)*F_eta(nlm-1) )   )

      ! Calculate theta tendency
      th_l2_f(i,j,nlm-1) = - fpt1 - (c1/m_l2(i,j,nlm-1)) *           &
           (G_flx_contrib_vert(nlm)-G_flx_contrib_vert(nlm-1)) *     &
                                                 inv_d_eta_l2(nlm-1)




   end do  !  i = 1, im
end do     !  j = 1, jm



end subroutine vertical_Takacs_theta_advec_2
!***************************************************************************






!***************************************************************************
subroutine vert_mass_flx_iterative_solver ( phi, th, w3, eta_dot_l2_prime )

! Iteratively solves for vertical mass flux which sets F(th,sgma) to its
! specified value of F_set

implicit none

!---------------------------------------------------------------------------
! INTENT IN 
!---------------------------------------------------------------------------
real (kind=dbl_kind), dimension(im,jm,nlm), intent(in) ::          &
          phi,  &        ! geopotential at layer centers (J/kg)
          th             ! potential temperature at layer centers (K)
real (kind=dbl_kind), intent(in) :: w3    ! Time stepping weight

!---------------------------------------------------------------------------
! INTENT OUT
!---------------------------------------------------------------------------
real (kind=dbl_kind), dimension(im,jm,1:nlm-1), intent(out) ::     &
          eta_dot_l2_prime    ! generalized vertical velocity (eta/s)

!---------------------------------------------------------------------------
! LOCAL
!---------------------------------------------------------------------------
integer :: i, j, k
integer :: iter                        ! iteration counter
integer, parameter :: max_iter = 100   ! max. number of iterations allowed
real (kind=dbl_kind) ::                                              &
           Delta_eta    ! difference between F_star_l2 and F_set, i.e.,
                        ! between the calculated and specified F(th,sgma)
real (kind=dbl_kind), parameter ::                                   &
      converge = 1.E-12_dbl_kind       ! convergence criterion
real (kind=dbl_kind) ::                                              &
      F_star_l2         ! functional evaluation of vertical coordinate eta (-)
real (kind=dbl_kind) ::                                              &
          dF_dsigma_l2,   & ! partial derivative of F (functional definition
                            ! of vertical coordinate eta) w.r.t. sigma holding
                            ! theta constant (-)
          dF_dtheta_l2      ! partial derivative of F w.r.t. theta holding
                            ! sigma constant (K^-1)
real (kind=dbl_kind) :: eta_dot_l2_incr   ! incremental eta_dot_l2_prime
real (kind=dbl_kind) :: fpt1, fpt2        ! working variables



do k = 1, nlm-1

   do j = 1, jm
      do i = 1, im

         ! Initialize eta_dot_l2_prime
         eta_dot_l2_prime(i,j,k) = c0
      
         iter = 0
         Delta_eta = 1.E+10_dbl_kind

         do while ( abs(Delta_eta) .ge. converge )

            ! Calculate F_star_l2 based on new values of sigma_l2_star
            ! and th_l2_star.         
            fpt1 = (c1-sigma_l2_star(i,j,k))**r_fac
            F_star_l2 = theta_min * fpt1 +                           &
                        dth_dsigma_min * (c1-sigma_l2_star(i,j,k)) * &
                        (c1-r_fac_param1*fpt1) +                     &
                        th_l2_star(i,j,k) * (c1-fpt1)

            ! Calculate Delta_eta
            Delta_eta = F_star_l2 - F_set(i,j,k)

            ! Calculate dF_dsigma_l2's and dF_dtheta_l2's
            dF_dsigma_l2 = (th_l2_star(i,j,k)-theta_min) *           &
                           r_fac * fpt1/(c1-sigma_l2_star(i,j,k)) -  &
                           dth_dsigma_min * (c1-fpt1)
            dF_dtheta_l2 = c1-fpt1


            !---------------------------------------------------------------
            ! Increment eta_dot_l2_prime based on new value of Delta_eta
            !---------------------------------------------------------------
            fpt1 = dF_dtheta_l2*(th(i,j,k+1)-th(i,j,k)) +            &
                   dF_dsigma_l2*(phi(i,j,k+1)-phi(i,j,k))/           &
                                          (h_1(i,j)*grav)
            fpt2 = Delta_eta * d_eta_l2(k) / (w3 * dt)
            eta_dot_l2_incr = fpt2 / fpt1
            eta_dot_l2_prime(i,j,k) = eta_dot_l2_prime(i,j,k) +      &
                                                   eta_dot_l2_incr


            !---------------------------------------------------------------
            ! Recalculate sigma_l2_star and th_l2_star based on incremental
            ! eta_dot_l2
            !---------------------------------------------------------------
            fpt1 = - eta_dot_l2_incr * (phi(i,j,k+1)-phi(i,j,k)) *   &
                                                     inv_d_eta_l2(k)
            sigma_l2_star(i,j,k) = sigma_l2_star(i,j,k) +            &
                                 dt * w3 * fpt1 / (h_1(i,j) * grav)

            fpt1 = - eta_dot_l2_incr * (th(i,j,k+1)-th(i,j,k)) *     &
                                                     inv_d_eta_l2(k)
            th_l2_star(i,j,k) = th_l2_star(i,j,k) + dt * w3 * fpt1

            iter = iter + 1
            if (iter .ge. max_iter) then
               print *
               print *,   &
               "Max. iteration count reached in eta_dot_l2_prime_diagnosis."
               print *, "Delta_eta = ", Delta_eta
               ! temporary lines
               print *
               print *, "i =", i, "   j =", j, "   k =", k
               print *, "In previous time step:"
               fpt2 = (phi_l2(i,j,k,n4)-phis(i,j))/(grav*h_1(i,j))
               print *, "   sigma_l2 =", fpt2
               print *, "   th_l2 =", th_l2(i,j,k,n4)
               fpt1 = (c1-fpt2)**r_fac
               F_star_l2 = theta_min*fpt1 + dth_dsigma_min*(c1-fpt2) * &
                               (c1-r_fac_param1*fpt1) +                &
                                  th_l2(i,j,k,n4) * (c1-fpt1)
               print *, "   F(th_l2,sigma_l2) =", F_star_l2
               print *, "   F_set =", F_set(i,j,k)
               print *, "   eta_l2 =", eta_l2(k)
               print *, "   eta_dot_l2_regrid =", eta_dot_l2_regrid(i,j,k)
               print *
               print *, "PROGRAM STOPPING "
               print *
               stop
            end if

   
         end do   !  while (Delta_eta .ge. converge)

      end do  !  i = 1, im
   end do   ! j = 1, jm

end do  !  k = 1, nlm-1


end subroutine vert_mass_flx_iterative_solver
!***************************************************************************





end module eta_dot_diagnosis

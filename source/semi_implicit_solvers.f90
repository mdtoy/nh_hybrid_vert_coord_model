module semi_implicit_solvers

!-----------------------------------------------------------------------
! PURPOSE: Time-steps variables using semi-implicit time differencing
!-----------------------------------------------------------------------

use kinds
use model_parameters
use physical_parameters
use eta_coordinate
use prognostics

implicit none
save


real (kind = dbl_kind), parameter ::                                 &
           beta_si =  0.55_dbl_kind,  & ! implicit weighting factor (0.5 for trapezoidal)
           alpha_si = 0.45_dbl_kind     ! explicit weighting factor (0.5 for trapezoidal)
                             ! for semi_implicit_solve_m_w_phi




contains



!======================================================================
! BEGINNING OF SEMI_IMPLICIT_SOLVE_M_W_PHI
!======================================================================

subroutine semi_implicit_solve_m_w_phi ( phi, exp_horiz_phi_advec,   &
                  eta_dot_sigma_coeff, moist_corr1, moist_corr2_l2 )

!---------------------------------------------------------------------------
! PURPOSE:
!   Computes final n+1 time level m and w_l2 variables using semi-
!   implicit (trapezoidal) time differencing to handle vertically
!   propagating acoustic waves.  The vertical pressure gradient
!   term and vertical mass flux divergence terms are integrated
!   in this manner.
!---------------------------------------------------------------------------

implicit none

!---------------------------------------------------------------------------
! INTENT IN
!---------------------------------------------------------------------------
real (kind=dbl_kind), dimension(im,jm,nlm), intent(in) ::            &
          phi            ! layer center geopotential (J/kg)

real (kind=dbl_kind), dimension(im,jm,1:nlm-1), intent(in) ::        &
          exp_horiz_phi_advec, &  ! explicitly weighted phi advection
          eta_dot_sigma_coeff     ! coefficient for calculating
                                  ! eta_dot_l2_sigma

real (kind=dbl_kind), dimension(im,jm,nlm), intent(in) ::            &
          moist_corr1     ! moisture correction for equation of
                          ! state at layer centers (-)

real (kind=dbl_kind), dimension(im,jm,0:nlm), intent(in) ::          &
          moist_corr2_l2  ! moisture correction for pressure gradient
                          ! force at layer edges (-)

!---------------------------------------------------------------------------
! LOCAL
!---------------------------------------------------------------------------
integer :: i,j,k, ii

integer, parameter :: jacdim = 3*nlm - 2, &    ! Jacobian matrix dimension
                      jacdim_red = 2*nlm - 1   ! Reduced Jacob. matrix dim.

real (kind = dbl_kind), dimension(nlm-1) ::                          &
           dz_deta_l2_n,      &   ! dz_deta_l2 at time step n
           G_phi_l2, G_w_l2,  &
           phi_l2_iter, w_l2_iter, m_l2_iter,  &
           detadot_dwl2, eta_dot_l2_sig_iter,  &
           u_elem, r_elem, b_elem, c_elem, d_elem, g_elem,   &
           a1_red_elem, F2_red_vec,  &
           F1_vec, F2_vec         ! vectors of residuals 1 & 2

real (kind = dbl_kind), dimension(0:nlm) ::                          &
           th_l2_np1,  &    ! potential temperature at time step n+1
           moist_corr2_l2_n ! moisture correction for pressure gradient
                            ! force at layer edges at time step n (-)

real (kind = dbl_kind), dimension(nlm) ::                            &
           G_m, th_np1,          &
           m_iter, P_exnr_iter,  &
           a3_elem,              &
           F3_vec,  &      ! vector of residuals 3
           moist_corr1_n   ! moisture correction for equation of
                           ! state at layer centers (-)

real (kind = dbl_kind), dimension(2:nlm-1) ::                        &
           t_elem, q_elem

real (kind = dbl_kind), dimension(nlm-2) ::                          &
           s_elem, p_elem

real (kind = dbl_kind), dimension(2:nlm) ::                          &
           e_elem, h_elem

real (kind = dbl_kind) ::                                            &
           invp0gamma_mns1,    &
           a1_elem, a2_elem, inv_a1_elem,  &
           rho_iter, pl_iter,  &
           fpt1, fpt2, fpt3, fpt4, fpt5, fpt6   ! working variables

real (kind = dbl_kind), dimension(jacdim_red) ::                     &
           AA_elem, FF_vec, XX_vec, YY_vec

real (kind = dbl_kind), dimension(jacdim_red-1) ::                   &
           CC_elem

real (kind = dbl_kind), dimension(jacdim_red-2) ::                   &
           EE_elem

real (kind = dbl_kind), dimension(2:jacdim_red) ::                   &
           BB_elem

real (kind = dbl_kind), dimension(3:jacdim_red) ::                   &
           DD_elem

real (kind = dbl_kind), dimension(2:jacdim_red,2) ::                 &
           LL_elem

real (kind = dbl_kind), dimension(1:jacdim_red,3) ::                 &
           UU_elem

real (kind = dbl_kind), dimension(jacdim) :: step_change



integer :: iter
integer, parameter :: max_iter = 100

real (kind = dbl_kind), parameter ::                                 &
           converge_cond = 1.0E-09

real (kind = dbl_kind) ::                                            &
           max_step_change

! Temporary line
integer :: max_number_iter





! Assign parameter
invp0gamma_mns1 = c1/(p0_ref**gamma_mns1)


! Assign constant elements of Jacobian matrix (do this once)
a1_elem = invdt
a2_elem = invdt
inv_a1_elem = c1/a1_elem



! Temporary line
max_number_iter = 0

do j = 1, jm
   do i = 1, im


      ! Calculate invariant functions and values and initialize
      ! variables and form invariant elements of Jacobian
      ! matrix
      th_l2_np1(0) = th_l2(i,j,0,n3)
      moist_corr2_l2_n(0) = moist_corr2_l2(i,j,0)
      do k = 1, nlm-1
         dz_deta_l2_n(k) = inv_d_eta_l2(k)*invgrav*                  &
                           (phi(i,j,k+1)-phi(i,j,k))
         G_phi_l2(k) = invdt*phi_l2(i,j,k,n3)
         G_w_l2(k) = invdt*w_l2(i,j,k,n3)
         G_m(k) = invdt*m(i,j,k,n3)
         th_l2_np1(k) = th_l2(i,j,k,n3)
         th_np1(k) = p5*(th_l2_np1(k-1)+th_l2_np1(k))
         phi_l2_iter(k) = phi_l2(i,j,k,n4)
         m_iter(k) = m(i,j,k,n4)
         w_l2_iter(k) = w_l2(i,j,k,n4)
         detadot_dwl2(k) = grav*eta_dot_sigma_coeff(i,j,k)
         u_elem(k) = beta_si*grav*(detadot_dwl2(k)*dz_deta_l2_n(k)-  &
                        c1)
         moist_corr2_l2_n(k) = moist_corr2_l2(i,j,k)
         moist_corr1_n(k) = moist_corr1(i,j,k)
      end do
      G_m(nlm) = invdt*m(i,j,nlm,n3)
      th_l2_np1(nlm) = th_l2(i,j,nlm,n3)
      th_np1(nlm) = p5*(th_l2_np1(nlm-1)+th_l2_np1(nlm))
      m_iter(nlm) = m(i,j,nlm,n4)
      moist_corr2_l2_n(nlm) = moist_corr2_l2(i,j,nlm)
      moist_corr1_n(nlm) = moist_corr1(i,j,nlm)


      !
      ! Perform iterations
      !
      iter = 0
      max_step_change = 1.E+10_dbl_kind
      
      do while ( max_step_change.ge.converge_cond )

      
         ! Initialize max_step_change
         max_step_change = c0


         !
         ! Calculate diagnostic variables and vectors of residuals
         ! Fn_vec
         !
         
         rho_iter = m_iter(1)*d_eta(1)*grav/(phi_l2_iter(1) -        &
                       phis(i,j))
         pl_iter = invp0gamma_mns1 * ( rho_iter*gas_const_R*         &
                      moist_corr1_n(1)*th_np1(1) ) ** gamma
         P_exnr_iter(1) = spec_heat_cp*(pl_iter*inv_p0_ref)**kappa
         m_l2_iter(1) = p5*inv_d_eta_l2(1)*( d_eta(1)*m_iter(1) +    &
                           d_eta(2)*m_iter(2) )
         eta_dot_l2_sig_iter(1) = eta_dot_sigma_coeff(i,j,1)*        &
            ( w_l2_iter(1)*grav + exp_horiz_phi_advec(i,j,1) )
         do k = 2, nlm-1
            rho_iter = m_iter(k)*d_eta(k)*grav/(phi_l2_iter(k) -     &
                          phi_l2_iter(k-1))
            pl_iter = invp0gamma_mns1 * ( rho_iter*gas_const_R*      &
                         moist_corr1_n(k)*th_np1(k) ) ** gamma
            P_exnr_iter(k) = spec_heat_cp*(pl_iter*inv_p0_ref)**kappa
            m_l2_iter(k) = p5*inv_d_eta_l2(k)*( d_eta(k)*m_iter(k) + &
                              d_eta(k+1)*m_iter(k+1) )
            eta_dot_l2_sig_iter(k) = eta_dot_sigma_coeff(i,j,k)*     &
               ( w_l2_iter(k)*grav + exp_horiz_phi_advec(i,j,k) )
         end do
         rho_iter = m_iter(nlm)*d_eta(nlm)*grav/(grav*z_top -        &
                       phi_l2_iter(nlm-1))
         pl_iter = invp0gamma_mns1 * ( rho_iter*gas_const_R*         &
                      moist_corr1_n(nlm)*th_np1(nlm) ) ** gamma
         P_exnr_iter(nlm) = spec_heat_cp*(pl_iter*inv_p0_ref)**kappa


         ! Note this is the "negative" residual
         F1_vec(1) = -invdt*phi_l2_iter(1) - beta_si*grav*           &
            (eta_dot_l2_sig_iter(1)*dz_deta_l2_n(1)-w_l2_iter(1)) +  &
            G_phi_l2(1)
         F2_vec(1) = -invdt*w_l2_iter(1) - beta_si*                  &
                      moist_corr2_l2_n(1)*th_l2_np1(1)*              &
                      2*grav*(P_exnr_iter(2)-P_exnr_iter(1))/        &
                       (phi_l2_iter(2)-phis(i,j)) + G_w_l2(1)
         F3_vec(1) = -invdt*m_iter(1) - beta_si*inv_d_eta(1)*        &
                        m_l2_iter(1)*eta_dot_l2_sig_iter(1) + G_m(1)
         do k = 2, nlm-2
            F1_vec(k) = -invdt*phi_l2_iter(k) - beta_si*grav*        &
             (eta_dot_l2_sig_iter(k)*dz_deta_l2_n(k)-w_l2_iter(k)) + &
              G_phi_l2(k)
            F2_vec(k) = -invdt*w_l2_iter(k)-beta_si*                 &
                   moist_corr2_l2_n(k)*th_l2_np1(k)*                 &
                   2*grav*(P_exnr_iter(k+1)-P_exnr_iter(k))/         &
                     (phi_l2_iter(k+1)-phi_l2_iter(k-1)) + G_w_l2(k)
            F3_vec(k) = -invdt*m_iter(k) - beta_si*inv_d_eta(k)*     &
               (m_l2_iter(k)*eta_dot_l2_sig_iter(k)-m_l2_iter(k-1)*  &
                eta_dot_l2_sig_iter(k-1)) + G_m(k)
         end do
         F1_vec(nlm-1) = -invdt*phi_l2_iter(nlm-1) - beta_si*grav*   &
             (eta_dot_l2_sig_iter(nlm-1)*dz_deta_l2_n(nlm-1)-        &
                 w_l2_iter(nlm-1)) + G_phi_l2(nlm-1)
         F2_vec(nlm-1) = -invdt*w_l2_iter(nlm-1) - beta_si*          &
                  moist_corr2_l2_n(nlm-1)*th_l2_np1(nlm-1)*2*grav*   &
                  (P_exnr_iter(nlm)-P_exnr_iter(nlm-1))/             &
                     (grav*z_top-phi_l2_iter(nlm-2)) + G_w_l2(nlm-1)
         F3_vec(nlm-1) = -invdt*m_iter(nlm-1) - beta_si*             &
                  inv_d_eta(nlm-1)*(m_l2_iter(nlm-1)*                &
                      eta_dot_l2_sig_iter(nlm-1)-m_l2_iter(nlm-2)*   &
                         eta_dot_l2_sig_iter(nlm-2)) + G_m(nlm-1)
         F3_vec(nlm) = -invdt*m_iter(nlm) + beta_si*inv_d_eta(nlm)*  &
              m_l2_iter(nlm-1)*eta_dot_l2_sig_iter(nlm-1) + G_m(nlm)



         !
         ! Calculate Jacobian matrix
         !

         ! D[F2]
         
         fpt1 = phi_l2_iter(2) - phis(i,j)
         fpt2 = c1 / fpt1
         fpt3 = c1 / (phi_l2_iter(2) - phi_l2_iter(1))
         fpt4 = c1 / (phi_l2_iter(1) - phis(i,j))
         fpt5 = 2*grav*beta_si*moist_corr2_l2_n(1)*th_l2_np1(1)
         s_elem(1) = fpt5*( -gamma_mns1*P_exnr_iter(2)*fpt1*fpt3 -         &
              P_exnr_iter(2) + P_exnr_iter(1) ) * fpt2**2
         r_elem(1) = fpt5*gamma_mns1*                                      &
            ( P_exnr_iter(2)*fpt3 + P_exnr_iter(1)*fpt4) * fpt2
         fpt6 = fpt5*gamma_mns1*fpt2
         b_elem(1) = - fpt6*P_exnr_iter(1)/m_iter(1)
         c_elem(1) = fpt6*P_exnr_iter(2)/m_iter(2)
         do k = 2, nlm-2
            fpt1 = phi_l2_iter(k+1) - phi_l2_iter(k-1)
            fpt2 = c1 / fpt1
            fpt3 = c1 / (phi_l2_iter(k+1) - phi_l2_iter(k))
            fpt4 = c1 / (phi_l2_iter(k) - phi_l2_iter(k-1))
            fpt5 = 2*grav*beta_si*moist_corr2_l2_n(k)*th_l2_np1(k)
            s_elem(k) = fpt5*( -gamma_mns1*P_exnr_iter(k+1)*fpt1*fpt3 -    &
                 P_exnr_iter(k+1) + P_exnr_iter(k) ) * fpt2**2
            r_elem(k) = fpt5*gamma_mns1*                                   &
               ( P_exnr_iter(k+1)*fpt3 + P_exnr_iter(k)*fpt4 ) * fpt2
            t_elem(k) = fpt5*( -gamma_mns1*P_exnr_iter(k)*fpt1*fpt4 +      &
                 P_exnr_iter(k+1) - P_exnr_iter(k) ) * fpt2**2
            fpt6 = fpt5*gamma_mns1*fpt2
            b_elem(k) = - fpt6*P_exnr_iter(k)/m_iter(k)
            c_elem(k) = fpt6*P_exnr_iter(k+1)/m_iter(k+1)
         end do
         fpt1 = grav*z_top - phi_l2_iter(nlm-2)
         fpt2 = c1 / fpt1
         fpt3 = c1 / (grav*z_top - phi_l2_iter(nlm-1))
         fpt4 = c1 / (phi_l2_iter(nlm-1) - phi_l2_iter(nlm-2))
         fpt5 = 2*grav*beta_si*moist_corr2_l2_n(nlm-1)*th_l2_np1(nlm-1)
         r_elem(nlm-1) = fpt5*gamma_mns1*                                  &
            ( P_exnr_iter(nlm)*fpt3 + P_exnr_iter(nlm-1)*fpt4 ) * fpt2
         t_elem(nlm-1) = fpt5*( -gamma_mns1*P_exnr_iter(nlm-1)*fpt1*fpt4 + &
              P_exnr_iter(nlm) - P_exnr_iter(nlm-1) ) * fpt2**2
         fpt6 = fpt5*gamma_mns1*fpt2
         b_elem(nlm-1) = - fpt6*P_exnr_iter(nlm-1)/m_iter(nlm-1)
         c_elem(nlm-1) = fpt6*P_exnr_iter(nlm)/m_iter(nlm)

         ! D[F3]

         a3_elem(1) = invdt + p5*beta_si*                                  &
                         inv_d_eta_l2(1)*eta_dot_l2_sig_iter(1)
         d_elem(1) = beta_si*inv_d_eta(1)*m_l2_iter(1)*detadot_dwl2(1)
         g_elem(1) = p5*beta_si*d_eta(2)*eta_dot_l2_sig_iter(1)*           &
                        inv_d_eta(1)*inv_d_eta_l2(1)

         do k = 2, nlm-1
            a3_elem(k) = invdt + p5*beta_si*                               &
                         ( inv_d_eta_l2(k)*eta_dot_l2_sig_iter(k) -        &
                           inv_d_eta_l2(k-1)*eta_dot_l2_sig_iter(k-1) )
            d_elem(k) = beta_si*inv_d_eta(k)*m_l2_iter(k)*detadot_dwl2(k)
            e_elem(k) = - beta_si*inv_d_eta(k)*m_l2_iter(k-1)*             &
                             detadot_dwl2(k-1)
            g_elem(k) = p5*beta_si*d_eta(k+1)*eta_dot_l2_sig_iter(k)*      &
                            inv_d_eta(k)*inv_d_eta_l2(k)
            h_elem(k) = - p5*beta_si*d_eta(k-1)*eta_dot_l2_sig_iter(k-1)*  &
                            inv_d_eta(k)*inv_d_eta_l2(k-1)
         end do
         a3_elem(nlm) = invdt - p5*beta_si*                                &
                         inv_d_eta_l2(nlm-1)*eta_dot_l2_sig_iter(nlm-1)
         e_elem(nlm) = - beta_si*inv_d_eta(nlm)*m_l2_iter(nlm-1)*          &
                             detadot_dwl2(nlm-1)
         h_elem(nlm) = - p5*beta_si*d_eta(nlm-1)*                          &
            eta_dot_l2_sig_iter(nlm-1)*inv_d_eta(nlm)*inv_d_eta_l2(nlm-1)


         !
         ! Form "reduced" Jacobian and F_vec
         !
         
         a1_red_elem(1) = a2_elem - inv_a1_elem*r_elem(1)*u_elem(1)
         p_elem(1) = - inv_a1_elem*s_elem(1)*u_elem(2)
         do k = 2, nlm-2
            a1_red_elem(k) = a2_elem - inv_a1_elem*r_elem(k)*u_elem(k)
            p_elem(k) = - inv_a1_elem*s_elem(k)*u_elem(k+1)
            q_elem(k) = - inv_a1_elem*t_elem(k)*u_elem(k-1)
         end do
         a1_red_elem(nlm-1) = a2_elem - inv_a1_elem*r_elem(nlm-1)*      &
                                 u_elem(nlm-1)
         q_elem(nlm-1) = - inv_a1_elem*t_elem(nlm-1)*u_elem(nlm-2)
         
         F2_red_vec(1) = F2_vec(1) - inv_a1_elem*                       &
                         ( r_elem(1)*F1_vec(1) + s_elem(1)*F1_vec(2) )
         do k = 2, nlm-2
            F2_red_vec(k) = F2_vec(k) - inv_a1_elem*                    &
               ( t_elem(k)*F1_vec(k-1) + r_elem(k)*F1_vec(k) +          &
                 s_elem(k)*F1_vec(k+1) )
         end do
         F2_red_vec(nlm-1) = F2_vec(nlm-1) - inv_a1_elem*               &
           ( t_elem(nlm-1)*F1_vec(nlm-2) + r_elem(nlm-1)*F1_vec(nlm-1) )


         !
         ! Rearrange rows and columns to give pentadiagonal Jacobian matrix
         ! and rearranged F_red_vec
         !
         
         AA_elem(1) = a3_elem(1)
         AA_elem(2) = a1_red_elem(1)
         CC_elem(1) = d_elem(1)
         CC_elem(2) = c_elem(1)
         EE_elem(1) = g_elem(1)
         EE_elem(2) = p_elem(1)
         BB_elem(2) = b_elem(1)
         FF_vec(1) = F3_vec(1)
         FF_vec(2) = F2_red_vec(1)
         do k = 2, nlm-2
            AA_elem(1+2*(k-1)) = a3_elem(k)
            AA_elem(2*k) = a1_red_elem(k)
            CC_elem(1+2*(k-1)) = d_elem(k)
            CC_elem(2*k) = c_elem(k)
            EE_elem(1+2*(k-1)) = g_elem(k)
            EE_elem(2*k) = p_elem(k)
            BB_elem(1+2*(k-1)) = e_elem(k)
            BB_elem(2*k) = b_elem(k)
            DD_elem(1+2*(k-1)) = h_elem(k)
            DD_elem(2*k) = q_elem(k)
            FF_vec(1+2*(k-1)) = F3_vec(k)
            FF_vec(2*k) = F2_red_vec(k)
         end do
         AA_elem(jacdim_red-2) = a3_elem(nlm-1)
         AA_elem(jacdim_red-1) = a1_red_elem(nlm-1)
         AA_elem(jacdim_red) = a3_elem(nlm)
         CC_elem(jacdim_red-2) = d_elem(nlm-1)
         CC_elem(jacdim_red-1) = c_elem(nlm-1)
         EE_elem(jacdim_red-2) = g_elem(nlm-1)
         BB_elem(jacdim_red-2) = e_elem(nlm-1)
         BB_elem(jacdim_red-1) = b_elem(nlm-1)
         BB_elem(jacdim_red) = e_elem(nlm)
         DD_elem(jacdim_red-2) = h_elem(nlm-1)
         DD_elem(jacdim_red-1) = q_elem(nlm-1)
         DD_elem(jacdim_red) = h_elem(nlm)
         FF_vec(jacdim_red-2) = F3_vec(nlm-1)
         FF_vec(jacdim_red-1) = F2_red_vec(nlm-1)
         FF_vec(jacdim_red) = F3_vec(nlm)


         !
         ! Solve linear system
         !
         
         ! Perform LU decomposition
         
         UU_elem(1,1) = AA_elem(1)
         UU_elem(1,2) = CC_elem(1)
         UU_elem(1,3) = EE_elem(1)

         LL_elem(2,1) = BB_elem(2)/UU_elem(1,1)
         UU_elem(2,1) = AA_elem(2) - LL_elem(2,1)*UU_elem(1,2)
         UU_elem(2,2) = CC_elem(2) - LL_elem(2,1)*UU_elem(1,3)
         UU_elem(2,3) = EE_elem(2)
         
         do ii = 3, jacdim_red-2
            LL_elem(ii,1) = DD_elem(ii)/UU_elem(ii-2,1)
            LL_elem(ii,2) = ( BB_elem(ii) -                          &
                 LL_elem(ii,1)*UU_elem(ii-2,2) ) / UU_elem(ii-1,1)
            UU_elem(ii,1) = AA_elem(ii) - LL_elem(ii,1)*             &
                 UU_elem(ii-2,3) - LL_elem(ii,2)*UU_elem(ii-1,2)
            UU_elem(ii,2) = CC_elem(ii) - LL_elem(ii,2)*             &
                 UU_elem(ii-1,3)
            UU_elem(ii,3) = EE_elem(ii)
         end do
         
         LL_elem(jacdim_red-1,1) = DD_elem(jacdim_red-1)/            &
              UU_elem(jacdim_red-3,1)
         LL_elem(jacdim_red-1,2) = ( BB_elem(jacdim_red-1) -         &
              LL_elem(jacdim_red-1,1)*UU_elem(jacdim_red-3,2) ) /    &
                 UU_elem(jacdim_red-2,1)
         UU_elem(jacdim_red-1,1) = AA_elem(jacdim_red-1) -           &
              LL_elem(jacdim_red-1,1)*UU_elem(jacdim_red-3,3) -      &
                 LL_elem(jacdim_red-1,2)*UU_elem(jacdim_red-2,2)
         UU_elem(jacdim_red-1,2) = CC_elem(jacdim_red-1) -           &
              LL_elem(jacdim_red-1,2)*UU_elem(jacdim_red-2,3)

         LL_elem(jacdim_red,1) = DD_elem(jacdim_red)/                &
              UU_elem(jacdim_red-2,1)
         LL_elem(jacdim_red,2) = ( BB_elem(jacdim_red) -             &
              LL_elem(jacdim_red,1)*UU_elem(jacdim_red-2,2) ) /      &
                 UU_elem(jacdim_red-1,1)
         UU_elem(jacdim_red,1) = AA_elem(jacdim_red) -               &
              LL_elem(jacdim_red,1)*UU_elem(jacdim_red-2,3) -        &
                 LL_elem(jacdim_red,2)*UU_elem(jacdim_red-1,2)


         ! Forward elimination
         XX_vec(1) = FF_vec(1)
         XX_vec(2) = FF_vec(2) - LL_elem(2,1)*XX_vec(1)
         do ii = 3, jacdim_red
            XX_vec(ii) = FF_vec(ii) - LL_elem(ii,1)*XX_vec(ii-2) -   &
                            LL_elem(ii,2)*XX_vec(ii-1)
         end do
         
         ! Back substitution
         YY_vec(jacdim_red) = XX_vec(jacdim_red) /                   &
                                 UU_elem(jacdim_red,1)
         YY_vec(jacdim_red-1) = ( XX_vec(jacdim_red-1) -             &
                 UU_elem(jacdim_red-1,2)*YY_vec(jacdim_red) ) /      &
                    UU_elem(jacdim_red-1,1)
         do ii = jacdim_red-2, 1, -1
            YY_vec(ii) = ( XX_vec(ii) - UU_elem(ii,2)*YY_vec(ii+1) - &
                     UU_elem(ii,3)*YY_vec(ii+2) ) / UU_elem(ii,1)
         end do


         ! Finally, calculate step_change
         do k = 1, nlm-1
            step_change(nlm-1+k) = YY_vec(2*k)  ! change in w_l2_iter
            step_change(jacdim_red-1+k) = YY_vec(1+2*(k-1))
                                                  ! change in m_iter
         end do
         step_change(jacdim) = YY_vec(jacdim_red) ! change in m_iter
         
         do k = 1, nlm-1
            step_change(k) = inv_a1_elem *                          &
                     ( F1_vec(k) - u_elem(k)*step_change(nlm-1+k) )
                         ! change in phi_l2_iter
         end do
         
         
         ! Increment iterative variables phi_l2_iter, w_l2_iter
         ! and m_iter and determine max step changes
         do k = 1, nlm-1
            phi_l2_iter(k) = phi_l2_iter(k) + step_change(k)
            w_l2_iter(k) = w_l2_iter(k) + step_change(nlm-1+k)
            m_iter(k) = m_iter(k) + step_change(jacdim_red-1+k)
            max_step_change = max ( max_step_change,                 &
               abs(step_change(k)),abs(step_change(nlm-1+k)),        &
               abs(step_change(jacdim_red-1+k)) )
         end do
         m_iter(nlm) = m_iter(nlm) + step_change(jacdim)
         max_step_change = max ( max_step_change,                    &
                                 abs(step_change(jacdim)) )




         iter = iter + 1
         if (iter .ge. max_iter) then
            print *
            print *, "Max. iteration count reached in subroutine &
                     &SEMI_IMPLICIT_SOLVE_M_W_PHI"
            print *, "max_step_change = ", max_step_change
            print *, "i =", i, "  j =", j
            print *
            print *, "PROGRAM STOPPING "
            print *
            stop
         end if
         
      end do    ! while ( max_step_change.ge.converge_cond )
      
      
      ! Finalize prognostic variables phi_l2, w_l2 and m
      do k = 1, nlm-1
         phi_l2(i,j,k,n3) = phi_l2_iter(k)
         w_l2(i,j,k,n3) = w_l2_iter(k)
         m(i,j,k,n3) = m_iter(k)
      end do
      m(i,j,nlm,n3) = m_iter(nlm)


      ! Temporary lines
      max_number_iter = max(max_number_iter, iter)

   end do    ! i = 1, im
end do       ! j = 1, jm


! Temporary line
print "(A61,I5)", "In this time step, semi-implicit solver &
                   &max_number_iter =", max_number_iter



end subroutine semi_implicit_solve_m_w_phi

!======================================================================
! END OF SEMI_IMPLICIT_SOLVE_M_W_PHI
!======================================================================




end module semi_implicit_solvers

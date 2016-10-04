module tracer_transport

!-----------------------------------------------------------------------
! PURPOSE: Calculates tracer tendencies.
!
!-----------------------------------------------------------------------

use kinds
use model_parameters
use physical_parameters
use eta_coordinate
use prognostics


implicit none
save


contains



!======================================================================
! BEGINNING OF GET_horiz_div_mvqc
!======================================================================

subroutine get_horiz_div_mvqc (F_x, F_y, mDqc_Dt)

!---------------------------------------------------------------------------
! PURPOSE:
!   Computes the horizontal divergence of flux of cloud water mixing ratio
!   at layer edges for calculation of diabatic heating.  In the
!   process, "third-order" horizontal mass fluxes  F_x and F_y are 
!   calculated.
!
!---------------------------------------------------------------------------

implicit none

!---------------------------------------------------------------------------
! INTENT OUT
!---------------------------------------------------------------------------
real (kind=dbl_kind), dimension(im,jm,nlm), intent(out) ::           &
           F_x,   &  ! "3rd-order" mass flux in x-direction
                     ! (kg/s/m/eta) (colocated with u-points)
           F_y       ! "3rd-order" mass flux in y-direction
                     ! (kg/s/m/eta) (colocated with v-points)

real (kind=dbl_kind), dimension(im,jm,0:nlm), intent(out) ::         &
           mDqc_Dt   ! mass (kg/m^2/eta) times material time derivative
                     ! of cloud-water mixing ratio (kg/kg/s)

!---------------------------------------------------------------------------
! LOCAL
!---------------------------------------------------------------------------
integer :: k

real (kind=dbl_kind), dimension(im,jm) ::                            &
          G_flx_contrib

real (kind=dbl_kind), dimension(im,jm) ::                            &
          F_tr_x, F_tr_y, F_horiz_l2



!---------------------------------------------------------------------------
! Calculate horizontal mass fluxes
!---------------------------------------------------------------------------

do k = 1, nlm

   ! x-contribution
   G_flx_contrib(:,:) = - inv12 * abs(u(:,:,k,n4)) *                 &
      ( 3*(m(:,:,k,n4)-m(im1(:),:,k,n4)) - m(ip1(:),:,k,n4) +        &
        m(im2(:),:,k,n4) )
   F_x(:,:,k) = inv12 * u(:,:,k,n4) *                                &
      ( 7*(m(:,:,k,n4)+m(im1(:),:,k,n4)) -                           &
        m(ip1(:),:,k,n4) - m(im2(:),:,k,n4) ) + G_flx_contrib(:,:)

   ! y-contribution
   G_flx_contrib(:,:) = - inv12 * abs(v(:,:,k,n4)) *                 &
      ( 3*(m(:,:,k,n4)-m(:,jm1(:),k,n4)) - m(:,jp1(:),k,n4) +        &
        m(:,jm2(:),k,n4) )
   F_y(:,:,k) = inv12 * v(:,:,k,n4) *                                &
      ( 7*(m(:,:,k,n4)+m(:,jm1(:),k,n4)) -                           &
        m(:,jp1(:),k,n4) - m(:,jm2(:),k,n4) ) + G_flx_contrib(:,:)

end do



!---------------------------------------------------------------------------
! Diagnose horizontal advection of qc -- cloud water mixing ratio
!---------------------------------------------------------------------------


! Horizontal advection

! Start at bottom level (surface)

k = 0

! x-contribution
F_horiz_l2(:,:) = F_x(:,:,k+1)
G_flx_contrib(:,:) = - inv12 * abs(F_horiz_l2(:,:)) *                &
   ( 3*(qc_l2(:,:,k)-qc_l2(im1(:),:,k)) - qc_l2(ip1(:),:,k) +        &
     qc_l2(im2(:),:,k) )
F_tr_x(:,:) = inv12 * F_horiz_l2(:,:) *                              &
      ( 7*(qc_l2(:,:,k)+qc_l2(im1(:),:,k)) -                         &
        qc_l2(ip1(:),:,k) - qc_l2(im2(:),:,k) ) + G_flx_contrib(:,:)

! y-contribution
F_horiz_l2(:,:) = F_y(:,:,k+1)
G_flx_contrib(:,:) = - inv12 * abs(F_horiz_l2(:,:)) *                &
   ( 3*(qc_l2(:,:,k)-qc_l2(:,jm1(:),k)) - qc_l2(:,jp1(:),k) +        &
     qc_l2(:,jm2(:),k) )
F_tr_y(:,:) = inv12 * F_horiz_l2(:,:) *                              &
      ( 7*(qc_l2(:,:,k)+qc_l2(:,jm1(:),k)) -                         &
        qc_l2(:,jp1(:),k) - qc_l2(:,jm2(:),k) ) + G_flx_contrib(:,:)

mDqc_Dt(:,:,k) =  invdx * ( F_tr_x(ip1(:),:) - F_tr_x(:,:) ) +       &
                  invdy * ( F_tr_y(:,jp1(:)) - F_tr_y(:,:) )


! Move upward

do k = 1, nlm-1

   ! x-contribution
   F_horiz_l2(:,:) = p5*inv_d_eta_l2(k)*                             &
                       (d_eta(k)*F_x(:,:,k)+d_eta(k+1)*F_x(:,:,k+1))
   G_flx_contrib(:,:) = - inv12 * abs(F_horiz_l2(:,:)) *             &
      ( 3*(qc_l2(:,:,k)-qc_l2(im1(:),:,k)) - qc_l2(ip1(:),:,k) +     &
        qc_l2(im2(:),:,k) )
   F_tr_x(:,:) = inv12 * F_horiz_l2(:,:) *                           &
      ( 7*(qc_l2(:,:,k)+qc_l2(im1(:),:,k)) -                         &
        qc_l2(ip1(:),:,k) - qc_l2(im2(:),:,k) ) + G_flx_contrib(:,:)

   ! y-contribution
   F_horiz_l2(:,:) = p5*inv_d_eta_l2(k)*                             &
                       (d_eta(k)*F_y(:,:,k)+d_eta(k+1)*F_y(:,:,k+1))
   G_flx_contrib(:,:) = - inv12 * abs(F_horiz_l2(:,:)) *             &
      ( 3*(qc_l2(:,:,k)-qc_l2(:,jm1(:),k)) - qc_l2(:,jp1(:),k) +     &
        qc_l2(:,jm2(:),k) )
   F_tr_y(:,:) = inv12 * F_horiz_l2(:,:) *                           &
      ( 7*(qc_l2(:,:,k)+qc_l2(:,jm1(:),k)) -                         &
        qc_l2(:,jp1(:),k) - qc_l2(:,jm2(:),k) ) + G_flx_contrib(:,:)

   mDqc_Dt(:,:,k) =  invdx * ( F_tr_x(ip1(:),:) - F_tr_x(:,:) ) +    &
                     invdy * ( F_tr_y(:,jp1(:)) - F_tr_y(:,:) )

end do


! Top level (model top)

k = nlm

! x-contribution
F_horiz_l2(:,:) = F_x(:,:,k)
G_flx_contrib(:,:) = - inv12 * abs(F_horiz_l2(:,:)) *                &
   ( 3*(qc_l2(:,:,k)-qc_l2(im1(:),:,k)) - qc_l2(ip1(:),:,k) +        &
     qc_l2(im2(:),:,k) )
F_tr_x(:,:) = inv12 * F_horiz_l2(:,:) *                              &
      ( 7*(qc_l2(:,:,k)+qc_l2(im1(:),:,k)) -                         &
        qc_l2(ip1(:),:,k) - qc_l2(im2(:),:,k) ) + G_flx_contrib(:,:)

! y-contribution
F_horiz_l2(:,:) = F_y(:,:,k)
G_flx_contrib(:,:) = - inv12 * abs(F_horiz_l2(:,:)) *                &
   ( 3*(qc_l2(:,:,k)-qc_l2(:,jm1(:),k)) - qc_l2(:,jp1(:),k) +        &
     qc_l2(:,jm2(:),k) )
F_tr_y(:,:) = inv12 * F_horiz_l2(:,:) *                              &
      ( 7*(qc_l2(:,:,k)+qc_l2(:,jm1(:),k)) -                         &
        qc_l2(:,jp1(:),k) - qc_l2(:,jm2(:),k) ) + G_flx_contrib(:,:)

mDqc_Dt(:,:,k) =  invdx * ( F_tr_x(ip1(:),:) - F_tr_x(:,:) ) +       &
                  invdy * ( F_tr_y(:,jp1(:)) - F_tr_y(:,:) )




! Vertical advection:  To be done later in module step



end subroutine get_horiz_div_mvqc

!======================================================================
! END OF GET_horiz_div_mvqc
!======================================================================





!======================================================================
! BEGINNING OF GET_TRACER_TENDENCIES
!======================================================================

subroutine get_tracer_tendencies (F_x, F_y, F_eta_l2_AB3,            &
                                     F_eta_l2_EF, eta_dot_l2)

!---------------------------------------------------------------------------
! PURPOSE:
!   Computes the time tendency of tracers --
!   water vapor + cloud mixing ratio and rain water mixing ratio
!   I.E., Total moisture at layer edges (upstream-weighted Wicker and
!   Skamarock (2002)) and cloud-water vertical advective tendency
!   (for latent heating) calculation (second-order centered).
!   
!---------------------------------------------------------------------------

implicit none

!---------------------------------------------------------------------------
! INTENT IN
!---------------------------------------------------------------------------
real (kind = dbl_kind), dimension(im,jm,nlm), intent(in) ::          &
           F_x,          &  ! "3rd-order" mass flux in x-direction
                            ! (kg/s/m/eta) (colocated with u-points)
           F_y              ! "3rd-order" mass flux in y-direction
                            ! (kg/s/m/eta) (colocated with v-points)

real (kind = dbl_kind), dimension(im,jm,1:nlm-1), intent(in) ::      &
           F_eta_l2_AB3, &  ! vertical mass flux (kg/s/m^2)
           F_eta_l2_EF,  &  ! (colocated with eta_dot_l2 points)
                            ! AB3 and Euler-forward, respectively
           eta_dot_l2       ! Total generalized vertical velocity (K/s)

!---------------------------------------------------------------------------
! LOCAL
!---------------------------------------------------------------------------
integer :: k

integer :: kp0, kp1     ! indices for k level (kp0 = k+0) and the k+1
                        ! level (kp1 = k+1)  (they toggle between 0 and 1)

real (kind=dbl_kind), dimension(im,jm,0:nlm) ::                      &
          div_horiz_tr1_m_flux, &  ! horiz. divergence of horiz. (mass
          div_horiz_tr2_m_flux     ! times tracer) flux (kg/m^2/eta/s)

real (kind=dbl_kind), dimension(im,jm,0:nlm) ::                      &
          div_vert_tr1_m_flux, &   ! vert. divergence of vert. mass flux
          div_vert_tr2_m_flux      ! (kg/m^2/eta/s)

real (kind=dbl_kind), dimension(im,jm) ::                            &
          G_flx_contrib, G_flx_contrib_2,   &
          F_tr1_x, F_tr1_y, F_tr2_x, F_tr2_y, F_horiz_l2

real (kind=dbl_kind), dimension(im,jm,0:1) ::                        &
          F_eta_l2    ! Vertical flux with total eta_dot_l2

real (kind=dbl_kind), dimension(im,jm) ::                            &
          F_eta_tot   ! Total vert. flux interp to layer centers

real (kind=dbl_kind), dimension(im,jm,0:1) ::                        &
          F_tr1_eta, F_tr2_eta, F_tr_eta_2

real (kind=dbl_kind), dimension(im,jm,nlm) ::                        &
          F_eta   ! vertical mass flux interpolated to layer centers






!---------------------------------------------------------------------------
! Advect qvpc -- total moisture mixing ratio
! and compute vertical advection of qc -- cloud water mixing ratio
!---------------------------------------------------------------------------


! Horizontal advection

! Start at bottom level (surface)

k = 0

! x-contribution
F_horiz_l2(:,:) = F_x(:,:,k+1)
G_flx_contrib(:,:) = - inv12 * abs(F_horiz_l2(:,:)) *                  &
   ( 3*(qvpc_l2(:,:,k)-qvpc_l2(im1(:),:,k)) - qvpc_l2(ip1(:),:,k) +    &
     qvpc_l2(im2(:),:,k) )
F_tr1_x(:,:) = inv12 * F_horiz_l2(:,:) *                               &
     ( 7*(qvpc_l2(:,:,k)+qvpc_l2(im1(:),:,k)) -                        &
       qvpc_l2(ip1(:),:,k) - qvpc_l2(im2(:),:,k) ) + G_flx_contrib(:,:)
G_flx_contrib(:,:) = - inv12 * abs(F_horiz_l2(:,:)) *                  &
   ( 3*(qr_l2(:,:,k)-qr_l2(im1(:),:,k)) - qr_l2(ip1(:),:,k) +          &
     qr_l2(im2(:),:,k) )
F_tr2_x(:,:) = inv12 * F_horiz_l2(:,:) *                               &
     ( 7*(qr_l2(:,:,k)+qr_l2(im1(:),:,k)) -                            &
       qr_l2(ip1(:),:,k) - qr_l2(im2(:),:,k) ) + G_flx_contrib(:,:)

! y-contribution
F_horiz_l2(:,:) = F_y(:,:,k+1)
G_flx_contrib(:,:) = - inv12 * abs(F_horiz_l2(:,:)) *                  &
   ( 3*(qvpc_l2(:,:,k)-qvpc_l2(:,jm1(:),k)) - qvpc_l2(:,jp1(:),k) +    &
     qvpc_l2(:,jm2(:),k) )
F_tr1_y(:,:) = inv12 * F_horiz_l2(:,:) *                               &
     ( 7*(qvpc_l2(:,:,k)+qvpc_l2(:,jm1(:),k)) -                        &
       qvpc_l2(:,jp1(:),k) - qvpc_l2(:,jm2(:),k) ) + G_flx_contrib(:,:)
G_flx_contrib(:,:) = - inv12 * abs(F_horiz_l2(:,:)) *                  &
   ( 3*(qr_l2(:,:,k)-qr_l2(:,jm1(:),k)) - qr_l2(:,jp1(:),k) +          &
     qr_l2(:,jm2(:),k) )
F_tr2_y(:,:) = inv12 * F_horiz_l2(:,:) *                               &
     ( 7*(qr_l2(:,:,k)+qr_l2(:,jm1(:),k)) -                            &
       qr_l2(:,jp1(:),k) - qr_l2(:,jm2(:),k) ) + G_flx_contrib(:,:)

div_horiz_tr1_m_flux(:,:,k) =                                          &
                 invdx * ( F_tr1_x(ip1(:),:) - F_tr1_x(:,:) ) +        &
                 invdy * ( F_tr1_y(:,jp1(:)) - F_tr1_y(:,:) )
div_horiz_tr2_m_flux(:,:,k) =                                          &
                 invdx * ( F_tr2_x(ip1(:),:) - F_tr2_x(:,:) ) +        &
                 invdy * ( F_tr2_y(:,jp1(:)) - F_tr2_y(:,:) )


! Move upward

do k = 1, nlm-1

   ! x-contribution
   F_horiz_l2(:,:) = p5*inv_d_eta_l2(k)*                             &
                       (d_eta(k)*F_x(:,:,k)+d_eta(k+1)*F_x(:,:,k+1))
   G_flx_contrib(:,:) = - inv12 * abs(F_horiz_l2(:,:)) *               &
      ( 3*(qvpc_l2(:,:,k)-qvpc_l2(im1(:),:,k)) - qvpc_l2(ip1(:),:,k) + &
        qvpc_l2(im2(:),:,k) )
   F_tr1_x(:,:) = inv12 * F_horiz_l2(:,:) *                            &
     ( 7*(qvpc_l2(:,:,k)+qvpc_l2(im1(:),:,k)) -                        &
       qvpc_l2(ip1(:),:,k) - qvpc_l2(im2(:),:,k) ) + G_flx_contrib(:,:)
   G_flx_contrib(:,:) = - inv12 * abs(F_horiz_l2(:,:)) *               &
      ( 3*(qr_l2(:,:,k)-qr_l2(im1(:),:,k)) - qr_l2(ip1(:),:,k) +       &
        qr_l2(im2(:),:,k) )
   F_tr2_x(:,:) = inv12 * F_horiz_l2(:,:) *                            &
     ( 7*(qr_l2(:,:,k)+qr_l2(im1(:),:,k)) -                            &
       qr_l2(ip1(:),:,k) - qr_l2(im2(:),:,k) ) + G_flx_contrib(:,:)

   ! y-contribution
   F_horiz_l2(:,:) = p5*inv_d_eta_l2(k)*                             &
                       (d_eta(k)*F_y(:,:,k)+d_eta(k+1)*F_y(:,:,k+1))
   G_flx_contrib(:,:) = - inv12 * abs(F_horiz_l2(:,:)) *               &
      ( 3*(qvpc_l2(:,:,k)-qvpc_l2(:,jm1(:),k)) - qvpc_l2(:,jp1(:),k) + &
        qvpc_l2(:,jm2(:),k) )
   F_tr1_y(:,:) = inv12 * F_horiz_l2(:,:) *                            &
     ( 7*(qvpc_l2(:,:,k)+qvpc_l2(:,jm1(:),k)) -                        &
       qvpc_l2(:,jp1(:),k) - qvpc_l2(:,jm2(:),k) ) + G_flx_contrib(:,:)
   G_flx_contrib(:,:) = - inv12 * abs(F_horiz_l2(:,:)) *               &
      ( 3*(qr_l2(:,:,k)-qr_l2(:,jm1(:),k)) - qr_l2(:,jp1(:),k) +       &
        qr_l2(:,jm2(:),k) )
   F_tr2_y(:,:) = inv12 * F_horiz_l2(:,:) *                            &
     ( 7*(qr_l2(:,:,k)+qr_l2(:,jm1(:),k)) -                            &
       qr_l2(:,jp1(:),k) - qr_l2(:,jm2(:),k) ) + G_flx_contrib(:,:)

   div_horiz_tr1_m_flux(:,:,k) =                                       &
                    invdx * ( F_tr1_x(ip1(:),:) - F_tr1_x(:,:) ) +     &
                    invdy * ( F_tr1_y(:,jp1(:)) - F_tr1_y(:,:) )
   div_horiz_tr2_m_flux(:,:,k) =                                       &
                    invdx * ( F_tr2_x(ip1(:),:) - F_tr2_x(:,:) ) +     &
                    invdy * ( F_tr2_y(:,jp1(:)) - F_tr2_y(:,:) )

end do


! Top level (model top)

k = nlm

! x-contribution
F_horiz_l2(:,:) = F_x(:,:,k)
G_flx_contrib(:,:) = - inv12 * abs(F_horiz_l2(:,:)) *                  &
   ( 3*(qvpc_l2(:,:,k)-qvpc_l2(im1(:),:,k)) - qvpc_l2(ip1(:),:,k) +    &
     qvpc_l2(im2(:),:,k) )
F_tr1_x(:,:) = inv12 * F_horiz_l2(:,:) *                               &
     ( 7*(qvpc_l2(:,:,k)+qvpc_l2(im1(:),:,k)) -                        &
       qvpc_l2(ip1(:),:,k) - qvpc_l2(im2(:),:,k) ) + G_flx_contrib(:,:)
G_flx_contrib(:,:) = - inv12 * abs(F_horiz_l2(:,:)) *                  &
   ( 3*(qr_l2(:,:,k)-qr_l2(im1(:),:,k)) - qr_l2(ip1(:),:,k) +          &
     qr_l2(im2(:),:,k) )
F_tr2_x(:,:) = inv12 * F_horiz_l2(:,:) *                               &
     ( 7*(qr_l2(:,:,k)+qr_l2(im1(:),:,k)) -                            &
       qr_l2(ip1(:),:,k) - qr_l2(im2(:),:,k) ) + G_flx_contrib(:,:)

! y-contribution
F_horiz_l2(:,:) = F_y(:,:,k)
G_flx_contrib(:,:) = - inv12 * abs(F_horiz_l2(:,:)) *                  &
   ( 3*(qvpc_l2(:,:,k)-qvpc_l2(:,jm1(:),k)) - qvpc_l2(:,jp1(:),k) +    &
     qvpc_l2(:,jm2(:),k) )
F_tr1_y(:,:) = inv12 * F_horiz_l2(:,:) *                               &
     ( 7*(qvpc_l2(:,:,k)+qvpc_l2(:,jm1(:),k)) -                        &
       qvpc_l2(:,jp1(:),k) - qvpc_l2(:,jm2(:),k) ) + G_flx_contrib(:,:)
G_flx_contrib(:,:) = - inv12 * abs(F_horiz_l2(:,:)) *                  &
   ( 3*(qr_l2(:,:,k)-qr_l2(:,jm1(:),k)) - qr_l2(:,jp1(:),k) +          &
     qr_l2(:,jm2(:),k) )
F_tr2_y(:,:) = inv12 * F_horiz_l2(:,:) *                               &
     ( 7*(qr_l2(:,:,k)+qr_l2(:,jm1(:),k)) -                            &
       qr_l2(:,jp1(:),k) - qr_l2(:,jm2(:),k) ) + G_flx_contrib(:,:)

div_horiz_tr1_m_flux(:,:,k) =                                          &
                 invdx * ( F_tr1_x(ip1(:),:) - F_tr1_x(:,:) ) +        &
                 invdy * ( F_tr1_y(:,jp1(:)) - F_tr1_y(:,:) )
div_horiz_tr2_m_flux(:,:,k) =                                          &
                 invdx * ( F_tr2_x(ip1(:),:) - F_tr2_x(:,:) ) +        &
                 invdy * ( F_tr2_y(:,jp1(:)) - F_tr2_y(:,:) )





! Vertical advection -- AB3 fluxes


! Calculate vertical mass flux interpolated to layer centers -- AB3 fluxes
F_eta(:,:,1) = p5*F_eta_l2_AB3(:,:,1)
do k = 2, nlm-1
   F_eta(:,:,k) = p5*(F_eta_l2_AB3(:,:,k)+F_eta_l2_AB3(:,:,k-1))
end do
F_eta(:,:,nlm) = p5*F_eta_l2_AB3(:,:,nlm-1)


! Start at bottom level (surface)

k = 0
kp0 = mod(k,2)
kp1 = mod(k+1,2)

! G_flx_contrib(:,:) = c0
! G_flx_contrib_2(:,:) = c0

F_tr1_eta(:,:,kp1) = p5 * F_eta(:,:,k+1) *                               &
         (qvpc_l2(:,:,k+1)+qvpc_l2(:,:,k))
F_tr_eta_2(:,:,kp1) = p5 * F_eta(:,:,k+1) *                              &
         (qc_l2(:,:,k+1)+qc_l2(:,:,k))

div_vert_tr1_m_flux(:,:,k) = inv_d_eta_l2(k) * F_tr1_eta(:,:,kp1)
qc_vert_adv_ab3(:,:,k,n3_f) = inv_d_eta_l2(k) * F_tr_eta_2(:,:,kp1)


! Move on up

do k = 1, nlm-2

   kp0 = mod(k,2)
   kp1 = mod(k+1,2)

   G_flx_contrib(:,:) = - inv12 * abs(F_eta(:,:,k+1)) *                 &
      ( 3*(qvpc_l2(:,:,k+1)-qvpc_l2(:,:,k))-                            &
           qvpc_l2(:,:,k+2)+qvpc_l2(:,:,k-1) )
   G_flx_contrib_2(:,:) = - inv12 * abs(F_eta(:,:,k+1)) *               &
      ( 3*(qc_l2(:,:,k+1)-qc_l2(:,:,k))-qc_l2(:,:,k+2)+qc_l2(:,:,k-1) )

   F_tr1_eta(:,:,kp1) = inv12 * F_eta(:,:,k+1) *                        &
       ( 7*(qvpc_l2(:,:,k+1)+qvpc_l2(:,:,k))-                           &
            qvpc_l2(:,:,k+2)-qvpc_l2(:,:,k-1) ) + G_flx_contrib(:,:)
   F_tr_eta_2(:,:,kp1) = inv12 * F_eta(:,:,k+1) *                       &
       ( 7*(qc_l2(:,:,k+1)+qc_l2(:,:,k))-                               &
            qc_l2(:,:,k+2)-qc_l2(:,:,k-1) ) + G_flx_contrib_2(:,:)


   div_vert_tr1_m_flux(:,:,k) = inv_d_eta_l2(k) *                       &
                            ( F_tr1_eta(:,:,kp1) - F_tr1_eta(:,:,kp0) )
   qc_vert_adv_ab3(:,:,k,n3_f) = inv_d_eta_l2(k) *                      &
                          ( F_tr_eta_2(:,:,kp1) - F_tr_eta_2(:,:,kp0) )
end do


! Almost there!

k = nlm-1
kp0 = mod(k,2)
kp1 = mod(k+1,2)

! G_flx_contrib(:,:) = c0
! G_flx_contrib_2(:,:) = c0

F_tr1_eta(:,:,kp1) = p5 * F_eta(:,:,k+1) *                              &
         (qvpc_l2(:,:,k+1)+qvpc_l2(:,:,k))
F_tr_eta_2(:,:,kp1) = p5 * F_eta(:,:,k+1) *                             &
         (qc_l2(:,:,k+1)+qc_l2(:,:,k))

div_vert_tr1_m_flux(:,:,k) = inv_d_eta_l2(k) *                          &
                           ( F_tr1_eta(:,:,kp1) - F_tr1_eta(:,:,kp0) )
qc_vert_adv_ab3(:,:,k,n3_f) = inv_d_eta_l2(k) *                         &
                       ( F_tr_eta_2(:,:,kp1) - F_tr_eta_2(:,:,kp0) )


! Top layer

k = nlm
kp0 = mod(k,2)
kp1 = mod(k+1,2)

div_vert_tr1_m_flux(:,:,k) = - inv_d_eta_l2(k) * F_tr1_eta(:,:,kp0)
qc_vert_adv_ab3(:,:,k,n3_f) = - inv_d_eta_l2(k) * F_tr_eta_2(:,:,kp0)



!---------------------------------------------------------------------------
! Finally, compute time tendency of m_l2*qvpc_l2 due to AB3 fluxes
!---------------------------------------------------------------------------

do k = 0, nlm
   qvpc_l2_m_l2_f(:,:,k,n3_f) = - div_horiz_tr1_m_flux(:,:,k) -      &
                                          div_vert_tr1_m_flux(:,:,k)
end do



! Vertical advection -- EF fluxes


! Calculate vertical mass flux interpolated to layer centers -- EF fluxes
F_eta(:,:,1) = p5*F_eta_l2_EF(:,:,1)
do k = 2, nlm-1
   F_eta(:,:,k) = p5*(F_eta_l2_EF(:,:,k)+F_eta_l2_EF(:,:,k-1))
end do
F_eta(:,:,nlm) = p5*F_eta_l2_EF(:,:,nlm-1)


! Start at bottom level (surface)

k = 0
kp0 = mod(k,2)
kp1 = mod(k+1,2)

! G_flx_contrib(:,:) = c0
! G_flx_contrib_2(:,:) = c0

F_tr1_eta(:,:,kp1) = p5 * F_eta(:,:,k+1) *                           &
         (qvpc_l2(:,:,k+1)+qvpc_l2(:,:,k))
F_tr_eta_2(:,:,kp1) = p5 * F_eta(:,:,k+1) *                          &
         (qc_l2(:,:,k+1)+qc_l2(:,:,k))

div_vert_tr1_m_flux(:,:,k) = inv_d_eta_l2(k) * F_tr1_eta(:,:,kp1)
qc_vert_adv_ef(:,:,k) = inv_d_eta_l2(k) * F_tr_eta_2(:,:,kp1)


! Move on up

do k = 1, nlm-2

   kp0 = mod(k,2)
   kp1 = mod(k+1,2)

   G_flx_contrib(:,:) = - inv12 * abs(F_eta(:,:,k+1)) *                 &
      ( 3*(qvpc_l2(:,:,k+1)-qvpc_l2(:,:,k))-                            &
           qvpc_l2(:,:,k+2)+qvpc_l2(:,:,k-1) )
   G_flx_contrib_2(:,:) = - inv12 * abs(F_eta(:,:,k+1)) *               &
      ( 3*(qc_l2(:,:,k+1)-qc_l2(:,:,k))-qc_l2(:,:,k+2)+qc_l2(:,:,k-1) )

   F_tr1_eta(:,:,kp1) = inv12 * F_eta(:,:,k+1) *                        &
    ( 7*(qvpc_l2(:,:,k+1)+qvpc_l2(:,:,k))-                              &
         qvpc_l2(:,:,k+2)-qvpc_l2(:,:,k-1) ) + G_flx_contrib(:,:)
   F_tr_eta_2(:,:,kp1) = inv12 * F_eta(:,:,k+1) *                       &
    ( 7*(qc_l2(:,:,k+1)+qc_l2(:,:,k))-qc_l2(:,:,k+2)-qc_l2(:,:,k-1) ) + &
       G_flx_contrib_2(:,:)

   div_vert_tr1_m_flux(:,:,k) = inv_d_eta_l2(k) *                       &
                          ( F_tr1_eta(:,:,kp1) - F_tr1_eta(:,:,kp0) )
   qc_vert_adv_ef(:,:,k) = inv_d_eta_l2(k) *                            &
                          ( F_tr_eta_2(:,:,kp1) - F_tr_eta_2(:,:,kp0) )

end do


! Almost there!

k = nlm-1
kp0 = mod(k,2)
kp1 = mod(k+1,2)

! G_flx_contrib(:,:) = c0
! G_flx_contrib_2(:,:) = c0

F_tr1_eta(:,:,kp1) = p5 * F_eta(:,:,k+1) *                              &
         (qvpc_l2(:,:,k+1)+qvpc_l2(:,:,k))
F_tr_eta_2(:,:,kp1) = p5 * F_eta(:,:,k+1) *                             &
         (qc_l2(:,:,k+1)+qc_l2(:,:,k))

div_vert_tr1_m_flux(:,:,k) = inv_d_eta_l2(k) *                          &
                       ( F_tr1_eta(:,:,kp1) - F_tr1_eta(:,:,kp0) )
qc_vert_adv_ef(:,:,k) = inv_d_eta_l2(k) *                               &
                       ( F_tr_eta_2(:,:,kp1) - F_tr_eta_2(:,:,kp0) )


! Top layer

k = nlm
kp0 = mod(k,2)
kp1 = mod(k+1,2)

div_vert_tr1_m_flux(:,:,k) = - inv_d_eta_l2(k) * F_tr1_eta(:,:,kp0)
qc_vert_adv_ef(:,:,k) = - inv_d_eta_l2(k) * F_tr_eta_2(:,:,kp0)



!---------------------------------------------------------------------------
! Finally, compute time tendency of m_l2*qvpc_l2 due to EF fluxes
!---------------------------------------------------------------------------

do k = 0, nlm
   qvpc_l2_m_l2_f_ef(:,:,k) = - div_vert_tr1_m_flux(:,:,k)
end do




!---------------------------------------------------------------------------
! Advect miscellaneous tracers with total eta_dot_l2, i.e,
!    1)  Rain water mixing ratio qr_l2
!---------------------------------------------------------------------------

! Start at bottom level (surface)

k = 0
kp0 = mod(k,2)
kp1 = mod(k+1,2)

! G_flx_contrib(:,:) = c0

F_eta_l2(:,:,kp1) = m_l2(:,:,1) * eta_dot_l2(:,:,1)
F_eta_tot(:,:) = p5*F_eta_l2(:,:,kp1)

F_tr2_eta(:,:,kp1) = p5 * F_eta_tot(:,:) *                           &
                                  (qr_l2(:,:,k+1)+qr_l2(:,:,k))

div_vert_tr2_m_flux(:,:,k) = inv_d_eta_l2(k) * F_tr2_eta(:,:,kp1)


! Move on up

do k = 1, nlm-2

   kp0 = mod(k,2)
   kp1 = mod(k+1,2)

   F_eta_l2(:,:,kp1) = m_l2(:,:,k+1) * eta_dot_l2(:,:,k+1)
   F_eta_tot(:,:) = p5*(F_eta_l2(:,:,kp1)+F_eta_l2(:,:,kp0))
   
   G_flx_contrib(:,:) = - inv12 * abs(F_eta_tot(:,:)) *                 &
      ( 3*(qr_l2(:,:,k+1)-qr_l2(:,:,k))-qr_l2(:,:,k+2)+qr_l2(:,:,k-1) )
   F_tr2_eta(:,:,kp1) = inv12 * F_eta_tot(:,:) *                        &
    ( 7*(qr_l2(:,:,k+1)+qr_l2(:,:,k))-qr_l2(:,:,k+2)-qr_l2(:,:,k-1) ) + &
       G_flx_contrib(:,:)

   div_vert_tr2_m_flux(:,:,k) = inv_d_eta_l2(k) *                       &
                         ( F_tr2_eta(:,:,kp1) - F_tr2_eta(:,:,kp0) )

end do


! Almost there

k = nlm-1
kp0 = mod(k,2)
kp1 = mod(k+1,2)

! G_flx_contrib(:,:) = c0

F_eta_tot(:,:) = p5*F_eta_l2(:,:,kp0)

F_tr2_eta(:,:,kp1) = p5 * F_eta_tot(:,:) *                           &
                             (qr_l2(:,:,k+1)+qr_l2(:,:,k))

div_vert_tr2_m_flux(:,:,k) = inv_d_eta_l2(k) *                       &
                      ( F_tr2_eta(:,:,kp1) - F_tr2_eta(:,:,kp0) )


! Top layer

k = nlm
kp0 = mod(k,2)
kp1 = mod(k+1,2)

div_vert_tr2_m_flux(:,:,k) = - inv_d_eta_l2(k) * F_tr2_eta(:,:,kp0)


!---------------------------------------------------------------------------
! Finally, compute time tendency of m_l2*qr_l2 due to total vertical fluxes
!---------------------------------------------------------------------------

do k = 0, nlm
   qr_l2_m_l2_f(:,:,k,n3_f) = - div_horiz_tr2_m_flux(:,:,k) -        &
                                          div_vert_tr2_m_flux(:,:,k)
end do




end subroutine get_tracer_tendencies

!======================================================================
! END OF GET_TRACER_TENDENCIES
!======================================================================



end module tracer_transport

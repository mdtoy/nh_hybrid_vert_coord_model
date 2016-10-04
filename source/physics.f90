module physics

!-----------------------------------------------------------------------
! PURPOSE: Calculates parameterized subgrid-scale turbulence fluxes
!          of momentum and heat.
!-----------------------------------------------------------------------

use kinds
use model_parameters
use physical_parameters
use eta_coordinate
use prognostics


implicit none
save



! Declare deformation tensor   (sec^-1)
real (kind = dbl_kind), dimension(im,jm,nlm) ::            &
           D_11, D_12, D_22, D_33
real (kind = dbl_kind), dimension(im,jm,0:nlm) ::          &
           D_13, D_23

! Declare Brunt-Vaisala frequency (squared) and Richardson number
real (kind = dbl_kind), dimension(im,jm,nlm) ::            &
           Nsq,    &  ! square of B-V frequency  (s^-2)
           Ri         ! Richardson number (-)

! Declare turbulent mixing coefficients
real (kind = dbl_kind), dimension(im,jm,nlm) ::            &
           Kmh, Kmv,  &  ! horizontal and vertical turb. mixing coeffs
                         ! for momentum (m^2/s)
           KHh, KHv      ! horizontal and vertical turb. mixing coeffs
                         ! potential temperature and tracers (m^2/s)

! Declare Reynolds stresses  (m^2/s^2)
real (kind = dbl_kind), dimension(im,jm,nlm) ::            &
           tau_11, tau_12, tau_21, tau_22, tau_33
real (kind = dbl_kind), dimension(im,jm,0:nlm) ::          &
           tau_13, tau_23, tau_31, tau_32

! Declare turbulent heat fluxes  (K m/s)
real (kind = dbl_kind), dimension(im,jm,nlm) ::            &
           H_th_3, H_qvpc_3, H_qr_3, H_qc_3
real (kind = dbl_kind), dimension(im,jm,0:nlm) ::          &
           H_th_1, H_th_2, H_qvpc_1, H_qvpc_2,   &
           H_qr_1, H_qr_2, H_qc_1, H_qc_2

! Declare contributions to momentum and potential temperature
! tendencies due to sub-grid scale turbulent flux divergence
real (kind = dbl_kind), dimension(im,jm,nlm) ::            &
           F_turb_u, F_turb_v  ! (m/s^2)
real (kind = dbl_kind), dimension(im,jm,1:nlm-1) ::        &
           F_turb_w_l2         ! (m/s^2)
real (kind = dbl_kind), dimension(im,jm,0:nlm) ::          &
           F_turb_th_l2,   &      ! (K/s)
           F_turb_qvpc_l2, F_turb_qr_l2, F_turb_qc_l2   ! (kg/kg/s)

! Various parameters needed for the parameterization
real (kind = dbl_kind), parameter ::                       &
           Pr = 0.7_dbl_kind,   &   ! turbulent Prandtl number
           invPr = c1/Pr,  &        ! inverse of Prandtl number
           dxdy = dx*dy,   &        ! horizontal grid area (m^2)
           k_mxng_const_sq = 0.21_dbl_kind**2 ! square of empirical const.






contains




!======================================================================
! BEGINNING OF CALC_SUBGRID_TURBULENT_FLX_DIV
!======================================================================

subroutine calc_subgrid_turbulent_flx_div( th, phi, rho, n4 )

!------------------------------------------------------------------------
! PURPOSE:
!   1) Computes the flow field's deformation tensor which is used to
!   calculate the Reynolds stresses and mixing coefficients in
!   the subgrid-scale turbulence parameterization.
!   Also, the square of the Brunt-Vaisala frequency and the Richardson
!   number are calculated here.
!   2) Computes the turbulent mixing coefficents for momentum, potential
!   temperature and tracers.
!   3) Computes subgrid-scale turbulent fluxes, i.e., Reynolds stresses
!   and turbulent heat fluxes.
!   4) Calculates the divergence of the subgrid-scale turbulent fluxes
!   which are components of the momentum and potential temperature
!   tendency equations.
!------------------------------------------------------------------------

implicit none

!------------------------------------------------------------------------
! INTENT IN
!------------------------------------------------------------------------
real (kind=dbl_kind), dimension(im,jm,nlm), intent(in) ::            &
               th,    &    ! layer center potential temperature (K)
               phi,   &    ! layer center geopotential (J/kg)
               rho         ! layer center density (kg/m^3)

integer, intent(in) :: n4                    ! current time index

!------------------------------------------------------------------------
! LOCAL
!------------------------------------------------------------------------
integer :: i,j,k

integer :: kp0, kp1     ! indices for k level (kp0 = k+0) and the k+1
                        ! level (kp1 = k+1)  (they toggle between 0 and 1)

real (kind=dbl_kind), parameter ::                                   &
               epsln = 1.E-10_dbl_kind  ! parameter to prevent denominator of
                                        ! Richardson number from becoming zero

real (kind=dbl_kind), dimension(im,jm) ::                            &
               f1, f2, f3, f4, f5,    &
               f6, f7, f8, f9, f10,   &   ! working variables
               u_l2, v_l2,   &   ! horizontal velocity averaged to lyr edges
               w                 ! vertical velocity averaged to lyr centers

real (kind=dbl_kind), dimension(im,jm,0:1) ::                        &
               phi_l2_i_j, phi_i, phi_j, th_i, th_j,      &
               qvpc_i, qvpc_j, qr_i, qr_j, qc_i, qc_j,  &
               u_l2_i, u_l2_j, v_l2_i, v_l2_j, w_i, w_j,  &  ! horiz. averages
               H_1_hacek_i, H_2_hacek_j, H_qvpc_1_hacek_i,  &
               H_qvpc_2_hacek_j, H_qr_1_hacek_i, H_qr_2_hacek_j, &
               H_qc_1_hacek_i, H_qc_2_hacek_j, tau_31_hacek_i,  &
               tau_32_hacek_j, tau_11_hat_i, tau_12_hat_j, &
               tau_21_hat_i, tau_22_hat_j   ! horiz. and vert. averages

real (kind = dbl_kind) ::                                            &
               Defsq,     &             ! "magnitude" of deformation
               fpt1, fpt2, fpt3, fpt4   ! working variables

real (kind=dbl_kind) ::                                              &
               KHh_hat_i, KHh_hat_j ! KH interpolated to heat flux points




!
! 1) Compute deformation tensor.
!    (Also, calculate the square of the Brunt-Vaisala frequency and the
!    Richardson number.)
!

D_13(:,:,0) = c0
D_13(:,:,nlm) = c0
D_23(:,:,0) = c0
D_23(:,:,nlm) = c0


!---------------------------------------------------------------------------
! start at k = 1
!---------------------------------------------------------------------------

kp0 = mod(0,2)
kp1 = mod(1,2)

w(:,:) = p5*(w_l2_0(:,:)+w_l2(:,:,1,n4))       !  w(k=1)
u_l2(:,:) = p5*(u(:,:,1,n4)+u(:,:,2,n4))       !  u_l2(k=1)
v_l2(:,:) = p5*(v(:,:,1,n4)+v(:,:,2,n4))       !  v_l2(k=1)
do j = 1, jm
   do i = 1, im
      w_i(i,j,kp0) = p5*(w(i,j)+w(im1(i),j))
      w_j(i,j,kp0) = p5*(w(i,j)+w(i,jm1(j)))
      phi_i(i,j,kp0) = p5*(phi(i,j,1)+phi(im1(i),j,1))
      phi_i(i,j,kp1) = p5*(phi(i,j,2)+phi(im1(i),j,2))
      phi_j(i,j,kp0) = p5*(phi(i,j,1)+phi(i,jm1(j),1))
      phi_j(i,j,kp1) = p5*(phi(i,j,2)+phi(i,jm1(j),2))
      phi_l2_i_j(i,j,kp1) = p25*(phi_l2(i,j,1,n4)+phi_l2(im1(i),j,1,n4)+ &
                       phi_l2(im1(i),jm1(j),1,n4)+phi_l2(i,jm1(j),1,n4))
      f1(i,j) = p5*(u(ip1(i),j,1,n4)+u(i,j,1,n4))    ! u_i(k=1)
      f2(i,j) = p5*(u(i,j,1,n4)+u(i,jm1(j),1,n4))    ! u_j(k=1)
      f3(i,j) = p5*(v(i,j,1,n4)+v(im1(i),j,1,n4))    ! v_i(k=1)
      f4(i,j) = p5*(v(i,jp1(j),1,n4)+v(i,j,1,n4))    ! v_j(k=1)
      u_l2_i(i,j,kp1) = p5*(u_l2(ip1(i),j)+u_l2(i,j))
      u_l2_j(i,j,kp1) = p5*(u_l2(i,j)+u_l2(i,jm1(j)))
      v_l2_i(i,j,kp1) = p5*(v_l2(i,j)+v_l2(im1(i),j))
      v_l2_j(i,j,kp1) = p5*(v_l2(i,jp1(j))+v_l2(i,j))
   end do
   do i = 1, im
      f5(i,j) = p5*(phi_j(i,j,kp0)+phi_j(im1(i),j,kp0))   ! phi_i_j(k=1)
   end do
end do

w(:,:) = p5*(w_l2(:,:,1,n4)+w_l2(:,:,2,n4))    !  w(k=2)
do j = 1, jm
   do i = 1, im
      w_i(i,j,kp1) = p5*(w(i,j)+w(im1(i),j))
      w_j(i,j,kp1) = p5*(w(i,j)+w(i,jm1(j)))
   end do
end do


do j = 1, jm
   do i = 1, im
      D_11(i,j,1) = 2*invdx*( (u(ip1(i),j,1,n4)-u(i,j,1,n4)) -            &
         (phi_i(ip1(i),j,kp0)-phi_i(i,j,kp0))*(u_l2_i(i,j,kp1)-f1(i,j)) / &
            (phi_l2(i,j,1,n4)-phi(i,j,1)) )
      D_12(i,j,1) = invdy*(u(i,j,1,n4)-u(i,jm1(j),1,n4)) + invdx*         &
         (v(i,j,1,n4)-v(im1(i),j,1,n4)) - ( invdy*(phi_i(i,j,kp0)-        &
          phi_i(i,jm1(j),kp0))*(u_l2_j(i,j,kp1)-f2(i,j)) + invdx*         &
          (phi_j(i,j,kp0)-phi_j(im1(i),j,kp0))*(v_l2_i(i,j,kp1)-          &
          f3(i,j)) ) / (phi_l2_i_j(i,j,kp1)-f5(i,j))
      D_13(i,j,1) = invdx*(w_l2(i,j,1,n4)-w_l2(im1(i),j,1,n4)) +          &
         ( (u(i,j,2,n4)-u(i,j,1,n4))*grav - invdx*(phi_l2(i,j,1,n4)-      &
           phi_l2(im1(i),j,1,n4))*(w_i(i,j,kp1)-w_i(i,j,kp0)) ) /         &
           (phi_i(i,j,kp1)-phi_i(i,j,kp0))
      D_22(i,j,1) = 2*invdy*( (v(i,jp1(j),1,n4)-v(i,j,1,n4)) -            &
         (phi_j(i,jp1(j),kp0)-phi_j(i,j,kp0))*(v_l2_j(i,j,kp1)-f4(i,j)) / &
            (phi_l2(i,j,1,n4)-phi(i,j,1)) )
      D_23(i,j,1) = invdy*(w_l2(i,j,1,n4)-w_l2(i,jm1(j),1,n4)) +          &
         ( (v(i,j,2,n4)-v(i,j,1,n4))*grav - invdy*(phi_l2(i,j,1,n4)-      &
           phi_l2(i,jm1(j),1,n4))*(w_j(i,j,kp1)-w_j(i,j,kp0)) ) /         &
           (phi_j(i,j,kp1)-phi_j(i,j,kp0))
   end do
end do
D_33(:,:,1) = 2*(w_l2(:,:,1,n4)-w_l2_0(:,:)) * grav /                &
                                   (phi_l2(:,:,1,n4)-phis(:,:))

Nsq(:,:,1) = (grav**2)*(th_l2(:,:,1,n4)-th_l2(:,:,0,n4)) / ( p5*     &
     (th_l2(:,:,1,n4)+th_l2(:,:,0,n4))*(phi_l2(:,:,1,n4)-phis(:,:)) )

Ri(:,:,1) = Nsq(:,:,1)*(invgrav**2)*(phi_l2(:,:,1,n4)-phi(:,:,1))**2 /    &
   ( (u_l2_i(:,:,kp1)-f1(:,:))**2 + (v_l2_j(:,:,kp1)-f4(:,:))**2 + epsln )



!---------------------------------------------------------------------------
! move upward
!---------------------------------------------------------------------------

do k = 2, nlm-2

   kp0 = mod(k-1,2)
   kp1 = mod(k,2)

   w(:,:) = p5*(w_l2(:,:,k,n4)+w_l2(:,:,k+1,n4))    ! w(k+1)
   u_l2(:,:) = p5*(u(:,:,k,n4)+u(:,:,k+1,n4))       ! u_l2(k)
   v_l2(:,:) = p5*(v(:,:,k,n4)+v(:,:,k+1,n4))       ! v_l2(k)
   do j = 1, jm
      do i = 1, im
         w_i(i,j,kp1) = p5*(w(i,j)+w(im1(i),j))
         w_j(i,j,kp1) = p5*(w(i,j)+w(i,jm1(j)))
         phi_i(i,j,kp1) = p5*(phi(i,j,k+1)+phi(im1(i),j,k+1))
         phi_j(i,j,kp1) = p5*(phi(i,j,k+1)+phi(i,jm1(j),k+1))
         phi_l2_i_j(i,j,kp1) = p25*(phi_l2(i,j,k,n4)+phi_l2(im1(i),j,k,n4)+ &
                       phi_l2(im1(i),jm1(j),k,n4)+phi_l2(i,jm1(j),k,n4))
         u_l2_i(i,j,kp1) = p5*(u_l2(ip1(i),j)+u_l2(i,j))
         u_l2_j(i,j,kp1) = p5*(u_l2(i,j)+u_l2(i,jm1(j)))
         v_l2_i(i,j,kp1) = p5*(v_l2(i,j)+v_l2(im1(i),j))
         v_l2_j(i,j,kp1) = p5*(v_l2(i,jp1(j))+v_l2(i,j))
      end do
   end do


   do j = 1, jm
      do i = 1, im
         D_11(i,j,k) = 2*invdx*( (u(ip1(i),j,k,n4)-u(i,j,k,n4)) -          &
            (phi_i(ip1(i),j,kp0)-phi_i(i,j,kp0))*(u_l2_i(i,j,kp1)-         &
               u_l2_i(i,j,kp0)) / (phi_l2(i,j,k,n4)-phi_l2(i,j,k-1,n4)) )
         D_12(i,j,k) = invdy*(u(i,j,k,n4)-u(i,jm1(j),k,n4)) + invdx*       &
            (v(i,j,k,n4)-v(im1(i),j,k,n4)) - ( invdy*(phi_i(i,j,kp0)-      &
             phi_i(i,jm1(j),kp0))*(u_l2_j(i,j,kp1)-u_l2_j(i,j,kp0)) +      &
             invdx* (phi_j(i,j,kp0)-phi_j(im1(i),j,kp0))*(v_l2_i(i,j,kp1)- &
             v_l2_i(i,j,kp0)) ) / (phi_l2_i_j(i,j,kp1)-phi_l2_i_j(i,j,kp0))
         D_13(i,j,k) = invdx*(w_l2(i,j,k,n4)-w_l2(im1(i),j,k,n4)) +        &
            ( (u(i,j,k+1,n4)-u(i,j,k,n4))*grav - invdx*(phi_l2(i,j,k,n4)-  &
               phi_l2(im1(i),j,k,n4))*(w_i(i,j,kp1)-w_i(i,j,kp0)) ) /      &
                (phi_i(i,j,kp1)-phi_i(i,j,kp0))
         D_22(i,j,k) = 2*invdy*( (v(i,jp1(j),k,n4)-v(i,j,k,n4)) -          &
            (phi_j(i,jp1(j),kp0)-phi_j(i,j,kp0))*(v_l2_j(i,j,kp1)-         &
               v_l2_j(i,j,kp0)) / (phi_l2(i,j,k,n4)-phi_l2(i,j,k-1,n4)) )
         D_23(i,j,k) = invdy*(w_l2(i,j,k,n4)-w_l2(i,jm1(j),k,n4)) +        &
            ( (v(i,j,k+1,n4)-v(i,j,k,n4))*grav - invdy*(phi_l2(i,j,k,n4)-  &
               phi_l2(i,jm1(j),k,n4))*(w_j(i,j,kp1)-w_j(i,j,kp0)) ) /      &
                (phi_j(i,j,kp1)-phi_j(i,j,kp0))
      end do
   end do
   D_33(:,:,k) = 2*(w_l2(:,:,k,n4)-w_l2(:,:,k-1,n4)) * grav /        &
                              (phi_l2(:,:,k,n4)-phi_l2(:,:,k-1,n4))

   Nsq(:,:,k) = (grav**2)*(th_l2(:,:,k,n4)-th_l2(:,:,k-1,n4)) /      &
                       ( p5*(th_l2(:,:,k,n4)+th_l2(:,:,k-1,n4))*     &
                          (phi_l2(:,:,k,n4)-phi_l2(:,:,k-1,n4)) )

   Ri(:,:,k) = Nsq(:,:,k)*(invgrav**2)*                              &
                  (phi_l2(:,:,k,n4)-phi_l2(:,:,k-1,n4))**2 /         &
                     ( (u_l2_i(:,:,kp1)-u_l2_i(:,:,kp0))**2 +        &
                       (v_l2_j(:,:,kp1)-v_l2_j(:,:,kp0))**2 + epsln )

end do


!---------------------------------------------------------------------------
! k = nlm-1
!---------------------------------------------------------------------------

k = nlm-1

kp0 = mod(k-1,2)
kp1 = mod(k,2)

w(:,:) = p5*w_l2(:,:,k,n4)                       ! w(nlm)
u_l2(:,:) = p5*(u(:,:,k,n4)+u(:,:,k+1,n4))       ! u_l2(nlm-1)
v_l2(:,:) = p5*(v(:,:,k,n4)+v(:,:,k+1,n4))       ! v_l2(nlm-1)
do j = 1, jm
   do i = 1, im
      w_i(i,j,kp1) = p5*(w(i,j)+w(im1(i),j))
      w_j(i,j,kp1) = p5*(w(i,j)+w(i,jm1(j)))
      phi_i(i,j,kp1) = p5*(phi(i,j,k+1)+phi(im1(i),j,k+1))
      phi_j(i,j,kp1) = p5*(phi(i,j,k+1)+phi(i,jm1(j),k+1))
      phi_l2_i_j(i,j,kp1) = p25*(phi_l2(i,j,k,n4)+phi_l2(im1(i),j,k,n4)+ &
                    phi_l2(im1(i),jm1(j),k,n4)+phi_l2(i,jm1(j),k,n4))
      u_l2_i(i,j,kp1) = p5*(u_l2(ip1(i),j)+u_l2(i,j))
      u_l2_j(i,j,kp1) = p5*(u_l2(i,j)+u_l2(i,jm1(j)))
      v_l2_i(i,j,kp1) = p5*(v_l2(i,j)+v_l2(im1(i),j))
      v_l2_j(i,j,kp1) = p5*(v_l2(i,jp1(j))+v_l2(i,j))
   end do
end do


do j = 1, jm
   do i = 1, im
      D_11(i,j,k) = 2*invdx*( (u(ip1(i),j,k,n4)-u(i,j,k,n4)) -          &
         (phi_i(ip1(i),j,kp0)-phi_i(i,j,kp0))*(u_l2_i(i,j,kp1)-         &
            u_l2_i(i,j,kp0)) / (phi_l2(i,j,k,n4)-phi_l2(i,j,k-1,n4)) )
      D_12(i,j,k) = invdy*(u(i,j,k,n4)-u(i,jm1(j),k,n4)) + invdx*       &
         (v(i,j,k,n4)-v(im1(i),j,k,n4)) - ( invdy*(phi_i(i,j,kp0)-      &
          phi_i(i,jm1(j),kp0))*(u_l2_j(i,j,kp1)-u_l2_j(i,j,kp0)) +      &
          invdx* (phi_j(i,j,kp0)-phi_j(im1(i),j,kp0))*(v_l2_i(i,j,kp1)- &
          v_l2_i(i,j,kp0)) ) / (phi_l2_i_j(i,j,kp1)-phi_l2_i_j(i,j,kp0))
      D_13(i,j,k) = invdx*(w_l2(i,j,k,n4)-w_l2(im1(i),j,k,n4)) +        &
         ( (u(i,j,k+1,n4)-u(i,j,k,n4))*grav - invdx*(phi_l2(i,j,k,n4)-  &
            phi_l2(im1(i),j,k,n4))*(w_i(i,j,kp1)-w_i(i,j,kp0)) ) /      &
             (phi_i(i,j,kp1)-phi_i(i,j,kp0))
      D_22(i,j,k) = 2*invdy*( (v(i,jp1(j),k,n4)-v(i,j,k,n4)) -          &
         (phi_j(i,jp1(j),kp0)-phi_j(i,j,kp0))*(v_l2_j(i,j,kp1)-         &
            v_l2_j(i,j,kp0)) / (phi_l2(i,j,k,n4)-phi_l2(i,j,k-1,n4)) )
      D_23(i,j,k) = invdy*(w_l2(i,j,k,n4)-w_l2(i,jm1(j),k,n4)) +        &
         ( (v(i,j,k+1,n4)-v(i,j,k,n4))*grav - invdy*(phi_l2(i,j,k,n4)-  &
            phi_l2(i,jm1(j),k,n4))*(w_j(i,j,kp1)-w_j(i,j,kp0)) ) /      &
             (phi_j(i,j,kp1)-phi_j(i,j,kp0))
   end do
end do
D_33(:,:,k) = 2*(w_l2(:,:,k,n4)-w_l2(:,:,k-1,n4)) * grav /           &
                           (phi_l2(:,:,k,n4)-phi_l2(:,:,k-1,n4))

Nsq(:,:,k) = (grav**2)*(th_l2(:,:,k,n4)-th_l2(:,:,k-1,n4)) /         &
                    ( p5*(th_l2(:,:,k,n4)+th_l2(:,:,k-1,n4))*        &
                       (phi_l2(:,:,k,n4)-phi_l2(:,:,k-1,n4)) )

Ri(:,:,k) = Nsq(:,:,k)*(invgrav**2)*                                 &
               (phi_l2(:,:,k,n4)-phi_l2(:,:,k-1,n4))**2 /            &
                  ( (u_l2_i(:,:,kp1)-u_l2_i(:,:,kp0))**2 +           &
                    (v_l2_j(:,:,kp1)-v_l2_j(:,:,kp0))**2 + epsln )


!---------------------------------------------------------------------------
! k = nlm (top layer)
!---------------------------------------------------------------------------

k = nlm

kp0 = mod(k-1,2)
kp1 = mod(k,2)

do j = 1, jm
   do i = 1, im
      f1(i,j) = p5*(u(ip1(i),j,k,n4)+u(i,j,k,n4))       ! u_i(k=nlm)
      f2(i,j) = p5*(u(i,j,k,n4)+u(i,jm1(j),k,n4))       ! u_j(k=nlm)
      f3(i,j) = p5*(v(i,j,k,n4)+v(im1(i),j,k,n4))       ! v_i(k=nlm)
      f4(i,j) = p5*(v(i,jp1(j),k,n4)+v(i,j,k,n4))       ! v_j(k=nlm)
      f5(i,j) = p5*(phi_j(i,j,kp0)+phi_j(im1(i),j,kp0)) ! phi_i_j(k=nlm)
   end do
end do


do j = 1, jm
   do i = 1, im
      D_11(i,j,k) = 2*invdx*( (u(ip1(i),j,k,n4)-u(i,j,k,n4)) -          &
         (phi_i(ip1(i),j,kp0)-phi_i(i,j,kp0))*(f1(i,j)-                 &
            u_l2_i(i,j,kp0)) / (phi(i,j,k)-phi_l2(i,j,k-1,n4)) )
      D_12(i,j,k) = invdy*(u(i,j,k,n4)-u(i,jm1(j),k,n4)) + invdx*       &
         (v(i,j,k,n4)-v(im1(i),j,k,n4)) - ( invdy*(phi_i(i,j,kp0)-      &
          phi_i(i,jm1(j),kp0))*(f2(i,j)-u_l2_j(i,j,kp0)) +              &
          invdx* (phi_j(i,j,kp0)-phi_j(im1(i),j,kp0))*(f3(i,j)-         &
          v_l2_i(i,j,kp0)) ) / (f5(i,j)-phi_l2_i_j(i,j,kp0))
      D_22(i,j,k) = 2*invdy*( (v(i,jp1(j),k,n4)-v(i,j,k,n4)) -          &
         (phi_j(i,jp1(j),kp0)-phi_j(i,j,kp0))*(f4(i,j)-                 &
            v_l2_j(i,j,kp0)) / (phi(i,j,k)-phi_l2(i,j,k-1,n4)) )
   end do
end do
D_33(:,:,k) = 2*(-w_l2(:,:,k-1,n4)) * grav /                         &
                           (grav*z_top-phi_l2(:,:,k-1,n4))


Nsq(:,:,k) = (grav**2)*(th_l2(:,:,k,n4)-th_l2(:,:,k-1,n4)) /         &
                    ( p5*(th_l2(:,:,k,n4)+th_l2(:,:,k-1,n4))*        &
                       (grav*z_top-phi_l2(:,:,k-1,n4)) )

Ri(:,:,k) = Nsq(:,:,k)*(invgrav**2)*                                 &
               (phi(:,:,k)-phi_l2(:,:,k-1,n4))**2 /                  &
                  ( (f1(:,:)-u_l2_i(:,:,kp0))**2 +                   &
                    (f4(:,:)-v_l2_j(:,:,kp0))**2 + epsln )





!
! 2) Compute the turbulent mixing coefficents for momentum, potential
!    temperature and tracers (if any).
!

! First, calculate Kmh and Kmv

! Bottom layer
k = 1
do j = 1, jm
   do i = 1, im
      fpt1 = p25*( D_12(i,j,k) + D_12(ip1(i),j,k) +                  &
                 D_12(ip1(i),jp1(j),k) + D_12(i,jp1(j),k) )
      fpt2 = p25*( D_13(i,j,k) + D_13(i,j,k-1) +                     &
                 D_13(ip1(i),j,k) + D_13(ip1(i),j,k-1) )
      fpt3 = p25*( D_23(i,j,k) + D_23(i,j,k-1) +                     &
                 D_23(i,jp1(j),k) + D_23(i,jp1(j),k-1) )
      Defsq = p5*( D_11(i,j,k)**2 + D_22(i,j,k)**2 +                 &
                   D_33(i,j,k)**2 ) + fpt1**2 + fpt2**2 + fpt3**2
      fpt1 = ( max(Defsq-Nsq(i,j,k)*invPr,c0) )**p5
      Kmh(i,j,k) = k_mxng_const_sq  * dxdy * fpt1
      Kmv(i,j,k) = k_mxng_const_sq *                                 &
          (invgrav*(phi_l2(i,j,k,n4)-phis(i,j)))**2 * fpt1
   end do
end do

! Move upward
do k = 2, nlm-1
   do j = 1, jm
      do i = 1, im
         fpt1 = p25*( D_12(i,j,k) + D_12(ip1(i),j,k) +               &
                    D_12(ip1(i),jp1(j),k) + D_12(i,jp1(j),k) )
         fpt2 = p25*( D_13(i,j,k) + D_13(i,j,k-1) +                  &
                    D_13(ip1(i),j,k) + D_13(ip1(i),j,k-1) )
         fpt3 = p25*( D_23(i,j,k) + D_23(i,j,k-1) +                  &
                    D_23(i,jp1(j),k) + D_23(i,jp1(j),k-1) )
         Defsq = p5*( D_11(i,j,k)**2 + D_22(i,j,k)**2 +              &
                      D_33(i,j,k)**2 ) + fpt1**2 + fpt2**2 + fpt3**2
         fpt1 = ( max(Defsq-Nsq(i,j,k)*invPr,c0) )**p5
         Kmh(i,j,k) = k_mxng_const_sq  * dxdy * fpt1
         Kmv(i,j,k) = k_mxng_const_sq *                              &
             (invgrav*(phi_l2(i,j,k,n4)-phi_l2(i,j,k-1,n4)))**2 * fpt1
      end do
   end do
end do

! Top layer
k = nlm
do j = 1, jm
   do i = 1, im
      fpt1 = p25*( D_12(i,j,k) + D_12(ip1(i),j,k) +                  &
                 D_12(ip1(i),jp1(j),k) + D_12(i,jp1(j),k) )
      fpt2 = p25*( D_13(i,j,k) + D_13(i,j,k-1) +                     &
                 D_13(ip1(i),j,k) + D_13(ip1(i),j,k-1) )
      fpt3 = p25*( D_23(i,j,k) + D_23(i,j,k-1) +                     &
                 D_23(i,jp1(j),k) + D_23(i,jp1(j),k-1) )
      Defsq = p5*( D_11(i,j,k)**2 + D_22(i,j,k)**2 +                 &
                   D_33(i,j,k)**2 ) + fpt1**2 + fpt2**2 + fpt3**2
      fpt1 = ( max(Defsq-Nsq(i,j,k)*invPr,c0) )**p5
      Kmh(i,j,k) = k_mxng_const_sq  * dxdy * fpt1
      Kmv(i,j,k) = k_mxng_const_sq *                                 &
          (invgrav*(grav*z_top-phi_l2(i,j,k-1,n4)))**2 * fpt1
   end do
end do


! Now, calculate KHh and KHv from Km's
KHh = invPr*Kmh
KHv = invPr*Kmv




!
! 3) Compute subgrid-scale turbulent fluxes, i.e., Reynolds stresses
!    and turbulent heat fluxes.
!


!---------------------------------------------------------------------------
! start at surface (k = 0)
!---------------------------------------------------------------------------

kp0 = mod(0,2)
kp1 = mod(1,2)

! "No slip" lower boundary condition
tau_13(:,:,0) = c0
tau_23(:,:,0) = c0
tau_31(:,:,0) = c0
tau_32(:,:,0) = c0

do j = 1, jm
   do i = 1, im
      f1(i,j) = p5*(th_l2(i,j,0,n4)+th_l2(im1(i),j,0,n4))   ! th_l2_i(k=0)
      f2(i,j) = p5*(th_l2(i,j,0,n4)+th_l2(i,jm1(j),0,n4))   ! th_l2_j(k=0)
      f3(i,j) = p5*(phis(i,j)+phis(im1(i),j))   ! phis_i
      f4(i,j) = p5*(phis(i,j)+phis(i,jm1(j)))   ! phis_j
      f5(i,j) = p5*(qvpc_l2(i,j,0)+qvpc_l2(im1(i),j,0))   ! qvpc_l2_i(k=0)
      f6(i,j) = p5*(qvpc_l2(i,j,0)+qvpc_l2(i,jm1(j),0))   ! qvpc_l2_j(k=0)
      f7(i,j) = p5*(qr_l2(i,j,0)+qr_l2(im1(i),j,0))       ! qr_l2_i(k=0)
      f8(i,j) = p5*(qr_l2(i,j,0)+qr_l2(i,jm1(j),0))       ! qr_l2_j(k=0)
      f9(i,j) = p5*(qc_l2(i,j,0)+qc_l2(im1(i),j,0))       ! qc_l2_i(k=0)
      f10(i,j) = p5*(qc_l2(i,j,0)+qc_l2(i,jm1(j),0))       ! qc_l2_j(k=0)
      th_i(i,j,kp1) = p5*(th(i,j,1)+th(im1(i),j,1))
      th_j(i,j,kp1) = p5*(th(i,j,1)+th(i,jm1(j),1))
      phi_i(i,j,kp1) = p5*(phi(i,j,1)+phi(im1(i),j,1))
      phi_j(i,j,kp1) = p5*(phi(i,j,1)+phi(i,jm1(j),1))
      qvpc_i(i,j,kp1) = p25*( qvpc_l2(i,j,1)+qvpc_l2(im1(i),j,1) +   &
                              qvpc_l2(i,j,0)+qvpc_l2(im1(i),j,0) )
      qvpc_j(i,j,kp1) = p25*( qvpc_l2(i,j,1)+qvpc_l2(i,jm1(j),1) +   &
                              qvpc_l2(i,j,0)+qvpc_l2(i,jm1(j),0) )
      qr_i(i,j,kp1) = p25*( qr_l2(i,j,1)+qr_l2(im1(i),j,1) +         &
                            qr_l2(i,j,0)+qr_l2(im1(i),j,0) )
      qr_j(i,j,kp1) = p25*( qr_l2(i,j,1)+qr_l2(i,jm1(j),1) +         &
                            qr_l2(i,j,0)+qr_l2(i,jm1(j),0) )
      qc_i(i,j,kp1) = p25*( qc_l2(i,j,1)+qc_l2(im1(i),j,1) +         &
                            qc_l2(i,j,0)+qc_l2(im1(i),j,0) )
      qc_j(i,j,kp1) = p25*( qc_l2(i,j,1)+qc_l2(i,jm1(j),1) +         &
                            qc_l2(i,j,0)+qc_l2(i,jm1(j),0) )
   end do
end do

do j = 1, jm
   do i = 1, im
      ! turbulent heat fluxes
      KHh_hat_i = p5*(KHh(i,j,1)+KHh(im1(i),j,1))
      KHh_hat_j = p5*(KHh(i,j,1)+KHh(i,jm1(j),1))
      fpt1 = p5*(rho(i,j,1)+rho(im1(i),j,1))  ! rho_l2_i(k=0)
      fpt2 = p5*(rho(i,j,1)+rho(i,jm1(j),1))  ! rho_l2_j(k=0)
      H_th_1(i,j,0) = fpt1*KHh_hat_i*invdx*(th_l2(i,j,0,n4)-         &
         th_l2(im1(i),j,0,n4)-(phis(i,j)-phis(im1(i),j))*            &
         (th_i(i,j,kp1)-f1(i,j))/(phi_i(i,j,kp1)-f3(i,j)))
      H_th_2(i,j,0) = fpt2*KHh_hat_j*invdy*(th_l2(i,j,0,n4)-         &
         th_l2(i,jm1(j),0,n4)-(phis(i,j)-phis(i,jm1(j)))*            &
         (th_j(i,j,kp1)-f2(i,j))/(phi_j(i,j,kp1)-f4(i,j)))
      H_qvpc_1(i,j,0) = fpt1*KHh_hat_i*invdx*(qvpc_l2(i,j,0)-        &
         qvpc_l2(im1(i),j,0)-(phis(i,j)-phis(im1(i),j))*             &
         (qvpc_i(i,j,kp1)-f5(i,j))/(phi_i(i,j,kp1)-f3(i,j)))
      H_qvpc_2(i,j,0) = fpt2*KHh_hat_j*invdy*(qvpc_l2(i,j,0)-        &
         qvpc_l2(i,jm1(j),0)-(phis(i,j)-phis(i,jm1(j)))*             &
         (qvpc_j(i,j,kp1)-f6(i,j))/(phi_j(i,j,kp1)-f4(i,j)))
      H_qr_1(i,j,0) = fpt1*KHh_hat_i*invdx*(qr_l2(i,j,0)-            &
         qr_l2(im1(i),j,0)-(phis(i,j)-phis(im1(i),j))*               &
         (qr_i(i,j,kp1)-f7(i,j))/(phi_i(i,j,kp1)-f3(i,j)))
      H_qr_2(i,j,0) = fpt2*KHh_hat_j*invdy*(qr_l2(i,j,0)-            &
         qr_l2(i,jm1(j),0)-(phis(i,j)-phis(i,jm1(j)))*               &
         (qr_j(i,j,kp1)-f8(i,j))/(phi_j(i,j,kp1)-f4(i,j)))
      H_qc_1(i,j,0) = fpt1*KHh_hat_i*invdx*(qc_l2(i,j,0)-            &
         qc_l2(im1(i),j,0)-(phis(i,j)-phis(im1(i),j))*               &
         (qc_i(i,j,kp1)-f9(i,j))/(phi_i(i,j,kp1)-f3(i,j)))
      H_qc_2(i,j,0) = fpt2*KHh_hat_j*invdy*(qc_l2(i,j,0)-            &
         qc_l2(i,jm1(j),0)-(phis(i,j)-phis(i,jm1(j)))*               &
         (qc_j(i,j,kp1)-f10(i,j))/(phi_j(i,j,kp1)-f4(i,j)))
   end do
end do


!---------------------------------------------------------------------------
! move upward
!---------------------------------------------------------------------------

do k = 1, nlm-1

   kp0 = mod(k,2)
   kp1 = mod(k+1,2)

   do j = 1, jm
      do i = 1, im
         th_i(i,j,kp1) = p5*(th(i,j,k+1)+th(im1(i),j,k+1))
         th_j(i,j,kp1) = p5*(th(i,j,k+1)+th(i,jm1(j),k+1))
         phi_i(i,j,kp1) = p5*(phi(i,j,k+1)+phi(im1(i),j,k+1))
         phi_j(i,j,kp1) = p5*(phi(i,j,k+1)+phi(i,jm1(j),k+1))
         qvpc_i(i,j,kp1) = p25*( qvpc_l2(i,j,k+1)+qvpc_l2(im1(i),j,k+1) +  &
                                 qvpc_l2(i,j,k)+qvpc_l2(im1(i),j,k) )
         qvpc_j(i,j,kp1) = p25*( qvpc_l2(i,j,k+1)+qvpc_l2(i,jm1(j),k+1) +  &
                                 qvpc_l2(i,j,k)+qvpc_l2(i,jm1(j),k) )
         qr_i(i,j,kp1) = p25*( qr_l2(i,j,k+1)+qr_l2(im1(i),j,k+1) +        &
                               qr_l2(i,j,k)+qr_l2(im1(i),j,k) )
         qr_j(i,j,kp1) = p25*( qr_l2(i,j,k+1)+qr_l2(i,jm1(j),k+1) +        &
                               qr_l2(i,j,k)+qr_l2(i,jm1(j),k) )
         qc_i(i,j,kp1) = p25*( qc_l2(i,j,k+1)+qc_l2(im1(i),j,k+1) +        &
                               qc_l2(i,j,k)+qc_l2(im1(i),j,k) )
         qc_j(i,j,kp1) = p25*( qc_l2(i,j,k+1)+qc_l2(i,jm1(j),k+1) +        &
                               qc_l2(i,j,k)+qc_l2(i,jm1(j),k) )
      end do
   end do

   do j = 1, jm
      do i = 1, im

         ! Reynolds stresses
         fpt1 = p25*(rho(i,j,k)+rho(im1(i),j,k)+rho(i,j,k+1)+        &
                     rho(im1(i),j,k+1))  ! rho_l2_i
         fpt2 = p25*(rho(i,j,k)+rho(i,jm1(j),k)+rho(i,j,k+1)+        &
                     rho(i,jm1(j),k+1))  ! rho_l2_j
         tau_13(i,j,k) = fpt1*p25*(Kmv(i,j,k)+Kmv(i,j,k+1)+          &
                   Kmv(im1(i),j,k)+Kmv(im1(i),j,k+1))*D_13(i,j,k)
         tau_23(i,j,k) = fpt2*p25*(Kmv(i,j,k)+Kmv(i,j,k+1)+          &
                   Kmv(i,jm1(j),k)+Kmv(i,jm1(j),k+1))*D_23(i,j,k)
         tau_31(i,j,k) = fpt1*p25*(Kmh(i,j,k)+Kmh(i,j,k+1)+          &
                   Kmh(im1(i),j,k)+Kmh(im1(i),j,k+1))*D_13(i,j,k)
         tau_32(i,j,k) = fpt2*p25*(Kmh(i,j,k)+Kmh(i,j,k+1)+          &
                   Kmh(i,jm1(j),k)+Kmh(i,jm1(j),k+1))*D_23(i,j,k)

         tau_11(i,j,k) = rho(i,j,k)*Kmh(i,j,k)*D_11(i,j,k)
         tau_22(i,j,k) = rho(i,j,k)*Kmh(i,j,k)*D_22(i,j,k)
         tau_33(i,j,k) = rho(i,j,k)*Kmv(i,j,k)*D_33(i,j,k)

         fpt4 = p25*(rho(i,j,k)+rho(im1(i),j,k)+                     &
                     rho(im1(i),jm1(j),k)+rho(i,jm1(j),k)) ! rho_i_j
         tau_12(i,j,k) = fpt4*p25*(Kmh(i,j,k)+Kmh(im1(i),j,k)+       &
                   Kmh(im1(i),jm1(j),k)+Kmh(i,jm1(j),k))*D_12(i,j,k)

         ! turbulent heat fluxes
         KHh_hat_i = p25*(KHh(i,j,k+1)+KHh(im1(i),j,k+1)+            &
                                KHh(i,j,k)+KHh(im1(i),j,k))
         KHh_hat_j = p25*(KHh(i,j,k+1)+KHh(i,jm1(j),k+1)+            &
                                KHh(i,j,k)+KHh(i,jm1(j),k))
         H_th_1(i,j,k) = fpt1*KHh_hat_i*invdx*(th_l2(i,j,k,n4)-      &
            th_l2(im1(i),j,k,n4)-(phi_l2(i,j,k,n4)-                  &
            phi_l2(im1(i),j,k,n4))*(th_i(i,j,kp1)-th_i(i,j,kp0))/    &
            (phi_i(i,j,kp1)-phi_i(i,j,kp0)))
         H_th_2(i,j,k) = fpt2*KHh_hat_j*invdy*(th_l2(i,j,k,n4)-      &
            th_l2(i,jm1(j),k,n4)-(phi_l2(i,j,k,n4)-                  &
            phi_l2(i,jm1(j),k,n4))*(th_j(i,j,kp1)-th_j(i,j,kp0))/    &
            (phi_j(i,j,kp1)-phi_j(i,j,kp0)))
         H_qvpc_1(i,j,k) = fpt1*KHh_hat_i*invdx*(qvpc_l2(i,j,k)-        &
            qvpc_l2(im1(i),j,k)-(phi_l2(i,j,k,n4)-                      &
            phi_l2(im1(i),j,k,n4))*(qvpc_i(i,j,kp1)-qvpc_i(i,j,kp0))/   &
            (phi_i(i,j,kp1)-phi_i(i,j,kp0)))
         H_qvpc_2(i,j,k) = fpt2*KHh_hat_j*invdy*(qvpc_l2(i,j,k)-        &
            qvpc_l2(i,jm1(j),k)-(phi_l2(i,j,k,n4)-                      &
            phi_l2(i,jm1(j),k,n4))*(qvpc_j(i,j,kp1)-qvpc_j(i,j,kp0))/   &
            (phi_j(i,j,kp1)-phi_j(i,j,kp0)))
         H_qr_1(i,j,k) = fpt1*KHh_hat_i*invdx*(qr_l2(i,j,k)-         &
            qr_l2(im1(i),j,k)-(phi_l2(i,j,k,n4)-                     &
            phi_l2(im1(i),j,k,n4))*(qr_i(i,j,kp1)-qr_i(i,j,kp0))/    &
            (phi_i(i,j,kp1)-phi_i(i,j,kp0)))
         H_qr_2(i,j,k) = fpt2*KHh_hat_j*invdy*(qr_l2(i,j,k)-         &
            qr_l2(i,jm1(j),k)-(phi_l2(i,j,k,n4)-                     &
            phi_l2(i,jm1(j),k,n4))*(qr_j(i,j,kp1)-qr_j(i,j,kp0))/    &
            (phi_j(i,j,kp1)-phi_j(i,j,kp0)))
         H_qc_1(i,j,k) = fpt1*KHh_hat_i*invdx*(qc_l2(i,j,k)-         &
            qc_l2(im1(i),j,k)-(phi_l2(i,j,k,n4)-                     &
            phi_l2(im1(i),j,k,n4))*(qc_i(i,j,kp1)-qc_i(i,j,kp0))/    &
            (phi_i(i,j,kp1)-phi_i(i,j,kp0)))
         H_qc_2(i,j,k) = fpt2*KHh_hat_j*invdy*(qc_l2(i,j,k)-         &
            qc_l2(i,jm1(j),k)-(phi_l2(i,j,k,n4)-                     &
            phi_l2(i,jm1(j),k,n4))*(qc_j(i,j,kp1)-qc_j(i,j,kp0))/    &
            (phi_j(i,j,kp1)-phi_j(i,j,kp0)))

      end do
   end do

end do


!---------------------------------------------------------------------------
! model top (k = nlm)
!---------------------------------------------------------------------------

k = nlm

! "No slip" upper boundary condition
tau_13(:,:,k) = c0
tau_23(:,:,k) = c0
tau_31(:,:,k) = c0
tau_32(:,:,k) = c0

do j = 1, jm
   do i = 1, im

      ! Reynolds stresses
      tau_11(i,j,k) = rho(i,j,k)*Kmh(i,j,k)*D_11(i,j,k)
      tau_22(i,j,k) = rho(i,j,k)*Kmh(i,j,k)*D_22(i,j,k)
      tau_33(i,j,k) = rho(i,j,k)*Kmv(i,j,k)*D_33(i,j,k)

      fpt4 = p25*(rho(i,j,k)+rho(im1(i),j,k)+rho(im1(i),jm1(j),k)+   &
                  rho(i,jm1(j),k))   ! rho_i_j
      tau_12(i,j,k) = fpt4*p25*(Kmh(i,j,k)+Kmh(im1(i),j,k)+          &
              Kmh(im1(i),jm1(j),k)+Kmh(i,jm1(j),k))*D_12(i,j,k)

      ! turbulent heat fluxes
      KHh_hat_i = p5*(KHh(i,j,k)+KHh(im1(i),j,k))
      KHh_hat_j = p5*(KHh(i,j,k)+KHh(i,jm1(j),k))
      fpt1 = p5*(rho(i,j,k)+rho(im1(i),j,k))  ! rho_l2_i(k=nlm)
      fpt2 = p5*(rho(i,j,k)+rho(i,jm1(j),k))  ! rho_l2_j(k=nlm)
      H_th_1(i,j,k) = fpt1*KHh_hat_i*invdx*                          &
                        (th_l2(i,j,k,n4)-th_l2(im1(i),j,k,n4))
      H_th_2(i,j,k) = fpt2*KHh_hat_j*invdy*                          &
                        (th_l2(i,j,k,n4)-th_l2(i,jm1(j),k,n4))
      H_qvpc_1(i,j,k) = fpt1*KHh_hat_i*invdx*                        &
                        (qvpc_l2(i,j,k)-qvpc_l2(im1(i),j,k))
      H_qvpc_2(i,j,k) = fpt2*KHh_hat_j*invdy*                        &
                        (qvpc_l2(i,j,k)-qvpc_l2(i,jm1(j),k))
      H_qr_1(i,j,k) = fpt1*KHh_hat_i*invdx*                          &
                        (qr_l2(i,j,k)-qr_l2(im1(i),j,k))
      H_qr_2(i,j,k) = fpt2*KHh_hat_j*invdy*                          &
                        (qr_l2(i,j,k)-qr_l2(i,jm1(j),k))
      H_qc_1(i,j,k) = fpt1*KHh_hat_i*invdx*                          &
                        (qc_l2(i,j,k)-qc_l2(im1(i),j,k))
      H_qc_2(i,j,k) = fpt2*KHh_hat_j*invdy*                          &
                        (qc_l2(i,j,k)-qc_l2(i,jm1(j),k))
   end do
end do



! tau_21 for all layers is easily calculated
tau_21 = tau_12


! Finally, calculate H_th_3 for all layers
H_th_3(:,:,1) = rho(:,:,1)*KHv(:,:,1)*(th_l2(:,:,1,n4)-              &
        th_l2(:,:,0,n4))*grav / (phi_l2(:,:,1,n4)-phis(:,:))
H_qvpc_3(:,:,1) = rho(:,:,1)*KHv(:,:,1)*(qvpc_l2(:,:,1)-             &
        qvpc_l2(:,:,0))*grav / (phi_l2(:,:,1,n4)-phis(:,:))
H_qr_3(:,:,1) = rho(:,:,1)*KHv(:,:,1)*(qr_l2(:,:,1)-                 &
        qr_l2(:,:,0))*grav / (phi_l2(:,:,1,n4)-phis(:,:))
H_qc_3(:,:,1) = rho(:,:,1)*KHv(:,:,1)*(qc_l2(:,:,1)-                 &
        qc_l2(:,:,0))*grav / (phi_l2(:,:,1,n4)-phis(:,:))
do k = 2, nlm-1
   H_th_3(:,:,k) = rho(:,:,k)*KHv(:,:,k)*(th_l2(:,:,k,n4)-           &
      th_l2(:,:,k-1,n4))*grav / (phi_l2(:,:,k,n4)-phi_l2(:,:,k-1,n4))
   H_qvpc_3(:,:,k) = rho(:,:,k)*KHv(:,:,k)*(qvpc_l2(:,:,k)-          &
      qvpc_l2(:,:,k-1))*grav / (phi_l2(:,:,k,n4)-phi_l2(:,:,k-1,n4))
   H_qr_3(:,:,k) = rho(:,:,k)*KHv(:,:,k)*(qr_l2(:,:,k)-              &
      qr_l2(:,:,k-1))*grav / (phi_l2(:,:,k,n4)-phi_l2(:,:,k-1,n4))
   H_qc_3(:,:,k) = rho(:,:,k)*KHv(:,:,k)*(qc_l2(:,:,k)-              &
      qc_l2(:,:,k-1))*grav / (phi_l2(:,:,k,n4)-phi_l2(:,:,k-1,n4))
end do
H_th_3(:,:,nlm) = rho(:,:,nlm)*KHv(:,:,nlm)*(th_l2(:,:,nlm,n4)-      &
        th_l2(:,:,nlm-1,n4))*grav / (grav*z_top-phi_l2(:,:,nlm-1,n4))
H_qvpc_3(:,:,nlm) = rho(:,:,nlm)*KHv(:,:,nlm)*(qvpc_l2(:,:,nlm)-     &
        qvpc_l2(:,:,nlm-1))*grav / (grav*z_top-phi_l2(:,:,nlm-1,n4))
H_qr_3(:,:,nlm) = rho(:,:,nlm)*KHv(:,:,nlm)*(qr_l2(:,:,nlm)-         &
        qr_l2(:,:,nlm-1))*grav / (grav*z_top-phi_l2(:,:,nlm-1,n4))
H_qc_3(:,:,nlm) = rho(:,:,nlm)*KHv(:,:,nlm)*(qc_l2(:,:,nlm)-         &
        qc_l2(:,:,nlm-1))*grav / (grav*z_top-phi_l2(:,:,nlm-1,n4))




!
! 4) Calculate the divergence of the subgrid-scale turbulent fluxes,
!    and their contributions to the momentum and potential temperature
!    tendencies.
!


!---------------------------------------------------------------------------
! start at surface (k = 0)
!---------------------------------------------------------------------------

kp0 = mod(0,2)
kp1 = mod(1,2)

do j = 1, jm
   do i = 1, im
      f1(i,j) = p5*(H_th_1(i,j,0)+H_th_1(ip1(i),j,0))  ! H_th_1_i(k=0)
      f2(i,j) = p5*(H_th_2(i,j,0)+H_th_2(i,jp1(j),0))  ! H_th_2_j(k=0)
      f3(i,j) = p5*(H_qvpc_1(i,j,0)+H_qvpc_1(ip1(i),j,0))  ! H_qvpc_1_i(k=0)
      f4(i,j) = p5*(H_qvpc_2(i,j,0)+H_qvpc_2(i,jp1(j),0))  ! H_qvpc_2_j(k=0)
      f5(i,j) = p5*(H_qr_1(i,j,0)+H_qr_1(ip1(i),j,0))  ! H_qr_1_i(k=0)
      f6(i,j) = p5*(H_qr_2(i,j,0)+H_qr_2(i,jp1(j),0))  ! H_qr_2_j(k=0)
      f7(i,j) = p5*(H_qc_1(i,j,0)+H_qc_1(ip1(i),j,0))  ! H_qc_1_i(k=0)
      f8(i,j) = p5*(H_qc_2(i,j,0)+H_qc_2(i,jp1(j),0))  ! H_qc_2_j(k=0)
      H_1_hacek_i(i,j,kp1) = p5*( f1(i,j) +                          &
                              p5*(H_th_1(i,j,1)+H_th_1(ip1(i),j,1)) )
      H_2_hacek_j(i,j,kp1) = p5*(f2(i,j) +                           &
                              p5*(H_th_2(i,j,1)+H_th_2(i,jp1(j),1)) )
      H_qvpc_1_hacek_i(i,j,kp1) = p5*( f3(i,j) +                     &
                         p5*(H_qvpc_1(i,j,1)+H_qvpc_1(ip1(i),j,1)) )
      H_qvpc_2_hacek_j(i,j,kp1) = p5*(f4(i,j) +                      &
                         p5*(H_qvpc_2(i,j,1)+H_qvpc_2(i,jp1(j),1)) )
      H_qr_1_hacek_i(i,j,kp1) = p5*( f5(i,j) +                       &
                         p5*(H_qr_1(i,j,1)+H_qr_1(ip1(i),j,1)) )
      H_qr_2_hacek_j(i,j,kp1) = p5*(f6(i,j) +                        &
                         p5*(H_qr_2(i,j,1)+H_qr_2(i,jp1(j),1)) )
      H_qc_1_hacek_i(i,j,kp1) = p5*( f7(i,j) +                       &
                         p5*(H_qc_1(i,j,1)+H_qc_1(ip1(i),j,1)) )
      H_qc_2_hacek_j(i,j,kp1) = p5*(f8(i,j) +                        &
                         p5*(H_qc_2(i,j,1)+H_qc_2(i,jp1(j),1)) )
      tau_31_hacek_i(i,j,kp1) = p25*( tau_31(i,j,0)+                 &
            tau_31(ip1(i),j,0)+tau_31(i,j,1)+tau_31(ip1(i),j,1) )
      tau_32_hacek_j(i,j,kp1) = p25*( tau_32(i,j,0)+                 &
            tau_32(i,jp1(j),0)+tau_32(i,j,1)+tau_32(i,jp1(j),1) )  
   end do
end do

do j = 1, jm
   do i = 1, im
      fpt1 = grav*H_th_3(i,j,1) - invdx*p5*(phis(ip1(i),j)-            &
         phis(im1(i),j))*(H_1_hacek_i(i,j,kp1) - f1(i,j)) - invdy*p5*  &
         (phis(i,jp1(j))-phis(i,jm1(j)))*(H_2_hacek_j(i,j,kp1)-f2(i,j))
      F_turb_th_l2(i,j,0) = ( invdx*(H_th_1(ip1(i),j,0)-H_th_1(i,j,0)) +  &
         invdy*(H_th_2(i,jp1(j),0)-H_th_2(i,j,0)) + fpt1 /                &
         (phi(i,j,1)-phis(i,j)) ) / rho(i,j,1)
      fpt1 = grav*H_qvpc_3(i,j,1) - invdx*p5*(phis(ip1(i),j)-             &
        phis(im1(i),j))*(H_qvpc_1_hacek_i(i,j,kp1) - f3(i,j)) - invdy*p5* &
        (phis(i,jp1(j))-phis(i,jm1(j)))*(H_qvpc_2_hacek_j(i,j,kp1)-f4(i,j))
      F_turb_qvpc_l2(i,j,0) = ( invdx*(H_qvpc_1(ip1(i),j,0)-              &
        H_qvpc_1(i,j,0)) + invdy*(H_qvpc_2(i,jp1(j),0)-H_qvpc_2(i,j,0)) + &
         fpt1 / (phi(i,j,1)-phis(i,j)) ) / rho(i,j,1)
      fpt1 = grav*H_qr_3(i,j,1) - invdx*p5*(phis(ip1(i),j)-             &
        phis(im1(i),j))*(H_qr_1_hacek_i(i,j,kp1) - f5(i,j)) - invdy*p5* &
        (phis(i,jp1(j))-phis(i,jm1(j)))*(H_qr_2_hacek_j(i,j,kp1)-f6(i,j))
      F_turb_qr_l2(i,j,0) = ( invdx*(H_qr_1(ip1(i),j,0)-              &
        H_qr_1(i,j,0)) + invdy*(H_qr_2(i,jp1(j),0)-H_qr_2(i,j,0)) + &
         fpt1 / (phi(i,j,1)-phis(i,j)) ) / rho(i,j,1)
      fpt1 = grav*H_qc_3(i,j,1) - invdx*p5*(phis(ip1(i),j)-             &
        phis(im1(i),j))*(H_qc_1_hacek_i(i,j,kp1) - f7(i,j)) - invdy*p5* &
        (phis(i,jp1(j))-phis(i,jm1(j)))*(H_qc_2_hacek_j(i,j,kp1)-f8(i,j))
      F_turb_qc_l2(i,j,0) = ( invdx*(H_qc_1(ip1(i),j,0)-              &
        H_qc_1(i,j,0)) + invdy*(H_qc_2(i,jp1(j),0)-H_qc_2(i,j,0)) + &
         fpt1 / (phi(i,j,1)-phis(i,j)) ) / rho(i,j,1)
   end do
end do


!---------------------------------------------------------------------------
! next level (k = 1)
!---------------------------------------------------------------------------

k = 1

kp0 = mod(k,2)
kp1 = mod(k+1,2)

do j = 1, jm
   do i = 1, im
      H_1_hacek_i(i,j,kp1) = p25*( H_th_1(i,j,k)+H_th_1(ip1(i),j,k)+   &
                               H_th_1(i,j,k+1)+H_th_1(ip1(i),j,k+1) )
      H_2_hacek_j(i,j,kp1) = p25*( H_th_2(i,j,k)+H_th_2(i,jp1(j),k)+   &
                               H_th_2(i,j,k+1)+H_th_2(i,jp1(j),k+1) )
      H_qvpc_1_hacek_i(i,j,kp1) = p25*( H_qvpc_1(i,j,k)+               &
         H_qvpc_1(ip1(i),j,k)+H_qvpc_1(i,j,k+1)+H_qvpc_1(ip1(i),j,k+1) )
      H_qvpc_2_hacek_j(i,j,kp1) = p25*( H_qvpc_2(i,j,k)+               &
         H_qvpc_2(i,jp1(j),k)+H_qvpc_2(i,j,k+1)+H_qvpc_2(i,jp1(j),k+1) )
      H_qr_1_hacek_i(i,j,kp1) = p25*( H_qr_1(i,j,k)+                   &
         H_qr_1(ip1(i),j,k)+H_qr_1(i,j,k+1)+H_qr_1(ip1(i),j,k+1) )
      H_qr_2_hacek_j(i,j,kp1) = p25*( H_qr_2(i,j,k)+                   &
         H_qr_2(i,jp1(j),k)+H_qr_2(i,j,k+1)+H_qr_2(i,jp1(j),k+1) )
      H_qc_1_hacek_i(i,j,kp1) = p25*( H_qc_1(i,j,k)+                   &
         H_qc_1(ip1(i),j,k)+H_qc_1(i,j,k+1)+H_qc_1(ip1(i),j,k+1) )
      H_qc_2_hacek_j(i,j,kp1) = p25*( H_qc_2(i,j,k)+                   &
         H_qc_2(i,jp1(j),k)+H_qc_2(i,j,k+1)+H_qc_2(i,jp1(j),k+1) )
      tau_31_hacek_i(i,j,kp1) = p25*( tau_31(i,j,k)+                   &
            tau_31(ip1(i),j,k)+tau_31(i,j,k+1)+tau_31(ip1(i),j,k+1) )
      tau_32_hacek_j(i,j,kp1) = p25*( tau_32(i,j,k)+                   &
            tau_32(i,jp1(j),k)+tau_32(i,j,k+1)+tau_32(i,jp1(j),k+1) )
      tau_11_hat_i(i,j,kp1) = p25*( tau_11(i,j,k)+tau_11(im1(i),j,k)+  &
            tau_11(i,j,k+1)+tau_11(im1(i),j,k+1) )
      tau_12_hat_j(i,j,kp1) = p25*( tau_12(i,j,k)+tau_12(i,jp1(j),k)+  &
            tau_12(i,j,k+1)+tau_12(i,jp1(j),k+1) )
      tau_21_hat_i(i,j,kp1) = p25*( tau_21(i,j,k)+tau_21(ip1(i),j,k)+  &
            tau_21(i,j,k+1)+tau_21(ip1(i),j,k+1) )
      tau_22_hat_j(i,j,kp1) = p25*( tau_22(i,j,k)+tau_22(i,jm1(j),k)+  &
            tau_22(i,j,k+1)+tau_22(i,jm1(j),k+1) )
      f1(i,j) = phi(i,j,k+1) - phi(i,j,k)              ! delta_phi
      f2(i,j) = p5*(phi_l2(ip1(i),j,k,n4)-phi_l2(im1(i),j,k,n4))
                                                       ! d_i_phi_l2_i
      f3(i,j) = p5*(phi_l2(i,jp1(j),k,n4)-phi_l2(i,jm1(j),k,n4))
                                                       ! d_j_phi_l2_j
      f4(i,j) = p25*(phi(i,j,k)+phi(im1(i),j,k)+phi(im1(i),jm1(j),k)+  &
                     phi(i,jm1(j),k))                  ! phi_i_j
      ! Note:  Here we borrow some tau_x_hat variables temporarily to
      !        represent layer-center ("non-hat") values.
      tau_11_hat_i(i,j,kp0) = p5*( tau_11(i,j,k)+tau_11(im1(i),j,k) )
      tau_12_hat_j(i,j,kp0) = p5*( tau_12(i,j,k)+tau_12(i,jp1(j),k) )
      tau_21_hat_i(i,j,kp0) = p5*( tau_21(i,j,k)+tau_21(ip1(i),j,k) )
      tau_22_hat_j(i,j,kp0) = p5*( tau_22(i,j,k)+tau_22(i,jm1(j),k) )      
   end do
end do

do j = 1, jm
   do i = 1, im

      fpt1 = grav*(H_th_3(i,j,k+1)-H_th_3(i,j,k)) - invdx*f2(i,j)*        &
        (H_1_hacek_i(i,j,kp1) - H_1_hacek_i(i,j,kp0)) - invdy*f3(i,j)*    &
        (H_2_hacek_j(i,j,kp1)-H_2_hacek_j(i,j,kp0))
      fpt4 = p5*(rho(i,j,k)+rho(i,j,k+1))  ! rho_l2
      F_turb_th_l2(i,j,k) = ( invdx*(H_th_1(ip1(i),j,k)-H_th_1(i,j,k)) +  &
        invdy*(H_th_2(i,jp1(j),k)-H_th_2(i,j,k)) + fpt1 / f1(i,j) ) / fpt4
      fpt1 = grav*(H_qvpc_3(i,j,k+1)-H_qvpc_3(i,j,k)) - invdx*f2(i,j)*    &
        (H_qvpc_1_hacek_i(i,j,kp1) - H_qvpc_1_hacek_i(i,j,kp0)) -         &
        invdy*f3(i,j)*(H_qvpc_2_hacek_j(i,j,kp1)-H_qvpc_2_hacek_j(i,j,kp0))
      F_turb_qvpc_l2(i,j,k) = ( invdx*(H_qvpc_1(ip1(i),j,k)-              &
        H_qvpc_1(i,j,k)) +invdy*(H_qvpc_2(i,jp1(j),k)-H_qvpc_2(i,j,k)) +  &
        fpt1 / f1(i,j) ) / fpt4
      fpt1 = grav*(H_qr_3(i,j,k+1)-H_qr_3(i,j,k)) - invdx*f2(i,j)*        &
        (H_qr_1_hacek_i(i,j,kp1) - H_qr_1_hacek_i(i,j,kp0)) -             &
        invdy*f3(i,j)*(H_qr_2_hacek_j(i,j,kp1)-H_qr_2_hacek_j(i,j,kp0))
      F_turb_qr_l2(i,j,k) = ( invdx*(H_qr_1(ip1(i),j,k)-                  &
        H_qr_1(i,j,k)) +invdy*(H_qr_2(i,jp1(j),k)-H_qr_2(i,j,k)) +        &
        fpt1 / f1(i,j) ) / fpt4
      fpt1 = grav*(H_qc_3(i,j,k+1)-H_qc_3(i,j,k)) - invdx*f2(i,j)*        &
        (H_qc_1_hacek_i(i,j,kp1) - H_qc_1_hacek_i(i,j,kp0)) -             &
        invdy*f3(i,j)*(H_qc_2_hacek_j(i,j,kp1)-H_qc_2_hacek_j(i,j,kp0))
      F_turb_qc_l2(i,j,k) = ( invdx*(H_qc_1(ip1(i),j,k)-                  &
        H_qc_1(i,j,k)) +invdy*(H_qc_2(i,jp1(j),k)-H_qc_2(i,j,k)) +        &
        fpt1 / f1(i,j) ) / fpt4

      fpt1 = grav*(tau_33(i,j,k+1)-tau_33(i,j,k)) - invdx*f2(i,j)*        &
        (tau_31_hacek_i(i,j,kp1)-tau_31_hacek_i(i,j,kp0)) - invdy*        &
        f3(i,j)*(tau_32_hacek_j(i,j,kp1)-tau_32_hacek_j(i,j,kp0))
      F_turb_w_l2(i,j,k) = ( invdx*(tau_31(ip1(i),j,k)-tau_31(i,j,k)) +   &
        invdy*(tau_32(i,jp1(j),k)-tau_32(i,j,k)) + fpt1 / f1(i,j) ) / fpt4

      fpt1 = invdx*(phi(i,j,k)-phi(im1(i),j,k))*                          &
                (tau_11_hat_i(i,j,kp1)-tau_11_hat_i(i,j,kp0)) +           &
             invdy*(f4(i,jp1(j))-f4(i,j))*                                &
                (tau_12_hat_j(i,j,kp1)-tau_12_hat_j(i,j,kp0))
      fpt2 = p5*(phi_l2(i,j,k,n4)+phi_l2(im1(i),j,k,n4)-                  &
                 phis(i,j)-phis(im1(i),j))
      fpt3 = p5*(phi_l2(i,j,k,n4)+phi_l2(im1(i),j,k,n4)-                  &
                 phi(i,j,k)-phi(im1(i),j,k))
      fpt4 = p5*(rho(i,j,k)+rho(im1(i),j,k))   ! rho_i
      F_turb_u(i,j,k) = ( invdx*(tau_11(i,j,k)-tau_11(im1(i),j,k)) +      &
        invdy*(tau_12(i,jp1(j),k)-tau_12(i,j,k)) + grav*(tau_13(i,j,k)-   &
        tau_13(i,j,k-1))/fpt2 - fpt1/fpt3 ) / fpt4

      fpt1 = invdx*(f4(ip1(i),j)-f4(i,j))*                                &
                (tau_21_hat_i(i,j,kp1)-tau_21_hat_i(i,j,kp0)) +           &
             invdy*(phi(i,j,k)-phi(i,jm1(j),k))*                          &
                (tau_22_hat_j(i,j,kp1)-tau_22_hat_j(i,j,kp0))
      fpt2 = p5*(phi_l2(i,j,k,n4)+phi_l2(i,jm1(j),k,n4)-                  &
                 phis(i,j)-phis(i,jm1(j)))
      fpt3 = p5*(phi_l2(i,j,k,n4)+phi_l2(i,jm1(j),k,n4)-                  &
                 phi(i,j,k)-phi(i,jm1(j),k))
      fpt4 = p5*(rho(i,j,k)+rho(i,jm1(j),k))   ! rho_j
      F_turb_v(i,j,k) = ( invdx*(tau_21(ip1(i),j,k)-tau_21(i,j,k)) +      &
        invdy*(tau_22(i,j,k)-tau_22(i,jm1(j),k)) + grav*(tau_23(i,j,k)-   &
        tau_23(i,j,k-1))/fpt2 - fpt1/fpt3 ) / fpt4

   end do
end do



!---------------------------------------------------------------------------
! move upward
!---------------------------------------------------------------------------

do k = 2, nlm-1

   kp0 = mod(k,2)
   kp1 = mod(k+1,2)
   
   do j = 1, jm
      do i = 1, im
         H_1_hacek_i(i,j,kp1) = p25*( H_th_1(i,j,k)+H_th_1(ip1(i),j,k)+   &
                                  H_th_1(i,j,k+1)+H_th_1(ip1(i),j,k+1) )
         H_2_hacek_j(i,j,kp1) = p25*( H_th_2(i,j,k)+H_th_2(i,jp1(j),k)+   &
                                  H_th_2(i,j,k+1)+H_th_2(i,jp1(j),k+1) )
         H_qvpc_1_hacek_i(i,j,kp1) = p25*( H_qvpc_1(i,j,k)+               &
            H_qvpc_1(ip1(i),j,k)+H_qvpc_1(i,j,k+1)+H_qvpc_1(ip1(i),j,k+1) )
         H_qvpc_2_hacek_j(i,j,kp1) = p25*( H_qvpc_2(i,j,k)+               &
            H_qvpc_2(i,jp1(j),k)+H_qvpc_2(i,j,k+1)+H_qvpc_2(i,jp1(j),k+1) )
         H_qr_1_hacek_i(i,j,kp1) = p25*( H_qr_1(i,j,k)+                   &
            H_qr_1(ip1(i),j,k)+H_qr_1(i,j,k+1)+H_qr_1(ip1(i),j,k+1) )
         H_qr_2_hacek_j(i,j,kp1) = p25*( H_qr_2(i,j,k)+                   &
            H_qr_2(i,jp1(j),k)+H_qr_2(i,j,k+1)+H_qr_2(i,jp1(j),k+1) )
         H_qc_1_hacek_i(i,j,kp1) = p25*( H_qc_1(i,j,k)+                   &
            H_qc_1(ip1(i),j,k)+H_qc_1(i,j,k+1)+H_qc_1(ip1(i),j,k+1) )
         H_qc_2_hacek_j(i,j,kp1) = p25*( H_qc_2(i,j,k)+                   &
            H_qc_2(i,jp1(j),k)+H_qc_2(i,j,k+1)+H_qc_2(i,jp1(j),k+1) )
         tau_31_hacek_i(i,j,kp1) = p25*( tau_31(i,j,k)+                   &
                 tau_31(ip1(i),j,k)+tau_31(i,j,k+1)+tau_31(ip1(i),j,k+1) )
         tau_32_hacek_j(i,j,kp1) = p25*( tau_32(i,j,k)+                   &
                 tau_32(i,jp1(j),k)+tau_32(i,j,k+1)+tau_32(i,jp1(j),k+1) )
         tau_11_hat_i(i,j,kp1) = p25*( tau_11(i,j,k)+tau_11(im1(i),j,k)+  &
               tau_11(i,j,k+1)+tau_11(im1(i),j,k+1) )
         tau_12_hat_j(i,j,kp1) = p25*( tau_12(i,j,k)+tau_12(i,jp1(j),k)+  &
               tau_12(i,j,k+1)+tau_12(i,jp1(j),k+1) )
         tau_21_hat_i(i,j,kp1) = p25*( tau_21(i,j,k)+tau_21(ip1(i),j,k)+  &
               tau_21(i,j,k+1)+tau_21(ip1(i),j,k+1) )
         tau_22_hat_j(i,j,kp1) = p25*( tau_22(i,j,k)+tau_22(i,jm1(j),k)+  &
               tau_22(i,j,k+1)+tau_22(i,jm1(j),k+1) )
         f1(i,j) = phi(i,j,k+1) - phi(i,j,k)  ! delta_phi
         f2(i,j) = p5*(phi_l2(ip1(i),j,k,n4)-phi_l2(im1(i),j,k,n4))
                                                           ! d_i_phi_l2_i
         f3(i,j) = p5*(phi_l2(i,jp1(j),k,n4)-phi_l2(i,jm1(j),k,n4))
                                                           ! d_j_phi_l2_j
         f4(i,j) = p25*(phi(i,j,k)+phi(im1(i),j,k)+phi(im1(i),jm1(j),k)+  &
                        phi(i,jm1(j),k))                  ! phi_i_j
      end do
   end do

   do j = 1, jm
      do i = 1, im

         fpt1 = grav*(H_th_3(i,j,k+1)-H_th_3(i,j,k)) - invdx*f2(i,j)*     &
           (H_1_hacek_i(i,j,kp1) - H_1_hacek_i(i,j,kp0)) -                &
           invdy*f3(i,j)*(H_2_hacek_j(i,j,kp1)-H_2_hacek_j(i,j,kp0))
         fpt4 = p5*(rho(i,j,k)+rho(i,j,k+1))  ! rho_l2
         F_turb_th_l2(i,j,k) = ( invdx*(H_th_1(ip1(i),j,k)-               &
           H_th_1(i,j,k)) + invdy*(H_th_2(i,jp1(j),k)-H_th_2(i,j,k)) +    &
           fpt1 / f1(i,j) ) / fpt4
         fpt1 = grav*(H_qvpc_3(i,j,k+1)-H_qvpc_3(i,j,k)) - invdx*f2(i,j)* &
           (H_qvpc_1_hacek_i(i,j,kp1) - H_qvpc_1_hacek_i(i,j,kp0)) -      &
           invdy*f3(i,j)*(H_qvpc_2_hacek_j(i,j,kp1)-                      &
              H_qvpc_2_hacek_j(i,j,kp0))
         F_turb_qvpc_l2(i,j,k) = ( invdx*(H_qvpc_1(ip1(i),j,k)-           &
            H_qvpc_1(i,j,k)) +invdy*(H_qvpc_2(i,jp1(j),k)-                &
            H_qvpc_2(i,j,k)) + fpt1 / f1(i,j) ) / fpt4
         fpt1 = grav*(H_qr_3(i,j,k+1)-H_qr_3(i,j,k)) - invdx*f2(i,j)*     &
            (H_qr_1_hacek_i(i,j,kp1) - H_qr_1_hacek_i(i,j,kp0)) -         &
            invdy*f3(i,j)*(H_qr_2_hacek_j(i,j,kp1)-H_qr_2_hacek_j(i,j,kp0))
         F_turb_qr_l2(i,j,k) = ( invdx*(H_qr_1(ip1(i),j,k)-               &
            H_qr_1(i,j,k)) +invdy*(H_qr_2(i,jp1(j),k)-H_qr_2(i,j,k)) +    &
            fpt1 / f1(i,j) ) / fpt4
         fpt1 = grav*(H_qc_3(i,j,k+1)-H_qc_3(i,j,k)) - invdx*f2(i,j)*     &
            (H_qc_1_hacek_i(i,j,kp1) - H_qc_1_hacek_i(i,j,kp0)) -         &
            invdy*f3(i,j)*(H_qc_2_hacek_j(i,j,kp1)-H_qc_2_hacek_j(i,j,kp0))
         F_turb_qc_l2(i,j,k) = ( invdx*(H_qc_1(ip1(i),j,k)-               &
            H_qc_1(i,j,k)) +invdy*(H_qc_2(i,jp1(j),k)-H_qc_2(i,j,k)) +    &
            fpt1 / f1(i,j) ) / fpt4

         fpt1 = grav*(tau_33(i,j,k+1)-tau_33(i,j,k)) - invdx*f2(i,j)*     &
           (tau_31_hacek_i(i,j,kp1)-tau_31_hacek_i(i,j,kp0)) - invdy*     &
           f3(i,j)*(tau_32_hacek_j(i,j,kp1)-tau_32_hacek_j(i,j,kp0))
         F_turb_w_l2(i,j,k) = ( invdx*(tau_31(ip1(i),j,k)-                &
           tau_31(i,j,k)) + invdy*(tau_32(i,jp1(j),k)-tau_32(i,j,k)) +    &
           fpt1 / f1(i,j) ) / fpt4

         fpt1 = grav*(tau_13(i,j,k)-tau_13(i,j,k-1)) -                    &
                invdx*(phi(i,j,k)-phi(im1(i),j,k))*                       &
                   (tau_11_hat_i(i,j,kp1)-tau_11_hat_i(i,j,kp0)) -        &
                invdy*(f4(i,jp1(j))-f4(i,j))*                             &
                   (tau_12_hat_j(i,j,kp1)-tau_12_hat_j(i,j,kp0))
         fpt2 = p5*(phi_l2(i,j,k,n4)+phi_l2(im1(i),j,k,n4)-               &
                    phi_l2(i,j,k-1,n4)-phi_l2(im1(i),j,k-1,n4))
         fpt4 = p5*(rho(i,j,k)+rho(im1(i),j,k))   ! rho_i
         F_turb_u(i,j,k) = ( invdx*(tau_11(i,j,k)-tau_11(im1(i),j,k)) +   &
           invdy*(tau_12(i,jp1(j),k)-tau_12(i,j,k)) + fpt1/fpt2 ) / fpt4

         fpt1 = grav*(tau_23(i,j,k)-tau_23(i,j,k-1))  -                   &
                invdx*(f4(ip1(i),j)-f4(i,j))*                             &
                   (tau_21_hat_i(i,j,kp1)-tau_21_hat_i(i,j,kp0)) -        &
                invdy*(phi(i,j,k)-phi(i,jm1(j),k))*                       &
                   (tau_22_hat_j(i,j,kp1)-tau_22_hat_j(i,j,kp0))
         fpt2 = p5*(phi_l2(i,j,k,n4)+phi_l2(i,jm1(j),k,n4)-               &
                    phi_l2(i,j,k-1,n4)-phi_l2(i,jm1(j),k-1,n4))
         fpt4 = p5*(rho(i,j,k)+rho(i,jm1(j),k))   ! rho_j
         F_turb_v(i,j,k) = ( invdx*(tau_21(ip1(i),j,k)-tau_21(i,j,k)) +   &
           invdy*(tau_22(i,j,k)-tau_22(i,jm1(j),k)) + fpt1/fpt2 ) / fpt4

      end do
   end do

end do



!---------------------------------------------------------------------------
! model top (k = nlm)
!---------------------------------------------------------------------------

k = nlm

kp0 = mod(k,2)
kp1 = mod(k+1,2)

do j = 1, jm
   do i = 1, im
      ! Note:  Here we borrow some tau_x_hat variables temporarily to
      !        represent layer-center ("non-hat") values.
      tau_11_hat_i(i,j,kp1) = p5*( tau_11(i,j,k)+tau_11(im1(i),j,k) )
      tau_12_hat_j(i,j,kp1) = p5*( tau_12(i,j,k)+tau_12(i,jp1(j),k) )
      tau_21_hat_i(i,j,kp1) = p5*( tau_21(i,j,k)+tau_21(ip1(i),j,k) )
      tau_22_hat_j(i,j,kp1) = p5*( tau_22(i,j,k)+tau_22(i,jm1(j),k) )
      f4(i,j) = p25*(phi(i,j,k)+phi(im1(i),j,k)+phi(im1(i),jm1(j),k)+  &
                     phi(i,jm1(j),k))                  ! phi_i_j
   end do
end do

do j = 1, jm
   do i = 1, im
      F_turb_th_l2(i,j,k) = ( invdx*(H_th_1(ip1(i),j,k)-H_th_1(i,j,k)) +  &
        invdy*(H_th_2(i,jp1(j),k)-H_th_2(i,j,k)) - grav*H_th_3(i,j,k) /   &
        (grav*z_top-phi(i,j,k)) ) / rho(i,j,k)
      F_turb_qvpc_l2(i,j,k) = ( invdx*(H_qvpc_1(ip1(i),j,k)-              &
        H_qvpc_1(i,j,k)) + invdy*(H_qvpc_2(i,jp1(j),k)-H_qvpc_2(i,j,k)) - &
        grav*H_qvpc_3(i,j,k) / (grav*z_top-phi(i,j,k)) ) / rho(i,j,k)
      F_turb_qr_l2(i,j,k) = ( invdx*(H_qr_1(ip1(i),j,k)-                  &
        H_qr_1(i,j,k)) + invdy*(H_qr_2(i,jp1(j),k)-H_qr_2(i,j,k)) -       &
        grav*H_qr_3(i,j,k) / (grav*z_top-phi(i,j,k)) ) / rho(i,j,k)
      F_turb_qc_l2(i,j,k) = ( invdx*(H_qc_1(ip1(i),j,k)-                  &
        H_qc_1(i,j,k)) + invdy*(H_qc_2(i,jp1(j),k)-H_qc_2(i,j,k)) -       &
        grav*H_qc_3(i,j,k) / (grav*z_top-phi(i,j,k)) ) / rho(i,j,k)

      fpt1 = invdx*(phi(i,j,k)-phi(im1(i),j,k))*                          &
                (tau_11_hat_i(i,j,kp1)-tau_11_hat_i(i,j,kp0)) +           &
             invdy*(f4(i,jp1(j))-f4(i,j))*                                &
                (tau_12_hat_j(i,j,kp1)-tau_12_hat_j(i,j,kp0))
      fpt2 = grav*z_top - p5*(phi_l2(i,j,k-1,n4)+phi_l2(im1(i),j,k-1,n4))
      fpt3 = p5*(phi(i,j,k)+phi(im1(i),j,k)-                              &
                 phi_l2(i,j,k-1,n4)-phi_l2(im1(i),j,k-1,n4))
      fpt4 = p5*(rho(i,j,k)+rho(im1(i),j,k))   ! rho_i
      F_turb_u(i,j,k) = ( invdx*(tau_11(i,j,k)-tau_11(im1(i),j,k)) +      &
        invdy*(tau_12(i,jp1(j),k)-tau_12(i,j,k)) + grav*(tau_13(i,j,k)-   &
        tau_13(i,j,k-1))/fpt2 - fpt1/fpt3 ) / fpt4

      fpt1 = invdx*(f4(ip1(i),j)-f4(i,j))*                                &
                (tau_21_hat_i(i,j,kp1)-tau_21_hat_i(i,j,kp0)) +           &
             invdy*(phi(i,j,k)-phi(i,jm1(j),k))*                          &
                (tau_22_hat_j(i,j,kp1)-tau_22_hat_j(i,j,kp0))
      fpt2 = grav*z_top - p5*(phi_l2(i,j,k-1,n4)+phi_l2(i,jm1(j),k-1,n4))
      fpt3 = p5*(phi(i,j,k)+phi(i,jm1(j),k)-                              &
                 phi_l2(i,j,k-1,n4)-phi_l2(i,jm1(j),k-1,n4))
      fpt4 = p5*(rho(i,j,k)+rho(i,jm1(j),k))   ! rho_j
      F_turb_v(i,j,k) = ( invdx*(tau_21(ip1(i),j,k)-tau_21(i,j,k)) +      &
        invdy*(tau_22(i,j,k)-tau_22(i,jm1(j),k)) + grav*(tau_23(i,j,k)-   &
        tau_23(i,j,k-1))/fpt2 - fpt1/fpt3 ) / fpt4

   end do
end do




end subroutine calc_subgrid_turbulent_flx_div

!======================================================================
! END OF CALC_SUBGRID_TURBULENT_FLX_DIV
!======================================================================




end module physics

module mass_tendency

!-----------------------------------------------------------------------
! PURPOSE: Calculates the tendency of pseuo-density at layer centers.
!          Takacs "third-order" scheme is used.
!
!-----------------------------------------------------------------------

use kinds
use model_parameters
use physical_parameters
use eta_coordinate
use semi_implicit_solvers


implicit none
save


contains


!======================================================================
! BEGINNING OF GET_MF_1
!======================================================================

subroutine get_mf_1 (eta_dot_l2_sigma, eta_dot_l2_AB3_part, eta_dot_l2_EF, &
                   m_l2, F_x, F_y, F_eta_l2_AB3, F_eta_l2_EF, m_f, m_f_ef, &
                   m_f_exp_trap)

!---------------------------------------------------------------------------
! PURPOSE:
!   Computes the time tendency of pseudo-density at layer centers.
!   Takacs "third-order" uncentered advection scheme used in horizontal.
!   Simple second-order centered scheme used in vertical.
!   Scheme is in flux-form and is mass-conserving.
!
!---------------------------------------------------------------------------

implicit none

!---------------------------------------------------------------------------
! INTENT IN
!---------------------------------------------------------------------------

real (kind=dbl_kind), dimension(im,jm,1:nlm-1), intent(in) ::        &
          eta_dot_l2_sigma,    &
          eta_dot_l2_AB3_part, & ! generalized vert. velocity (eta/s)
          eta_dot_l2_EF          ! AB3 and Euler-forward, respectively

real (kind=dbl_kind), dimension(im,jm,0:nlm), intent(in) ::          &
          m_l2              ! layer-edge pseudo-density (kg/m^2/eta)

real (kind = dbl_kind), dimension(im,jm,nlm), intent(in) ::          &
          F_x,            & ! "3rd-order" mass flux in x-direction
                            ! (kg/s/m/eta) (colocated with u-points)
          F_y               ! "3rd-order" mass flux in y-direction
                            ! (kg/s/m/eta) (colocated with v-points)

!---------------------------------------------------------------------------
! INTENT OUT
!---------------------------------------------------------------------------
real (kind = dbl_kind), dimension(im,jm,1:nlm-1), intent(out) ::     &
           F_eta_l2_AB3, &  ! vertical mass flux (kg/s/m^2)
           F_eta_l2_EF      ! (colocated with eta_dot_l2 points)
                            ! AB3 and Euler-forward, respectively

real (kind=dbl_kind), dimension(im,jm,nlm), intent(out) ::           &
           m_f, m_f_ef,  &  ! time tendency of m (kg/m^2/eta/s)
           m_f_exp_trap     ! AB3, Euler-forward, and explicit part
                            ! of trapezoidal scheme respectively

!---------------------------------------------------------------------------
! LOCAL
!---------------------------------------------------------------------------
integer :: k

real (kind=dbl_kind), dimension(im,jm,nlm) ::                        &
          div_horiz_m_flux,  &  ! horiz. divergence of horiz. mass flux
                                ! (kg/m^2/eta/s)

          div_vert_m_flux       ! vert. divergence of vert. mass flux
                                ! (kg/m^2/eta/s)

real (kind=dbl_kind), dimension(im,jm,0:1) ::                        &
          f1     ! working variable

integer :: km0, km1     ! indices for k level (km0 = k-0) and the k-1
                        ! level (km1 = k-1)  (they toggle between 0 and 1)



!---------------------------------------------------------------------------
! calculate horizontal divergence of horizontal mass flux
!---------------------------------------------------------------------------

do k = 1, nlm

   div_horiz_m_flux(:,:,k) = invdx*(F_x(ip1(:),:,k)-F_x(:,:,k)) +    &
                             invdy*(F_y(:,jp1(:),k)-F_y(:,:,k))

end do



!---------------------------------------------------------------------------
! calculate vertical divergence of vertical mass flux -- AB3 fluxes
!---------------------------------------------------------------------------

F_eta_l2_AB3(:,:,1) = m_l2(:,:,1) * eta_dot_l2_AB3_part(:,:,1)
div_vert_m_flux(:,:,1) = inv_d_eta(1) * F_eta_l2_AB3(:,:,1)

do k = 2, nlm-2
   F_eta_l2_AB3(:,:,k) = m_l2(:,:,k) * eta_dot_l2_AB3_part(:,:,k)
   div_vert_m_flux(:,:,k) = inv_d_eta(k) *                           &
                      ( F_eta_l2_AB3(:,:,k) - F_eta_l2_AB3(:,:,k-1) )
end do

F_eta_l2_AB3(:,:,nlm-1) = m_l2(:,:,nlm-1) * eta_dot_l2_AB3_part(:,:,nlm-1)
div_vert_m_flux(:,:,nlm-1) = inv_d_eta(nlm-1) *                      &
                ( F_eta_l2_AB3(:,:,nlm-1) - F_eta_l2_AB3(:,:,nlm-2) )

div_vert_m_flux(:,:,nlm) = - inv_d_eta(nlm) * F_eta_l2_AB3(:,:,nlm-1)


!---------------------------------------------------------------------------
! compute time tendency of m -- due to AB3 fluxes
!---------------------------------------------------------------------------

do k = 1, nlm
   m_f(:,:,k) = -div_horiz_m_flux(:,:,k) - div_vert_m_flux(:,:,k)
end do


!---------------------------------------------------------------------------
! compute time tendency of m -- due to explicit part of trapezoidal scheme
! and update F_eta_l2_AB3 to include eta_dot_l2_sigma for later tracer
! transport calculations.
!---------------------------------------------------------------------------

km0 = mod(1,2)
km1 = mod(0,2)
f1(:,:,km0) = m_l2(:,:,1) * eta_dot_l2_sigma(:,:,1)
div_vert_m_flux(:,:,1) = inv_d_eta(1) * f1(:,:,km0)
F_eta_l2_AB3(:,:,1) = F_eta_l2_AB3(:,:,1) + f1(:,:,km0)

do k = 2, nlm-2
   km0 = mod(k,2)
   km1 = mod(k-1,2)
   f1(:,:,km0) = m_l2(:,:,k) * eta_dot_l2_sigma(:,:,k)
   div_vert_m_flux(:,:,k) = inv_d_eta(k) *                           &
                      ( f1(:,:,km0) - f1(:,:,km1) )
   F_eta_l2_AB3(:,:,k) = F_eta_l2_AB3(:,:,k) + f1(:,:,km0)
end do

km0 = mod(nlm-1,2)
km1 = mod(nlm-2,2)
f1(:,:,km0) = m_l2(:,:,nlm-1) * eta_dot_l2_sigma(:,:,nlm-1)
div_vert_m_flux(:,:,nlm-1) = inv_d_eta(nlm-1) *                      &
                      ( f1(:,:,km0) - f1(:,:,km1) )
F_eta_l2_AB3(:,:,nlm-1) = F_eta_l2_AB3(:,:,nlm-1) + f1(:,:,km0)

km0 = mod(nlm,2)
km1 = mod(nlm-1,2)
div_vert_m_flux(:,:,nlm) = - inv_d_eta(nlm) * f1(:,:,km1)


!---------------------------------------------------------------------------
! compute time tendency of m -- due to explicit part of trapezoidal scheme
!---------------------------------------------------------------------------

do k = 1, nlm
   m_f_exp_trap(:,:,k) = - alpha_si*div_vert_m_flux(:,:,k)
end do



!---------------------------------------------------------------------------
! calculate vertical divergence of vertical mass flux -- EF fluxes
!---------------------------------------------------------------------------

F_eta_l2_EF(:,:,1) = m_l2(:,:,1) * eta_dot_l2_EF(:,:,1)
div_vert_m_flux(:,:,1) = inv_d_eta(1) * F_eta_l2_EF(:,:,1)

do k = 2, nlm-2
   F_eta_l2_EF(:,:,k) = m_l2(:,:,k) * eta_dot_l2_EF(:,:,k)
   div_vert_m_flux(:,:,k) = inv_d_eta(k) *                           &
                      ( F_eta_l2_EF(:,:,k) - F_eta_l2_EF(:,:,k-1) )
end do

F_eta_l2_EF(:,:,nlm-1) = m_l2(:,:,nlm-1) * eta_dot_l2_EF(:,:,nlm-1)
div_vert_m_flux(:,:,nlm-1) = inv_d_eta(nlm-1) *                      &
                ( F_eta_l2_EF(:,:,nlm-1) - F_eta_l2_EF(:,:,nlm-2) )

div_vert_m_flux(:,:,nlm) = - inv_d_eta(nlm) * F_eta_l2_EF(:,:,nlm-1)


!---------------------------------------------------------------------------
! compute time tendency of m -- due to EF fluxes
!---------------------------------------------------------------------------

do k = 1, nlm
   m_f_ef(:,:,k) = - div_vert_m_flux(:,:,k)
end do



end subroutine get_mf_1

!======================================================================
! END OF GET_MF_1
!======================================================================





!======================================================================
! BEGINNING OF GET_MF_2
!======================================================================

subroutine get_mf_2 (eta_dot_l2_AB3, eta_dot_l2_EF, m_l2, w3,        &
                     F_eta_l2_AB3, F_eta_l2_EF, m_f, m_f_temp, m_f_ef)

!---------------------------------------------------------------------------
! PURPOSE:
!   Computes the time tendency of pseudo-density at layer centers.
!   Continuation of semi-implicit algorithm.
!
!---------------------------------------------------------------------------

implicit none

!---------------------------------------------------------------------------
! INTENT IN
!---------------------------------------------------------------------------

real (kind=dbl_kind), dimension(im,jm,1:nlm-1), intent(in) ::        &
          eta_dot_l2_AB3,  &  ! generalized vert. velocity (eta/s)
          eta_dot_l2_EF       ! (partial) AB3 and Euler-forward, respectively

real (kind=dbl_kind), dimension(im,jm,0:nlm), intent(in) ::          &
          m_l2              ! layer-edge pseudo-density (kg/m^2/eta)

real (kind=dbl_kind), intent(in) :: w3    ! Time stepping weight

!---------------------------------------------------------------------------
! INTENT INOUT
!---------------------------------------------------------------------------
real (kind=dbl_kind), dimension(im,jm,1:nlm-1), intent(inout) ::     &
           F_eta_l2_AB3, &  ! vertical mass flux (kg/s/m^2)
           F_eta_l2_EF      ! (colocated with eta_dot_l2 points)
                            ! AB3 and Euler-forward, respectively

real (kind=dbl_kind), dimension(im,jm,nlm), intent(inout) ::         &
           m_f              ! time tendency of m (kg/m^2/eta/s)

!---------------------------------------------------------------------------
! INTENT OUT
!---------------------------------------------------------------------------
real (kind=dbl_kind), dimension(im,jm,nlm), intent(out) ::           &
           m_f_temp,  &   ! additional AB3 time tendency of m (kg/m^2/eta/s)
           m_f_ef         ! additional EF time tendency of m (kg/m^2/eta/s)

!---------------------------------------------------------------------------
! LOCAL
!---------------------------------------------------------------------------
integer :: k

real (kind=dbl_kind), dimension(im,jm,nlm) ::                        &
          div_vert_m_flux       ! vert. divergence of vert. mass flux
                                ! (kg/m^2/eta/s)

real (kind=dbl_kind), dimension(im,jm,0:1) ::                        &
          f1     ! working variable

integer :: km0, km1     ! indices for k level (km0 = k-0) and the k-1
                        ! level (km1 = k-1)  (they toggle between 0 and 1)



!---------------------------------------------------------------------------
! update vertical divergence of vertical mass flux -- AB3 fluxes
!---------------------------------------------------------------------------

km0 = mod(1,2)
km1 = mod(0,2)
f1(:,:,km0) = m_l2(:,:,1) * eta_dot_l2_AB3(:,:,1)
div_vert_m_flux(:,:,1) = inv_d_eta(1) * f1(:,:,km0)
F_eta_l2_AB3(:,:,1) = F_eta_l2_AB3(:,:,1) + f1(:,:,km0)

do k = 2, nlm-2
   km0 = mod(k,2)
   km1 = mod(k-1,2)
   f1(:,:,km0) = m_l2(:,:,k) * eta_dot_l2_AB3(:,:,k)
   div_vert_m_flux(:,:,k) = inv_d_eta(k) *                           &
                      ( f1(:,:,km0) - f1(:,:,km1) )
   F_eta_l2_AB3(:,:,k) = F_eta_l2_AB3(:,:,k) + f1(:,:,km0)
end do

km0 = mod(nlm-1,2)
km1 = mod(nlm-2,2)
f1(:,:,km0) = m_l2(:,:,nlm-1) * eta_dot_l2_AB3(:,:,nlm-1)
div_vert_m_flux(:,:,nlm-1) = inv_d_eta(nlm-1) *                      &
                      ( f1(:,:,km0) - f1(:,:,km1) )
F_eta_l2_AB3(:,:,nlm-1) = F_eta_l2_AB3(:,:,nlm-1) + f1(:,:,km0)

km0 = mod(nlm,2)
km1 = mod(nlm-1,2)
div_vert_m_flux(:,:,nlm) = - inv_d_eta(nlm) * f1(:,:,km1)


!---------------------------------------------------------------------------
! update AB3 time tendencies of m
!---------------------------------------------------------------------------

m_f_temp(:,:,:) = - div_vert_m_flux(:,:,:)
m_f(:,:,:) = m_f(:,:,:) - div_vert_m_flux(:,:,:)




!---------------------------------------------------------------------------
! update vertical divergence of vertical mass flux -- EF fluxes
!---------------------------------------------------------------------------

km0 = mod(1,2)
km1 = mod(0,2)
f1(:,:,km0) = m_l2(:,:,1) * eta_dot_l2_EF(:,:,1)
div_vert_m_flux(:,:,1) = inv_d_eta(1) * f1(:,:,km0)
F_eta_l2_EF(:,:,1) = F_eta_l2_EF(:,:,1) + f1(:,:,km0)

do k = 2, nlm-2
   km0 = mod(k,2)
   km1 = mod(k-1,2)
   f1(:,:,km0) = m_l2(:,:,k) * eta_dot_l2_EF(:,:,k)
   div_vert_m_flux(:,:,k) = inv_d_eta(k) *                           &
                      ( f1(:,:,km0) - f1(:,:,km1) )
   F_eta_l2_EF(:,:,k) = F_eta_l2_EF(:,:,k) + f1(:,:,km0)
end do

km0 = mod(nlm-1,2)
km1 = mod(nlm-2,2)
f1(:,:,km0) = m_l2(:,:,nlm-1) * eta_dot_l2_EF(:,:,nlm-1)
div_vert_m_flux(:,:,nlm-1) = inv_d_eta(nlm-1) *                      &
                      ( f1(:,:,km0) - f1(:,:,km1) )
F_eta_l2_EF(:,:,nlm-1) = F_eta_l2_EF(:,:,nlm-1) + f1(:,:,km0)

km0 = mod(nlm,2)
km1 = mod(nlm-1,2)
div_vert_m_flux(:,:,nlm) = - inv_d_eta(nlm) * f1(:,:,km1)


!---------------------------------------------------------------------------
! calculate m tendency
!---------------------------------------------------------------------------

m_f_ef(:,:,:) = - div_vert_m_flux(:,:,:)




end subroutine get_mf_2

!======================================================================
! END OF GET_MF_2
!======================================================================



end module mass_tendency
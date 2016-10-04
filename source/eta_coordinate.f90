module eta_coordinate

!--------------------------------------------------------------------------
!
!  this module defines the number of vertical levels and includes
!     a subroutine that defines the coordinate levels.
!--------------------------------------------------------------------------

use kinds
use model_parameters
use physical_parameters

implicit none
save


integer, parameter :: nlm = 40                              ! number of layers
real (kind=dbl_kind), parameter :: z_top = 20000._dbl_kind    ! model top height (m)

real (kind=dbl_kind), parameter ::       &
         theta_min = 280._dbl_kind,      &  ! "theta_min" in coord. definition (K)
         dth_dsigma_min = 0._dbl_kind, &  ! "dth_dsigma_min" in coord. definition (K)
         r_fac = 16._dbl_kind,          &  ! for hybrid coordinate
         dF_deta_cutoff = 0.7_dbl_kind, & ! min. static stab. cutoff
         r_fac_sgma = 28._dbl_kind,    &  ! for new sigma_z coordinate
         m_fac_sgma = 1.5_dbl_kind,    &  ! for new sigma_z coordinate
         z0 = 19.80198019802_dbl_kind, &  ! for new sigma_z coordinate
         H0 = z_top - z0                  ! for new sigma_z coordinate

real (kind=dbl_kind), dimension(0:nlm) :: eta_l2          ! vert. coord. (at lyr. edges) (-)
real (kind=dbl_kind), dimension(nlm)   :: eta             ! vert. coord. (at lyr. ctrs.) (-)
real (kind=dbl_kind), dimension(nlm)   :: d_eta           ! delta eta (at layer centers) (-)
real (kind=dbl_kind), dimension(0:nlm) :: d_eta_l2        ! delta eta (at layer edges)   (-)
real (kind=dbl_kind), dimension(nlm)   :: inv_d_eta       ! inverse of d_eta             (-)
real (kind=dbl_kind), dimension(0:nlm) :: inv_d_eta_l2    ! inverse of d_eta_l2          (-)


contains


!======================================================================
! BEGINNING OF INIT_VERTICAL_COORDINATE
!======================================================================

subroutine init_vertical_coordinate

!---------------------------------------------------------------------------
!  This subroutine defines the vertical eta levels at the layer edges and
!  centers and computes the eta thicknesses of each layer.
!---------------------------------------------------------------------------

implicit none

!-------------------------------------------------------------------
! LOCAL
!-------------------------------------------------------------------
integer :: i, k


open (unit = 22, file = "./data/eta.defined", action = "read", form = "formatted")

print *
print *
print *, " ETA IS DEFINED FROM DATA FILE "
do k = 0,nlm
   read (22,*) i, eta_l2(i)
   print "(I6,F25.15)", i, eta_l2(i)
end do
print *
print *

close (22)

!---------------------------------------------------------------------------
! compute eta at layer centers
!---------------------------------------------------------------------------
do k = 1, nlm
   eta(k) = p5 * ( eta_l2(k) + eta_l2(k-1) )
end  do

!---------------------------------------------------------------------------
! compute delta eta for layer centers
!---------------------------------------------------------------------------
do k = 1,nlm
   d_eta(k) = eta_l2(k) - eta_l2(k-1)
end do

!---------------------------------------------------------------------------
! compute delta eta for layer edges
!---------------------------------------------------------------------------
d_eta_l2(0) = p5 * d_eta(1)
do k = 1, nlm-1
   d_eta_l2(k) = p5 * ( d_eta(k) + d_eta(k+1) )
end do
d_eta_l2(nlm) = p5 * d_eta(nlm)

!---------------------------------------------------------------------------
! compute inverses of d_eta and d_eta_l2
!---------------------------------------------------------------------------
inv_d_eta = c1/d_eta
inv_d_eta_l2 = c1/d_eta_l2



end subroutine init_vertical_coordinate

!===========================================================================
! END OF INIT_VERTICAL_COORDINATE
!===========================================================================



end module eta_coordinate


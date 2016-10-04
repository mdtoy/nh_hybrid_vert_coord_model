module initialize

!-----------------------------------------------------------------
!   This module sets the initial conditions of the model and
!   and sets the values of the eta levels.
!-----------------------------------------------------------------

use kinds
use model_parameters
use physical_parameters
use eta_coordinate
use prognostics
use step

implicit none
save


contains


!===================================================================
! BEGINNING OF INITIALIZE_MODEL
!===================================================================

subroutine initialize_model(step_count,tau_initial)

implicit none

!-------------------------------------------------------------------
! INTENT IN 
!-------------------------------------------------------------------
integer (kind = int_kind), intent(in) :: step_count

!-------------------------------------------------------------------
! INTENT INOUT
!-------------------------------------------------------------------
real (kind = dbl_kind), intent(inout) :: tau_initial

!-------------------------------------------------------------------
! LOCAL
!-------------------------------------------------------------------
integer :: i, j, k

real (kind = dbl_kind) :: fpt1, fpt2, fpt3   ! working variables



call init_prognostics

call init_vertical_coordinate

open (unit = 40, file = "./junk_time.txt", form = "formatted")
write (40,*) "hello"
close (40)

open (unit = 20, file = "./data/ic_prog", action = "read", form = "unformatted")


! Initialize prognostic and tendency time step indices
n4 = mod(step_count+1,2) + 1
n3 = mod(step_count,2) + 1
n3_f = mod(step_count+2,3) + 1
n2_f = mod(step_count+1,3) + 1
n1_f = mod(step_count,3) + 1



! Read in initial value of time (tau)
read (20) tau_initial  ! comment out if 1st line of ic_prog doesn't contain initial time
! tau_initial = c0   ! comment out if 1st line of ic_prog contains initial time


! Read initial conditions of layer-center prognostic variables
do k = 1,nlm
   do j = 1,jm
      read (20) (u(i,j,k,n4),  i = 1,im)
      read (20) (v(i,j,k,n4),  i = 1,im)
      read (20) (m(i,j,k,n4),  i = 1,im)
   end do
end do


! Read initial conditions of layer-center prognostic variables
do k = 1,nlm-1
   do j = 1,jm
      read (20) (w_l2(i,j,k,n4),   i = 1,im)
      read (20) (phi_l2(i,j,k,n4), i = 1,im)
   end do
end do
do k = 0,nlm
   do j = 1,jm
      read (20) (th_l2(i,j,k,n4),  i = 1,im)
   end do
end do


! Initialize pseudo-density at layer edges (m_l2)
m_l2(:,:,0) = m(:,:,1,n4)
do k = 1, nlm-1
   m_l2(:,:,k) = p5 * ( m(:,:,k,n4)*d_eta(k) +                       &
                        m(:,:,k+1,n4)*d_eta(k+1) ) * inv_d_eta_l2(k)
end do
m_l2(:,:,nlm) = m(:,:,nlm,n4)


! Read initial conditions of tracers if from restart (or ic_prog) file
do k = 0,nlm
   do j = 1,jm
      read (20) (qvpc_l2(i,j,k),  i = 1,im)   ! Water vapor + cloud
                                              ! mixing ratios
   end do
end do
do k = 0,nlm
   do j = 1,jm
      read (20) (qc_l2(i,j,k),  i = 1,im)   ! Cloud mixing ratio
   end do
end do
! Temporary changes
qr_l2 = c0
! qr_l2(:,:,11) = 0.0075_dbl_kind
! do k = 12,20
!    do j = 19, 24 
!       do i = 19, 24
!          qr_l2(i,j,k) = 0.015_dbl_kind
!          qr_l2(:,:,k) = 0.015_dbl_kind
!       end do
!    end do
! end do
! qr_l2(:,:,21) = 0.0075_dbl_kind
! do k = 0,nlm
!    do j = 1,jm
!       read (20) (qr_l2(i,j,k),  i = 1,im)   ! Rain mixing ratio
!    end do
! end do
! Now calculate initial value of the prognostic variables qvpc_l2_m_l2
! and qr_l2_m_l2
qvpc_l2_m_l2(:,:,:,n4) = m_l2(:,:,:) * qvpc_l2(:,:,:)
qr_l2_m_l2(:,:,:,n4) = m_l2(:,:,:) * qr_l2(:,:,:)
print *
print *, "Tracers initialized from restart file."
print *

close (20)



! Read lower boundary info. (i.e. surface geopotential)
open (unit = 21, file = "./data/bc_surf", action = "read", form = "unformatted")

do j = 1,jm
   read (21) (phis(i,j), i = 1,im)
end do

close (21)



! Calculate model "thickness" (z_top - z_surface)
h_1(:,:) = z_top - phis(:,:) / grav
inv_h_1(:,:) = c1/h_1(:,:)



! Initialize F_set
! Report if starting from a restart file and points are
! regridded points -- this is true if F_th_sgma is not
! equal to target value eta_l2
print *
do k = 1, nlm-1
   do j = 1, jm
      do i = 1, im
         ! calculate F_th_sgma (fpt3)
         fpt1 = invgrav*(phi_l2(i,j,k,n4)-phis(i,j))/h_1(i,j)  ! sigma_l2
         fpt2 = (c1-fpt1)**r_fac
         fpt3 = theta_min * fpt2 + dth_dsigma_min * (c1-fpt1) *      &
              (c1-(c1/(r_fac+c1))*fpt2) + th_l2(i,j,k,n4) * (c1-fpt2)
         F_set(i,j,k) = fpt3    ! F_set initialized
         if ( abs(fpt3-eta_l2(k)) .gt. 1.0E-09 ) then
            print "(A34,I6,A10,I6,A10,I6,A31,ES28.15E3)",            &
               "Initial regridded point at   k =", k, "j =", j,      &
               "i =", i, "F_th_sgma - eta_l2(k) =", eta_l2(k)-fpt3
         end if
      end do
   end do
end do
print *




! Initialize in_cloud, i.e., if there is a cloud present at initial time-step
in_cloud(:,:,:) = .false.
do j = 1, jm
   do i = 1, im
      do k = 1, nlm-1
         if ( qc_l2(i,j,k).gt.1.E-12_dbl_kind )                      &
            in_cloud(i,j,k) = .true.
      end do
   end do
end do


! Initialize accumulated rain
RAIN_SUM(:,:) = c0



! Set horizontal grid point "+1/-1" indices
do i = 1,im
   ip1(i) = mod(i+im,im) + 1
   im1(i) = mod(i+im-2,im) + 1
   ip2(i) = mod(i+im+1,im) + 1
   im2(i) = mod(i+im-3,im) + 1
end do
do j = 1,jm
   jp1(j) = mod(j+jm,jm) + 1
   jm1(j) = mod(j+jm-2,jm) + 1
   jp2(j) = mod(j+jm+1,jm) + 1
   jm2(j) = mod(j+jm-3,jm) + 1
end do


end subroutine initialize_model

!===================================================================
! END OF INITIALIZE_MODEL
!===================================================================


end module initialize

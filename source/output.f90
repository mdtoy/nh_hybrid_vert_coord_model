module output

!-----------------------------------------------------------------
!   This module outputs the model data in ascii and/or binary 
!   form.
!-----------------------------------------------------------------

use kinds
use model_parameters
use physical_parameters
use eta_coordinate
use prognostics
use physics
use step


implicit none


contains



subroutine initial_output (tau)

implicit none

real (kind = dbl_kind), intent(in) :: tau    ! Time in hours

integer :: nt,i,j,k
integer :: ntm

ntm = int(tau_duration*3600._dbl_kind/(dt*out_freq))+1

write (31) ntm, im, jm, nlm
write (32) ntm, im, jm, nlm
write (33) ntm, im, jm, nlm
write (34) ntm, im, jm, nlm
write (35) ntm, im, jm, nlm
write (36) ntm, im, jm, nlm
write (37) ntm, im, jm, nlm
write (38) ntm, im, jm, nlm
write (48) ntm, im, jm, nlm
write (39) ntm, im, jm
write (71) ntm, im, jm, nlm
write (72) ntm, im, jm, nlm
write (73) ntm, im, jm, nlm
write (74) ntm, im, jm, nlm
write (31) (tau + nt*tau_duration/(ntm-1),    nt = 0,ntm-1)
write (32) (tau + nt*tau_duration/(ntm-1),    nt = 0,ntm-1)
write (33) (tau + nt*tau_duration/(ntm-1),    nt = 0,ntm-1)
write (34) (tau + nt*tau_duration/(ntm-1),    nt = 0,ntm-1)
write (35) (tau + nt*tau_duration/(ntm-1),    nt = 0,ntm-1)
write (36) (tau + nt*tau_duration/(ntm-1),    nt = 0,ntm-1)
write (37) (tau + nt*tau_duration/(ntm-1),    nt = 0,ntm-1)
write (38) (tau + nt*tau_duration/(ntm-1),    nt = 0,ntm-1)
write (48) (tau + nt*tau_duration/(ntm-1),    nt = 0,ntm-1)
write (39) (tau + nt*tau_duration/(ntm-1),    nt = 0,ntm-1)
write (71) (tau + nt*tau_duration/(ntm-1),    nt = 0,ntm-1)
write (72) (tau + nt*tau_duration/(ntm-1),    nt = 0,ntm-1)
write (73) (tau + nt*tau_duration/(ntm-1),    nt = 0,ntm-1)
write (74) (tau + nt*tau_duration/(ntm-1),    nt = 0,ntm-1)
write (31) (i*dx,                  i  = 1,im   )
write (32) (i*dx,                  i  = 1,im   )
write (33) (i*dx,                  i  = 1,im   )
write (34) (i*dx,                  i  = 1,im   )
write (35) (i*dx,                  i  = 1,im   )
write (36) (i*dx,                  i  = 1,im   )
write (37) (i*dx,                  i  = 1,im   )
write (38) (i*dx,                  i  = 1,im   )
write (48) (i*dx,                  i  = 1,im   )
write (71) (i*dx,                  i  = 1,im   )
write (72) (i*dx,                  i  = 1,im   )
write (73) (i*dx,                  i  = 1,im   )
write (74) (i*dx,                  i  = 1,im   )
write (39) (i*dx,                  i  = 1,im   )
write (31) (j*dy,                  j  = 1,jm   )
write (32) (j*dy,                  j  = 1,jm   )
write (33) (j*dy,                  j  = 1,jm   )
write (34) (j*dy,                  j  = 1,jm   )
write (35) (j*dy,                  j  = 1,jm   )
write (36) (j*dy,                  j  = 1,jm   )
write (37) (j*dy,                  j  = 1,jm   )
write (38) (j*dy,                  j  = 1,jm   )
write (48) (j*dy,                  j  = 1,jm   )
write (39) (j*dy,                  j  = 1,jm   )
write (71) (j*dy,                  j  = 1,jm   )
write (72) (j*dy,                  j  = 1,jm   )
write (73) (j*dy,                  j  = 1,jm   )
write (74) (j*dy,                  j  = 1,jm   )
write (31) (eta_l2(k),             k  = 0,nlm-1)
write (32) (eta_l2(k),             k  = 0,nlm-1)
write (33) (eta_l2(k),             k  = 0,nlm-1)
write (34) (eta_l2(k),             k  = 0,nlm-1)
write (35) (eta_l2(k),             k  = 0,nlm-1)
write (36) (eta_l2(k),             k  = 0,nlm-1)
write (37) (eta_l2(k),             k  = 0,nlm-1)
write (38) (eta_l2(k),             k  = 0,nlm-1)
write (48) (eta_l2(k),             k  = 0,nlm-1)
write (71) (eta_l2(k),             k  = 0,nlm-1)
write (72) (eta_l2(k),             k  = 0,nlm-1)
write (73) (eta_l2(k),             k  = 0,nlm-1)
write (74) (eta_l2(k),             k  = 0,nlm-1)
write (31) z_top
write (32) z_top
write (33) z_top
write (34) z_top
write (35) z_top
write (36) z_top
write (37) z_top
write (38) z_top
write (48) z_top
write (71) z_top
write (72) z_top
write (73) z_top
write (74) z_top


do j = 1,jm
   write (31) (phis(i,j)/grav,        i = 1,im)
   write (32) (phis(i,j)/grav,        i = 1,im)
   write (33) (phis(i,j)/grav,        i = 1,im)
   write (34) (phis(i,j)/grav,        i = 1,im)
   write (35) (phis(i,j)/grav,        i = 1,im)
   write (36) (phis(i,j)/grav,        i = 1,im)
   write (37) (phis(i,j)/grav,        i = 1,im)
   write (38) (phis(i,j)/grav,        i = 1,im)
   write (48) (phis(i,j)/grav,        i = 1,im)
   write (39) (phis(i,j)/grav,        i = 1,im)
   write (71) (phis(i,j)/grav,        i = 1,im)
   write (72) (phis(i,j)/grav,        i = 1,im)
   write (73) (phis(i,j)/grav,        i = 1,im)
   write (74) (phis(i,j)/grav,        i = 1,im)
end do


end subroutine initial_output



subroutine output_data(stp_cnt,tau)

implicit none


!----------------------------------------------------------------------------
! INTENT IN 
!----------------------------------------------------------------------------
integer (kind = int_kind) :: stp_cnt         ! Number of time step
real (kind = dbl_kind), intent(in) :: tau    ! Time in hours

!----------------------------------------------------------------------------
! LOCAL
!----------------------------------------------------------------------------
integer :: i,j,k

real (kind = dbl_kind), dimension (im,jm,0:nlm) ::    &
           vort_x,  &     ! x-component of relative vorticity
           vort_y         ! y-component of relative vorticity

real (kind = dbl_kind), dimension(im,jm,0:nlm) ::     &
           m_l2_n4,  &    ! current time-step value of m_l2
           RH_l2_n4       ! current time-step value of relative humidity
                          ! (for output purposes only)

! Declare moisture-related layer-edge diagnostic variables
real (kind = dbl_kind), dimension(im,jm) ::                          &
           rho_l2,    &        ! density (kg/m^3)
           T_temp_l2, &        ! temperature (K)
           esat,      &        ! saturation vapor pressure (Pa)
           qv_l2               ! vapor mixing ratio (kg/kg)




! Calculate m_l2_n4 and RH_l2_n4 for output purposes
! NOTE:  qc_l2_n4 was stored in module step

m_l2_n4(:,:,0) = m(:,:,1,n4)
do k = 1, nlm-1
   m_l2_n4(:,:,k) = p5 * ( m(:,:,k,n4)*d_eta(k) +                    &
                        m(:,:,k+1,n4)*d_eta(k+1) ) * inv_d_eta_l2(k)
end do
m_l2_n4(:,:,nlm) = m(:,:,nlm,n4)

k = 0   ! bottom layer-edge
qv_l2(:,:) = qvpc_l2_m_l2(:,:,k,n4)/m_l2_n4(:,:,k) - qc_l2_n4(:,:,k)
rho_l2(:,:) = p0_ref*inv_gas_const_R*(P_exnr_l2_0(:,:)*              &
                 inv_spec_heat_cp)**inv_kappa_mns1/((c1+invepsilon*  &
                    qv_l2(:,:))*th_l2(:,:,0,n4))
T_temp_l2(:,:) = th_l2(:,:,k,n4)**gamma*(rho_l2(:,:)*gas_const_R*    &
                   (c1+invepsilon*qv_l2(:,:))/p0_ref)**gamma_mns1
esat(:,:) = 611.2_dbl_kind*exp(17.67_dbl_kind*                       &
   (T_temp_l2(:,:)-273.15_dbl_kind)/(T_temp_l2(:,:)-29.65_dbl_kind))
RH_l2_n4(:,:,k) = qv_l2(:,:)*rho_l2(:,:)*gas_const_Rv*               &
                     T_temp_l2(:,:)/esat(:,:)

k = 1   ! next layer-edge
rho_l2(:,:) = m_l2_n4(:,:,k)*grav*d_eta_l2(k)*                       &
                 c2/(phi_l2(:,:,k+1,n4)-phis(:,:))
qv_l2(:,:) = qvpc_l2_m_l2(:,:,k,n4)/m_l2_n4(:,:,k) - qc_l2_n4(:,:,k)
T_temp_l2(:,:) = th_l2(:,:,k,n4)**gamma*(rho_l2(:,:)*gas_const_R*    &
                   (c1+invepsilon*qv_l2(:,:))/p0_ref)**gamma_mns1
esat(:,:) = 611.2_dbl_kind*exp(17.67_dbl_kind*                       &
   (T_temp_l2(:,:)-273.15_dbl_kind)/(T_temp_l2(:,:)-29.65_dbl_kind))
RH_l2_n4(:,:,k) = qv_l2(:,:)*rho_l2(:,:)*gas_const_Rv*               &
                     T_temp_l2(:,:)/esat(:,:)

do k = 2, nlm-2   ! move upward
   rho_l2(:,:) = m_l2_n4(:,:,k)*grav*d_eta_l2(k)*                    &
                    c2/(phi_l2(:,:,k+1,n4)-phi_l2(:,:,k-1,n4))
   qv_l2(:,:) = qvpc_l2_m_l2(:,:,k,n4)/m_l2_n4(:,:,k) - qc_l2_n4(:,:,k)
   T_temp_l2(:,:) = th_l2(:,:,k,n4)**gamma*(rho_l2(:,:)*gas_const_R* &
                      (c1+invepsilon*qv_l2(:,:))/p0_ref)**gamma_mns1
   esat(:,:) = 611.2_dbl_kind*exp(17.67_dbl_kind*                    &
     (T_temp_l2(:,:)-273.15_dbl_kind)/(T_temp_l2(:,:)-29.65_dbl_kind))
   RH_l2_n4(:,:,k) = qv_l2(:,:)*rho_l2(:,:)*gas_const_Rv*            &
                        T_temp_l2(:,:)/esat(:,:)
end do

k = nlm-1   ! second from top layer-edge
rho_l2(:,:) = m_l2_n4(:,:,k)*grav*d_eta_l2(k)*                       &
                 c2/(grav*z_top-phi_l2(:,:,k-1,n4))
qv_l2(:,:) = qvpc_l2_m_l2(:,:,k,n4)/m_l2_n4(:,:,k) - qc_l2_n4(:,:,k)
T_temp_l2(:,:) = th_l2(:,:,k,n4)**gamma*(rho_l2(:,:)*gas_const_R*    &
                   (c1+invepsilon*qv_l2(:,:))/p0_ref)**gamma_mns1
esat(:,:) = 611.2_dbl_kind*exp(17.67_dbl_kind*                       &
   (T_temp_l2(:,:)-273.15_dbl_kind)/(T_temp_l2(:,:)-29.65_dbl_kind))
RH_l2_n4(:,:,k) = qv_l2(:,:)*rho_l2(:,:)*gas_const_Rv*               &
                     T_temp_l2(:,:)/esat(:,:)





! Calculate x-component of relative vorticity
vort_x(:,:,0) = c0
!  Note: (p5*(phi_l2(:,:,k+1,n4)-phi_l2(:,:,k-1,n4))/grav) equals delta z
!        at layer edges
vort_x(:,:,1) = (w_l2(:,:,1,n4)-w_l2(:,jm1(:),1,n4))/dy -            &
                (v(:,:,2,n4)-v(:,:,1,n4))/                           &
                   (p5*(phi_l2(:,:,2,n4)-phis(:,:))/grav)
do k = 2, nlm-2
   vort_x(:,:,k) = (w_l2(:,:,k,n4)-w_l2(:,jm1(:),k,n4))/dy -         &
                   (v(:,:,k+1,n4)-v(:,:,k,n4))/                      &
                   (p5*(phi_l2(:,:,k+1,n4)-phi_l2(:,:,k-1,n4))/grav)
end do
vort_x(:,:,nlm-1) = (w_l2(:,:,nlm-1,n4)-w_l2(:,jm1(:),nlm-1,n4))/dy - &
                    (v(:,:,nlm,n4)-v(:,:,nlm-1,n4))/                  &
                    (p5*(grav*z_top-phi_l2(:,:,nlm-2,n4))/grav)
vort_x(:,:,nlm) = c0



! Calculate y-component of relative vorticity
vort_y(:,:,0) = c0
!  Note: (p5*(phi_l2(:,:,k+1,n4)-phi_l2(:,:,k-1,n4))/grav) equals delta z
!        at layer edges
vort_y(:,:,1) = (u(:,:,2,n4)-u(:,:,1,n4))/                           &
                   (p5*(phi_l2(:,:,2,n4)-phis(:,:))/grav) -          &
                   (w_l2(:,:,1,n4)-w_l2(im1(:),:,1,n4))/dx
do k = 2, nlm-2
   vort_y(:,:,k) = (u(:,:,k+1,n4)-u(:,:,k,n4))/                      &
                 (p5*(phi_l2(:,:,k+1,n4)-phi_l2(:,:,k-1,n4))/grav) - &
                   (w_l2(:,:,k,n4)-w_l2(im1(:),:,k,n4))/dx
end do
vort_y(:,:,nlm-1) = (u(:,:,nlm,n4)-u(:,:,nlm-1,n4))/                 &
                       (p5*(grav*z_top-phi_l2(:,:,nlm-2,n4))/grav) - &
                    (w_l2(:,:,nlm-1,n4)-w_l2(im1(:),:,nlm-1,n4))/dx
vort_y(:,:,nlm) = c0




do j = 1, jm
   write (31) (u(i,j,1,n4),           i = 1,im)
   write (31) (v(i,j,1,n4),           i = 1,im)
   write (31) (w_l2_0(i,j),           i = 1,im)
end do

do k = 1, nlm-1
   do j = 1, jm
      write (31) (u(i,j,k+1,n4),      i = 1,im)
      write (31) (v(i,j,k+1,n4),      i = 1,im)
      write (31) (w_l2(i,j,k,n4),     i = 1,im)
   end do
end do


do j = 1, jm
   write (32) (th_l2(i,j,0,n4),           i = 1,im)
   write (32) (phis(i,j),                 i = 1,im)
end do

do k = 1, nlm-1
   do j = 1,jm
      write (32) (th_l2(i,j,k,n4),           i = 1,im)
      write (32) (phi_l2(i,j,k,n4),          i = 1,im)
   end do
end do


do k = 0, nlm-1
   do j = 1,jm
      write (33) (pv(i,j,k+1),            i = 1,im)
      write (33) (vort_x(i,j,k),          i = 1,im)
      write (33) (vort_y(i,j,k),          i = 1,im)
      write (34) (pl(i,j,k+1),            i = 1,im)
      write (34) (pl(i,j,k+1)/                                       &
           (gas_const_R*T_temp(i,j,k+1)), i = 1,im)   ! rho
      write (34) (T_temp(i,j,k+1),        i = 1,im)
   end do
end do


do j = 1, jm
   write (35) (c0,                        i = 1,im)
   write (35) (c0,                        i = 1,im)
   write (35) (c0,                        i = 1,im)
   write (35) (c0,                        i = 1,im)
   write (36) (c0,                        i = 1,im)
   write (36) (c0,                        i = 1,im)
   write (36) (c0,                        i = 1,im)
end do

do k = 1, nlm-1
   do j = 1, jm
      write (35) (eta_dot_l2(i,j,k),                  i = 1,im)
      write (35) (eta_dot_l2_sigma(i,j,k),            i = 1,im)
      write (35) (eta_dot_l2_theta(i,j,k),            i = 1,im)
      write (35) (eta_dot_l2_Q_diab(i,j,k),           i = 1,im)
      write (36) (eta_dot_l2_tau(i,j,k),              i = 1,im)
      write (36) (eta_dot_l2_regrid(i,j,k),           i = 1,im)
      write (36) (eta_dot_l2_prime(i,j,k),            i = 1,im)
   end do
end do


do k = 1, nlm
   do j = 1, jm
      write (37) (m(i,j,k,n4),      i = 1,im)
   end do
end do


do k = 0, nlm-1
   do j = 1, jm
      write (38) (qT_l2(i,j,k),                           i = 1,im)
      write (38) (qc_l2_n4(i,j,k),                        i = 1,im)
      write (38) (RH_l2_n4(i,j,k),                        i = 1,im)
      write (38) (Q_diabatic(i,j,k),                      i = 1,im)
   end do
end do

do j = 1,jm
   write (48) (RAIN_RATE(i,j),                            i = 1,im)
   write (48) (RAIN_SUM(i,j),                             i = 1,im)
end do
do k = 0, nlm-1
   do j = 1, jm
      write (48) (qr_l2_m_l2(i,j,k,n4)/m_l2_n4(i,j,k),    i = 1,im)
      write (48) (Q_diab_rain_evap(i,j,k),                i = 1,im)
      write (48) (Vt(i,j,k),                              i = 1,im)
   end do
end do


do j = 1,jm
   write (39) (p0_ref*(P_exnr_l2_0(i,j)/spec_heat_cp)**(c1/kappa),   &
                                          i = 1,im)   ! surface pressure
end do


do k = 1, nlm
   do j = 1, jm
      write (71) (D_11(i,j,k),      i = 1, im)
      write (71) (D_12(i,j,k),      i = 1, im)
      write (71) (D_13(i,j,k-1),    i = 1, im)
      write (71) (D_22(i,j,k),      i = 1, im)
      write (71) (D_23(i,j,k-1),    i = 1, im)
      write (71) (D_33(i,j,k),      i = 1, im)
      write (72) (Nsq(i,j,k),       i = 1, im)
      write (72) (Ri(i,j,k),        i = 1, im)
   end do
end do


do k = 1, nlm
   do j = 1, jm
      write (73) (Kmh(i,j,k),      i = 1, im)
      write (73) (Kmv(i,j,k),      i = 1, im)
   end do
end do


do j = 1, jm
   write (74) (F_turb_u(i,j,1),           i = 1, im)
   write (74) (F_turb_v(i,j,1),           i = 1, im)
   write (74) (c0,                        i = 1, im)
   write (74) (F_turb_th_l2(i,j,0),       i = 1, im)
end do
do k = 2, nlm
   do j = 1, jm
      write (74) (F_turb_u(i,j,k),           i = 1, im)
      write (74) (F_turb_v(i,j,k),           i = 1, im)
      write (74) (F_turb_w_l2(i,j,k-1),      i = 1, im)
      write (74) (F_turb_th_l2(i,j,k-1),     i = 1, im)
   end do
end do




call calc_gmeans(tau,m_l2_n4,qc_l2_n4)

print *, "Output data has been written."


! Test for restart file write
if ( (stp_cnt .ne. 0) .and. (restart_freq .ne. 0) .and.              &
     (mod(stp_cnt,restart_freq) .eq. 0) )  then
   call write_restart(tau,m_l2_n4)
end if



end subroutine output_data





subroutine calc_gmeans (tau, m_l2_n4, qc_l2_n4)
! Calculates and outputs global mass-weighted means of various variables

implicit none

!----------------------------------------------------------------------------
! INTENT IN 
!----------------------------------------------------------------------------
real (kind = dbl_kind), intent(in) :: tau    ! Time in hours

real (kind = dbl_kind), dimension(im,jm,0:nlm), intent(in) ::        &
      m_l2_n4,    &       ! current time-step value of m_l2
      qc_l2_n4            ! current time-step value of qc_l2
                          ! (for output purposes only)

!---------------------------------------------------------------------------
! LOCAL
!---------------------------------------------------------------------------
integer :: i,j,k

real (kind = dbl_kind) :: m_bar_temp  ! temp. storage for interpolated mass

real (kind = dbl_kind) ::                                            &
      m_bar,            & ! area-weighted mean mass (kg/m^2)
      theta_bar,        & ! mass-weighted mean potential temperature (K)
      theta_vnc_bar,    & ! mass-weighted mean pot. temp. variance (K^2)
      ta_bar,           & ! mass-weighted mean temperature (K)
      phi_bar,          & ! mass-weighted mean geopotential (J/kg)
      phi_vnc_bar,      & ! mass-weighted mean geopotential variance ((J/kg)^2)
      ke_bar,           & ! mass-weighted mean kinetic energy (J/kg)
      pv_bar,           & ! mass-weighted mean potential vort. (Pa^-1 s^-1)
      pot_enstr_bar,    & ! mass-weighted mean potential enstrophy (Pa^-2 s^-2)
      total_energy_bar, & ! mass-weighted mean total energy (J/kg)
      Q_diab_htg_bar,   & ! mass-weighted mean diabatic heating rate (W/kg)
      qT_bar,           & ! mass-weighted mean total moisture mix. ratio (-)
      qc_bar,           & ! mass-weighted mean cloud-water mix. ratio (-)
      qr_bar,           & ! mass-weighted mean rain mix. ratio (-)
      Q_diab_rain_evap_bar ! mass-weighted mean diabatic heating rate (W/kg)


! Initialize
m_bar = c0
theta_bar = c0
theta_vnc_bar = c0
ta_bar = c0
phi_bar = c0
phi_vnc_bar = c0
ke_bar = c0
pv_bar = c0
pot_enstr_bar = c0
total_energy_bar = c0
Q_diab_htg_bar = c0
qT_bar = c0
qc_bar = c0
qr_bar = c0
Q_diab_rain_evap_bar = c0


! Calculate area-weighted mass and temperature
do j = 1,jm
   do i = 1,im
      do k = 1,nlm
         m_bar = m_bar + m(i,j,k,n4) * d_eta(k)
         ta_bar = ta_bar + m(i,j,k,n4) * d_eta(k) *                  &
                  p5 * ( th_l2(i,j,k,n4) + th_l2(i,j,k-1,n4) ) *     &
                  ( pl(i,j,k)/p0_ref )**kappa
      end do
   end do
end do
m_bar = m_bar / (im*jm)
ta_bar = ta_bar / ( m_bar*im*jm )


! Calculate global mass-weighted potential temperature, 
! geopotential, qT, qc, qr and kinetic energy
do j = 1,jm
   do i = 1,im
      do k = 0,nlm
         theta_bar = theta_bar + m_l2_n4(i,j,k)*d_eta_l2(k)*th_l2(i,j,k,n4)
         qT_bar = qT_bar + m_l2_n4(i,j,k)*d_eta_l2(k)*qT_l2(i,j,k)
         qc_bar = qc_bar + m_l2_n4(i,j,k)*d_eta_l2(k)*qc_l2_n4(i,j,k)
         qr_bar = qr_bar + qr_l2_m_l2(i,j,k,n4)*d_eta_l2(k)
      end do
      phi_bar = phi_bar + m_l2_n4(i,j,0)*d_eta_l2(0)*phis(i,j)
      do k = 1,nlm-1
         phi_bar = phi_bar + m_l2_n4(i,j,k)*d_eta_l2(k)*phi_l2(i,j,k,n4)
      end do
      phi_bar = phi_bar + m_l2_n4(i,j,nlm)*d_eta_l2(nlm)*grav*z_top
   end do
end do
do j = 1,jm
   do i = 1,im
      ke_bar = ke_bar + m(i,j,1,n4)*d_eta(1)*( ke_horiz(i,j,1) +        &
                  p5 * (p5*w_l2_0(i,j)**2 + p5*w_l2(i,j,1,n4)**2) )
      do k = 2,nlm-1
         ke_bar = ke_bar + m(i,j,k,n4)*d_eta(k)*( ke_horiz(i,j,k) +     &
                  p5 * (p5*w_l2(i,j,k-1,n4)**2 + p5*w_l2(i,j,k,n4)**2) )
      end do
      ke_bar = ke_bar + m(i,j,nlm,n4)*d_eta(nlm)*( ke_horiz(i,j,nlm) +  &
                  p5 * (p5*w_l2(i,j,nlm-1,n4)**2) )
   end do
end do
theta_bar = theta_bar / ( m_bar*im*jm )
phi_bar = phi_bar / ( m_bar*im*jm )
ke_bar = ke_bar / ( m_bar*im*jm )
qT_bar = qT_bar / ( m_bar*im*jm )
qc_bar = qc_bar / ( m_bar*im*jm )
qr_bar = qr_bar / ( m_bar*im*jm )


! Calculate global mass-weighted variance of potential temperature
do j = 1,jm
   do i = 1,im
      do k = 0,nlm
         theta_vnc_bar = theta_vnc_bar + m_l2_n4(i,j,k)*d_eta_l2(k)* &
                         ( th_l2(i,j,k,n4) - theta_bar )**2
      end do
   end do
end do
theta_vnc_bar = theta_vnc_bar / ( m_bar*im*jm )


! Calculate global mass-weighted variance of geopotential
do j = 1,jm
   do i = 1,im
      phi_vnc_bar = phi_vnc_bar + m_l2_n4(i,j,0)*d_eta_l2(0)*        &
                         ( phis(i,j) - phi_bar )**2
      do k = 1,nlm-1
         phi_vnc_bar = phi_vnc_bar + m_l2_n4(i,j,k)*d_eta_l2(k)*     &
                         ( phi_l2(i,j,k,n4) - phi_bar )**2
      end do
      phi_vnc_bar = phi_vnc_bar + m_l2_n4(i,j,nlm)*d_eta_l2(nlm)*    &
                         ( grav*z_top - phi_bar )**2
   end do
end do
phi_vnc_bar = phi_vnc_bar / ( m_bar*im*jm )


! Calculate global mass-weighted potential vorticity and enstrophy
do j = 1,jm
   do i = 1,im
      do k = 1,nlm
         m_bar_temp = p25*( m(i,j,k,n4) + m(im1(i),j,k,n4) +          &
                            m(im1(i),jm1(j),k,n4) + m(i,jm1(j),k,n4) )
         pv_bar = pv_bar + m_bar_temp*d_eta(k)*pv(i,j,k)
         pot_enstr_bar = pot_enstr_bar + m_bar_temp*d_eta(k)*         &
                                         p5*pv(i,j,k)**2
      end do
   end do
end do
pv_bar = pv_bar / ( m_bar*im*jm )
pot_enstr_bar = pot_enstr_bar / ( m_bar*im*jm )


! Calculate global mass-weighted mean total energy (i.e. sum of internal,
! geopotential, and kinetic energy)
total_energy_bar = spec_heat_cv*ta_bar + phi_bar + ke_bar


! Calculate global mass-weighted mean diabatic heating rate
do j = 1,jm
   do i = 1,im
      do k = 0,nlm-1
         Q_diab_htg_bar = Q_diab_htg_bar + m_l2_n4(i,j,k)*d_eta_l2(k)*  &
                             Q_diabatic(i,j,k)
         Q_diab_rain_evap_bar = Q_diab_rain_evap_bar +                  &
                  m_l2_n4(i,j,k)*d_eta_l2(k)*Q_diab_rain_evap(i,j,k)
      end do
   end do
end do
Q_diab_htg_bar = Q_diab_htg_bar / ( m_bar*im*jm )
Q_diab_rain_evap_bar = Q_diab_rain_evap_bar / ( m_bar*im*jm )


! Write global means to file
write (45, "(F24.10,F24.6,15(ES24.14E3))" )                            &
             tau, tau*3600._dbl_kind, m_bar, theta_bar, theta_vnc_bar, &
             ta_bar, phi_bar, phi_vnc_bar, ke_bar, total_energy_bar,   &
             Q_diab_htg_bar, pv_bar, pot_enstr_bar,                    &
             qT_bar, qc_bar, qr_bar, Q_diab_rain_evap_bar

print "(A9,F24.13,A26,F24.13,A15,F24.13)", "ke_bar =", ke_bar,         &
         "total_energy_bar =", total_energy_bar, "m_bar =", m_bar

end subroutine calc_gmeans




subroutine write_restart (tau, m_l2_n4)
! Writes restart files

implicit none

!----------------------------------------------------------------------------
! INTENT IN 
!----------------------------------------------------------------------------
real (kind = dbl_kind), intent(in) :: tau  ! Time in hours

real (kind = dbl_kind), dimension(im,jm,0:nlm), intent(in) ::        &
        m_l2_n4           ! current time-step value of m_l2
                          ! (for output purposes only)

!---------------------------------------------------------------------------
! LOCAL
!---------------------------------------------------------------------------
character (len=5) :: dys_str
character (len=3) :: hrs_str,mts_str,sec_str
character (len=13) :: time_label
integer :: dys, hrs, mts, sec
integer :: i,j,k


! Create time_label

dys = int(tau/24._dbl_kind)                                 ! Days
hrs = int(tau - dys*24._dbl_kind)                           ! Hours
mts = int( (tau - dys*24._dbl_kind - hrs)*60._dbl_kind )    ! Minutes
sec = nint( (tau - dys*24._dbl_kind - hrs -  &
               mts/60._dbl_kind)*3600._dbl_kind )           ! Seconds

write (unit=dys_str, fmt='(i5)') dys
write (unit=hrs_str, fmt='(i3)') hrs
write (unit=mts_str, fmt='(i3)') mts
write (unit=sec_str, fmt='(i3)') sec

dys_str = '00'//trim(adjustl(dys_str))
hrs_str = '0'//trim(adjustl(hrs_str))
mts_str = '0'//trim(adjustl(mts_str))
sec_str = '0'//trim(adjustl(sec_str))

dys_str = adjustr(dys_str)
hrs_str = adjustr(hrs_str)
mts_str = adjustr(mts_str)
sec_str = adjustr(sec_str)

time_label = dys_str(3:5)//'d'//hrs_str(2:3)//'h'//                  &
             mts_str(2:3)//'m'//sec_str(2:3)//'s'



open (unit = 20, file = './restarts/model_restart.w_tracers.'//time_label,  &
                           action = "write", form = "unformatted")


! Write value of time at restart
write (20) tau

do k = 1,nlm
   do j = 1,jm
      write (20) (u(i,j,k,n4),  i = 1,im)
      write (20) (v(i,j,k,n4),  i = 1,im)
      write (20) (m(i,j,k,n4),  i = 1,im)
   end do
end do

do k = 1,nlm-1
   do j = 1,jm
      write (20) (w_l2(i,j,k,n4),   i = 1,im)
      write (20) (phi_l2(i,j,k,n4), i = 1,im)
   end do
end do

do k = 0,nlm
   do j = 1,jm
      write (20) (th_l2(i,j,k,n4),  i = 1,im)
   end do
end do

do k = 0,nlm
   do j = 1,jm
      write (20) (qvpc_l2_m_l2(i,j,k,n4)/m_l2_n4(i,j,k),  i = 1,im)
   end do
end do

do k = 0,nlm
   do j = 1,jm
      write (20) (qc_l2_n4(i,j,k),  i = 1,im)
   end do
end do

do k = 0,nlm
   do j = 1,jm
      write (20) (qr_l2_m_l2(i,j,k,n4)/m_l2_n4(i,j,k),  i = 1,im)
   end do
end do


close (20)


print *, "Restart file has been written -- time_label = ", time_label


end subroutine write_restart




end module output

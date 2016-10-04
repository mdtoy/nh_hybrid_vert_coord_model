!=====================================================================
! PROGRAM HYB_THZ_COORD_NONHYDRO_MODEL_II
!=====================================================================

program HYB_THZ_COORD_NONHYDRO_MODEL_II

!---------------------------------------------------------------------
! PURPOSE:  Driver for a non-hydrostatic atmospheric model using a
!           hybrid theta-terrainfollowing (height-based) vertical
!           coordinate.   
!
!           This is stage II -- non-theta-conserving (but an
!           initially isentropic atmosphere remains isentropic).
!           Total energy is conserved and the "mountain torque"
!           constraint is met in the limit of pure height
!           coordinates.  Calculation of density (rho) at layer
!           edges is based on the definition of the pseudo-density
!           (m) at layer edges.    
!---------------------------------------------------------------------

!---------------------------------------------------------------------
! DESCRIPTION:  A finite-difference dry dynamical core based on the 
! 3-dimensional, non-hydrostatic, equations of motion on an f-plane
! (or beta-plane).
!
! The vertical coordinate is a hybrid theta - terrain-following,
! normalized-height coordinate.
!
! Horizontal grid is quadrilateral using cartesian coordinates.
! Periodic horizontal domain.  Staggering that of Arakawa 'C-grid'.
! 
! Adams-Bashforth 3rd order time stepping scheme is used.
!
!
! AUTHOR:  Michael Toy,  11/2006
!---------------------------------------------------------------------


use kinds
use model_parameters
use physical_parameters
use initialize
use step
use output


implicit none


integer (kind=int_kind) :: step_count

real (kind=dbl_kind) ::   &
                tau,      &           ! time in hours
                tau_end               ! run termination time (hours)
                      

!
! Initialize step_count
!
step_count = 0


!
! Initialize the model (i.e. prognostic variables, initial and
!                     boundary conditions, eta_coordinate and time)
!
call initialize_model(step_count,tau)


print "(I10,A10,F18.10,A6)", step_count, "Time =", tau, "hours"


!
! Calculate run termination time
!
tau_end = tau + tau_duration


!
! Open output files to be written in subroutine output_data.
!
open (unit = 31, file = "./output/velocity.out", form = "unformatted")
open (unit = 32, file = "./output/pot_temp_geopotential.out", form = "unformatted")
open (unit = 33, file = "./output/vorticity.out", form = "unformatted")
open (unit = 34, file = "./output/press_rho_temp.out", form = "unformatted")
open (unit = 35, file = "./output/eta_dot_l2s_1.out", form = "unformatted")
open (unit = 36, file = "./output/eta_dot_l2s_2.out", form = "unformatted")
open (unit = 37, file = "./output/pseudo_density.out", form = "unformatted")
open (unit = 38, file = "./output/moisture.out", form = "unformatted")
open (unit = 48, file = "./output/moisture.precip.out", form = "unformatted")
open (unit = 39, file = "./output/surface_pressure.out", form = "unformatted")
open (unit = 71, file = "./output/deformation_tensor.out", form = "unformatted")
open (unit = 72, file = "./output/b-v_freq_Ri_number.out", form = "unformatted")
open (unit = 73, file = "./output/Kmh_Kmv.out", form = "unformatted")
open (unit = 74, file = "./output/turb_flux_divergence.out", form = "unformatted")
open (unit = 45, file = "./output/gmeans.out")


call initial_output (tau)

write (45, "(A16,16(A24))" ) "Tau(hr)", "Tau(sec)", "m", "theta",          &
                     "theta_variance", "ta", "phi", "phi_variance",        &
                     "ke", "total_energy", "diabatic_htg_rate",            &
                     "pv", "pot_enstrophy", "qT_bar", "qc_bar", "qr_bar",  &
                     "diab_rain_evap"



!
!  Take first timestep (Euler forward)
!
call update_diagnostics(tau,w1_ef,w2_ef,w3_ef)
call output_data(step_count,tau)        ! Output initial conditions
step_count = step_count + 1
tau = tau + dt/3600._dbl_kind
call step_dynamics(step_count,w1_ef,w2_ef,w3_ef)
print "(I10,A10,F18.10,A6)", step_count, "Time =", tau, "hours"
call update_diagnostics(tau,w1_ab2,w2_ab2,w3_ab2)
if (mod(step_count,out_freq)==0) then
   call output_data(step_count,tau)
end if

!
!  Take second (AB2 forward) timestep
!
step_count = step_count + 1
tau = tau + dt/3600._dbl_kind
call step_dynamics(step_count,w1_ab2,w2_ab2,w3_ab2)
print "(I10,A10,F18.10,A6)", step_count, "Time =", tau, "hours"
call update_diagnostics(tau,w1_ab3,w2_ab3,w3_ab3)
if (mod(step_count,out_freq)==0) then
   call output_data(step_count,tau)
end if


do while (tau .lt. tau_end)
   step_count = step_count + 1
   tau = tau + dt/3600._dbl_kind
   call step_dynamics(step_count,w1_ab3,w2_ab3,w3_ab3)
   print "(I10,A10,F18.10,A6)", step_count, "Time =", tau, "hours"
   call update_diagnostics(tau,w1_ab3,w2_ab3,w3_ab3)
   if (mod(step_count,out_freq)==0) then
      call output_data(step_count,tau)
   end if
end do


print *
print *, "*******  END PROGRAM  ******* "
print *


close (31)
close (32)
close (33)
close (34)
close (35)
close (36)
close (37)
close (38)
close (48)
close (39)
close (71)
close (72)
close (73)
close (74)
close (45)


end program HYB_THZ_COORD_NONHYDRO_MODEL_II

!=====================================================================
! END OF HYB_THZ_COORD_NONHYDRO_MODEL_II
!=====================================================================

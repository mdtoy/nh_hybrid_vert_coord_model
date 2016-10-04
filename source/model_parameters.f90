module model_parameters

!-----------------------------------------------------------------
!    This module specifies the grid size, time step length, and
!    run duration.
!-----------------------------------------------------------------


use kinds


implicit none
save

!
! Set grid dimensions
!
integer, parameter :: im = 80,    &    ! number of x-direction grid points
                      jm = 80            ! number of y-direction grid points

real (kind=dbl_kind), parameter :: &
                      dx = 1000.E+00_dbl_kind,     &   ! x-direction grid spacing (m)
                      dy = 1000.E+00_dbl_kind          ! y-direction grid spacing (m)


!
! Set time step length
!
real (kind=dbl_kind), parameter :: dt = 0.5_dbl_kind      ! time step length (s)


!
! Set termination time
!
real (kind=dbl_kind), parameter ::  &
                       tau_duration = 24._dbl_kind  !  2._dbl_kind        ! run duration (hours)


!
! Set frequency of output and restarts 
! (i.e. number of timesteps per output calls and restart file creations)
!
integer, parameter :: out_freq = 600  ! 120
integer, parameter :: restart_freq = 0  ! Note:  restart_freq must 
                                        ! be a multiple of out_freq
                                        ! (0 for no restart files)


real (kind=dbl_kind), parameter :: &
                      invdx  = 1.0_dbl_kind/dx,     &   ! inverse dx
                      invdx2 = 1.0_dbl_kind/dx**2,  &   ! inverse dx**2
                      dx2 = dx**2,                  &   ! square of dx
                      invdx4 = 1.0_dbl_kind/dx**4,  &   ! inverse dx**4
                      dx4 = dx**4,                  &   ! 4th power of dx
                      invdy  = 1.0_dbl_kind/dy,     &   ! inverse dy
                      invdy2 = 1.0_dbl_kind/dy**2,  &   ! inverse dy**2
                      dy2 = dy**2,                  &   ! square of dy
                      invdy4 = 1.0_dbl_kind/dy**4,  &   ! inverse dy**4
                      dy4 = dy**4,                  &   ! 4th power of dy
                      invdx2dy2 = invdx2*invdy2,    &   ! dx**-2*dy**-2
                      invdt  = 1.0_dbl_kind/dt          ! inverse dt


integer, parameter ::   &
           ntprog = 2,  &   ! no. of time slots needed for prognostic variables
           nttend = 3,  &   ! no. of time slots needed for tendency variables
           ntmisc = 2       ! no. of time slots needed for dqc_dt calc.


integer :: n3,  &                  ! index for time step n values
                                   ! i.e. previous time step
                                   
           n4                      ! index for time step n+1 values
                                   ! i.e. current time step

integer :: n3_f, n2_f, n1_f        ! index for {n, n-1, n-2} tendencies



! Declare horizontal grid point "+n/-n" indices:
!  ip1 represents "i + 1", im1 represents "i - 1"
!  ip2 represents "i + 2", im2 represents "i - 2"
!  Periodic boundary conditions are accounted for in that grid point
!  index zero has the value im, and index im+1 gets the value 1, etc.
!  Ditto the above for "j".
!  Indices are set in subroutine initialize_model.
integer, dimension (1:im) :: ip1, im1, ip2, im2
integer, dimension (1:jm) :: jp1, jm1, jp2, jm2



!
! Specify weighting factors for time-stepping
!

! Euler forward factors
real (kind=dbl_kind), parameter ::    &
          w3_ef = 1.0_dbl_kind,       &
          w2_ef = 0.0_dbl_kind,       &
          w1_ef = 0.0_dbl_kind

! Adams-Bashworth 2nd. order factors
real (kind=dbl_kind), parameter ::    &
          w3_ab2 =  1.5_dbl_kind,     &
          w2_ab2 = -0.5_dbl_kind,     &
          w1_ab2 =  0.0_dbl_kind

! Adams-Bashworth 3rd. order factors
real (kind=dbl_kind), parameter ::    &
          w3_ab3 =  1.91666666666666666666666666666666_dbl_kind,     &
          w2_ab3 = -1.33333333333333333333333333333333_dbl_kind,     &
          w1_ab3 =  0.41666666666666666666666666666666_dbl_kind




end module model_parameters

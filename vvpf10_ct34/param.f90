! --------------------------------------------------------------------------------- !
!                                                                                   !
!     Copyright (C) 2002, Jon E. Wallevik, The Norwegian University of              !
!                         Science and Technology (NTNU).                            !
!                                                                                   !
!     This file is part of Viscometric-ViscoPlastic-Flow (VVPF).                    !
!                                                                                   !
!     Viscometric-ViscoPlastic-Flow, is free software; you can redistribute it      !
!     and/or modify it under the terms of the GNU General Public License as         !
!     published by the Free Software Foundation; either version 2 of the            !
!     License, or (at your option) any later version.                               !
!                                                                                   !
!     Viscometric-ViscoPlastic-Flow, is distributed in the hope that it will be     !
!     useful, but WITHOUT ANY WARRANTY; without even the implied warranty of        !
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General      !
!     Public License for more details.                                              !
!                                                                                   !
!     You should have received a copy of the GNU General Public License             !
!     along with Viscometric-ViscoPlastic-Flow; if not, write to the Free           !
!     Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,                !
!     MA  02111-1307  USA                                                           !
!                                                                                   !
! --------------------------------------------------------------------------------- !
! File name: param.f90 (MODULE)                                                     !
! This code defines and sets all variables of relevance, like R_i, R_o, h, dr, dz,  !
! dt, tol, tol_RMS, f_min, f_max and so forth.                                      !
! The specific values for the different parameters shown here, corresponds to the   !
! case of VHMW Na, at t = 72 min and t = 102 min, in Section 9.5 (and Section 9.4). !
! --------------------------------------------------------------------------------- !
! rho => density of the test material in kg/m^3 
! ------------
! ZERO_TIME (Section 7.11): maximum iteration time in seconds when solving the 
! elliptical problem: Pseudotransient method is used. It will either be this 
! time or the term "tol_Plastic" that will terminate the elliptical iteration.
! ------------
! "REAL_TIME": real time duration of iteration when solving the parabolic problem:
! ------------
! CALCULATE_TIME_DEPENDENT_PROBL: Inform the main routine if parabolic problem 
! should be solved.
! ------------
! tol_Newton, tol_Plastic and tol_RMS: See Equation 7.73.
! tol => Tolerance to tell when to quit the successive substitution iteration and
! when, in the end, to quit iteration for the complete elliptical problem: 
! Pseudotransient Method is used.
! ------------
! dt_Newton, dt_Plastic, dt => Time step. When solving for viscoplastic fluid, 
! then much smaller time step is required compared to when solving for 
! Newtonian fluid.
! ------------
! k_max = IDNINT(ZERO_TIME/dt_Plastic): Maximum amount of time steps 
! when solving the elliptical problem.   
! ------------
! NUMBER_OF_TIME_ITERATIONS = IDNINT(REAL_TIME/dt_Plastic)
! Number of time steps when solving the parabolic problem.
! ------------
! count_max => Maximum number of iterations for each successive substitution.
! --------------------------------------------------------------------------------- !
MODULE CONSTANTS_AND_PARAMETERS
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: WHAT_TYPE_OF_VISCOMETER,ConTec_CONSTANTS,VELOCITY_AND_TIME_ConTec,&
            BML_CONSTANTS,VELOCITY_AND_TIME_BML
CONTAINS
! ================================================================================= !
SUBROUTINE WHAT_TYPE_OF_VISCOMETER(ConTec_v4,ConTec_BML_v3)
! ------------
LOGICAL,INTENT(OUT) :: ConTec_v4,ConTec_BML_v3
! ------------
ConTec_v4     = .TRUE.
ConTec_BML_v3 = .FALSE.   
! --------------------------------------------------------------------------------- !
RETURN 
END SUBROUTINE WHAT_TYPE_OF_VISCOMETER
! ================================================================================= !
SUBROUTINE ConTec_CONSTANTS(rho,REAL_TIME,CALCULATE_TIME_DEPENDENT_PROBL,&
                            ZERO_TIME,tol_Newton,tol_Plastic,tol_RMS,&
                            dt_Plastic,dt_Newton,count_max,R_i,R_o,h1,H2,H3)

DOUBLE PRECISION,INTENT(OUT) :: rho,REAL_TIME,ZERO_TIME,tol_Newton,&
                                tol_Plastic,tol_RMS,&
                                dt_Plastic,dt_Newton,R_i,R_o,h1,H2,H3
LOGICAL,INTENT(OUT)          :: CALCULATE_TIME_DEPENDENT_PROBL
INTEGER,INTENT(OUT)          :: count_max

DOUBLE PRECISION             :: TIME_INTERVAL,f,f_min,f_max,PERC
INTEGER                      :: NUMBER_OF_POINTS
LOGICAL                      :: SMOOTH
! --------------------------------------------------------------------------------- !
! Only ZERO_TIME and REAL_TIME are used from the 
! call statement below:
CALL VELOCITY_AND_TIME_ConTec(ZERO_TIME,REAL_TIME,&
     TIME_INTERVAL,f,f_min,f_max,PERC,NUMBER_OF_POINTS,SMOOTH)
! ---------------------------------------------------------
rho = 2090D0 ! Density of the test material in kg/m^3.
! ---------------------------------------------------------
! CALCULATE_TIME_DEPENDENT_PROBL = .FALSE.
CALCULATE_TIME_DEPENDENT_PROBL = .TRUE.
! ---------------------------------------------------------
tol_Newton  = 0.5D-4  ! Used as condition for time independence in the Newtonian case. 
                      ! Also used as tolerance for the successive substitution 
                      ! (in this case it acts as a dummy variable since always two
                      ! successive steps are made for the Newtonian case).
tol_Plastic = 0.5D-5  ! For the successive substitution tolerance (Equation 7.73).
tol_RMS     = 1.0D-15 ! Condition for time independence (Equation 7.75).
! ---------------------------------------------------------
dt_Newton   = 1.0D-3
dt_Plastic  = 0.1D-4  
count_max   = 15 ! Maximum number of successive (substitution) iterations, for each 
                 ! time step k.
R_i = 0.085D0    ! =>  8.5 cm  = Inner radius of viscometer.
R_o = 0.101D0    ! => 10.1 cm  = Outer radius of viscometer.
h1  = 0.002D0    ! =>  0.2 cm  = Distance between bottom plate of viscometer and the 
                 !               lowest part of cone.
H2  = 0.130D0    ! => 13.0 cm  = Total height of inner cylinder.
H3  = 0.116D0    ! => 11.6 cm  = Height where torque is measured (from top and downward).
! --------------------------------------------------------------------------------- !
RETURN 
END SUBROUTINE ConTec_CONSTANTS
! ================================================================================= !
SUBROUTINE VELOCITY_AND_TIME_ConTec(ZERO_TIME,REAL_TIME,&
           TIME_INTERVAL,f,f_min,f_max,PERC,NUMBER_OF_POINTS,SMOOTH)

DOUBLE PRECISION,INTENT(OUT) :: ZERO_TIME,REAL_TIME,&
                                TIME_INTERVAL,f,f_min,f_max,PERC
INTEGER,INTENT(OUT)          :: NUMBER_OF_POINTS
LOGICAL,INTENT(OUT)          :: SMOOTH
! --------------------------------------------------------------------------------- !
! If the condition "(DABS(VELOCITY_kp1(i,j) - VELOCITY_k(i,j)) > tol)", in 
! "main.f90", never gets fulfilled, then it will be ZERO_TIME 
! that will terminate the pseudotransient iteration (i.e. elliptical iteration).
ZERO_TIME = 0.5D0
! ---------------------------------------------------------
! Time interval for each constant angular velocity (in seconds).
! The total time is then TIME_INTERVAL*(NUMBER_OF_POINTS + 1)
!(Only used if CALCULATE_TIME_DEPENDENT_PROBL = .TRUE.)
TIME_INTERVAL = 5.0D0
! ---------------------------------------------------------
! Total number of measuring points (up and down).
! This number must be an odd number, beginning with 3: 3,5,7,9,...
!(Only used if CALCULATE_TIME_DEPENDENT_PROBL = .TRUE.)
NUMBER_OF_POINTS = 9
! ---------------------------------------------------------
! Total time of measurements (simulation), in seconds:
!(Only used if CALCULATE_TIME_DEPENDENT_PROBL = .TRUE.)
REAL_TIME = TIME_INTERVAL*DBLE(NUMBER_OF_POINTS+1) 
! ---------------------------------------------------------
! Only used if CALCULATE_TIME_DEPENDENT_PROBL = .FALSE.
f      = 0.10D0 ! -> Rotational frequency of the outer cylinder in 1/s (or rps).
! ---------------------------------------------------------
! Only used if CALCULATE_TIME_DEPENDENT_PROBL = .TRUE.
f_min  = 0.10D0 ! -> Minimum rotational frequency of the outer cylinder in 1/s (or rps).
f_max  = 0.65D0 ! -> Maximum rotational frequency of the outer cylinder in 1/s (or rps).
PERC   = 0.18D0 !
SMOOTH = .TRUE. ! SMOOTH = .FALSE.
! --------------------------------------------------------------------------------- !
RETURN
END SUBROUTINE VELOCITY_AND_TIME_ConTec
! ================================================================================= !
SUBROUTINE BML_CONSTANTS(rho,REAL_TIME,CALCULATE_TIME_DEPENDENT_PROBL,&
                         ZERO_TIME,tol_Newton,tol_Plastic,tol_RMS,&
                         dt_Plastic,dt_Newton,count_max,R_i,R_o,h1,H2,H3)

DOUBLE PRECISION,INTENT(OUT) :: rho,REAL_TIME,ZERO_TIME,tol_Newton,&
                                tol_Plastic,tol_RMS,&
                                dt_Plastic,dt_Newton,R_i,R_o,h1,H2,H3
LOGICAL,INTENT(OUT)          :: CALCULATE_TIME_DEPENDENT_PROBL
INTEGER,INTENT(OUT)          :: count_max

DOUBLE PRECISION             :: TIME_INTERVAL,f,f_min,f_max,PERC
INTEGER                      :: NUMBER_OF_POINTS
LOGICAL                      :: SMOOTH
! --------------------------------------------------------------------------------- !
! Only ZERO_TIME and REAL_TIME are used from the 
! call statement below:
CALL VELOCITY_AND_TIME_BML(ZERO_TIME,REAL_TIME,&
     TIME_INTERVAL,f,f_min,f_max,PERC,NUMBER_OF_POINTS,SMOOTH)
! ---------------------------------------------------------
rho = 2354D0 ! kg/m^3 
! ---------------------------------------------------------
CALCULATE_TIME_DEPENDENT_PROBL = .FALSE.
! CALCULATE_TIME_DEPENDENT_PROBL = .TRUE.
! ---------------------------------------------------------
tol_Newton  = 1.0D-3  ! Used as condition for time independence in the Newtonian case. 
                      ! Also used as tolerance for the successive substitution 
                      ! (in this case it acts as a dummy variable since always two 
                      ! successive steps are made for the Newtonian case). 
tol_Plastic = 1.0D-10 ! For the successive substitution tolerance (Equation 7.73).
tol_RMS     = 1.0D-30 ! Condition for time independence (Equation 7.75).
! ---------------------------------------------------------
dt_Newton   = 1.0D-1
dt_Plastic  = 0.1D-6
count_max   = 15 ! Maximum number of successive (substitution) iterations, for each
                 ! time step k.
R_i = 0.100D0    ! => 10.0 cm  = Inner radius of viscometer.
R_o = 0.145D0    ! => 14.5 cm  = Outer radius of viscometer.
h1  = 0.020D0    ! =>  2.0 cm  = Distance between bottom plate of viscometer and the 
                 !               lowest part of cone.
H2  = 0.245D0    ! => 24.5 cm  = Total height of inner cylinder 
H3  = 0.199D0    ! => 19.9 cm  = Height where torque is measured (from top and downward).
! --------------------------------------------------------------------------------- !
RETURN 
END SUBROUTINE BML_CONSTANTS
! ================================================================================= !
SUBROUTINE VELOCITY_AND_TIME_BML(ZERO_TIME,REAL_TIME,&
           TIME_INTERVAL,f,f_min,f_max,PERC,NUMBER_OF_POINTS,SMOOTH)

DOUBLE PRECISION,INTENT(OUT) :: ZERO_TIME,REAL_TIME,&
                                TIME_INTERVAL,f,f_min,f_max,PERC
INTEGER,INTENT(OUT)          :: NUMBER_OF_POINTS
LOGICAL,INTENT(OUT)          :: SMOOTH
! --------------------------------------------------------------------------------- !
! If the condition "(DABS(VELOCITY_kp1(i,j) - VELOCITY_k(i,j)) > tol)", in 
! "main.f90", never gets fulfilled, then it will be ZERO_TIME 
! that will terminate the pseudotransient iteration (i.e. elliptical iteration).
ZERO_TIME = 1.0D0  
! ---------------------------------------------------------
! Time interval for each constant angular velocity (in seconds).
! The total time is then TIME_INTERVAL*(NUMBER_OF_POINTS + 1)
!(Only used if CALCULATE_TIME_DEPENDENT_PROBL = .TRUE.)
TIME_INTERVAL = 5.0D0 
! ---------------------------------------------------------
! Total number of measuring points (up and down).
! This number must be an odd number, beginning with 3: 3,5,7,9,...
!(Only used if CALCULATE_TIME_DEPENDENT_PROBL = .TRUE.)
NUMBER_OF_POINTS = 9 
! ---------------------------------------------------------
! Total time of measurements (simulation), in seconds:
!(Only used if CALCULATE_TIME_DEPENDENT_PROBL = .TRUE.)
REAL_TIME = TIME_INTERVAL*DBLE(NUMBER_OF_POINTS+1) 
! ---------------------------------------------------------
! Only used if CALCULATE_TIME_DEPENDENT_PROBL = .FALSE.
! f    = 0.50D0 ! -> Rotational frequency of the outer cylinder in 1/s (or rps).
f      = 3.0D0/(2.0D0*DACOS(-1.0D0)) ! -> omega = 3 rad/s (PI = DACOS(-1.0D0))
! ---------------------------------------------------------
! Only used if CALCULATE_TIME_DEPENDENT_PROBL = .TRUE.
f_min  = 0.10D0 ! -> Minimum rotational frequency of the outer cylinder in 1/s (or rps).
f_max  = 0.65D0 ! -> Maximum rotational frequency of the outer cylinder in 1/s (or rps).
PERC   = 0.18D0 ! 0.30D0       
SMOOTH = .TRUE. ! SMOOTH= .FALSE.
! --------------------------------------------------------------------------------- !
RETURN
END SUBROUTINE VELOCITY_AND_TIME_BML
! ================================================================================= !
END MODULE CONSTANTS_AND_PARAMETERS
! --------------------------------------------------------------------------------- !
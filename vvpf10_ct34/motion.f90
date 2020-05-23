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
! File name: motion.f90 (MODULE)                                                    !
! This file reads the basic information from param.f90 to produce the angular       !
! velocity "omega". The information about the angular velocity is requested by the  !
! routine main.f90.                                                                 !
! --------------------------------------------------------------------------------- !
MODULE ROTATION
  USE CONSTANTS_AND_PARAMETERS
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ANGULAR_VELOCITY
CONTAINS
! ================================================================================= !
SUBROUTINE ANGULAR_VELOCITY(double_prec_k,dt,omega)

DOUBLE PRECISION,INTENT(IN)  :: double_prec_k,dt
DOUBLE PRECISION,INTENT(OUT) :: omega

DOUBLE PRECISION             :: rho,REAL_TIME,ZERO_TIME,tol_Newton,&
                                tol_Plastic,tol_RMS,dt_Plastic,dt_Newton,&
                                R_i,R_o,h1,H2,H3
INTEGER                      :: count_max
LOGICAL                      :: ConTec_v4,ConTec_BML_v3,&
                                CALCULATE_TIME_DEPENDENT_PROBL
! --------------------------------------------------------------------------------- !
CALL WHAT_TYPE_OF_VISCOMETER(ConTec_v4,ConTec_BML_v3)
! ---------------------------------------------------------
IF ((ConTec_v4).AND.(ConTec_BML_v3)) THEN
  PRINT *, " ERROR: Both ConTec_v4 AND ConTec_BML_v3 = .TRUE.  "
  STOP
ELSE IF ((.NOT.ConTec_v4).AND.(.NOT.ConTec_BML_v3)) THEN
  PRINT *, " ERROR: Both ConTec_v4 AND ConTec_BML_v3 = .FALSE. "
  STOP
END IF 
! ---------------------------------------------------------
IF (ConTec_v4) THEN
  CALL ConTec_CONSTANTS(rho,REAL_TIME,CALCULATE_TIME_DEPENDENT_PROBL,&
                        ZERO_TIME,tol_Newton,tol_Plastic,tol_RMS,&
                        dt_Plastic,dt_Newton,count_max,R_i,R_o,h1,H2,H3)
END IF

IF (ConTec_BML_v3) THEN
  CALL BML_CONSTANTS(rho,REAL_TIME,CALCULATE_TIME_DEPENDENT_PROBL,&
                     ZERO_TIME,tol_Newton,tol_Plastic,tol_RMS,&
                     dt_Plastic,dt_Newton,count_max,R_i,R_o,h1,H2,H3)
END IF
! ---------------------------------------------------------
IF (CALCULATE_TIME_DEPENDENT_PROBL) THEN
  CALL STANDARD_OMEGA_UP_DOWN(double_prec_k,dt,omega)
ELSE
  CALL CONSTANT_OMEGA(omega)
END IF
! --------------------------------------------------------------------------------- !
RETURN
END SUBROUTINE ANGULAR_VELOCITY
! ================================================================================= !
SUBROUTINE CONSTANT_OMEGA(omega)

DOUBLE PRECISION,INTENT(OUT) :: omega

DOUBLE PRECISION             :: PI,ZERO_TIME,REAL_TIME,&
                                TIME_INTERVAL,f,f_min,f_max,PERC
INTEGER                      :: NUMBER_OF_POINTS
LOGICAL                      :: ConTec_v4,ConTec_BML_v3,SMOOTH
! --------------------------------------------------------------------------------- !
PI = DACOS(-1.0D0)
! ---------------------------------------------------------
CALL WHAT_TYPE_OF_VISCOMETER(ConTec_v4,ConTec_BML_v3)
! ---------------------------------------------------------
IF ((ConTec_v4).AND.(ConTec_BML_v3)) THEN
  PRINT *, " ERROR: Both ConTec_v4 AND ConTec_BML_v3 = .TRUE.  "
  STOP
ELSE IF ((.NOT.ConTec_v4).AND.(.NOT.ConTec_BML_v3)) THEN
  PRINT *, " ERROR: Both ConTec_v4 AND ConTec_BML_v3 = .FALSE. "
  STOP
END IF 
! ---------------------------------------------------------
IF (ConTec_v4) THEN
CALL VELOCITY_AND_TIME_ConTec(ZERO_TIME,REAL_TIME,&
     TIME_INTERVAL,f,f_min,f_max,PERC,NUMBER_OF_POINTS,SMOOTH)
END IF

IF (ConTec_BML_v3) THEN
CALL VELOCITY_AND_TIME_BML(ZERO_TIME,REAL_TIME,&
     TIME_INTERVAL,f,f_min,f_max,PERC,NUMBER_OF_POINTS,SMOOTH)
END IF
! ---------------------------------------------------------
omega = 2.0D0*PI*f
! --------------------------------------------------------------------------------- !
RETURN
END SUBROUTINE CONSTANT_OMEGA
! ================================================================================= !
SUBROUTINE STANDARD_OMEGA_UP_DOWN(double_prec_k,dt,omega)

DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: omega_POINT,TIME_POINT,&
                                             omega_POINT_correction
INTEGER,ALLOCATABLE,DIMENSION(:)          :: OMEGA_STEP_VECTOR

DOUBLE PRECISION,INTENT(IN)  :: double_prec_k,dt
DOUBLE PRECISION,INTENT(OUT) :: omega

DOUBLE PRECISION             :: PI,f,f_min,f_max,OMEGA_MIN,OMEGA_MAX,&
                                TIME_INTERVAL,REAL_TIME,OMEGA_STEP,&
                                ZERO_TIME,PERC,EPS,time,SLOPE,step,&
                                value,value_1,value_2,value_3,value_4,&
                                f_start,f_end,OMEGA_start,OMEGA_end,&
                                PERC_static_begin,PERC_static_end

INTEGER                      :: NUMBER_OF_POINTS,NT_LOCK,problem,I,J,L

LOGICAL                      :: ConTec_v4,ConTec_BML_v3,SMOOTH
! --------------------------------------------------------------------------------- !
CALL WHAT_TYPE_OF_VISCOMETER(ConTec_v4,ConTec_BML_v3)
! ---------------------------------------------------------
IF ((ConTec_v4).AND.(ConTec_BML_v3)) THEN
  PRINT *, " ERROR: Both ConTec_v4 AND ConTec_BML_v3 = .TRUE.  "
  STOP
ELSE IF ((.NOT.ConTec_v4).AND.(.NOT.ConTec_BML_v3)) THEN
  PRINT *, " ERROR: Both ConTec_v4 AND ConTec_BML_v3 = .FALSE. "
  STOP
END IF 
! ---------------------------------------------------------
IF (ConTec_v4) THEN
CALL VELOCITY_AND_TIME_ConTec(ZERO_TIME,REAL_TIME,&
     TIME_INTERVAL,f,f_min,f_max,PERC,NUMBER_OF_POINTS,SMOOTH)
END IF

IF (ConTec_BML_v3) THEN
CALL VELOCITY_AND_TIME_BML(ZERO_TIME,REAL_TIME,&
     TIME_INTERVAL,f,f_min,f_max,PERC,NUMBER_OF_POINTS,SMOOTH)
END IF
! --------------------------------------------------------------------------------- !
time   = double_prec_k*dt
PI     = DACOS(-1.0D0)
EPS    = 1.0D-9
! ---------------------------------------------------------
! Parameters that should be ported into param.f90 in future:
f_start = 0.02D0               ! rps
f_end   = (7.8D0/11.7D0)*f_max ! rps

PERC_static_begin = 0.12D0
PERC_static_end   = 0.21D0
! --------------------------------------------------------------------------------- !
ALLOCATE(omega_POINT_correction(NUMBER_OF_POINTS+1),stat=problem)
IF (problem/=0) THEN 
  PRINT *, " PROGRAM COULD NOT ALLOCATE SPACE FOR THE VECTOR "
  PRINT *, " omega_POINT_correction!                         " 
  PRINT *, " EXECUTION TERMINATED!                           "
  STOP
END IF

! Here, omega_POINT_correction is in rps.
omega_POINT_correction = (/ -0.00061363240021D0,0.00012914263445D0,&
                             0.00138781741912D0,0.00249033048894D0,&
                             0.00323310552361D0,0.00266225011894D0,&
                             0.00138781741912D0,0.00030113269445D0,&
                            -0.00061363240021D0,0.00170215702079D0/)

! omega_POINT_correction = 0.0D0
! ---------------------------------------------------------
IF (PERC_static_begin.LT.EPS) THEN
  f_start = f_min
END IF

! Converting omega_POINT_correction from rps to rad/s:
omega_POINT_correction = 2.0D0*PI*omega_POINT_correction

OMEGA_start = 2.0D0*PI*f_start ! rad/s
OMEGA_end   = 2.0D0*PI*f_end   ! rad/s
OMEGA_MIN   = 2.0D0*PI*f_min   ! rad/s
OMEGA_MAX   = 2.0D0*PI*f_max   ! rad/s
OMEGA_STEP  = (OMEGA_MAX - OMEGA_MIN)/(DBLE(NUMBER_OF_POINTS - 1)/2.0D0)
! --------------------------------------------------------------------------------- !
ALLOCATE(OMEGA_STEP_VECTOR(NUMBER_OF_POINTS),stat=problem)
IF (problem/=0) THEN 
  PRINT *, " PROGRAM COULD NOT ALLOCATE SPACE FOR THE VECTOR "
  PRINT *, " OMEGA_STEP_VECTOR!                              " 
  PRINT *, " EXECUTION TERMINATED!                           "
  STOP
END IF

OMEGA_STEP_VECTOR = 0

DO I=0,((NUMBER_OF_POINTS-1)/2),1
  OMEGA_STEP_VECTOR(I+1) = I
END DO

DO I=(((NUMBER_OF_POINTS-1)/2)-1),0,-1
  OMEGA_STEP_VECTOR(NUMBER_OF_POINTS-I) = I
END DO
! ---------------------------------------------------------
ALLOCATE(omega_POINT(2*(NUMBER_OF_POINTS+1)+1),&
         TIME_POINT(2*(NUMBER_OF_POINTS+1)+1),stat=problem)
IF (problem/=0) THEN 
  PRINT *, " PROGRAM COULD NOT ALLOCATE SPACE FOR THE VECTOR "
  PRINT *, " omega_POINT OR TIME_POINT!                      "
  PRINT *, " EXECUTION TERMINATED!                           "
  STOP
END IF

omega_POINT = 0.0D0
TIME_POINT  = 0.0D0
! --------------------------------------------------------------------------------- !
DO J=3,2*NUMBER_OF_POINTS-1,2 ! J=1,19,2
  TIME_POINT(J)   = DBLE((J+1)/2-1)*TIME_INTERVAL
  TIME_POINT(J+1) = TIME_POINT(J)+PERC*TIME_INTERVAL
  IF (J <= (NUMBER_OF_POINTS+1)) THEN
    omega_POINT(J)   = (OMEGA_MIN+OMEGA_STEP*&
                       DBLE(OMEGA_STEP_VECTOR((J+1)/2)-1))
    omega_POINT(J+1) = (OMEGA_MIN+OMEGA_STEP*&
                       DBLE(OMEGA_STEP_VECTOR((J+1)/2)))
  ELSE
    omega_POINT(J)   = (2.0D0*OMEGA_STEP+OMEGA_MIN+OMEGA_STEP*&
                       DBLE(OMEGA_STEP_VECTOR((J+1)/2)-1))
    omega_POINT(J+1) = (OMEGA_MIN+OMEGA_STEP*&
                       DBLE(OMEGA_STEP_VECTOR((J+1)/2)))
  END IF
END DO

J = 1
TIME_POINT(J)      = DBLE((J+1)/2-1)*TIME_INTERVAL
TIME_POINT(J+1)    = TIME_POINT(J) + PERC_static_begin*TIME_INTERVAL
omega_POINT(J)     = OMEGA_start
omega_POINT(J+1)   = OMEGA_MIN

J = 2*NUMBER_OF_POINTS + 1
TIME_POINT(J)      = DBLE((J+1)/2-1)*TIME_INTERVAL ! 45 sec
TIME_POINT(J+1)    = TIME_POINT(J) + PERC_static_end*TIME_INTERVAL
omega_POINT(J)     = OMEGA_MIN
omega_POINT(J+1)   = OMEGA_end

J = 2*NUMBER_OF_POINTS + 3
TIME_POINT(J)      = DBLE((J+1)/2-1)*TIME_INTERVAL
omega_POINT(J)     = OMEGA_end

DO J=2,2*NUMBER_OF_POINTS+2,2
  omega_POINT(J)   = omega_POINT(J)   + omega_POINT_correction(J/2)
  omega_POINT(J+1) = omega_POINT(J+1) + omega_POINT_correction(J/2)
END DO

! Legacy from debugging period (might come in use again):
! DO J=1,2*NUMBER_OF_POINTS+2
!   PRINT "( ' ', (F8.5,5X,F8.5) )  ",TIME_POINT(J),omega_POINT(J)
! END DO
! PRINT *, "time"
! READ *,time

DO L=1,2*NUMBER_OF_POINTS+2
  IF ((TIME_POINT(L).LE.time).AND.(TIME_POINT(L+1).GE.time)) THEN
    NT_LOCK=L
  END IF
END DO

L = NT_LOCK
! --------------------------------------------------------------------------------- !
! The following is to be used if smoothing is to be applied in only "one directions".
! EXP() - function is used.
! ---------------------------------------------------------
! When PERC = 0.17D0
value_1 = 5.8D0
value_2 = 7.0D0
value_3 = 6.5D0 
value_4 = 6.0D0
! ---------------------------------------------------------
SLOPE = (omega_POINT(L+1) - omega_POINT(L))/(TIME_POINT(L+1) - TIME_POINT(L) + EPS)

OMEGA_STEP = DABS(omega_POINT(L+1) - omega_POINT(L))

IF (SMOOTH) THEN
! ---------------------------------------------------------
  IF (DABS(SLOPE).LT.1D-5) THEN
    omega = omega_POINT(L)
  ELSE
    IF (L.LT.2) THEN 
      step  = (time - TIME_POINT(L))/(TIME_POINT(L+1) - TIME_POINT(L)) ! \in [0,1]
      value = (1.0D0 - DEXP(-step*value_1))                            ! \in [0,1]
      omega = OMEGA_STEP*value + omega_POINT(L)
    ELSE IF (L.LT.(NUMBER_OF_POINTS+1)) THEN !10
      step  = (time - TIME_POINT(L))/(TIME_POINT(L+1) - TIME_POINT(L)) ! \in [0,1]
      value = (1.0D0 - DEXP(-step*value_2))                            ! \in [0,1]
      omega = OMEGA_STEP*value + omega_POINT(L)
    ELSE IF (L.LT.(2*NUMBER_OF_POINTS)) THEN !18
      step  = (time - TIME_POINT(L))/(TIME_POINT(L+1) - TIME_POINT(L)) 
      value = - (1.0D0 - DEXP(-step*value_3))
      omega = OMEGA_STEP*value + omega_POINT(L)
    ELSE 
      step  = (time - TIME_POINT(L))/(TIME_POINT(L+1) - TIME_POINT(L)) 
      value = (1.0D0 - DEXP(-step*value_4))     
      omega = OMEGA_STEP*value + omega_POINT(L)
    END IF
  END IF
! ---------------------------------------------------------
ELSE 
  omega = omega_POINT(L) + SLOPE*(time - TIME_POINT(L))
END IF
! --------------------------------------------------------------------------------- !
DEALLOCATE(OMEGA_STEP_VECTOR,stat=problem)
IF (problem/=0) THEN 
  PRINT *, " PROGRAM COULD NOT DEALLOCATE SPACE FOR THE VECTOR "
  PRINT *, " OMEGA_STEP_VECTOR!                                " 
END IF

DEALLOCATE(omega_POINT,omega_POINT_correction,TIME_POINT,stat=problem)
IF (problem/=0) THEN 
  PRINT *, " PROGRAM COULD NOT DEALLOCATE SPACE FOR THE VECTOR "
  PRINT *, " omega_POINT OR TIME_POINT!                        " 
END IF
! --------------------------------------------------------------------------------- !
RETURN
END SUBROUTINE STANDARD_OMEGA_UP_DOWN
! ================================================================================= !
END MODULE ROTATION
! --------------------------------------------------------------------------------- !
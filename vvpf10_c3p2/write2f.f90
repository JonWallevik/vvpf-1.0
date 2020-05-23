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
! File name: write2f.f90 (MODULE)                                                   !
! This file takes care of writing all computed data into the different files.       !
! It is only the source main.f90 that makes such request.                           !
! --------------------------------------------------------------------------------- !
MODULE WRITE_INFORMATION
  USE SHEAR_VISCOSITY  
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: WARNING_FOR_WRITING,WRITE2FILE_k,WRITE2FILE_kp1,&
            WRITE2FILE_debug,WRITE2FILE_rms
CONTAINS
! ================================================================================= !
SUBROUTINE WARNING_FOR_WRITING(NY)
INTEGER,INTENT(IN) :: NY

IF (NY > 500) THEN
  PRINT *, " ERROR: NY2 > 500 ( NY2 = ",NY,")"
  PRINT *, " FORMAT STATEMENT IN THE FILE 'write2f.f90' IS TO SHORT:     "
  PRINT *, " ERROR -> 10 FORMAT(1X,500(F7.4,1X))                         "
  PRINT *, " PLEASE MAKE THE NECESSARY ADJUSTMENT IN ALL THE SUBROUTINES "
  PRINT *, " OF THIS FILE. TERMINAL ERROR!                               "
  STOP
END IF

RETURN
END SUBROUTINE WARNING_FOR_WRITING
! ================================================================================= !
SUBROUTINE WRITE2FILE_k(VELOCITY,NX)

DOUBLE PRECISION,DIMENSION(:,:),INTENT(IN) :: VELOCITY
INTEGER,INTENT(IN)                         :: NX
INTEGER                                    :: problem,i
! --------------------------------------------------------------------------------- !
10 FORMAT(1X,500(F7.4,1X))
! --------------------------------------------------------------------------------- !
OPEN(unit=8,file="vel_testing.dat",status="replace",action="write",&
     position="rewind",iostat=problem)
IF (problem/=0) THEN
  PRINT *," Could not create the file: vel_testing.dat! "
  RETURN
ELSE
  DO i = 1,NX
    WRITE (unit=8,fmt=10) VELOCITY(i,:)
  END DO  
END IF 

CLOSE (UNIT=8)
! --------------------------------------------------------------------------------- !
RETURN
END SUBROUTINE WRITE2FILE_k
! ================================================================================= !
SUBROUTINE WRITE2FILE_kp1(V,NX1,NX1_corner,NX2,NY1,NY1_corner,NY2,dr,dz)

DOUBLE PRECISION,DIMENSION(:,:),INTENT(IN)  :: V
INTEGER,INTENT(IN)                          :: NX1,NX1_corner,NX2,&
                                               NY1,NY1_corner,NY2
DOUBLE PRECISION,INTENT(IN)                 :: dr,dz

DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: SR,ETA,von_Mises
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)   :: TORQUE_ON_DISK

DOUBLE PRECISION  :: ETA_tmp,SR_tmp,TIME,Lambda,r,PI,Measured_Torque,&
                     shear_stress_theta_r
INTEGER           :: problem,i,j,k
! --------------------------------------------------------------------------------- !
PI = DACOS(-1.0D0)
! --------------------------------------------------------------------------------- !
ALLOCATE(SR(NX2,NY2),ETA(NX2,NY2),von_Mises(NX2,NY2),stat=problem)
IF (problem/=0) THEN 
  PRINT *," WRITE2FILE_kp1: The program could not allocate space! "
  PRINT *," Error code 1 in write2f and execution terminated!     "
  STOP
END IF

ALLOCATE(TORQUE_ON_DISK(NX1),stat=problem)
IF (problem/=0) THEN 
  PRINT *," WRITE2FILE_kp1: The program could not allocate space! "
  PRINT *," Error code 2 in write2f and execution terminated!     "
  STOP
END IF

SR             = 0.0D0
ETA            = 0.0D0
von_Mises      = 0.0D0
TORQUE_ON_DISK = 0.0D0

TIME   = 0.0D0
Lambda = 1.0D0
! --------------------------------------------------------------------------------- !
! CALCULATING THE SHEAR RATE AT ALL POINTS: SR 
CALL SR_PROFILE(V,NX1,NX1_corner,NX2,NY1,NY1_corner,NY2,dr,dz,SR)
! --------------------------------------------------------------------------------- !
! CALCULATING THE SHEAR VISCOSITY AT ALL POINTS:
ETA_tmp = 0.0D0
DO i=1,NX2
  DO j=1,NY2
    SR_tmp   = SR(i,j)
    CALL VISCOSITY(TIME,Lambda,SR_tmp,ETA_tmp)
    ETA(i,j) = ETA_tmp
  END DO
END DO
! --------------------------------------------------------------------------------- !
! CALCULATING THE von Mises SHEAR STRESS AT ALL POINTS:
DO i=1,NX2
  DO j=1,NY2
    von_Mises(i,j) = SR(i,j)*ETA(i,j)
  END DO
END DO
! --------------------------------------------------------------------------------- !
! CALCULATING THE TORQUE ON THE TOP PLATE (z = h_gap):
j = NY1_corner
DO i=2,NX1-1
  r = DBLE(i-1)*dr
  shear_stress_theta_r = - ETA(i,j)*((-4.0D0*V(i,j-1) + V(i,j-2) &
                         + 3.0D0*V(i,j))/(2.0D0*dz))
  TORQUE_ON_DISK(i)    = r*(shear_stress_theta_r*dr*(2*PI*r))
END DO
! ---------------------------------------------------------
i = NX1
r = DBLE(i-1)*dr
shear_stress_theta_r = - ETA(i,j)*((-4.0D0*V(i,j-1) + V(i,j-2) &
                       + 3.0D0*V(i,j))/(2.0D0*dz))
TORQUE_ON_DISK(NX1)  = r*(shear_stress_theta_r*(dr/2)*(2*PI*r)) 
! --------------------------------------------------------------------------------- !
Measured_Torque = 0.0D0
DO i=1,NX1
  Measured_Torque = TORQUE_ON_DISK(i) + Measured_Torque
END DO
! --------------------------------------------------------------------------------- !
10 FORMAT(1X,500(F10.4,1X))
12 FORMAT(1X,500(F14.8,1X))
16 FORMAT(1X,500(F14.4,1X))
! --------------------------------------------------------------------------------- !
OPEN(unit=8,file="vel_t0.dat",status="replace",action="write",&
     position="rewind",iostat=problem)
IF (problem/=0) THEN
  PRINT *," Could not create the file: vel_t0.dat! "
  RETURN
ELSE
  DO i = 1,NX2
    WRITE (unit=8,fmt=10) V(i,:)
  END DO  
END IF 

CLOSE (UNIT=8)
! --------------------------------------------------------------------------------- !
OPEN(unit=8,file="SR_t0.dat",status="replace",action="write",&
     position="rewind",iostat=problem)
IF (problem/=0) THEN
  PRINT *," Could not create the file: SR_t0.dat! "
  RETURN
ELSE
  DO i = 1,NX2
    WRITE (unit=8,fmt=10) SR(i,:)
  END DO  
END IF 

CLOSE (UNIT=8)
! --------------------------------------------------------------------------------- !
OPEN(unit=8,file="ETA_t0.dat",status="replace",action="write",&
     position="rewind",iostat=problem)
IF (problem/=0) THEN
  PRINT *," Could not create the file: ETA_t0.dat! "
  RETURN
ELSE
  DO i = 1,NX2
    WRITE (unit=8,fmt=16) ETA(i,:)
  END DO  
END IF 

CLOSE (UNIT=8)
! --------------------------------------------------------------------------------- !
OPEN(unit=8,file="vonMises_t0.dat",status="replace",action="write",&
     position="rewind",iostat=problem)
IF (problem/=0) THEN
  PRINT *," Could not create the file: vonMises_t0.dat! "
  RETURN
ELSE
  DO i = 1,NX2
    WRITE (unit=8,fmt=10) von_Mises(i,:)
  END DO  
END IF 

CLOSE (UNIT=8)
! --------------------------------------------------------------------------------- !
PRINT *, " --------------------------------------------------------- "
PRINT *, " Measured torque (on the bottom disk plate):               "
PRINT *,   Measured_Torque, " Nm i_z                                 "
PRINT *, " --------------------------------------------------------- "
PRINT *, "                                                           "
! --------------------------------------------------------------------------------- !
24 FORMAT(1X,"4) Measured torque (on the bottom disk plate) = ",(F14.10,1X),"Nm")
! --------------------------------------------------------------------------------- !
OPEN(unit=8,file="TORQUE_t0.dat",status="replace",action="write",&
     position="rewind",iostat=problem)
IF (problem/=0) THEN
  PRINT *," Could not create the file: TORQUE_t0.dat! "
  RETURN
ELSE
  WRITE (unit=8,fmt=24) Measured_Torque
END IF 

CLOSE (UNIT=8)
! --------------------------------------------------------------------------------- !
OPEN(unit=8,file="TORQUE_disc.dat",status="replace",action="write",&
     position="rewind",iostat=problem)
IF (problem/=0) THEN
  PRINT *," Could not create the file: TORQUE_disc.dat! "
  RETURN
ELSE
  WRITE (unit=8,fmt=12) TORQUE_ON_DISK(1:NX1)
END IF 

CLOSE (UNIT=8)
! --------------------------------------------------------------------------------- !
DEALLOCATE(SR,ETA,von_Mises,stat=problem)
IF (problem/=0) THEN 
  PRINT *," WRITE2FILE_kp1: The program could not deallocate space! "
  PRINT *," Error code 3 in write2f!                                "
END IF

DEALLOCATE(TORQUE_ON_DISK,stat=problem)
IF (problem/=0) THEN 
  PRINT *," WRITE2FILE_kp1: The program could not deallocate space! "
  PRINT *," Error code 4 in write2f!                                "
END IF
! --------------------------------------------------------------------------------- !
RETURN
END SUBROUTINE WRITE2FILE_kp1
! ================================================================================= !
SUBROUTINE WRITE2FILE_rms(k,rms) 

DOUBLE PRECISION,INTENT(IN) :: rms
INTEGER,INTENT(IN)          :: k
INTEGER                     :: problem
! --------------------------------------------------------------------------------- !
OPEN(unit=8,file="log.dat",status="old",action="write",&
     position="append",iostat=problem)
IF (problem/=0) THEN
  PRINT *," Could not write into the existing file: log.dat! "
  RETURN
ELSE
  WRITE (unit=8,fmt=*) k,rms
END IF 

CLOSE (UNIT=8)
! --------------------------------------------------------------------------------- !
RETURN
END SUBROUTINE WRITE2FILE_rms
! ================================================================================= !
SUBROUTINE WRITE2FILE_debug(M,K_M,v_new,N)

DOUBLE PRECISION,DIMENSION(:,:),INTENT(IN) :: M
DOUBLE PRECISION,DIMENSION(:),INTENT(IN)   :: K_M,v_new
INTEGER,INTENT(IN)                         :: N
INTEGER                                    :: problem,i
! --------------------------------------------------------------------------------- !
10 FORMAT(1X,201(F10.4,1X))
11 FORMAT(1X,F10.4)
! --------------------------------------------------------------------------------- !
OPEN(unit=8,file="MM_debug.dat",status="replace",action="write",&
     position="rewind",iostat=problem)
IF (problem/=0) THEN
  PRINT *," Could not create the file: MM_debug.dat! "
  RETURN
ELSE
  DO i = 1,N
    WRITE (unit=8,fmt=10) M(i,:)
  END DO
END IF 

CLOSE (UNIT=8)
! ---------------------------------------------------------
OPEN(unit=8,file="KK_debug.dat",status="replace",action="write",&
     position="rewind",iostat=problem)
IF (problem/=0) THEN
  PRINT *," Could not create the file: KK_debug.dat! "
  RETURN
ELSE
  DO i = 1,N
    WRITE (unit=8,fmt=11) K_M(i)
  END DO  
END IF 

CLOSE (UNIT=8)
! ---------------------------------------------------------
OPEN(unit=8,file="vel_debug.dat",status="replace",action="write",&
     position="rewind",iostat=problem)
IF (problem/=0) THEN
  PRINT *," Could not create the file: vel_debug.dat! "
  RETURN
ELSE
  DO i = 1,N
    WRITE (unit=8,fmt=11) v_new(i)
  END DO  
END IF 

CLOSE (UNIT=8)
! --------------------------------------------------------------------------------- !
RETURN
END SUBROUTINE WRITE2FILE_debug
! ================================================================================= !
SUBROUTINE SR_PROFILE(V,NX1,NX1_corner,NX2,NY1,NY1_corner,NY2,dr,dz,SR)

DOUBLE PRECISION,DIMENSION(:,:),INTENT(IN)  :: V
DOUBLE PRECISION,DIMENSION(:,:),INTENT(OUT) :: SR

INTEGER,INTENT(IN)                          :: NX1,NX1_corner,NX2,&
                                               NY1,NY1_corner,NY2
DOUBLE PRECISION,INTENT(IN)                 :: dr,dz

INTEGER                                     :: i,j,k,problem
DOUBLE PRECISION                            :: r,SR1_ij,SR2_ij,EPS
! --------------------------------------------------------------------------------- !
EPS = 1.0D-15
! --------------------------------------------------------------------------------- !
! CALCULATING THE SHEAR RATE IN THE BULK:
! See Section 7.5 about the formulas for the shear rate (SR). Note that ROS and SR
! means the same thing: ROS = rate of shear = SR = shear rate.
! --------------------------------------------------------------------------------- !
DO i=2,NX2-1
  r = DBLE(i-1)*dr
  DO j=2,NY2-1
    SR1_ij  = (V(i+1,j) - V(i-1,j))/(2.0D0*dr) - V(i,j)/r
    SR2_ij  = (V(i,j+1) - V(i,j-1))/(2.0D0*dz)
    SR(i,j) = DSQRT(SR1_ij**2.0D0 + SR2_ij**2.0D0)
  END DO
END DO
! --------------------------------------------------------------------------------- !
! CALCULATING THE SHEAR RATE ON THE BOTTOM PLATE (z=0):
j = 1
DO i=2,NX2-1
  r = DBLE(i-1)*dr
  SR1_ij  = (V(i+1,j) - V(i-1,j))/(2.0D0*dr) - V(i,j)/r
  SR2_ij  = (4.0D0*V(i,j+1) - V(i,j+2) - 3.0D0*V(i,j))/(2.0D0*dz)
  SR(i,j) = DSQRT(SR1_ij**2.0D0 + SR2_ij**2.0D0)
END DO
! ---------------------------------------------------------
! CALCULATING THE SHEAR RATE ON THE "LEFT" WALL (r=R_i):
! SR2_ij is actually zero since the Dirichlet boundary 
! condition is not chancing with z!
i = NX1_corner
r = DBLE(i-1)*dr
DO j=NY1+1,NY2-1
  SR1_ij  = (4.0D0*V(i+1,j) - V(i+2,j) - 3.0D0*V(i,j))/(2.0D0*dr) - V(i,j)/r
  SR2_ij  = (V(i,j+1) - V(i,j-1))/(2.0D0*dz)
  SR(i,j) = DSQRT(SR1_ij**2.0D0 + SR2_ij**2.0D0)
END DO
! ---------------------------------------------------------
! CALCULATING THE SHEAR RATE ON THE "RIGHT" WALL (r=R_o):
! SR2_ij is actually zero since the Dirichlet boundary 
! condition is not chancing with z!
i = NX2
r = DBLE(i-1)*dr
DO j=2,NY2-1
  SR1_ij  = (-4.0D0*V(i-1,j) + V(i-2,j) + 3.0D0*V(i,j))/(2.0D0*dr) - V(i,j)/r
  SR2_ij  = (V(i,j+1) - V(i,j-1))/(2.0D0*dz)
  SR(i,j) = DSQRT(SR1_ij**2.0D0 + SR2_ij**2.0D0)
END DO
! ---------------------------------------------------------
! CALCULATING THE SHEAR RATE ON TOP PLATE (z=h_gap):
j = NY1_corner
DO i=2,NX1-1
  r = DBLE(i-1)*dr
  SR1_ij  = (V(i+1,j) - V(i-1,j))/(2.0D0*dr) - V(i,j)/r
  SR2_ij  = (-4.0D0*V(i,j-1) + V(i,j-2) + 3.0D0*V(i,j))/(2.0D0*dz)
  SR(i,j) = DSQRT(SR1_ij**2.0D0 + SR2_ij**2.0D0)
END DO
! ---------------------------------------------------------
! CALCULATING THE SHEAR RATE AT THE OPEN BOUNDARY (z=H):
j = NY2
DO i=NX1_corner+1,NX2-1
  r = DBLE(i-1)*dr
  SR1_ij  = (V(i+1,j) - V(i-1,j))/(2.0D0*dr) - V(i,j)/r
  SR2_ij  = (-4.0D0*V(i,j-1) + V(i,j-2) + 3.0D0*V(i,j))/(2.0D0*dz)
  SR(i,j) = DSQRT(SR1_ij**2.0D0 + SR2_ij**2.0D0)
END DO
! --------------------------------------------------------------------------------- !
! CALCULATING THE SHEAR RATE AT THE CENTER OF THE BOTTOM PLATE (r=0,z=0):
! This calculation is redundant since the rate of shear at the center line
! ($r=0 \forall z \in [0,H]$) is zero due to symmetry in the r-direction and 
! due to the Dirichlet boundary condition $v_{\rm i,j}=0$ at the center line.
! i = 1
! j = 1
! r = DBLE(i-1)*dr
! SR1_ij  = (V(i+1,j) - V(i+1,j))/(2.0D0*dr) - V(i,j)/r ! due to symmetry
! SR2_ij  = ( 4.0D0*V(i,j+1) - V(i,j+2) - 3.0D0*V(i,j))/(2.0D0*dz)
! SR(i,j) = DSQRT(SR1_ij**2.0D0 + SR2_ij**2.0D0)
! ------------
! Rather enforcing a zero rate of shear at the center line:
SR(1,1:NY2) = 0.0D0
! --------------------------------------------------------------------------------- !
! CALCULATING THE SHEAR RATE ON THE LOWER RIGHT CORNER (r=R_o,z=0):
i = NX2
j = 1
r = DBLE(i-1)*dr
SR1_ij  = (-4.0D0*V(i-1,j) + V(i-2,j) + 3.0D0*V(i,j))/(2.0D0*dr) - V(i,j)/r
SR2_ij  = ( 4.0D0*V(i,j+1) - V(i,j+2) - 3.0D0*V(i,j))/(2.0D0*dz)
SR(i,j) = DSQRT(SR1_ij**2.0D0 + SR2_ij**2.0D0)
! ---------------------------------------------------------
! CALCULATING THE SHEAR RATE ON THE TOP RIGHT CORNER (r=R_o,z=H):
i = NX2
j = NY2
r = DBLE(i-1)*dr
SR1_ij  = (-4.0D0*V(i-1,j) + V(i-2,j) + 3.0D0*V(i,j))/(2.0D0*dr) - V(i,j)/r
SR2_ij  = (-4.0D0*V(i,j-1) + V(i,j-2) + 3.0D0*V(i,j))/(2.0D0*dz)
SR(i,j) = DSQRT(SR1_ij**2.0D0 + SR2_ij**2.0D0)
! ---------------------------------------------------------
! CALCULATING THE SHEAR RATE ON THE TOP LEFT CORNER (r=R_i,z=H):
i = NX1_corner
j = NY2
r = DBLE(i-1)*dr
SR1_ij  = ( 4.0D0*V(i+1,j) - V(i+2,j) - 3.0D0*V(i,j))/(2.0D0*dr) - V(i,j)/r
SR2_ij  = (-4.0D0*V(i,j-1) + V(i,j-2) + 3.0D0*V(i,j))/(2.0D0*dz)
SR(i,j) = DSQRT(SR1_ij**2.0D0 + SR2_ij**2.0D0)
! ---------------------------------------------------------
! CALCULATING THE SHEAR RATE AT SMOOTH CORNER (r=[NX1,NX1_corner];
!                                          and z=[NY1_corner,NY1]):
DO k = 1,NY1-NY1_corner+1 ! k=1,11
  i = NX1 + (k-1)
  j = NY1_corner + (k-1)
  r = DBLE(i-1)*dr
  SR1_ij  = ( 4.0D0*V(i+1,j) - V(i+2,j) - 3.0D0*V(i,j))/(2.0D0*dr) - V(i,j)/r
  SR2_ij  = (-4.0D0*V(i,j-1) + V(i,j-2) + 3.0D0*V(i,j))/(2.0D0*dz)
  SR(i,j) = DSQRT(SR1_ij**2.0D0 + SR2_ij**2.0D0)
END DO
! --------------------------------------------------------------------------------- !
RETURN
END SUBROUTINE SR_PROFILE
! ================================================================================= !
END MODULE WRITE_INFORMATION
! --------------------------------------------------------------------------------- !
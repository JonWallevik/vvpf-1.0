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
  PUBLIC :: WARNING_FOR_WRITING,WRITE2FILE_k,WRITE2FILE_kp1,WRITE2FILE_time,&
            WRITE2FILE_torque_ZERO,WRITE2FILE_torque,WRITE2FILE_debug,&
            WRITE2FILE_rms,ROS_PROFILE
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
OPEN(unit=8,file="vel_time.dat",status="replace",action="write",&
     position="rewind",iostat=problem)
IF (problem/=0) THEN
  PRINT *," Could not create the file: vel_time.dat! "
  RETURN
ELSE
  DO i = 1,NX
    WRITE (unit=8,fmt=10) VELOCITY(i,:)
  END DO  
END IF 

CLOSE (UNIT=8)
! --------------------------------------------------------------------------------- !
! Creating file to log time, for the transient calculation:
OPEN(unit=8,file="time_k.dat",status="replace",action="write",&
     position="rewind",iostat=problem)
IF (problem/=0) THEN
  PRINT *," Could not create the file: time_k.dat! "
  RETURN
ELSE
  WRITE (unit=8,fmt=10) 0.0D0
END IF 

CLOSE (UNIT=8)
! --------------------------------------------------------------------------------- !
RETURN
END SUBROUTINE WRITE2FILE_k
! ================================================================================= !
SUBROUTINE WRITE2FILE_kp1(V,NX1,NY1,NX2,NY2,dr,dz)

DOUBLE PRECISION,DIMENSION(:,:),INTENT(IN)  :: V
INTEGER,INTENT(IN)                          :: NX1,NY1,NX2,NY2
DOUBLE PRECISION,INTENT(IN)                 :: dr,dz

DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: SR,ETA,von_Mises
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)   :: SR_r

DOUBLE PRECISION  :: ETA_tmp,SR_tmp,TIME,Lambda
INTEGER           :: problem,i,j,NX,NYt08
! --------------------------------------------------------------------------------- !
NYt08    = IDNINT(0.8D0*DBLE(NY2))
NX       = NX2 - NX1 + 1 
! --------------------------------------------------------------------------------- !
ALLOCATE(SR(NX2,NY2),SR_r(NX),ETA(NX2,NY2),von_Mises(NX2,NY2),stat=problem)
IF (problem/=0) THEN 
  PRINT *," WRITE2FILE_kp1: The program could not allocate space! "
  PRINT *," Error code 1 in write2f and execution terminated!     "
  STOP
END IF

SR        = 0.0D0
ETA       = 0.0D0
von_Mises = 0.0D0

TIME      = 0.0D0
Lambda    = 1.0D0
! --------------------------------------------------------------------------------- !
! CALCULATING THE SHEAR RATE AT ALL POINTS: SR 
! Rate of shear at z = NYt08 as function of r: SR_r (not used here!)
CALL ROS_PROFILE(V,NX1,NX2,NY1,NY2,NYt08,dr,dz,SR,SR_r)
! --------------------------------------------------------------------------------- !
! CALCULATING THE SHEAR VISCOSITY AT ALL POINTS:
ETA_tmp = 0.0D0

DO i=1,NX2
  DO j=1,NY2
    SR_tmp   = SR(i,j) 
    ! Using dt = 0.0D0 here, because there is no iteration forward in time. The objectives
    ! here is simply to write the last data (of time = k*dt, where k=0) to the hard drive. 
    ! I.e. we want to have Gamma=ALPHA_II*(FMSR+ALPHA_I*SR*dt)=ALPHA_II*(0+ALPHA_I*SR*0)=0, 
    ! because no shear rate history exists at k = 0 (FMSR is already fully updated at k=0).
    ! [ The same consideration applies for FMCR, BETA_I, BETA_II and H(SR) ]. 
    CALL VISCOSITY(0.0D0,TIME,Lambda,SR_tmp,0.0D0,0.0D0,ETA_tmp) 
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
10 FORMAT(1X,500(F10.4,1X))
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
OPEN(unit=8,file="ROS_t0.dat",status="replace",action="write",&
     position="rewind",iostat=problem)
IF (problem/=0) THEN
  PRINT *," Could not create the file: ROS_t0.dat! "
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
DEALLOCATE(SR,SR_r,ETA,von_Mises,stat=problem)
IF (problem/=0) THEN 
  PRINT *," WRITE2FILE_kp1: The program could not deallocate space! "
  PRINT *," Error code 2 in write2f!                                "
END IF
! --------------------------------------------------------------------------------- !
RETURN
END SUBROUTINE WRITE2FILE_kp1
! ================================================================================= !
SUBROUTINE WRITE2FILE_time(VELOCITY,NX,kp1,dt)

DOUBLE PRECISION,DIMENSION(:,:),INTENT(IN) :: VELOCITY
DOUBLE PRECISION,INTENT(IN)                :: dt
INTEGER,INTENT(IN)                         :: NX,kp1
INTEGER                                    :: problem,i
! --------------------------------------------------------------------------------- !
10 FORMAT(1X,500(F7.4,1X))
! --------------------------------------------------------------------------------- !
OPEN(unit=8,file="vel_time.dat",status="old",action="write",&
     position="append",iostat=problem)
IF (problem/=0) THEN
  PRINT *," Could not write into the existing file: vel_time.dat! "
  RETURN
ELSE
  DO i = 1,NX
    WRITE (unit=8,fmt=10) VELOCITY(i,:)
  END DO  
END IF 

CLOSE (UNIT=8)
! --------------------------------------------------------------------------------- !
OPEN(unit=8,file="time_k.dat",status="old",action="write",&
     position="append",iostat=problem)
IF (problem/=0) THEN
  PRINT *," Could not write into the existing file: time_k.dat! "
  RETURN
ELSE
  WRITE (unit=8,fmt=10) DBLE(kp1)*dt
END IF 

CLOSE (UNIT=8)
! --------------------------------------------------------------------------------- !
RETURN
END SUBROUTINE WRITE2FILE_time
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
SUBROUTINE WRITE2FILE_torque_ZERO(V,NX1,NX2,NY1,NY2,kp1,dt,Lambda,dr,dz,H3,omega)

DOUBLE PRECISION,DIMENSION(:,:),INTENT(IN)  :: V
DOUBLE PRECISION,INTENT(IN)                 :: dt,Lambda,dr,dz,H3,omega
INTEGER,INTENT(IN)                          :: NX1,NX2,NY1,NY2,kp1

DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: SR,ETA
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)   :: SR_r,ETA_r,ROS_torque,ETA_torque,&
                                               TORQUE_ON_INNER_CYLINDER,&
                                               TORQUE_ON_OUTER_CYLINDER,&
                                               TORQUE_ON_BOTTOM_PLATE
DOUBLE PRECISION :: ETA_tmp,SR_tmp,&
                    TIME,TORQUE_sum_Ri,&
                    TORQUE_sum_work
INTEGER          :: NY2mH3,NYt08,problem,i,j,NX,NY,&
                    count_i,count_j,NY_write
! --------------------------------------------------------------------------------- !
NYt08    = IDNINT(0.8D0*DBLE(NY2))
NY2mH3   = NY2 - IDNINT(H3/dz)
NY_write = 50
NX       = NX2 - NX1 + 1    ! -> Same as in main.f90 
NY       = NY2 - NY2mH3 + 1 ! -> Or equally, NY = IDNINT(H3/dz) + 1
TIME     = DBLE(kp1)*dt

IF (NY2mH3.LT.NY1+1) THEN
  PRINT *, "-----------------------------------------------"
  PRINT *, "  ERROR:                                       "
  PRINT *, "    write2f.f90 says: Something is wrong!      "
  PRINT *, "    NY2-H3/dz < NY1+1. Check out if H3         "
  PRINT *, "    (in the file param.f90) is correct!        "
  PRINT *, "-----------------------------------------------"
  STOP
END IF

ALLOCATE(SR(NX2,NY2),ETA(NX2,NY2),SR_r(NX),ETA_r(NX),ROS_torque(NY),&
         ETA_torque(NY),TORQUE_ON_INNER_CYLINDER(NY),&
         TORQUE_ON_OUTER_CYLINDER(NY2),TORQUE_ON_BOTTOM_PLATE(NX2),&
         stat=problem)
IF (problem/=0) THEN 
  PRINT *," WRITE2FILE_torque_ZERO: The program could not allocate space! "
  PRINT *," Error code 3 in write2f and execution terminated!             "
  STOP
END IF

SR         = 0.0D0
ETA        = 0.0D0
SR_r       = 0.0D0
ETA_r      = 0.0D0
ROS_torque = 0.0D0
ETA_torque = 0.0D0

TORQUE_ON_INNER_CYLINDER = 0.0D0
TORQUE_ON_OUTER_CYLINDER = 0.0D0
TORQUE_ON_BOTTOM_PLATE   = 0.0D0
! --------------------------------------------------------------------------------- !
! CALCULATING THE SHEAR RATE AT ALL POINTS: SR 
! Rate of shear at z = NYt08 as function of r: SR_r
CALL ROS_PROFILE(V,NX1,NX2,NY1,NY2,NYt08,dr,dz,SR,SR_r)
! --------------------------------------------------------------------------------- !
! CALCULATING THE SHEAR VISCOSITY AT ALL POINTS: ETA
! Shear viscosity at z = NYt08 as function of r: ETA_r
ETA_tmp = 0.0D0
DO i=1,NX2
  DO j=1,NY2
    SR_tmp   = SR(i,j)
    ! Using dt = 0.0D0 here, because there is no iteration forward in time. The objectives
    ! here is simply to write the last data (of time = k*dt, where k=0) to the hard drive. 
    ! I.e. we want to have Gamma=ALPHA_II*(FMSR+ALPHA_I*SR*dt)=ALPHA_II*(0+ALPHA_I*SR*0)=0, 
    ! because no shear rate history exists at k = 0 (FMSR is already fully updated at k=0).
    ! [ The same consideration applies for FMCR, BETA_I, BETA_II and H(SR) ]. 
    CALL VISCOSITY(0.0D0,TIME,Lambda,SR_tmp,0.0D0,0.0D0,ETA_tmp) 
    ETA(i,j) = ETA_tmp
  END DO
END DO
ETA_tmp = 0.0D0
DO i=1,NX
  SR_tmp   = SR_r(i)
  ! Using dt = 0.0D0 here, because there is no iteration forward in time. The objectives
  ! here is simply to write the last data (of time = k*dt, where k=0) to the hard drive. 
  ! I.e. we want to have Gamma=ALPHA_II*(FMSR+ALPHA_I*SR*dt)=ALPHA_II*(0+ALPHA_I*SR*0)=0, 
  ! because no shear rate history exists at k = 0 (FMSR is already fully updated at k=0).
  ! [ The same consideration applies for FMCR, BETA_I, BETA_II and H(SR) ]. 
  CALL VISCOSITY(0.0D0,TIME,Lambda,SR_tmp,0.0D0,0.0D0,ETA_tmp) 
  ETA_r(i) = ETA_tmp
END DO
! --------------------------------------------------------------------------------- !
! CALCULATING TORQUE APPLIED ON THE INNER CYLINDER, ON THE OUTER CYLINDER
! AND ON THE BOTTOM PLATE, FROM THE TEST MATERIAL:
CALL TORQUE(V,ETA,NX1,NX2,NY1,NY2,NY2mH3,dr,dz,TORQUE_ON_INNER_CYLINDER,&
            TORQUE_ON_OUTER_CYLINDER,TORQUE_ON_BOTTOM_PLATE)
! ---------------------------------------------------------
! ... AND THEN SUMMARIZING ALL THE TORQUE ELEMENTS
! for all z in H3 (i.e. at where torque is measured),...
TORQUE_sum_Ri = 0.0D0
DO j=1,NY
  TORQUE_sum_Ri = TORQUE_ON_INNER_CYLINDER(j) + TORQUE_sum_Ri
END DO
! ---------------------------------------------------------
! ... AND FOR WORK CALCULATIONS:
TORQUE_sum_work = 0.0D0
DO j=1,NY2
  TORQUE_sum_work = TORQUE_ON_OUTER_CYLINDER(j) + TORQUE_sum_work
END DO
DO i=1,NX2
  TORQUE_sum_work = TORQUE_ON_BOTTOM_PLATE(i) + TORQUE_sum_work
END DO
! --------------------------------------------------------------------------------- !
! WRITING INFORMATION TO FILES:
! --------------------------------------------------------------------------------- !
10 FORMAT(1X,500(F9.5,1X))
12 FORMAT(1X,2(F9.5,1X))
14 FORMAT(1X,500(F5.3,1X))
16 FORMAT(1X,(F6.2,1X),(F8.5,1X),(F10.5,1X),(F11.5,1X))
18 FORMAT(1X,500(F7.3,1X))
! --------------------------------------------------------------------------------- !
! WRITING INNER AND OUTER CYLINDER TO FILE:
OPEN(unit=8,file="RiRo_dr.dat",status="replace",action="write",&
     position="rewind",iostat=problem)
IF (problem/=0) THEN
  PRINT *," Could not create the file: RiRo_dr.dat! "
ELSE
  WRITE (unit=8,fmt=12) DBLE(NX1-1)*dr,DBLE(NX2-1)*dr 
  WRITE (unit=8,fmt=12) dr,dz
END IF 

CLOSE (UNIT=8)
! --------------------------------------------------------------------------------- !
! WRITING ANGULAR VELOCITY TO FILE:
OPEN(unit=8,file="time_omega.dat",status="replace",action="write",&
     position="rewind",iostat=problem)
IF (problem/=0) THEN
  PRINT *," Could not create the file: time_omega.dat! "
ELSE
  WRITE (unit=8,fmt=12) DBLE(kp1)*dt,omega
END IF 

CLOSE (UNIT=8)
! --------------------------------------------------------------------------------- !
! WRITING VELOCITY TO FILE:
OPEN(unit=8,file="velocity_r.dat",status="replace",action="write",&
     position="rewind",iostat=problem)
IF (problem/=0) THEN
  PRINT *," Could not create the file: velocity_r.dat! "
ELSE
  WRITE (unit=8,fmt=10) V(NX1:NX2,NYt08)
END IF 

CLOSE (UNIT=8)
! --------------------------------------------------------------------------------- !
! WRITING SHEAR RATE TO FILE:
OPEN(unit=8,file="ros_r.dat",status="replace",action="write",&
     position="rewind",iostat=problem)
IF (problem/=0) THEN
  PRINT *," Could not create the file: ros_r.dat! "
ELSE
  WRITE (unit=8,fmt=10) SR_r(:)
END IF 

CLOSE (UNIT=8)
! --------------------------------------------------------------------------------- !
! WRITING SHEAR VISCOSITY TO FILE:
OPEN(unit=8,file="eta_r.dat",status="replace",action="write",&
     position="rewind",iostat=problem)
IF (problem/=0) THEN
  PRINT *," Could not create the file: eta_r.dat! "
ELSE
  WRITE (unit=8,fmt=10) ETA_r(:)
END IF 

CLOSE (UNIT=8)
! --------------------------------------------------------------------------------- !
! WRITING TORQUE TO FILE:
OPEN(unit=8,file="torque_z.dat",status="replace",action="write",&
     position="rewind",iostat=problem)
IF (problem/=0) THEN
  PRINT *," Could not create the file: torque_z.dat! "
ELSE
  WRITE (unit=8,fmt=10) TORQUE_ON_INNER_CYLINDER(:)
END IF 

CLOSE (UNIT=8)
! --------------------------------------------------------------------------------- !
! WRITING TORQUE TO FILE:
OPEN(unit=8,file="torque.dat",status="replace",action="write",&
     position="rewind",iostat=problem)
IF (problem/=0) THEN
  PRINT *," Could not create the file: torque.dat! "
ELSE
  WRITE (unit=8,fmt=12) DBLE(kp1)*dt,TORQUE_sum_Ri
END IF 

CLOSE (UNIT=8)
! --------------------------------------------------------------------------------- !
! WRITING WORK DONE BY THE OUTER CYLINDER AND BOTTOM PLATE, TO FILE:
OPEN(unit=8,file="work.dat",status="replace",action="write",&
     position="rewind",iostat=problem)
IF (problem/=0) THEN
  PRINT *," Could not create the file: work.dat! "
ELSE
  WRITE (unit=8,fmt=16) DBLE(kp1)*dt,omega,-TORQUE_sum_work,-TORQUE_sum_work*omega
  ! See the comment in the SUBROUTINE TORQUE about why the term "TORQUE_sum_work"
  ! is written with a negative value into the file "work.dat".
END IF 

CLOSE (UNIT=8)
! --------------------------------------------------------------------------------- !
! WRITING VELOCITY OF CORNER TO FILE:
OPEN(unit=8,file="vel_corner.dat",status="replace",action="write",&
     position="rewind",iostat=problem)
IF (problem/=0) THEN
  PRINT *," Could not create the file: vel_corner.dat! "
ELSE
!NX1_begin = NX1/2 + DNINT(0.5D0*DBLE(NX1/2))
  DO i = NX2-50,NX2
    WRITE (unit=8,fmt=14) V(i,1:NY_write+1)
  END DO  
END IF 

CLOSE (UNIT=8)
! --------------------------------------------------------------------------------- !
! WRITING VELOCITY OF UPPER PART TO FILE:
OPEN(unit=8,file="vel_upper.dat",status="replace",action="write",&
     position="rewind",iostat=problem)
IF (problem/=0) THEN
  PRINT *," Could not create the file: vel_upper.dat! "
ELSE
  DO i = NX2-50,NX2
    WRITE (unit=8,fmt=14) V(i,NY2-NY_write:NY2)
  END DO  
END IF 

CLOSE (UNIT=8)
! --------------------------------------------------------------------------------- !
! WRITING SHEAR RATE (SR) OF CORNER TO FILE:
OPEN(unit=8,file="ROS_corner.dat",status="replace",action="write",&
     position="rewind",iostat=problem)
IF (problem/=0) THEN
  PRINT *," Could not create the file: ROS_corner.dat! "
ELSE
!NX1_begin = NX1/2 + DNINT(0.5D0*DBLE(NX1/2))
  DO i = NX2-50,NX2
    WRITE (unit=8,fmt=18) SR(i,1:NY_write+1)
  END DO
END IF

CLOSE (UNIT=8)
! --------------------------------------------------------------------------------- !
! WRITING SHEAR RATE (SR) OF UPPER PART TO FILE:
OPEN(unit=8,file="ROS_upper.dat",status="replace",action="write",&
     position="rewind",iostat=problem)
IF (problem/=0) THEN
  PRINT *," Could not create the file: ROS_upper.dat! "
ELSE
  DO i = NX2-50,NX2
    WRITE (unit=8,fmt=18) SR(i,NY2-NY_write:NY2)
  END DO
END IF

CLOSE (UNIT=8)
! --------------------------------------------------------------------------------- !
DEALLOCATE(SR,ETA,SR_r,ETA_r,ROS_torque,ETA_torque,&
           TORQUE_ON_INNER_CYLINDER,TORQUE_ON_OUTER_CYLINDER,&
           TORQUE_ON_BOTTOM_PLATE,stat=problem)
IF (problem/=0) THEN 
  PRINT *," WRITE2FILE_torque_ZERO: The program could not deallocate space! "
  PRINT *," Error code 4 in write2f!                                        "
END IF
! --------------------------------------------------------------------------------- !
RETURN
END SUBROUTINE WRITE2FILE_torque_ZERO
! ================================================================================= !
SUBROUTINE WRITE2FILE_torque(V,FMSR,FMCR,NX1,NX2,NY1,NY2,kp1,dt,Lambda,dr,dz,H3,omega)

DOUBLE PRECISION,DIMENSION(:,:),INTENT(IN)  :: V,FMSR,FMCR
DOUBLE PRECISION,INTENT(IN)                 :: dt,Lambda,dr,dz,H3,omega
INTEGER,INTENT(IN)                          :: NX1,NX2,NY1,NY2,kp1

DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: SR,ETA
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)   :: SR_r,ETA_r,ROS_torque,ETA_torque,&
                                               TORQUE_ON_INNER_CYLINDER,&
                                               TORQUE_ON_OUTER_CYLINDER,&
                                               TORQUE_ON_BOTTOM_PLATE
DOUBLE PRECISION :: ETA_tmp,SR_tmp,TIME,&
                    TORQUE_sum_Ri,FMSR_tmp,&
                    FMCR_tmp,TORQUE_sum_work

INTEGER          :: NY2mH3,NYt08,problem,i,j,NX,NY,&
                    count_i,count_j,NY_write
! --------------------------------------------------------------------------------- !
NYt08    = IDNINT(0.8D0*DBLE(NY2))
NY2mH3   = NY2 - IDNINT(H3/dz)
NY_write = 50
NX       = NX2 - NX1 + 1    ! -> Same as in main.f90
NY       = NY2 - NY2mH3 + 1 ! -> Or equally, NY = IDNINT(H3/dz) + 1
TIME     = DBLE(kp1)*dt

ALLOCATE(SR(NX2,NY2),ETA(NX2,NY2),SR_r(NX),ETA_r(NX),ROS_torque(NY),&
         ETA_torque(NY),TORQUE_ON_INNER_CYLINDER(NY),&
         TORQUE_ON_OUTER_CYLINDER(NY2),TORQUE_ON_BOTTOM_PLATE(NX2),&
         stat=problem)
IF (problem/=0) THEN 
  PRINT *," WRITE2FILE_torque: The program could not allocate space! "
  PRINT *," Error code 5 in write2f and execution terminated!        "
  STOP
END IF

SR         = 0.0D0
ETA        = 0.0D0
SR_r       = 0.0D0
ETA_r      = 0.0D0
ROS_torque = 0.0D0
ETA_torque = 0.0D0

TORQUE_ON_INNER_CYLINDER = 0.0D0
TORQUE_ON_OUTER_CYLINDER = 0.0D0
TORQUE_ON_BOTTOM_PLATE   = 0.0D0
! --------------------------------------------------------------------------------- !
! CALCULATING THE SHEAR RATE AT ALL POINTS: SR 
! Rate of shear at z = NYt08 as function of r: SR_r
CALL ROS_PROFILE(V,NX1,NX2,NY1,NY2,NYt08,dr,dz,SR,SR_r)
! --------------------------------------------------------------------------------- !
! CALCULATING THE SHEAR VISCOSITY AT ALL POINTS: ETA
! Shear viscosity at z = NYt08 as function of r: ETA_r
ETA_tmp = 0.0D0
DO i=1,NX2
  DO j=1,NY2
    SR_tmp   = SR(i,j)
    FMSR_tmp = FMSR(i,j)
    FMCR_tmp = FMCR(i,j) 
    ! Using dt = 0.0D0 here, because there is no iterating forward in time. The objectives
    ! here is simply to write the last data (of time = (k+1)*dt) to the hard drive. 
    ! I.e. we want to have Gamma=ALPHA_II*FMSR and not Gamma=ALPHA_II*(FMSR+ALPHA_I*SR*dt), 
    ! because FMSR is already fully updated relative to k+1.
    ! [ The same consideration applies for FMCR, BETA_I, BETA_II and H(SR) ]. 
    CALL VISCOSITY(0.0D0,TIME,Lambda,SR_tmp,FMSR_tmp,FMCR_tmp,ETA_tmp)
    ETA(i,j) = ETA_tmp
  END DO
END DO
ETA_tmp = 0.0D0
DO i=1,NX
  SR_tmp   = SR_r(i)
  FMSR_tmp = FMSR((NX1-1)+i,NY2mH3) 
  FMCR_tmp = FMCR((NX1-1)+i,NY2mH3) 
  ! Using dt = 0.0D0 here, because there is no iterating forward in time. The objectives
  ! here is simply to write the last data (of time = (k+1)*dt) to the hard drive.
  ! I.e. we want to have Gamma=ALPHA_II*FMSR and not Gamma=ALPHA_II*(FMSR+ALPHA_I*SR*dt), 
  ! because FMSR is already fully updated relative to k+1.
  ! [ The same consideration applies for FMCR, BETA_I, BETA_II and H(SR) ]. 
  CALL VISCOSITY(0.0D0,TIME,Lambda,SR_tmp,FMSR_tmp,FMCR_tmp,ETA_tmp) 
  ETA_r(i) = ETA_tmp
END DO
! --------------------------------------------------------------------------------- !
! CALCULATING TORQUE APPLIED ON THE INNER CYLINDER, ON THE OUTER CYLINDER
! AND ON THE BOTTOM PLATE, FROM THE TEST MATERIAL:
CALL TORQUE(V,ETA,NX1,NX2,NY1,NY2,NY2mH3,dr,dz,TORQUE_ON_INNER_CYLINDER,&
            TORQUE_ON_OUTER_CYLINDER,TORQUE_ON_BOTTOM_PLATE)
! ---------------------------------------------------------
! ... AND THEN SUMMARIZING ALL THE TORQUE ELEMENTS
! for all z in H3 (i.e. at where torque is measured),...
TORQUE_sum_Ri = 0.0D0
DO j=1,NY
  TORQUE_sum_Ri = TORQUE_ON_INNER_CYLINDER(j) + TORQUE_sum_Ri
END DO
! ---------------------------------------------------------
! ... AND FOR WORK CALCULATIONS:
TORQUE_sum_work = 0.0D0
DO j=1,NY2
  TORQUE_sum_work = TORQUE_ON_OUTER_CYLINDER(j) + TORQUE_sum_work
END DO
DO i=1,NX2
  TORQUE_sum_work = TORQUE_ON_BOTTOM_PLATE(i) + TORQUE_sum_work
END DO
! --------------------------------------------------------------------------------- !
! WRITING INFORMATION TO FILES:
! --------------------------------------------------------------------------------- !
10 FORMAT(1X,500(F9.5,1X))
12 FORMAT(1X,2(F9.5,1X))
14 FORMAT(1X,500(F5.3,1X))
16 FORMAT(1X,(F6.2,1X),(F8.5,1X),(F10.5,1X),(F11.5,1X))
18 FORMAT(1X,500(F7.3,1X))
! --------------------------------------------------------------------------------- !
! WRITING ANGULAR VELOCITY TO FILE:
OPEN(unit=8,file="time_omega.dat",status="old",action="write",&
     position="append",iostat=problem)
IF (problem/=0) THEN
  PRINT *," Could not write into the existing file: time_omega.dat! "
ELSE
  WRITE (unit=8,fmt=12) DBLE(kp1)*dt,omega
END IF 

CLOSE (UNIT=8)
! --------------------------------------------------------------------------------- !
! WRITING VELOCITY TO FILE:
OPEN(unit=8,file="velocity_r.dat",status="old",action="write",&
     position="append",iostat=problem)
IF (problem/=0) THEN
  PRINT *," Could not write into the existing file: velocity_r.dat! "
ELSE
  WRITE (unit=8,fmt=10) V(NX1:NX2,NYt08)
END IF 

CLOSE (UNIT=8)
! --------------------------------------------------------------------------------- !
! WRITING SHEAR RATE TO FILE:
OPEN(unit=8,file="ros_r.dat",status="old",action="write",&
     position="append",iostat=problem)
IF (problem/=0) THEN
  PRINT *," Could not write into the existing file: ros_r.dat! "
ELSE
  WRITE (unit=8,fmt=10) SR_r(:)
END IF 

CLOSE (UNIT=8)
! --------------------------------------------------------------------------------- !
! WRITING SHEAR VISCOSITY TO FILE:
OPEN(unit=8,file="eta_r.dat",status="old",action="write",&
     position="append",iostat=problem)
IF (problem/=0) THEN
  PRINT *," Could not write into the existing file: eta_r.dat! "
ELSE
  WRITE (unit=8,fmt=10) ETA_r(:)
END IF 

CLOSE (UNIT=8)
! --------------------------------------------------------------------------------- !
! WRITING TORQUE TO FILE:
OPEN(unit=8,file="torque_z.dat",status="old",action="write",&
     position="append",iostat=problem)
IF (problem/=0) THEN
  PRINT *," Could not write into the existing file: torque_z.dat! "
ELSE
  WRITE (unit=8,fmt=10) TORQUE_ON_INNER_CYLINDER(:)
END IF 

CLOSE (UNIT=8)
! --------------------------------------------------------------------------------- !
! WRITING TORQUE TO FILE:
OPEN(unit=8,file="torque.dat",status="old",action="write",&
     position="append",iostat=problem)
IF (problem/=0) THEN
  PRINT *," Could not write into the existing file: torque.dat! "
ELSE
  WRITE (unit=8,fmt=12) DBLE(kp1)*dt,TORQUE_sum_Ri
END IF 

CLOSE (UNIT=8)
! --------------------------------------------------------------------------------- !
! WRITING WORK DONE BY THE OUTER CYLINDER AND BOTTOM PLATE, TO FILE:
OPEN(unit=8,file="work.dat",status="old",action="write",&
     position="append",iostat=problem)
IF (problem/=0) THEN
  PRINT *," Could not write into the existing file: work.dat! "
ELSE
  WRITE (unit=8,fmt=16) DBLE(kp1)*dt,omega,-TORQUE_sum_work,-TORQUE_sum_work*omega
  ! See the comment in the SUBROUTINE TORQUE about why the term "TORQUE_sum_work"
  ! is written with a negative value into the file "work.dat".
END IF 

CLOSE (UNIT=8)
! --------------------------------------------------------------------------------- !
! WRITING VELOCITY OF CORNER TO FILE:
OPEN(unit=8,file="vel_corner.dat",status="old",action="write",&
     position="append",iostat=problem)
IF (problem/=0) THEN
  PRINT *," Could not write into the existing file: vel_corner.dat! "
ELSE
  DO i = NX2-50,NX2
    WRITE (unit=8,fmt=14) V(i,1:NY_write+1)
  END DO  
END IF 

CLOSE (UNIT=8)
! --------------------------------------------------------------------------------- !
! WRITING VELOCITY OF UPPER PART TO FILE:
OPEN(unit=8,file="vel_upper.dat",status="old",action="write",&
     position="append",iostat=problem)
IF (problem/=0) THEN
  PRINT *," Could not write into the existing file: vel_upper.dat! "
ELSE
  DO i = NX2-50,NX2
    WRITE (unit=8,fmt=14) V(i,NY2-NY_write:NY2)
  END DO  
END IF 

CLOSE (UNIT=8)
! --------------------------------------------------------------------------------- !
! WRITING SHEAR RATE (SR) OF CORNER TO FILE:
OPEN(unit=8,file="ROS_corner.dat",status="old",action="write",&
     position="append",iostat=problem)
IF (problem/=0) THEN
  PRINT *," Could not write into the existing file: ROS_corner.dat! "
ELSE
!NX1_begin = NX1/2 + DNINT(0.5D0*DBLE(NX1/2))
  DO i = NX2-50,NX2
    WRITE (unit=8,fmt=18) SR(i,1:NY_write+1)
  END DO  
END IF 

CLOSE (UNIT=8)
! --------------------------------------------------------------------------------- !
! WRITING SHEAR RATE (SR) OF UPPER PART TO FILE:
OPEN(unit=8,file="ROS_upper.dat",status="old",action="write",&
     position="append",iostat=problem)
IF (problem/=0) THEN
  PRINT *," Could not write into the existing file: ROS_upper.dat! "
ELSE
  DO i = NX2-50,NX2
    WRITE (unit=8,fmt=18) SR(i,NY2-NY_write:NY2)
  END DO  
END IF 

CLOSE (UNIT=8)
! --------------------------------------------------------------------------------- !
DEALLOCATE(SR,ETA,SR_r,ETA_r,ROS_torque,ETA_torque,&
           TORQUE_ON_INNER_CYLINDER,TORQUE_ON_OUTER_CYLINDER,&
           TORQUE_ON_BOTTOM_PLATE,stat=problem)
IF (problem/=0) THEN 
  PRINT *," WRITE2FILE_torque: The program could not deallocate space! "
  PRINT *," Error code 6 in write2f!                                   "
END IF
! --------------------------------------------------------------------------------- !
RETURN
END SUBROUTINE WRITE2FILE_torque
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
! --------------------------------------------------------------------------------- !
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
! --------------------------------------------------------------------------------- !
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
SUBROUTINE ROS_PROFILE(V,NX1,NX2,NY1,NY2,NYt08,dr,dz,SR,SR_r)

DOUBLE PRECISION,DIMENSION(:,:),INTENT(IN)  :: V
DOUBLE PRECISION,DIMENSION(:,:),INTENT(OUT) :: SR
DOUBLE PRECISION,DIMENSION(:),INTENT(OUT)   :: SR_r
INTEGER,INTENT(IN)                          :: NX1,NY1,NX2,NY2,NYt08
DOUBLE PRECISION,INTENT(IN)                 :: dr,dz

INTEGER                                     :: i,j
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
! ---------------------------------------------------------
! CALCULATING THE SHEAR RATE ON THE "LEFT" WALL (r=R_i):
! SR2_ij is actually zero since the Dirichlet boundary 
! condition is not chancing with z!
i = NX1
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
! CALCULATING THE SHEAR RATE ON THE BOTTOM PLATE (z=0):
j = 1
DO i=2,NX2-1
  r = DBLE(i-1)*dr
  SR1_ij  = (V(i+1,j) - V(i-1,j))/(2.0D0*dr) - V(i,j)/r
  SR2_ij  = (4.0D0*V(i,j+1) - V(i,j+2) - 3.0D0*V(i,j))/(2.0D0*dz)
  SR(i,j) = DSQRT(SR1_ij**2.0D0 + SR2_ij**2.0D0)
END DO
! ---------------------------------------------------------
! CALCULATING THE SHEAR RATE BELOW THE INNER CYLINDER (z=H3-H2):
j = NY1
DO i=NX1-50,NX1-1
  r = DBLE(i-1)*dr
  SR1_ij  = (V(i+1,j) - V(i-1,j))/(2.0D0*dr) - V(i,j)/r
  SR2_ij  = (-4.0D0*V(i,j-1) + V(i,j-2) + 3.0D0*V(i,j))/(2.0D0*dz)
  SR(i,j) = DSQRT(SR1_ij**2.0D0 + SR2_ij**2.0D0)
END DO
! ---------------------------------------------------------
! CALCULATING THE SHEAR RATE AT THE OPEN BOUNDARY (z=H):
j = NY2
DO i=NX1+1,NX2-1
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
i = NX1
j = NY2
r = DBLE(i-1)*dr
SR1_ij  = ( 4.0D0*V(i+1,j) - V(i+2,j) - 3.0D0*V(i,j))/(2.0D0*dr) - V(i,j)/r
SR2_ij  = (-4.0D0*V(i,j-1) + V(i,j-2) + 3.0D0*V(i,j))/(2.0D0*dz)
SR(i,j) = DSQRT(SR1_ij**2.0D0 + SR2_ij**2.0D0)
! ---------------------------------------------------------
! CALCULATING THE SHEAR RATE AT CORNER OF INNER CYLINDER (r=R_i,z=h_1+h_cone):
i = NX1
j = NY1
r = DBLE(i-1)*dr
SR1_ij  = ( 4.0D0*V(i+1,j) - V(i+2,j) - 3.0D0*V(i,j))/(2.0D0*dr) - V(i,j)/r
SR2_ij  = (-4.0D0*V(i,j-1) + V(i,j-2) + 3.0D0*V(i,j))/(2.0D0*dz)
SR(i,j) = DSQRT(SR1_ij**2.0D0 + SR2_ij**2.0D0)
! --------------------------------------------------------------------------------- !
SR_r = SR(NX1:NX2,NYt08)
! --------------------------------------------------------------------------------- !
RETURN
END SUBROUTINE ROS_PROFILE
! ================================================================================= !
SUBROUTINE TORQUE(V,ETA,NX1,NX2,NY1,NY2,NY2mH3,dr,dz,TORQUE_ON_INNER_CYLINDER,&
                  TORQUE_ON_OUTER_CYLINDER,TORQUE_ON_BOTTOM_PLATE)

DOUBLE PRECISION,DIMENSION(:,:),INTENT(IN) :: V,ETA
DOUBLE PRECISION,DIMENSION(:),INTENT(OUT)  :: TORQUE_ON_INNER_CYLINDER,&
                                              TORQUE_ON_OUTER_CYLINDER,&
                                              TORQUE_ON_BOTTOM_PLATE
INTEGER,INTENT(IN)                         :: NX1,NY1,NX2,NY2,NY2mH3
DOUBLE PRECISION,INTENT(IN)                :: dr,dz

INTEGER          :: i,j,count_j
DOUBLE PRECISION :: shear_stress_theta_r_i,shear_stress_theta_r_o,&
                    shear_stress_theta_r,r,PI
! --------------------------------------------------------------------------------- !
PI = DACOS(-1.0D0)
! --------------------------------------------------------------------------------- !
! The elements "TORQUE_ON_OUTER_CYLINDER", "TORQUE_ON_BOTTOM_PLATE" 
! and "TORQUE_ON_INNER_CYLINDER" are the torque applied ON the 
! corresponding wall boundary FROM the test material. The first two 
! elements "TORQUE_ON_OUTER_CYLINDER" and "TORQUE_ON_BOTTOM_PLATE" are
! used when generating the file "work.dat". There it is desired to 
! gain the torque applied FROM the wall boundary ON the test material.
! As such, when writing data into the file "work.dat", both the torque 
! elements will be written with an opposite sign (i.e. "+" -> "-").
! No corresponding considerations are needed for the toque element
! "TORQUE_ON_INNER_CYLINDER".
! --------------------------------------------------------------------------------- !
! CALCULATING THE TORQUE ON THE RIGHT WALL, FROM THE TEST MATERIAL (r=R_o):
i = NX2
r = DBLE(i-1)*dr
DO j=2,NY2-1
  shear_stress_theta_r_o      = - ETA(i,j)*((-4.0D0*V(i-1,j) + V(i-2,j) &
                                + 3.0D0*V(i,j))/(2.0D0*dr) - V(i,j)/r)
  TORQUE_ON_OUTER_CYLINDER(j) = r*(shear_stress_theta_r_o*dz*(2*PI*r))
END DO
! --------------------------------------------------------------------------------- !
! CALCULATING THE TORQUE ON TOP RIGHT WALL, FROM THE TEST MATERIAL (r=R_o;z=H):
i = NX2
j = NY2
r = DBLE(i-1)*dr
shear_stress_theta_r_o      = - ETA(i,j)*((-4.0D0*V(i-1,j) + V(i-2,j) &
                              + 3.0D0*V(i,j))/(2.0D0*dr) - V(i,j)/r)
TORQUE_ON_OUTER_CYLINDER(j) = r*(shear_stress_theta_r_o*(dz/2)*(2*PI*r))
! ---------------------------------------------------------
! CALCULATING THE TORQUE ON BOTTOM RIGHT WALL, FROM THE TEST MATERIAL (r=R_o;z=0):
i = NX2
j = 1
r = DBLE(i-1)*dr
shear_stress_theta_r_o      = - ETA(i,j)*((-4.0D0*V(i-1,j) + V(i-2,j) &
                              + 3.0D0*V(i,j))/(2.0D0*dr) - V(i,j)/r)
TORQUE_ON_OUTER_CYLINDER(j) = r*(shear_stress_theta_r_o*(dz/2)*(2*PI*r))
! ---------------------------------------------------------
! CALCULATING THE TORQUE ON THE BOTTOM PLATE, FROM THE TEST MATERIAL (z=0):
j = 1
DO i=2,NX2-1
  r = DBLE(i-1)*dr
  shear_stress_theta_r      = ETA(i,j)*((4.0D0*V(i,j+1) - V(i,j+2) &
                              - 3.0D0*V(i,j))/(2.0D0*dz))
  TORQUE_ON_BOTTOM_PLATE(i) = r*(shear_stress_theta_r*dr*(2*PI*r))
END DO
! ---------------------------------------------------------
! CALCULATING THE TORQUE ON THE BOTTOM PLATE - RIGHT CORNER, FROM THE 
! TEST MATERIAL (r=R_o;z=0):
j = 1
i = NX2
r = DBLE(i-1)*dr
shear_stress_theta_r      = ETA(i,j)*((4.0D0*V(i,j+1) - V(i,j+2) &
                            - 3.0D0*V(i,j))/(2.0D0*dz))
TORQUE_ON_BOTTOM_PLATE(i) = r*(shear_stress_theta_r*(dr/2)*(2*PI*r))
! --------------------------------------------------------------------------------- !
! CALCULATING THE TORQUE ON THE LEFT WALL, FROM THE TEST MATERIAL (r=R_i):
! for all z in H3 (i.e. at where torque is measured):
i = NX1
r = DBLE(i-1)*dr
count_j = 0
DO j=NY2mH3,NY2
  count_j = count_j + 1 
  shear_stress_theta_r_i = ETA(i,j)*((4.0D0*V(i+1,j) - V(i+2,j) & 
                           - 3.0D0*V(i,j))/(2.0D0*dr) - V(i,j)/r)
  IF ((j.EQ.NY2mH3).OR.(j.EQ.NY2)) THEN
    TORQUE_ON_INNER_CYLINDER(count_j) = r*(shear_stress_theta_r_i*(dz/2)*(2*PI*r))
  ELSE
    TORQUE_ON_INNER_CYLINDER(count_j) = r*(shear_stress_theta_r_i*dz*(2*PI*r))
  END IF 
END DO
! --------------------------------------------------------------------------------- !
RETURN
END SUBROUTINE TORQUE
! ================================================================================= !
END MODULE WRITE_INFORMATION
! --------------------------------------------------------------------------------- !
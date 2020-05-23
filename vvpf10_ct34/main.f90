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
! File name: main.f90 (PROGRAM)                                                     !
! This is the center of the whole software, holding and passing information to and  !
! from the different subroutines. Some subroutines interact directly with each      !
! other without going through the channels defined by main.f90 (this applies mostly !
! for the subroutines in the files update.f90, shear.f90 and viscous.f90).          !
! The geometry of the viscometer, including the bottom cone, is defined in this     !
! part of the software.                                                             !
! --------------------------------------------------------------------------------- !
PROGRAM MAIN_ROUTINE

USE CONSTANTS_AND_PARAMETERS
USE SHEAR_VISCOSITY
USE ROTATION
USE MATRIX
USE WRITE_INFORMATION

IMPLICIT NONE

DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: MX1,MX2,MY1,MY2
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)   :: K1,K2,L1,L2,DUMMY_2

DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: VELOCITY_k,VELOCITY_kp12,VELOCITY_kp12_new,&
                                               VELOCITY_kp1_new,VELOCITY_kp1,SR,H,FMSR,FMCR

DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)   :: v1r_ijp1k,v1r_ijk,v1r_ijm1k,&
                                               v1r_ijp1kp12,v1r_ijkp12,v1r_ijm1kp12,&
                                               v1r_c_ijp1k,v1r_c_ijk,v1r_c_ijm1k,&
                                               v1r_c_ijp1kp12,v1r_c_ijkp12,v1r_c_ijm1kp12,&
                                               v2r_ijp1k,v2r_ijk,v2r_ijm1k,&
                                               v2r_ijp1kp12,v2r_ijkp12,v2r_ijm1kp12

DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)   :: FMSR1r_ijp1,FMSR1r_ij,FMSR1r_ijm1,&
                                               FMSR1r_c_ijp1,FMSR1r_c_ij,FMSR1r_c_ijm1,&
                                               FMSR2r_ijp1,FMSR2r_ij,FMSR2r_ijm1

DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)   :: FMCR1r_ijp1,FMCR1r_ij,FMCR1r_ijm1,&
                                               FMCR1r_c_ijp1,FMCR1r_c_ij,FMCR1r_c_ijm1,&
                                               FMCR2r_ijp1,FMCR2r_ij,FMCR2r_ijm1

DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)   :: v1r_ijkp12_new,v2r_ijkp12_new

DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)   :: v1z_ip1jkp12,v1z_ijkp12,v1z_im1jkp12,&
                                               v1z_ip1jkp1,v1z_ijkp1,v1z_im1jkp1,&
                                               v1z_c_ip1jkp12,v1z_c_ijkp12,v1z_c_im1jkp12,&
                                               v1z_c_ip1jkp1,v1z_c_ijkp1,v1z_c_im1jkp1,&
                                               v2z_ip1jkp12,v2z_ijkp12,v2z_im1jkp12,&
                                               v2z_ip1jkp1,v2z_ijkp1,v2z_im1jkp1

DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)   :: FMSR2z_ip1j,FMSR2z_ij,FMSR2z_im1j,&
                                               FMSR1z_ip1j,FMSR1z_ij,FMSR1z_im1j,&
                                               FMSR1z_c_ip1j,FMSR1z_c_ij,FMSR1z_c_im1j

DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)   :: FMCR2z_ip1j,FMCR2z_ij,FMCR2z_im1j,&
                                               FMCR1z_ip1j,FMCR1z_ij,FMCR1z_im1j,&
                                               FMCR1z_c_ip1j,FMCR1z_c_ij,FMCR1z_c_im1j

DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)   :: v1z_ijkp1_new,v2z_ijkp1_new
INTEGER,ALLOCATABLE,DIMENSION(:)            :: x_cone,y_cone,y_cone_dynamic


DOUBLE PRECISION :: dr,dz,dt,rho,omega,Lambda,R_i,r_i_cone,R_o,h1,H2,H3,tol,a,b,&
                    dt_Newton,dt_Plastic,tol_Newton,tol_Plastic,tol_RMS,&
                    tol_RMS_active,ZERO_TIME,REAL_TIME,H_cone,h1_static,dt_OUTPUT,&
                    dt_OUTPUT_torque,RMS,vel_norm,small_zero,EPS,ZERO_TIME_tmp,&
                    REAL_TIME_tmp,TIME_INTERVAL,f,f_min,f_max,PERC,time,&
                    memory_beta,memory_alpha,BETA_I,ALPHA_I,SR_tmp,H_tmp,U_o

INTEGER          :: i,j,k,error,problem,NX,NX1,NX2,NY1,NY2,N_Lambda,N_Lambda_MAX,&
                    count,MAX_NUMBER_OF_ITERATIONS,count_max,NX1_Cone,&
                    NX2_Cone,NX_cone,NY1_Cone,NY2_cone,NX2mNX_cone,&
                    NUMBER_OF_TIME_ITERATIONS,N_cone_x,N_cone_y,&
                    NUMBER_OF_POINTS,N_dt,count_rms,NX1_cone_rms,k_OUTPUT_rms,&
                    k_OUTPUT_torque,k_OUTPUT,NY2mH3,DUMMY_1

LOGICAL          :: CONVERGENCE,TIME_INDEPENDENCE,WARNING_SIGN,ConTec_v4,&
                    ConTec_BML_v3,CALCULATE_TIME_DEPENDENT_PROBL,SMOOTH,&
                    FALSE_CONVERGENCE
                               
CHARACTER        :: IGNORED_INPUT
! --------------------------------------------------------------------------------- !
PRINT *, "      _______________________________________________       "
PRINT *, "      Viscometric-ViscoPlastic-Flow v1.0  (CT3 & CT4)       "
PRINT *, "                                                            "
PRINT *, "            Copyright (C) 2002, Jon E. Wallevik,            "
PRINT *, "                (jon.wallevik@bygg.ntnu.no)                 "
PRINT *, " The Norwegian University of Science and Technology (NTNU)  "
PRINT *, "            Department of Structural Engineering            "
PRINT *, " __________________________________________________________ "
PRINT *, "                                                            "
PRINT *, " This software is free software; you can redistribute it    "
PRINT *, " and/or modify it under the terms of the GNU General Public "
PRINT *, " License as published by the Free Software Foundation;      "
PRINT *, " either version 2 of the License, or (at your option) any   "
PRINT *, " later version. This software is distributed in the hope    "
PRINT *, " that it will be useful, but WITHOUT ANY WARRANTY; without  "
PRINT *, " even the implied warranty of MERCHANTABILITY or FITNESS    "
PRINT *, " FOR A PARTICULAR PURPOSE.                                  "
PRINT *, " See the GNU General Public License for more details.       "
PRINT *, " __________________________________________________________ "
PRINT *, "                                                            "
WRITE(*,"(A)",ADVANCE="NO") "                PRESS 'ENTER' TO CONTINUE"
PRINT *, " "
READ (*,"(A)") IGNORED_INPUT
PRINT *, " "
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
! --------------------------------------------------------------------------------- !
IF (ConTec_v4) THEN
  PRINT '(7X,A43)', "==========================================="
  PRINT '(7X,A43)', "      Solving for ConTec Viscometer 4      "
  PRINT '(7X,A43)', "==========================================="
  PRINT '(7X,A43)', "                                           "

  CALL ConTec_CONSTANTS(rho,REAL_TIME,CALCULATE_TIME_DEPENDENT_PROBL,&
                        ZERO_TIME,tol_Newton,tol_Plastic,tol_RMS,&
                        dt_Plastic,dt_Newton,count_max,R_i,R_o,h1,H2,H3)

  CALL VELOCITY_AND_TIME_ConTec(ZERO_TIME_tmp,REAL_TIME_tmp,&
       TIME_INTERVAL,f,f_min,f_max,PERC,NUMBER_OF_POINTS,SMOOTH)

  dr = 0.5D-3 ! => 0.5 mm = Spacing between grid points in r-direction.
  dz = 0.5D-3 ! => 0.5 mm = Spacing between grid points in z-direction.

  H_cone    = DBLE(21 - 5)*dz ! -> Static (non changeable) variable:: 0.8 cm.
  h1_static = DBLE(5 - 1)*dz  ! -> Static (non changeable) variable:: 0.2 cm.

  N_cone_x  = 16              ! -> Static (non changeable) variable 
  N_cone_y  = 39              ! -> Static (non changeable) variable 
  
  ALLOCATE(x_cone(N_cone_x),y_cone(N_cone_y),y_cone_dynamic(N_cone_y),&
           stat=problem)
  IF (problem/=0) THEN 
    PRINT *," MAIN_ROUTINE says: The program could not allocate space! "
    PRINT *," Error code 0 in main and execution terminated!           "
    STOP
  END IF

  ! NX1=171;NX2=203;NY1=21;NY2=301;NX=NX2-NX1+1=203-171+1=33;
  x_cone   = (/9,12,14,17,19,22,24,27,29,32,34,37,39,42,44,47/)
  NX1_cone = x_cone(1)
  NX2_cone = x_cone(16)

  y_cone = (/5,6,6,6,7,7,8,8,8,9,9,10,10,10,11,11,12,12,12,13,13,&
           14,14,14,15,15,16,16,16,17,17,18,18,18,19,19,20,20,20/)
  y_cone_dynamic = y_cone - IDNINT((h1_static - h1)/dz)
  y_cone         = y_cone_dynamic
  NY1_cone       = y_cone(1)
! NY2_cone       = is defined elsewhere!
! --------------------------------------------------------------------------------- !
ELSE IF (ConTec_BML_v3) THEN
  PRINT '(4X,A49)', "================================================="
  PRINT '(4X,A49)', "       Solving for ConTec BML Viscometer 3       "
  PRINT '(4X,A49)', "================================================="
  PRINT '(4X,A49)', "                                                 "
  
  CALL BML_CONSTANTS(rho,REAL_TIME,CALCULATE_TIME_DEPENDENT_PROBL,&
                     ZERO_TIME,tol_Newton,tol_Plastic,tol_RMS,&
                     dt_Plastic,dt_Newton,count_max,R_i,R_o,h1,H2,H3)

  CALL VELOCITY_AND_TIME_BML(ZERO_TIME_tmp,REAL_TIME_tmp,&
       TIME_INTERVAL,f,f_min,f_max,PERC,NUMBER_OF_POINTS,SMOOTH)

  dr = 1.0D-3 ! => 1.0 mm = Spacing between grid points in r-direction.
  dz = 1.0D-3 ! => 1.0 mm = Spacing between grid points in z-direction.

  H_cone    = DBLE(61 - 21)*dz ! -> Static (non changeable) variable:: 4.0 cm.
  h1_static = DBLE(21 - 1)*dz  ! -> Static (non changeable) variable:: 2.0 cm.

  N_cone_x = 40
  N_cone_y = 27

  ALLOCATE(x_cone(N_cone_x),y_cone(N_cone_y),y_cone_dynamic(N_cone_y),&
           stat=problem)
  IF (problem/=0) THEN 
    PRINT *," MAIN_ROUTINE says: The program could not allocate space! "
    PRINT *," Error code 0 in main and execution terminated!           "
    STOP
  END IF

  ! NX1=101;NX2=146;NY1=61;NY2=306;NX=NX2-NX1+1=146-101+1=46;
  x_cone   = (/1,2,2,3,4,4,5,6,6,7,8,8,9,10,10,11,12,12,13,14,14,15,&
             16,16,17,18,18,19,20,20,21,22,22,23,24,24,25,26,26,27/)
  NX1_cone = x_cone(1)
  NX2_cone = x_cone(40)

  y_cone   = (/21,22,24,25,27,28,30,31,33,34,36,37,39,40,42,43,&
             45,46,48,49,51,52,54,55,57,58,60/)
  y_cone_dynamic = y_cone - IDNINT((h1_static - h1)/dz)
  y_cone         = y_cone_dynamic
  NY1_cone       = y_cone(1)
! NY2_cone       = is defined elsewhere!

END IF
! --------------------------------------------------------------------------------- !
NX1    = IDNINT(R_i/dr) + 1 ! 171 for CT and 101 for BML
NX2    = IDNINT(R_o/dr) + 1 ! 203 for CT and 146 for BML
NX     = NX2 - NX1 + 1      ! NX=203-171+1=33 for CT and NX=146-101+1=46 for BML
NY1    = IDNINT((h1 + H_cone)/dz) + 1      ! = 21  for CT and 61  for BML
NY2    = IDNINT((h1 + H_cone + H2)/dz) + 1 ! = 301 for CT and 306 for BML 
NY2mH3 = NY2 - IDNINT(H3/dz)

R_i = dr*DBLE(NX1-1)
R_o = dr*DBLE(NX2-1)
! --------------------------------------------------------------------------------- !
CALL ANGULAR_VELOCITY(0.0D0,dt,omega)
! ---------------------------------------------------------
14 FORMAT(10X,"NX1 = ",(I3,2X),"; NX2 = ",(I3,2X),"; dr = ",F6.4,"m")
15 FORMAT(10X,"R_i = ",F6.4,"m  ; R_o = ",F6.4,"m")
19 FORMAT(10X,"NY2mH3 = ",(I3,2X)," ; NX = ",(I3,2X))
16 FORMAT(10X,"NY1 = ",(I3,2X),"; NY2 = ",(I3,2X),"; dz = ",F6.4,"m")
17 FORMAT(10X,"dt_Plastic = ",E9.3,"s; f_o = ",F6.4,"rps")
PRINT '(8X,A26)',"Geometric and time values:"
PRINT 14, NX1,NX2,dr
PRINT 16, NY1,NY2,dz
PRINT 19, NY2mH3,NX
PRINT 15, R_i,R_o
PRINT 17, dt_Plastic,omega/(2.0D0*ACOS(-1.0D0))
PRINT *, "    "
! --------------------------------------------------------------------------------- !
IF (NX1.GE.NX2) THEN
  PRINT *, " Inner radius 'R_i' is larger than the outer radius 'R_o'! "
  PRINT *, "              TERMINAL ERROR!                              "
  STOP
END IF

IF ((NX2_cone+2).GE.NX1) THEN
  PRINT *, " Inner radius 'R_i' is smaller than the largest radius "
  PRINT *, " of the bottom cone r_cone^max=((NX2_cone+2)-1)*dr!    "
  PRINT *, "              TERMINAL ERROR!                          "
  STOP
END IF

IF (h1.LT.2.0D0*dz) THEN
  PRINT *, " The variable 'h1' in 'param.f90' is too small! "
  PRINT *, " Minimum height of h1 must be 2*dz,             "
  PRINT *, " otherwise a logical error will occur.          "
  PRINT *, "  [h1(ConTec) = 0.001 (0.1 cm) ; dz=0.5 mm];    "
  PRINT *, "  [h1(BML)    = 0.002 (0.2 cm) ; dz=1.0 mm];    " 
  PRINT *, "              TERMINAL ERROR!                   "
  STOP
END IF

IF (H2.LT.2.0D0*dz) THEN
  PRINT *, " The variable 'H2' in 'param.f90' is too small! "
  PRINT *, " Minimum height of h1 must be 2*dz,             "
  PRINT *, " otherwise a logical error will occur.          "
  PRINT *, "  [h1(ConTec) = 0.001 (0.1 cm) ; dz=0.5 mm];    "
  PRINT *, "  [h1(BML)    = 0.002 (0.2 cm) ; dz=1.0 mm];    " 
  PRINT *, "              TERMINAL ERROR!                   "
  STOP
END IF
! --------------------------------------------------------------------------------- !
CALL WARNING_FOR_WRITING(NY2)

MAX_NUMBER_OF_ITERATIONS  = IDNINT(ZERO_TIME/dt_Plastic)   
NUMBER_OF_TIME_ITERATIONS = DINT(REAL_TIME/dt_Plastic) 

PRINT *, "      -------------------------------------------  "
PRINT "( '         MAX_NUMBER_OF_ITERATIONS: ', I10 )        ",&
                   MAX_NUMBER_OF_ITERATIONS
PRINT *, "      -------------------------------------------  "
PRINT *, "                                                   "

IF (CALCULATE_TIME_DEPENDENT_PROBL) THEN
  PRINT *, "      -------------------------------------------  "
  PRINT "( '         NUMBER_OF_TIME_ITERATIONS:       ', I10 ) ",&
                     NUMBER_OF_TIME_ITERATIONS                
  PRINT *, "      -------------------------------------------  "
  PRINT *, "                                                   "
END IF

PRINT *, "              ___________________________                   "
WRITE(*,"(A)",ADVANCE="NO") "                PRESS 'ENTER' TO CONTINUE"
PRINT *, " "
READ (*,"(A)") IGNORED_INPUT
PRINT *, " "
! --------------------------------------------------------------------------------- !
k_OUTPUT_rms = 25     ! -> Information output every dt_OUTPUT_rms 
                      !    times (to console and file).
! --------------------------------------------------------------------------------- !
! The term "small_zero" does usually not have to be changed.
small_zero = 0.1D-7  ! -> Used in relation to screen and file output.
EPS        = 1.0D-15 ! -> Used in relation to vel_norm.
! --------------------------------------------------------------------------------- !
DUMMY_1    = IDNINT(0.8D0*DBLE(NY2))
! --------------------------------------------------------------------------------- !
! Creating log file and making the first entry:
OPEN(unit=8,file="log.dat",status="replace",action="write",&
     position="rewind",iostat=problem)
IF (problem/=0) THEN
  PRINT *," Could not create the file: log.dat! "
  STOP
ELSE
  WRITE (unit=8,fmt=*) 0,0.0D0 
END IF 

CLOSE (UNIT=8)
! --------------------------------------------------------------------------------- !
ALLOCATE(MX1(NX2-2,NX2-2),MX2(NX-2,NX-2),MY1(NY1-2,NY1-2),MY2(NY2-2,NY2-2),&
         VELOCITY_k(NX2,NY2),VELOCITY_kp12(NX2,NY2),VELOCITY_kp12_new(NX2,NY2),&
         VELOCITY_kp1(NX2,NY2),VELOCITY_kp1_new(NX2,NY2),SR(NX2,NY2),H(NX2,NY2),&
         FMSR(NX2,NY2),FMCR(NX2,NY2),stat=problem)
IF (problem/=0) THEN 
  PRINT *," MAIN_ROUTINE says: The program could not allocate space! "
  PRINT *," Error code 1 in main and execution terminated!           "
  STOP
END IF

ALLOCATE(K1(NX2-2),K2(NX-2),L1(NY1-2),L2(NY2-2),DUMMY_2(NX),stat=problem)
IF (problem/=0) THEN 
  PRINT *," MAIN_ROUTINE says: The program could not allocate space! "
  PRINT *," Error code 2 in main and execution terminated!           "
  STOP
END IF

ALLOCATE(v1r_ijp1k(NX2),v1r_ijk(NX2),v1r_ijm1k(NX2),v1r_ijp1kp12(NX2),&
         v1r_ijkp12(NX2),v1r_ijm1kp12(NX2),v1r_c_ijp1k(NX2),v1r_c_ijk(NX2),&
         v1r_c_ijm1k(NX2),v1r_c_ijp1kp12(NX2),v1r_c_ijkp12(NX2),&
         v1r_c_ijm1kp12(NX2),v2r_ijp1k(NX),v2r_ijk(NX),v2r_ijm1k(NX),&
         v2r_ijp1kp12(NX),v2r_ijkp12(NX),v2r_ijm1kp12(NX),stat=problem)
IF (problem/=0) THEN 
  PRINT *," MAIN_ROUTINE says: The program could not allocate space! "
  PRINT *," Error code 3a in main and execution terminated!          "
  STOP
END IF

ALLOCATE(FMSR1r_ijp1(NX2),FMSR1r_ij(NX2),FMSR1r_ijm1(NX2),&
         FMSR1r_c_ijp1(NX2),FMSR1r_c_ij(NX2),FMSR1r_c_ijm1(NX2),&
         FMSR2r_ijp1(NX),FMSR2r_ij(NX),FMSR2r_ijm1(NX),stat=problem)
IF (problem/=0) THEN 
  PRINT *," MAIN_ROUTINE says: The program could not allocate space! "
  PRINT *," Error code 3b in main and execution terminated!          "
  STOP
END IF

ALLOCATE(FMCR1r_ijp1(NX2),FMCR1r_ij(NX2),FMCR1r_ijm1(NX2),&
         FMCR1r_c_ijp1(NX2),FMCR1r_c_ij(NX2),FMCR1r_c_ijm1(NX2),&
         FMCR2r_ijp1(NX),FMCR2r_ij(NX),FMCR2r_ijm1(NX),stat=problem)
IF (problem/=0) THEN 
  PRINT *," MAIN_ROUTINE says: The program could not allocate space! "
  PRINT *," Error code 3c in main and execution terminated!          "
  STOP
END IF

ALLOCATE(v1r_ijkp12_new(NX2-2),v2r_ijkp12_new(NX-2),stat=problem)
IF (problem/=0) THEN 
  PRINT *," MAIN_ROUTINE says: The program could not allocate space! "
  PRINT *," Error code 4 in main and execution terminated!           "
  STOP
END IF

ALLOCATE(v1z_ip1jkp12(NY1),v1z_ijkp12(NY1),v1z_im1jkp12(NY1),v1z_ip1jkp1(NY1),&
         v1z_ijkp1(NY1),v1z_im1jkp1(NY1),v1z_c_ip1jkp12(NY1),v1z_c_ijkp12(NY1),&
         v1z_c_im1jkp12(NY1),v1z_c_ip1jkp1(NY1),v1z_c_ijkp1(NY1),&
         v1z_c_im1jkp1(NY1),v2z_ip1jkp12(NY2),v2z_ijkp12(NY2),v2z_im1jkp12(NY2),&
         v2z_ip1jkp1(NY2),v2z_ijkp1(NY2),v2z_im1jkp1(NY2),stat=problem)
IF (problem/=0) THEN 
  PRINT *," MAIN_ROUTINE says: The program could not allocate space! "
  PRINT *," Error code 5a in main and execution terminated!          "
  STOP
END IF

ALLOCATE(FMSR1z_ip1j(NY1),FMSR1z_ij(NY1),FMSR1z_im1j(NY1),&
         FMSR1z_c_ip1j(NY1),FMSR1z_c_ij(NY1),FMSR1z_c_im1j(NY1),&
         FMSR2z_ip1j(NY2),FMSR2z_ij(NY2),FMSR2z_im1j(NY2),stat=problem)
IF (problem/=0) THEN 
  PRINT *," MAIN_ROUTINE says: The program could not allocate space! "
  PRINT *," Error code 5b in main and execution terminated!          "
  STOP
END IF

ALLOCATE(FMCR1z_ip1j(NY1),FMCR1z_ij(NY1),FMCR1z_im1j(NY1),&
         FMCR1z_c_ip1j(NY1),FMCR1z_c_ij(NY1),FMCR1z_c_im1j(NY1),&
         FMCR2z_ip1j(NY2),FMCR2z_ij(NY2),FMCR2z_im1j(NY2),stat=problem)
IF (problem/=0) THEN 
  PRINT *," MAIN_ROUTINE says: The program could not allocate space! "
  PRINT *," Error code 5c in main and execution terminated!          "
  STOP
END IF

ALLOCATE(v1z_ijkp1_new(NY1-2),v2z_ijkp1_new(NY2-2),stat=problem)
IF (problem/=0) THEN 
  PRINT *," MAIN_ROUTINE says: The program could not allocate space! "
  PRINT *," Error code 5d in main and execution terminated!          "
  STOP
END IF
! --------------------------------------------------------------------------------- !
! Initialization: 
! ================================================================================= !
SR                = 0.0D0
H                 = 0.0D0
FMSR              = 0.0D0
FMCR              = 0.0D0
! ================================================================================= !
MX1               = 0.0D0
MX2               = 0.0D0
K1                = 0.0D0
K2                = 0.0D0

VELOCITY_k        = 0.0D0
VELOCITY_kp12     = 0.0D0
VELOCITY_kp12_new = 0.0D0
! --------------------------------------------------------------------------------- !
v1r_ijkp12_new    = 0.0D0
v2r_ijkp12_new    = 0.0D0

! v1r is used in K1 and MX1 -> K1=K1(v1r) and MX1=MX1(v1r) to solve the system 
! MX1*v1r_new=K1. "v1r" could be called "v1r_old" since it is the velocity 
! from the previous iteration.
v1r_ijp1k         = 0.0D0  ! v1r -> K1 & MX1*v1r_new=K1
v1r_ijk           = 0.0D0
v1r_ijm1k         = 0.0D0
v1r_ijp1kp12      = 0.0D0
v1r_ijkp12        = 0.0D0
v1r_ijm1kp12      = 0.0D0
! --------------------------------------------------------------------------------- !
FMSR1r_ijp1       = 0.0D0
FMSR1r_ij         = 0.0D0
FMSR1r_ijm1       = 0.0D0

FMCR1r_ijp1       = 0.0D0
FMCR1r_ij         = 0.0D0
FMCR1r_ijm1       = 0.0D0
! --------------------------------------------------------------------------------- !
v1r_c_ijp1k       = 0.0D0  ! v1r_c -> K1(1:x) & MX1(1:x)*v1r_c_new(1:x)=K1(1:x)
v1r_c_ijk         = 0.0D0  ! ..._c -> ..._cone
v1r_c_ijm1k       = 0.0D0
v1r_c_ijp1kp12    = 0.0D0
v1r_c_ijkp12      = 0.0D0
v1r_c_ijm1kp12    = 0.0D0
! --------------------------------------------------------------------------------- !
FMSR1r_c_ijp1     = 0.0D0
FMSR1r_c_ij       = 0.0D0
FMSR1r_c_ijm1     = 0.0D0

FMCR1r_c_ijp1     = 0.0D0
FMCR1r_c_ij       = 0.0D0
FMCR1r_c_ijm1     = 0.0D0
! --------------------------------------------------------------------------------- !
v2r_ijp1k         = 0.0D0  ! v2r -> K2 & MX2*v2r_new=K2  
v2r_ijk           = 0.0D0
v2r_ijm1k         = 0.0D0
v2r_ijp1kp12      = 0.0D0
v2r_ijkp12        = 0.0D0
v2r_ijm1kp12      = 0.0D0
! --------------------------------------------------------------------------------- !
FMSR2r_ijp1       = 0.0D0
FMSR2r_ij         = 0.0D0
FMSR2r_ijm1       = 0.0D0

FMCR2r_ijp1       = 0.0D0
FMCR2r_ij         = 0.0D0
FMCR2r_ijm1       = 0.0D0
! ================================================================================= !
MY1               = 0.0D0
MY2               = 0.0D0
L1                = 0.0D0
L2                = 0.0D0

VELOCITY_kp1      = 0.0D0
VELOCITY_kp1_new  = 0.0D0
! --------------------------------------------------------------------------------- !
v1z_ijkp1_new     = 0.0D0
v2z_ijkp1_new     = 0.0D0

! v1z is used in L1 and MY1 -> L1=L1(v1z) and MY1=MY1(v1z) to solve the system 
! MY1*v1z_new=L1. "v1z" could be called "v1z_old" since it is the velocity 
! from the previous iteration.
v1z_ip1jkp12      = 0.0D0  ! v1z -> L1 & MY1*v1z_new=L1  
v1z_ijkp12        = 0.0D0
v1z_im1jkp12      = 0.0D0
v1z_ip1jkp1       = 0.0D0
v1z_ijkp1         = 0.0D0
v1z_im1jkp1       = 0.0D0
! --------------------------------------------------------------------------------- !
FMSR1z_ip1j       = 0.0D0
FMSR1z_ij         = 0.0D0
FMSR1z_im1j       = 0.0D0

FMCR1z_ip1j       = 0.0D0
FMCR1z_ij         = 0.0D0
FMCR1z_im1j       = 0.0D0
! --------------------------------------------------------------------------------- !
v1z_c_ip1jkp12    = 0.0D0  ! v1z_c -> L1(1:y) & MY1(1:y)*v1z_c_new(1:y)=L1(1:y)
v1z_c_ijkp12      = 0.0D0  ! ..._c -> ..._cone
v1z_c_im1jkp12    = 0.0D0
v1z_c_ip1jkp1     = 0.0D0
v1z_c_ijkp1       = 0.0D0
v1z_c_im1jkp1     = 0.0D0
! --------------------------------------------------------------------------------- !
FMSR1z_c_ip1j     = 0.0D0
FMSR1z_c_ij       = 0.0D0
FMSR1z_c_im1j     = 0.0D0

FMCR1z_c_ip1j     = 0.0D0
FMCR1z_c_ij       = 0.0D0
FMCR1z_c_im1j     = 0.0D0
! --------------------------------------------------------------------------------- !
v2z_ip1jkp12      = 0.0D0  ! v2z -> L2 & MY2*v2z_new=L2
v2z_ijkp12        = 0.0D0
v2z_im1jkp12      = 0.0D0
v2z_ip1jkp1       = 0.0D0
v2z_ijkp1         = 0.0D0
v2z_im1jkp1       = 0.0D0
! --------------------------------------------------------------------------------- !
FMSR2z_ip1j       = 0.0D0
FMSR2z_ij         = 0.0D0
FMSR2z_im1j       = 0.0D0

FMCR2z_ip1j       = 0.0D0
FMCR2z_ij         = 0.0D0
FMCR2z_im1j       = 0.0D0
! ================================================================================= !
WARNING_SIGN      = .FALSE. 
FALSE_CONVERGENCE = .FALSE.
! --------------------------------------------------------------------------------- !
! Initialization of boundary condition at t = 0.0 sec:
CALL ANGULAR_VELOCITY(0.0D0,dt,omega)
! Dirichlet boundary condition:
DO i = 1,NX2
  VELOCITY_k(i,1) = omega*DBLE(i-1)*dr
END DO
VELOCITY_k(1:NX1,NY1)       = 0.0D0
VELOCITY_k(1,1:NY1)         = 0.0D0
VELOCITY_k(NX1,NY1:NY2)     = 0.0D0
VELOCITY_k(NX2,:)           = R_o*omega

! ################################################################################# !
! In Section 7.11.1 is a detailed description of the algorithm, which is used in    !
! the following.                                                                    !
! --------------------------------------------------------------------------------- !
! Linear approximation to speed up convergence:
DO i = 2,NX1_cone+1
  a = VELOCITY_k(i,1)
  b = VELOCITY_k(i,NY1_cone)
  DO j = 2,NY1_cone-1
    VELOCITY_k(i,j) = a - (a - b)*DBLE(j-1)/DBLE(NY1_cone-1)
  END DO
END DO

DO i = NX1_cone+2,NX2_cone+1
  a = VELOCITY_k(i,1)
  b = VELOCITY_k(i,y_cone(i-NX1_cone))
  DO j = 2,y_cone(i-NX1_cone)-1
    VELOCITY_k(i,j) = a - (a - b)*DBLE(j-1)/DBLE(y_cone(i-NX1_cone)-1)
  END DO
END DO

DO i = NX2_cone+2,NX1
  a = VELOCITY_k(i,1)
  b = VELOCITY_k(i,NY1)
  DO j = 2,NY1-1
    VELOCITY_k(i,j) = a - (a - b)*DBLE(j-1)/DBLE(NY1-1)
  END DO
END DO

DO j = 2,NY2-1
  a = VELOCITY_k(NX1,j)
  b = VELOCITY_k(NX2,j)
  DO i = NX1+1,NX2-1
    VELOCITY_k(i,j) = a - (a - b)*DBLE(i-NX1)/DBLE(NX2-NX1)
  END DO
END DO

! Neumann boundary condition:
VELOCITY_k(NX1+1:NX2-1,NY2) = (4.0D0 * &
     VELOCITY_k(NX1+1:NX2-1,NY2-1) - &
     VELOCITY_k(NX1+1:NX2-1,NY2-2))/3.0D0

! CHECK OUT IF VELOCITY_k IS OK:
! CALL WRITE2FILE_k(VELOCITY_k,NX2)
! STOP
! ================================================================================= !
! ============================ Begin of CONTINUATION ============================== !
! ================================================================================= !
N_Lambda_MAX = 1
CONTINUATION: DO N_Lambda = 0,N_Lambda_MAX
! Lambda => The Continuation Method (see Section 7.8).
Lambda = DBLE(N_Lambda)/DBLE(N_Lambda_MAX)
PRINT *,"________________________________________________________"
PRINT *,"CONTINUATION:",Lambda

IF (N_Lambda == 0) THEN
  dt  = dt_Newton
  tol = tol_Newton
  tol_RMS_active = tol_Newton  ! -> See Equation 7.75.
ELSE
  dt  = dt_Plastic
  tol = tol_Plastic
  tol_RMS_active = tol_RMS     ! -> See Equation 7.75.
END IF

TIME_INDEPENDENCE = .FALSE.
! Initializing time for each CONTINUATION step:
k = 0
! ================================================================================= !
! =========================== Begin of the time loop ============================== !
! ================================================================================= !
ZERO_TIME_LOOP: DO WHILE (.NOT.TIME_INDEPENDENCE)
CONVERGENCE = .FALSE.
IF (ABS(MOD((k+1),k_OUTPUT_rms)).LT.small_zero) THEN
  PRINT *,"______________________________________________________"
  PRINT *,"                                                      "
  PRINT *,"  PSEUDO-TRANSIENT time step: k+1 = ",k+1
  PRINT *,"------------------------------------------------------"
END IF
! --------------------------------------------------------------------------------- !
! VELOCITY_kp12     = VELOCITY_k
! VELOCITY_kp1      = VELOCITY_k
! VELOCITY_kp1_new  = VELOCITY_kp1
! The following routines are to update the Dirichlet and Neumann boundary 
! conditions for the time step k+1/2 and k+1. Most of these routines are 
! redundant since the boundary conditions are not changing with time. However 
! it is a good practice to include them, if by some unfortunate accident 
! some of the boundary values are overwritten.
! --------------------------------------------------------------------------------- !
! A guess for the time step k+1/2
VELOCITY_kp12 = VELOCITY_k

! Updating boundary condition at k+1/2:
CALL ANGULAR_VELOCITY(0.0D0,dt,omega)

! Dirichlet boundary condition:
DO i = 1,NX2
  VELOCITY_kp12(i,1) = omega*DBLE(i-1)*dr
END DO
VELOCITY_kp12(1:NX1,NY1)       = 0.0D0
VELOCITY_kp12(1,1:NY1)         = 0.0D0
VELOCITY_kp12(NX1,NY1:NY2)     = 0.0D0
VELOCITY_kp12(NX2,:)           = R_o*omega

! Neumann boundary condition:
VELOCITY_kp12(NX1+1:NX2-1,NY2) = (4.0D0 * &
     VELOCITY_kp12(NX1+1:NX2-1,NY2-1) - &
     VELOCITY_kp12(NX1+1:NX2-1,NY2-2))/3.0D0

! Also updating boundary condition for "..._new":
VELOCITY_kp12_new = VELOCITY_kp12
! --------------------------------------------------------------------------------- !
! A guess for the time step k+1
VELOCITY_kp1 = VELOCITY_k

! Updating boundary condition at k+1:
CALL ANGULAR_VELOCITY(0.0D0,dt,omega)

! Dirichlet boundary condition:
DO i = 1,NX2
  VELOCITY_kp1(i,1) = omega*DBLE(i-1)*dr
END DO
VELOCITY_kp1(1:NX1,NY1)       = 0.0D0
VELOCITY_kp1(1,1:NY1)         = 0.0D0
VELOCITY_kp1(NX1,NY1:NY2)     = 0.0D0
VELOCITY_kp1(NX2,:)           = R_o*omega

! Neumann boundary condition:
VELOCITY_kp1(NX1+1:NX2-1,NY2) = (4.0D0 * &
     VELOCITY_kp1(NX1+1:NX2-1,NY2-1) - &
     VELOCITY_kp1(NX1+1:NX2-1,NY2-2))/3.0D0

! Also updating boundary condition for "..._new":
VELOCITY_kp1_new = VELOCITY_kp1
! --------------------------------------------------------------------------------- !
count = 0
! ================================================================================= !
! ====================== BEGIN OF SUCCESSIVE SUBSTITUTION ========================= !
! ================================================================================= !
! The iteration loop here is because of the non-linearity of the governing
! Equations 7.22 and 7.23. To come around this problem, the successive substitution
! approach is used (see Section 7.8).
CONVERGE: DO WHILE (.NOT.CONVERGENCE)

! If convergence is a problem, then this might help:
! VELOCITY_kp12 = (VELOCITY_k + VELOCITY_kp1_new)/2
count = count + 1
10 FORMAT(4X,"Successive substitution number = ",1(I3,1X))

IF (ABS(MOD((k+1),k_OUTPUT_rms)).LT.small_zero) THEN
  PRINT 10, count
END IF

!======================== BEGIN OF X-ITERATION =========================
! Iteration is made along r-direction (i.e. along the i-direction as 
! in A(i,j)). It starts at the bottom of the viscometer i=(2:NX2-1) at 
! j = 2 and then move upward with increasing j (see Figures 8.1 and 8.2).

DO j = 2,NY1_cone-1
  v1r_ijp1k    = VELOCITY_k(:,j+1)
  v1r_ijk      = VELOCITY_k(:,j)
  v1r_ijm1k    = VELOCITY_k(:,j-1)
  v1r_ijp1kp12 = VELOCITY_kp12(:,j+1)
  v1r_ijkp12   = VELOCITY_kp12(:,j)
  v1r_ijm1kp12 = VELOCITY_kp12(:,j-1)

  FMSR1r_ijp1  = FMSR(:,j+1)
  FMSR1r_ij    = FMSR(:,j)
  FMSR1r_ijm1  = FMSR(:,j-1)

  FMCR1r_ijp1  = FMCR(:,j+1)
  FMCR1r_ij    = FMCR(:,j)
  FMCR1r_ijm1  = FMCR(:,j-1)

  CALL MATRIX_UPDATE_X(rho,k,dt,Lambda,dr,0.0D0,dz,NX2,v1r_ijp1k,v1r_ijk,&
                       v1r_ijm1k,v1r_ijp1kp12,v1r_ijkp12,v1r_ijm1kp12,&
                       FMSR1r_ijp1,FMSR1r_ij,FMSR1r_ijm1,&
                       FMCR1r_ijp1,FMCR1r_ij,FMCR1r_ijm1,MX1,K1)

  CALL MATRIX_SOLVER(MX1,K1,v1r_ijkp12_new,NX2-2)
  VELOCITY_kp12_new(2:NX2-1,j) = v1r_ijkp12_new
END DO

! -------------------------- BEGIN OF X-CONE ------------------------
DO j = NY1_cone,NY1-1
  NX_cone                       = x_cone(j + 1 - NY1_cone)
  r_i_cone                      = NX_cone*dr
  NX2mNX_cone                   = NX2 - NX_cone
  v1r_c_ijp1k(1:NX2mNX_cone)    = VELOCITY_k(NX_cone+1:NX2,j+1)
  v1r_c_ijk(1:NX2mNX_cone)      = VELOCITY_k(NX_cone+1:NX2,j)
  v1r_c_ijm1k(1:NX2mNX_cone)    = VELOCITY_k(NX_cone+1:NX2,j-1)
  v1r_c_ijp1kp12(1:NX2mNX_cone) = VELOCITY_kp12(NX_cone+1:NX2,j+1)
  v1r_c_ijkp12(1:NX2mNX_cone)   = VELOCITY_kp12(NX_cone+1:NX2,j)
  v1r_c_ijm1kp12(1:NX2mNX_cone) = VELOCITY_kp12(NX_cone+1:NX2,j-1)

  FMSR1r_c_ijp1(1:NX2mNX_cone)  = FMSR(NX_cone+1:NX2,j+1)
  FMSR1r_c_ij(1:NX2mNX_cone)    = FMSR(NX_cone+1:NX2,j)
  FMSR1r_c_ijm1(1:NX2mNX_cone)  = FMSR(NX_cone+1:NX2,j-1)

  FMCR1r_c_ijp1(1:NX2mNX_cone)  = FMCR(NX_cone+1:NX2,j+1)
  FMCR1r_c_ij(1:NX2mNX_cone)    = FMCR(NX_cone+1:NX2,j)
  FMCR1r_c_ijm1(1:NX2mNX_cone)  = FMCR(NX_cone+1:NX2,j-1)

  CALL MATRIX_UPDATE_X(rho,k,dt,Lambda,dr,r_i_cone,dz,NX2mNX_cone,v1r_c_ijp1k,v1r_c_ijk,&
                       v1r_c_ijm1k,v1r_c_ijp1kp12,v1r_c_ijkp12,v1r_c_ijm1kp12,&
                       FMSR1r_c_ijp1,FMSR1r_c_ij,FMSR1r_c_ijm1,&
                       FMCR1r_c_ijp1,FMCR1r_c_ij,FMCR1r_c_ijm1,MX1,K1)

  CALL MATRIX_SOLVER(MX1,K1,v1r_ijkp12_new,NX2mNX_cone-2)
! CALL WRITE2FILE_debug(MX1,K1,v1r_ijkp12_new,NX2mNX_cone-2)
! STOP
  VELOCITY_kp12_new(NX_cone+2:NX2-1,j) = v1r_ijkp12_new(1:NX2mNX_cone-2)
END DO
! --------------------------- END OF X-CONE -------------------------

DO j = NY1,NY2-1
  v2r_ijp1k    = VELOCITY_k(NX1:NX2,j+1)
  v2r_ijk      = VELOCITY_k(NX1:NX2,j)
  v2r_ijm1k    = VELOCITY_k(NX1:NX2,j-1)
  v2r_ijp1kp12 = VELOCITY_kp12(NX1:NX2,j+1)
  v2r_ijkp12   = VELOCITY_kp12(NX1:NX2,j)
  v2r_ijm1kp12 = VELOCITY_kp12(NX1:NX2,j-1)

  FMSR2r_ijp1  = FMSR(NX1:NX2,j+1)
  FMSR2r_ij    = FMSR(NX1:NX2,j)
  FMSR2r_ijm1  = FMSR(NX1:NX2,j-1)

  FMCR2r_ijp1  = FMCR(NX1:NX2,j+1)
  FMCR2r_ij    = FMCR(NX1:NX2,j)
  FMCR2r_ijm1  = FMCR(NX1:NX2,j-1)

  CALL MATRIX_UPDATE_X(rho,k,dt,Lambda,dr,R_i,dz,NX,v2r_ijp1k,v2r_ijk,&
                       v2r_ijm1k,v2r_ijp1kp12,v2r_ijkp12,v2r_ijm1kp12,&
                       FMSR2r_ijp1,FMSR2r_ij,FMSR2r_ijm1,&
                       FMCR2r_ijp1,FMCR2r_ij,FMCR2r_ijm1,MX2,K2)

  CALL MATRIX_SOLVER(MX2,K2,v2r_ijkp12_new,NX-2)
  VELOCITY_kp12_new(NX1+1:NX2-1,j) = v2r_ijkp12_new
END DO

j = NY2
v2r_ijp1k    = VELOCITY_k(NX1:NX2,j-1)    ! => v(i,j+1) = v(i,j-1)
v2r_ijk      = VELOCITY_k(NX1:NX2,j)
v2r_ijm1k    = VELOCITY_k(NX1:NX2,j-1)
v2r_ijp1kp12 = VELOCITY_kp12(NX1:NX2,j-1) ! => v(i,j+1) = v(i,j-1)
v2r_ijkp12   = VELOCITY_kp12(NX1:NX2,j)
v2r_ijm1kp12 = VELOCITY_kp12(NX1:NX2,j-1)

FMSR2r_ijp1  = FMSR(NX1:NX2,j-1)          ! => v(i,j+1) = v(i,j-1) 
FMSR2r_ij    = FMSR(NX1:NX2,j)
FMSR2r_ijm1  = FMSR(NX1:NX2,j-1)

FMCR2r_ijp1  = FMCR(NX1:NX2,j-1)          ! => v(i,j+1) = v(i,j-1) 
FMCR2r_ij    = FMCR(NX1:NX2,j)
FMCR2r_ijm1  = FMCR(NX1:NX2,j-1)

CALL MATRIX_UPDATE_X(rho,k,dt,Lambda,dr,R_i,dz,NX,v2r_ijp1k,v2r_ijk,&
                     v2r_ijm1k,v2r_ijp1kp12,v2r_ijkp12,v2r_ijm1kp12,&
                     FMSR2r_ijp1,FMSR2r_ij,FMSR2r_ijm1,&
                     FMCR2r_ijp1,FMCR2r_ij,FMCR2r_ijm1,MX2,K2)

CALL MATRIX_SOLVER(MX2,K2,v2r_ijkp12_new,NX-2)
VELOCITY_kp12_new(NX1+1:NX2-1,j) = v2r_ijkp12_new

! --------- PAUSE FOR DEBUGGING --------- 
! CALL WRITE2FILE_k(VELOCITY_kp12_new,NX2)
! WRITE(*,"(A)",ADVANCE="NO") " PRESS 'ENTER' TO CONTINUE "
! PRINT *, " "
! READ (*,"(A)") IGNORED_INPUT
! PRINT *, " "

!========================= END OF X-ITERATION =========================

! Updating ..._kp12:
VELOCITY_kp12 = VELOCITY_kp12_new

!======================== BEGIN OF Y-ITERATION ========================
! Iteration is made along z-direction (i.e. along the j-direction as 
! in A(i,j)). It starts at the right side of the viscometer j=(2:NY2-1) at 
! i = NX2-1 and then moves to the left with decreasing i 
! (see Figures 8.1 and 8.2).
DO i = NX2-1,NX1+1,-1
  v2z_ip1jkp12 = VELOCITY_kp12(i+1,:)
  v2z_ijkp12   = VELOCITY_kp12(i,:)
  v2z_im1jkp12 = VELOCITY_kp12(i-1,:)
  v2z_ip1jkp1  = VELOCITY_kp1(i+1,:)
  v2z_ijkp1    = VELOCITY_kp1(i,:)
  v2z_im1jkp1  = VELOCITY_kp1(i-1,:)

  FMSR2z_ip1j  = FMSR(i+1,:)
  FMSR2z_ij    = FMSR(i,:)
  FMSR2z_im1j  = FMSR(i-1,:)
  
  FMCR2z_ip1j  = FMCR(i+1,:)
  FMCR2z_ij    = FMCR(i,:)
  FMCR2z_im1j  = FMCR(i-1,:)

  CALL MATRIX_UPDATE_Y(rho,k,dt,Lambda,dr,i,dz,NY2,v2z_ip1jkp12,v2z_ijkp12,&
                       v2z_im1jkp12,v2z_ip1jkp1,v2z_ijkp1,v2z_im1jkp1,&
                       FMSR2z_ip1j,FMSR2z_ij,FMSR2z_im1j,&
                       FMCR2z_ip1j,FMCR2z_ij,FMCR2z_im1j,MY2,L2,.TRUE.)

  CALL MATRIX_SOLVER(MY2,L2,v2z_ijkp1_new,NY2-2)
  VELOCITY_kp1_new(i,2:NY2-1) = v2z_ijkp1_new
END DO

! v(i,j+1) = v(i,j-1) => 
VELOCITY_kp1_new(NX1+1:NX2-1,NY2) = VELOCITY_kp1_new(NX1+1:NX2-1,NY2-2)
  
DO i = NX1,NX2_cone+2,-1
  v1z_ip1jkp12 = VELOCITY_kp12(i+1,1:NY1)
  v1z_ijkp12   = VELOCITY_kp12(i,1:NY1)
  v1z_im1jkp12 = VELOCITY_kp12(i-1,1:NY1)
  v1z_ip1jkp1  = VELOCITY_kp1(i+1,1:NY1)
  v1z_ijkp1    = VELOCITY_kp1(i,1:NY1)
  v1z_im1jkp1  = VELOCITY_kp1(i-1,1:NY1)

  FMSR1z_ip1j  = FMSR(i+1,1:NY1)
  FMSR1z_ij    = FMSR(i,1:NY1)
  FMSR1z_im1j  = FMSR(i-1,1:NY1)

  FMCR1z_ip1j  = FMCR(i+1,1:NY1)
  FMCR1z_ij    = FMCR(i,1:NY1)
  FMCR1z_im1j  = FMCR(i-1,1:NY1)

  CALL MATRIX_UPDATE_Y(rho,k,dt,Lambda,dr,i,dz,NY1,v1z_ip1jkp12,v1z_ijkp12,&
                       v1z_im1jkp12,v1z_ip1jkp1,v1z_ijkp1,v1z_im1jkp1,&
                       FMSR1z_ip1j,FMSR1z_ij,FMSR1z_im1j,&
                       FMCR1z_ip1j,FMCR1z_ij,FMCR1z_im1j,MY1,L1,.FALSE.)

  CALL MATRIX_SOLVER(MY1,L1,v1z_ijkp1_new,NY1-2)
  VELOCITY_kp1_new(i,2:NY1-1) = v1z_ijkp1_new
END DO

! -------------------------- BEGIN OF Y-CONE ------------------------
DO i = NX2_cone+1,NX1_cone+2,-1
  NY2_cone                    = y_cone(i-NX1_cone)
  v1z_c_ip1jkp12(1:NY2_cone)  = VELOCITY_kp12(i+1,1:NY2_cone)
  v1z_c_ijkp12(1:NY2_cone)    = VELOCITY_kp12(i,1:NY2_cone)
  v1z_c_im1jkp12(1:NY2_cone)  = VELOCITY_kp12(i-1,1:NY2_cone)
  v1z_c_ip1jkp1(1:NY2_cone)   = VELOCITY_kp1(i+1,1:NY2_cone)
  v1z_c_ijkp1(1:NY2_cone)     = VELOCITY_kp1(i,1:NY2_cone)
  v1z_c_im1jkp1(1:NY2_cone)   = VELOCITY_kp1(i-1,1:NY2_cone)

  FMSR1z_c_ip1j(1:NY2_cone)   = FMSR(i+1,1:NY2_cone)
  FMSR1z_c_ij(1:NY2_cone)     = FMSR(i,1:NY2_cone)
  FMSR1z_c_im1j(1:NY2_cone)   = FMSR(i-1,1:NY2_cone)

  FMCR1z_c_ip1j(1:NY2_cone)   = FMCR(i+1,1:NY2_cone)
  FMCR1z_c_ij(1:NY2_cone)     = FMCR(i,1:NY2_cone)
  FMCR1z_c_im1j(1:NY2_cone)   = FMCR(i-1,1:NY2_cone)

  CALL MATRIX_UPDATE_Y(rho,k,dt,Lambda,dr,i,dz,NY2_cone,v1z_c_ip1jkp12,v1z_c_ijkp12,&
                       v1z_c_im1jkp12,v1z_c_ip1jkp1,v1z_c_ijkp1,v1z_c_im1jkp1,&
                       FMSR1z_c_ip1j,FMSR1z_c_ij,FMSR1z_c_im1j,&
                       FMCR1z_c_ip1j,FMCR1z_c_ij,FMCR1z_c_im1j,MY1,L1,.FALSE.)

  CALL MATRIX_SOLVER(MY1,L1,v1z_ijkp1_new,NY2_cone-2)
  VELOCITY_kp1_new(i,2:NY2_cone-1) = v1z_ijkp1_new(1:NY2_cone-2)
END DO

DO i = NX1_cone+1,2,-1
  v1z_c_ip1jkp12(1:NY1_cone)  = VELOCITY_kp12(i+1,1:NY1_cone)
  v1z_c_ijkp12(1:NY1_cone)    = VELOCITY_kp12(i,1:NY1_cone)
  v1z_c_im1jkp12(1:NY1_cone)  = VELOCITY_kp12(i-1,1:NY1_cone)
  v1z_c_ip1jkp1(1:NY1_cone)   = VELOCITY_kp1(i+1,1:NY1_cone)
  v1z_c_ijkp1(1:NY1_cone)     = VELOCITY_kp1(i,1:NY1_cone)
  v1z_c_im1jkp1(1:NY1_cone)   = VELOCITY_kp1(i-1,1:NY1_cone)

  FMSR1z_c_ip1j(1:NY1_cone)   = FMSR(i+1,1:NY1_cone)
  FMSR1z_c_ij(1:NY1_cone)     = FMSR(i,1:NY1_cone)
  FMSR1z_c_im1j(1:NY1_cone)   = FMSR(i-1,1:NY1_cone)

  FMCR1z_c_ip1j(1:NY1_cone)   = FMCR(i+1,1:NY1_cone)
  FMCR1z_c_ij(1:NY1_cone)     = FMCR(i,1:NY1_cone)
  FMCR1z_c_im1j(1:NY1_cone)   = FMCR(i-1,1:NY1_cone)

  CALL MATRIX_UPDATE_Y(rho,k,dt,Lambda,dr,i,dz,NY1_cone,v1z_c_ip1jkp12,v1z_c_ijkp12,&
                       v1z_c_im1jkp12,v1z_c_ip1jkp1,v1z_c_ijkp1,v1z_c_im1jkp1,&
                       FMSR1z_c_ip1j,FMSR1z_c_ij,FMSR1z_c_im1j,&
                       FMCR1z_c_ip1j,FMCR1z_c_ij,FMCR1z_c_im1j,MY1,L1,.FALSE.)

  CALL MATRIX_SOLVER(MY1,L1,v1z_ijkp1_new,NY1_cone-2)
  VELOCITY_kp1_new(i,2:NY1_cone-1) = v1z_ijkp1_new(1:NY1_cone-2)
END DO
! --------------------------- END OF Y-CONE -------------------------

!========================= END OF Y-ITERATION =========================

CONVERGENCE = .TRUE.

! Settings for testing of convergence (or rather stability):
RMS      = 0.0D0
vel_norm = 1.0D0

S1: DO j = 2,NY1_cone-1   
      DO i = 2,NX2_cone+2
        vel_norm =  (VELOCITY_kp1_new(i,j) + VELOCITY_kp1(i,j))/2.0D0 + EPS
        RMS      = ((VELOCITY_kp1_new(i,j) - VELOCITY_kp1(i,j))/vel_norm)**2.0D0
        IF (RMS > tol) THEN
          CONVERGENCE = .FALSE.
          EXIT S1
        END IF
      END DO
    END DO S1

IF (CONVERGENCE) THEN
S2: DO j = NY1_cone,NY1-1
      NX1_cone_rms = x_cone(j-NY1_cone+1) + 2
      DO i = NX1_cone_rms,NX2_cone+2
        vel_norm =  (VELOCITY_kp1_new(i,j) + VELOCITY_kp1(i,j))/2.0D0 + EPS
        RMS      = ((VELOCITY_kp1_new(i,j) - VELOCITY_kp1(i,j))/vel_norm)**2.0D0
        IF (RMS > tol) THEN
          CONVERGENCE = .FALSE.
          EXIT S2
        END IF
      END DO
    END DO S2
END IF

IF (CONVERGENCE) THEN
S3: DO j = 2,NY1-1
      DO i = NX2_cone+3,NX2-1
        vel_norm =  (VELOCITY_kp1_new(i,j) + VELOCITY_kp1(i,j))/2.0D0 + EPS
        RMS      = ((VELOCITY_kp1_new(i,j) - VELOCITY_kp1(i,j))/vel_norm)**2.0D0
        IF (RMS > tol) THEN
          CONVERGENCE = .FALSE.
          EXIT S3
        END IF
      END DO
    END DO S3
END IF

IF (CONVERGENCE) THEN
S4: DO j = NY1,NY2
      DO i = NX1+1,NX2-1
        vel_norm =  (VELOCITY_kp1_new(i,j) + VELOCITY_kp1(i,j))/2.0D0 + EPS
        RMS      = ((VELOCITY_kp1_new(i,j) - VELOCITY_kp1(i,j))/vel_norm)**2.0D0
        IF (RMS > tol) THEN
          CONVERGENCE = .FALSE.
          EXIT S4
        END IF
      END DO
    END DO S4
END IF

! Updating ..._kp1
VELOCITY_kp1 = VELOCITY_kp1_new

! CALL WRITE2FILE_k(VELOCITY_k,NX2)

IF (count == count_max) THEN
  CONVERGENCE       = .TRUE.
  FALSE_CONVERGENCE = .TRUE.
  PRINT *, " WARNING: FALSE CONVERGENCE! TIME STEP k = ", k
  PRINT *, " Maximum amount of successive substitutions is = ", count_max
  PRINT *, " ----------------------------------------------------------- "
  PRINT *, " RECOMMENDATION: Kill this application and reduce the time   "
  PRINT *, " step by an order of magnitude: dt -> dt/10                  "
  PRINT *, " ----------------------------------------------------------- "
END IF

END DO CONVERGE 
! ================================================================================= !
! ======================= END OF SUCCESSIVE SUBSTITUTION ========================== !
! ================================================================================= !
! Checking for time independence:
count_rms = 0
RMS       = 0.0D0
vel_norm  = 1.0D0

DO j = 2,NY1_cone-1   
  DO i = 2,NX2_cone+2
    count_rms = count_rms + 1
    vel_norm  =  (VELOCITY_kp1(i,j) + VELOCITY_k(i,j))/2.0D0 + EPS
    RMS       = ((VELOCITY_kp1(i,j) - VELOCITY_k(i,j))/vel_norm)**2.0D0 + RMS
  END DO
END DO

DO j = NY1_cone,NY1-1
  NX1_cone_rms = x_cone(j-NY1_cone+1) + 2
  DO i = NX1_cone_rms,NX2_cone+2
    count_rms = count_rms + 1
    vel_norm  =  (VELOCITY_kp1(i,j) + VELOCITY_k(i,j))/2.0D0 + EPS
    RMS       = ((VELOCITY_kp1(i,j) - VELOCITY_k(i,j))/vel_norm)**2.0D0 + RMS
  END DO
END DO

DO j = 2,NY1-1
  DO i = NX2_cone+3,NX2-1
    count_rms = count_rms + 1
    vel_norm  =  (VELOCITY_kp1(i,j) + VELOCITY_k(i,j))/2.0D0 + EPS
    RMS       = ((VELOCITY_kp1(i,j) - VELOCITY_k(i,j))/vel_norm)**2.0D0 + RMS
  END DO
END DO

DO j = NY1,NY2
  DO i = NX1+1,NX2-1
    count_rms = count_rms + 1
    vel_norm  =  (VELOCITY_kp1(i,j) + VELOCITY_k(i,j))/2.0D0 + EPS
    RMS       = ((VELOCITY_kp1(i,j) - VELOCITY_k(i,j))/vel_norm)**2.0D0 + RMS
  END DO
END DO

RMS = DSQRT(RMS/count_rms)

IF (ABS(MOD((k+1),k_OUTPUT_rms)).LT.small_zero) THEN
  PRINT *, "   RMS   =", RMS
  CALL WRITE2FILE_rms(k+1,RMS)
END IF

IF (RMS.LT.tol_RMS_active) THEN
  TIME_INDEPENDENCE = .TRUE.
ELSE 
  TIME_INDEPENDENCE = .FALSE.
END IF

IF (k == MAX_NUMBER_OF_ITERATIONS) THEN
  TIME_INDEPENDENCE = .TRUE.
  WARNING_SIGN      = .TRUE. 
END IF

VELOCITY_k = VELOCITY_kp1_new

k = k + 1

END DO ZERO_TIME_LOOP
! ================================================================================= !
! ============================ End of the time loop =============================== !
! ================================================================================= !
END DO CONTINUATION
! ================================================================================= !
! ============================= End of CONTINUATION =============================== !
! ================================================================================= !
IF (WARNING_SIGN) THEN 
  PRINT *, " --------------------------------------------------------- "
  PRINT *, " WARNING: k = MAX_NUMBER_OF_ITERATIONS; See log.dat        "
  PRINT *, " RECOMMENDATIONS:                                          "
  PRINT *, " I)  Rerun this application with reduced time step.        "
  PRINT *, "     Try order of magnitude less: dt -> dt/10.             "
  PRINT *, " II) Increase ZERO_TIME in the file param.f90.             "
  PRINT *, " --------------------------------------------------------- "
END IF

IF (FALSE_CONVERGENCE) THEN
  PRINT *, " --------------------------------------------------------- "
  PRINT *, " WARNING: FALSE CONVERGENCE WAS ACHIEVED. Reduce the time  "
  PRINT *, " step by order of magnitude: dt -> dt/10 and then rerun    "
  PRINT *, " the application.                                          "
  PRINT *, " --------------------------------------------------------- "
END IF

PRINT *, " PSEUDO-TRANSIENT CALCULATION FINISHED! "
PRINT *, " -------------------------------------- "
! --------------------------------------------------------------------------------- !
IF (.NOT.CALCULATE_TIME_DEPENDENT_PROBL) THEN
  PRINT *, " Number of grid points (not including Dirichlet " 
  PRINT *, " boundary points) = ", count_rms
  PRINT *, "                                                "
  PRINT *, " NO TIME DEPENDENT CALCULATION IS DONE SINCE    "
  PRINT *, " CALCULATE_TIME_DEPENDENT_PROBL = .FALSE.       "
  PRINT *, "                                                "
  PRINT *, " WRITING INFORMATION TO FILE, STAND BY...       "
  CALL WRITE2FILE_k(VELOCITY_k,NX2)
  CALL WRITE2FILE_kp1(VELOCITY_kp1_new,NX1,NY1,NX2,NY2,dr,dz)
  CALL WRITE2FILE_torque_ZERO(VELOCITY_k,NX1,NX2,NY1,NY2,0,dt,Lambda,dr,dz,H3,omega)
  PRINT *, " ...DONE! "

  ! CLEARING SOME MAJOR VARIABLES FROM THE RANDOM ACCESS MEMORY:
  DEALLOCATE(x_cone,y_cone,stat=problem)
  IF (problem/=0) THEN 
    PRINT *," MAIN_ROUTINE says: The program could not deallocate space! "
    PRINT *," Error code 6 in main and execution terminated!             "
  END IF

  DEALLOCATE(MX1,MX2,MY1,MY2,VELOCITY_k,VELOCITY_kp12,VELOCITY_kp12_new,&
             VELOCITY_kp1,VELOCITY_kp1_new,stat=problem)
  IF (problem/=0) THEN 
    PRINT *," MAIN_ROUTINE says: The program could not deallocate space! "
    PRINT *," Error code 7 in main!                                      "
  END IF

  DEALLOCATE(K1,K2,L1,L2,stat=problem)
  IF (problem/=0) THEN 
    PRINT *," MAIN_ROUTINE says: The program could not deallocate space! "
    PRINT *," Error code 8 in main!                                      "
  END IF

  DEALLOCATE(v1r_ijp1k,v1r_ijk,v1r_ijm1k,v1r_ijp1kp12,v1r_ijkp12,&
             v1r_ijm1kp12,v1r_c_ijp1k,v1r_c_ijk,v1r_c_ijm1k,&
             v1r_c_ijp1kp12,v1r_c_ijkp12,v1r_c_ijm1kp12,&
             v2r_ijp1k,v2r_ijk,v2r_ijm1k,v2r_ijp1kp12,&
             v2r_ijkp12,v2r_ijm1kp12,stat=problem)
  IF (problem/=0) THEN 
    PRINT *," MAIN_ROUTINE says: The program could not deallocate space! "
    PRINT *," Error code 9 in main!                                      "
  END IF

  DEALLOCATE(v1r_ijkp12_new,v2r_ijkp12_new,stat=problem)
  IF (problem/=0) THEN 
    PRINT *," MAIN_ROUTINE says: The program could not deallocate space! "
    PRINT *," Error code 10 in main!                                     "
  END IF

  DEALLOCATE(v1z_ip1jkp12,v1z_ijkp12,v1z_im1jkp12,v1z_ip1jkp1,v1z_ijkp1,&
             v1z_im1jkp1,v1z_c_ip1jkp12,v1z_c_ijkp12,v1z_c_im1jkp12,&
             v1z_c_ip1jkp1,v1z_c_ijkp1,v1z_c_im1jkp1,&
             v2z_ip1jkp12,v2z_ijkp12,v2z_im1jkp12,&
             v2z_ip1jkp1,v2z_ijkp1,v2z_im1jkp1,stat=problem)
  IF (problem/=0) THEN 
    PRINT *," MAIN_ROUTINE says: The program could not deallocate space! "
    PRINT *," Error code 11 in main!                                     "
  END IF

  DEALLOCATE(v1z_ijkp1_new,v2z_ijkp1_new,stat=problem)
  IF (problem/=0) THEN 
    PRINT *," MAIN_ROUTINE says: The program could not deallocate space! "
    PRINT *," Error code 12 in main!                                     "
  END IF

  PRINT *, " EXECUTION FINISHED! "
  STOP

END IF
! ================================================================================= !
! ================================================================================= !
! ================================================================================= !
dt  = dt_Plastic
tol = tol_Plastic

! "Lambda = -1" means that time dependent calculations (basically thixotropic)
!  have begun with time dependent shear viscosity ETA = ETA(SR,t,...).
Lambda = - 1.0D0

! Note that at the moment, then VELOCITY_k = VELOCITY_kp1_new. 
PRINT *, "                                                "
PRINT *, " WRITING 't=0' INFORMATION TO FILE, STAND BY... "
CALL WRITE2FILE_k(VELOCITY_k,NX2)
CALL WRITE2FILE_kp1(VELOCITY_kp1_new,NX1,NY1,NX2,NY2,dr,dz)
CALL WRITE2FILE_torque_ZERO(VELOCITY_k,NX1,NX2,NY1,NY2,0,dt,Lambda,dr,dz,H3,omega)
PRINT *, " ...DONE! "

PRINT *, " NOW FINALLY, BEGINNING WITH THE TIME DEPENDENT PROBLEM... "

! ################################################################################# !
! In Section 7.11.2 is a detailed description of the algorithm, which is used in    !
! the following.                                                                    !
! --------------------------------------------------------------------------------- !
! In this software, 'k' corresponds to t=(k+1/2)*dt and t=(k+1)*dt, and hence 'k = 0' 
! means t = (1/2)*dt and t = 1*dt.
! This is due to the calling routine of the angular velocity:
! 'CALL ANGULAR_VELOCITY(DBLE(k)+1.0D0/2.0D0,dt,omega)' and 
! 'CALL ANGULAR_VELOCITY(DBLE(k)+1.0D0,dt,omega)'.

! The fading memory modules are zero at k=0 (see Section 9.3.1):
FMSR = 0.0D0 
FMCR = 0.0D0

N_dt = (NUMBER_OF_TIME_ITERATIONS - 1)  
! ================================================================================= !
! ======================= Begin of the time loop (TIME) =========================== !
! ================================================================================= !
TIME_LOOP: DO k=0,N_dt

CONVERGENCE = .FALSE.

IF (ABS(MOD((k+1),k_OUTPUT_rms)).LT.small_zero) THEN
  PRINT *,"______________________________________________________________________"
  PRINT *,"______________________________________________________________________"
  PRINT *,"                                                                      "
  PRINT *,"  Calculating now for the time step: k+1/2 and k+1 = ",k+1," -> ->    "
  PRINT *,"  time = (k+1/2)dt = ",(DBLE(k)+1.0D0/2.0D0)*dt," SEC"
  PRINT *,"  time = (k+1)dt   = ",(DBLE(k)+1.0D0)*dt," SEC"
  PRINT *,"----------------------------------------------------------------------"
END IF
! --------------------------------------------------------------------------------- !
! A guess for the time step k+1/2
VELOCITY_kp12 = VELOCITY_k

! Updating boundary condition at k+1/2:
! CALL ANGULAR_VELOCITY(k+1/2,dt,omega)
CALL ANGULAR_VELOCITY(DBLE(k)+1.0D0/2.0D0,dt,omega)

! Dirichlet boundary condition:
DO i = 1,NX2
  VELOCITY_kp12(i,1) = omega*DBLE(i-1)*dr
END DO
VELOCITY_kp12(1:NX1,NY1)       = 0.0D0
VELOCITY_kp12(1,1:NY1)         = 0.0D0
VELOCITY_kp12(NX1,NY1:NY2)     = 0.0D0
VELOCITY_kp12(NX2,:)           = R_o*omega

! Neumann boundary condition (this is actually redundant since the new 
! information of "omega" cannot reach "VELOCITY_kp12(NX1+1:NX2-1,NY2-1)" 
! or "VELOCITY_kp12(NX1+1:NX2-1,NY2-2)":
VELOCITY_kp12(NX1+1:NX2-1,NY2) = (4.0D0 * &
     VELOCITY_kp12(NX1+1:NX2-1,NY2-1) - &
     VELOCITY_kp12(NX1+1:NX2-1,NY2-2))/3.0D0

! Also updating boundary condition for "..._new":
VELOCITY_kp12_new = VELOCITY_kp12
! --------------------------------------------------------------------------------- !
! A guess for the time step k+1
VELOCITY_kp1 = VELOCITY_k

! Updating boundary condition at k+1:
! CALL ANGULAR_VELOCITY(k+1,dt,omega)
CALL ANGULAR_VELOCITY(DBLE(k)+1.0D0,dt,omega)

! Dirichlet boundary condition:
DO i = 1,NX2
  VELOCITY_kp1(i,1) = omega*DBLE(i-1)*dr
END DO
VELOCITY_kp1(1:NX1,NY1)       = 0.0D0
VELOCITY_kp1(1,1:NY1)         = 0.0D0
VELOCITY_kp1(NX1,NY1:NY2)     = 0.0D0
VELOCITY_kp1(NX2,:)           = R_o*omega

! Neumann boundary condition (this is actually redundant since the new 
! information of "omega" cannot reach "VELOCITY_kp12(NX1+1:NX2-1,NY2-1)" 
! or "VELOCITY_kp12(NX1+1:NX2-1,NY2-2)":
VELOCITY_kp1(NX1+1:NX2-1,NY2) = (4.0D0 * &
     VELOCITY_kp1(NX1+1:NX2-1,NY2-1) - &
     VELOCITY_kp1(NX1+1:NX2-1,NY2-2))/3.0D0

! Also updating boundary condition for "..._new":
VELOCITY_kp1_new = VELOCITY_kp1
! --------------------------------------------------------------------------------- !
count = 0
! ================================================================================= !
! =================== BEGIN OF SUCCESSIVE SUBSTITUTION (TIME) ===================== !
! ================================================================================= !
! The iteration loop here is because of the non-linearity of the governing
! Equations 7.22 and 7.23. To come around this problem, the successive substitution
! approach is used (see Section 7.8).
TIME_CONVERGE: DO WHILE (.NOT.CONVERGENCE)

! If convergence is a problem, then this might help:
! VELOCITY_kp12 = (VELOCITY_k + VELOCITY_kp1_new)/2
count = count + 1
12 FORMAT(4X,"Successive substitution (time) number = ",1(I3,1X))

IF (ABS(MOD((k+1),k_OUTPUT_rms)).LT.small_zero) THEN
  PRINT 12, count
END IF

!======================== BEGIN OF X-ITERATION (TIME) =========================
! Iteration is made along r-direction (i.e. along the i-direction as 
! in A(i,j)). It starts at the bottom of the viscometer i=(2:NX2-1) at 
! j = 2 and then move upward with increasing j (see Figures 8.1 and 8.2).

DO j = 2,NY1_cone-1
  v1r_ijp1k    = VELOCITY_k(:,j+1)
  v1r_ijk      = VELOCITY_k(:,j)
  v1r_ijm1k    = VELOCITY_k(:,j-1)
  v1r_ijp1kp12 = VELOCITY_kp12(:,j+1)
  v1r_ijkp12   = VELOCITY_kp12(:,j)
  v1r_ijm1kp12 = VELOCITY_kp12(:,j-1)

  FMSR1r_ijp1  = FMSR(:,j+1)
  FMSR1r_ij    = FMSR(:,j)
  FMSR1r_ijm1  = FMSR(:,j-1)

  FMCR1r_ijp1  = FMCR(:,j+1)
  FMCR1r_ij    = FMCR(:,j)
  FMCR1r_ijm1  = FMCR(:,j-1)

  CALL MATRIX_UPDATE_X(rho,k,dt,Lambda,dr,0.0D0,dz,NX2,v1r_ijp1k,v1r_ijk,&
                       v1r_ijm1k,v1r_ijp1kp12,v1r_ijkp12,v1r_ijm1kp12,&
                       FMSR1r_ijp1,FMSR1r_ij,FMSR1r_ijm1,&
                       FMCR1r_ijp1,FMCR1r_ij,FMCR1r_ijm1,MX1,K1)

  CALL MATRIX_SOLVER(MX1,K1,v1r_ijkp12_new,NX2-2)
  VELOCITY_kp12_new(2:NX2-1,j) = v1r_ijkp12_new
END DO

! -------------------------- BEGIN OF X-CONE ------------------------
DO j = NY1_cone,NY1-1
  NX_cone                       = x_cone(j + 1 - NY1_cone)
  r_i_cone                      = NX_cone*dr
  NX2mNX_cone                   = NX2 - NX_cone
  v1r_c_ijp1k(1:NX2mNX_cone)    = VELOCITY_k(NX_cone+1:NX2,j+1)
  v1r_c_ijk(1:NX2mNX_cone)      = VELOCITY_k(NX_cone+1:NX2,j)
  v1r_c_ijm1k(1:NX2mNX_cone)    = VELOCITY_k(NX_cone+1:NX2,j-1)
  v1r_c_ijp1kp12(1:NX2mNX_cone) = VELOCITY_kp12(NX_cone+1:NX2,j+1)
  v1r_c_ijkp12(1:NX2mNX_cone)   = VELOCITY_kp12(NX_cone+1:NX2,j)
  v1r_c_ijm1kp12(1:NX2mNX_cone) = VELOCITY_kp12(NX_cone+1:NX2,j-1)

  FMSR1r_c_ijp1(1:NX2mNX_cone)  = FMSR(NX_cone+1:NX2,j+1)
  FMSR1r_c_ij(1:NX2mNX_cone)    = FMSR(NX_cone+1:NX2,j)
  FMSR1r_c_ijm1(1:NX2mNX_cone)  = FMSR(NX_cone+1:NX2,j-1)

  FMCR1r_c_ijp1(1:NX2mNX_cone)  = FMCR(NX_cone+1:NX2,j+1)
  FMCR1r_c_ij(1:NX2mNX_cone)    = FMCR(NX_cone+1:NX2,j)
  FMCR1r_c_ijm1(1:NX2mNX_cone)  = FMCR(NX_cone+1:NX2,j-1)

  CALL MATRIX_UPDATE_X(rho,k,dt,Lambda,dr,r_i_cone,dz,NX2mNX_cone,v1r_c_ijp1k,v1r_c_ijk,&
                       v1r_c_ijm1k,v1r_c_ijp1kp12,v1r_c_ijkp12,v1r_c_ijm1kp12,&
                       FMSR1r_c_ijp1,FMSR1r_c_ij,FMSR1r_c_ijm1,&
                       FMCR1r_c_ijp1,FMCR1r_c_ij,FMCR1r_c_ijm1,MX1,K1)

  CALL MATRIX_SOLVER(MX1,K1,v1r_ijkp12_new,NX2mNX_cone-2)
  ! CALL WRITE2FILE_debug(MX1,K1,v1r_ijkp12_new,NX2mNX_cone-2)
  ! STOP
  VELOCITY_kp12_new(NX_cone+2:NX2-1,j) = v1r_ijkp12_new(1:NX2mNX_cone-2)
END DO
! --------------------------- END OF X-CONE -------------------------

DO j = NY1,NY2-1
  v2r_ijp1k    = VELOCITY_k(NX1:NX2,j+1)
  v2r_ijk      = VELOCITY_k(NX1:NX2,j)
  v2r_ijm1k    = VELOCITY_k(NX1:NX2,j-1)
  v2r_ijp1kp12 = VELOCITY_kp12(NX1:NX2,j+1)
  v2r_ijkp12   = VELOCITY_kp12(NX1:NX2,j)
  v2r_ijm1kp12 = VELOCITY_kp12(NX1:NX2,j-1)

  FMSR2r_ijp1  = FMSR(NX1:NX2,j+1)
  FMSR2r_ij    = FMSR(NX1:NX2,j)
  FMSR2r_ijm1  = FMSR(NX1:NX2,j-1)

  FMCR2r_ijp1  = FMCR(NX1:NX2,j+1)
  FMCR2r_ij    = FMCR(NX1:NX2,j)
  FMCR2r_ijm1  = FMCR(NX1:NX2,j-1)

  CALL MATRIX_UPDATE_X(rho,k,dt,Lambda,dr,R_i,dz,NX,v2r_ijp1k,v2r_ijk,&
                       v2r_ijm1k,v2r_ijp1kp12,v2r_ijkp12,v2r_ijm1kp12,&
                       FMSR2r_ijp1,FMSR2r_ij,FMSR2r_ijm1,&
                       FMCR2r_ijp1,FMCR2r_ij,FMCR2r_ijm1,MX2,K2)

  CALL MATRIX_SOLVER(MX2,K2,v2r_ijkp12_new,NX-2)
  VELOCITY_kp12_new(NX1+1:NX2-1,j) = v2r_ijkp12_new
END DO

j = NY2
v2r_ijp1k    = VELOCITY_k(NX1:NX2,j-1)    ! => v(i,j+1) = v(i,j-1)
v2r_ijk      = VELOCITY_k(NX1:NX2,j)
v2r_ijm1k    = VELOCITY_k(NX1:NX2,j-1)
v2r_ijp1kp12 = VELOCITY_kp12(NX1:NX2,j-1) ! => v(i,j+1) = v(i,j-1)
v2r_ijkp12   = VELOCITY_kp12(NX1:NX2,j)
v2r_ijm1kp12 = VELOCITY_kp12(NX1:NX2,j-1)

FMSR2r_ijp1  = FMSR(NX1:NX2,j-1)          ! => v(i,j+1) = v(i,j-1) 
FMSR2r_ij    = FMSR(NX1:NX2,j)
FMSR2r_ijm1  = FMSR(NX1:NX2,j-1)

FMCR2r_ijp1  = FMCR(NX1:NX2,j-1)          ! => v(i,j+1) = v(i,j-1) 
FMCR2r_ij    = FMCR(NX1:NX2,j)
FMCR2r_ijm1  = FMCR(NX1:NX2,j-1)

CALL MATRIX_UPDATE_X(rho,k,dt,Lambda,dr,R_i,dz,NX,v2r_ijp1k,v2r_ijk,&
                     v2r_ijm1k,v2r_ijp1kp12,v2r_ijkp12,v2r_ijm1kp12,&
                     FMSR2r_ijp1,FMSR2r_ij,FMSR2r_ijm1,&
                     FMCR2r_ijp1,FMCR2r_ij,FMCR2r_ijm1,MX2,K2)

CALL MATRIX_SOLVER(MX2,K2,v2r_ijkp12_new,NX-2)
VELOCITY_kp12_new(NX1+1:NX2-1,j) = v2r_ijkp12_new

! --------- PAUSE FOR DEBUGGING --------- 
! CALL WRITE2FILE_k(VELOCITY_kp12_new,NX2)
! WRITE(*,"(A)",ADVANCE="NO") " PRESS 'ENTER' TO CONTINUE "
! PRINT *, " "
! READ (*,"(A)") IGNORED_INPUT
! PRINT *, " "

!========================= END OF X-ITERATION (TIME) =========================

! Updating ..._kp12:
VELOCITY_kp12 = VELOCITY_kp12_new

!======================== BEGIN OF Y-ITERATION (TIME) ========================
! Iteration is made along z-direction (i.e. along the j-direction as 
! in A(i,j)). It starts at the right side of the viscometer j=(2:NY2-1) at 
! i = NX2-1 and then moves to the left with decreasing i
! (see Figures 8.1 and 8.2).

DO i = NX2-1,NX1+1,-1
  v2z_ip1jkp12 = VELOCITY_kp12(i+1,:)
  v2z_ijkp12   = VELOCITY_kp12(i,:)
  v2z_im1jkp12 = VELOCITY_kp12(i-1,:)
  v2z_ip1jkp1  = VELOCITY_kp1(i+1,:)
  v2z_ijkp1    = VELOCITY_kp1(i,:)
  v2z_im1jkp1  = VELOCITY_kp1(i-1,:)

  FMSR2z_ip1j  = FMSR(i+1,:)
  FMSR2z_ij    = FMSR(i,:)
  FMSR2z_im1j  = FMSR(i-1,:)
  
  FMCR2z_ip1j  = FMCR(i+1,:)
  FMCR2z_ij    = FMCR(i,:)
  FMCR2z_im1j  = FMCR(i-1,:)

  CALL MATRIX_UPDATE_Y(rho,k,dt,Lambda,dr,i,dz,NY2,v2z_ip1jkp12,v2z_ijkp12,&
                       v2z_im1jkp12,v2z_ip1jkp1,v2z_ijkp1,v2z_im1jkp1,&
                       FMSR2z_ip1j,FMSR2z_ij,FMSR2z_im1j,&
                       FMCR2z_ip1j,FMCR2z_ij,FMCR2z_im1j,MY2,L2,.TRUE.)

  CALL MATRIX_SOLVER(MY2,L2,v2z_ijkp1_new,NY2-2)
  VELOCITY_kp1_new(i,2:NY2-1) = v2z_ijkp1_new
END DO

! v(i,j+1) = v(i,j-1) => 
VELOCITY_kp1_new(NX1+1:NX2-1,NY2) = VELOCITY_kp1_new(NX1+1:NX2-1,NY2-2)
  
DO i = NX1,NX2_cone+2,-1
  v1z_ip1jkp12 = VELOCITY_kp12(i+1,1:NY1)
  v1z_ijkp12   = VELOCITY_kp12(i,1:NY1)
  v1z_im1jkp12 = VELOCITY_kp12(i-1,1:NY1)
  v1z_ip1jkp1  = VELOCITY_kp1(i+1,1:NY1)
  v1z_ijkp1    = VELOCITY_kp1(i,1:NY1)
  v1z_im1jkp1  = VELOCITY_kp1(i-1,1:NY1)

  FMSR1z_ip1j  = FMSR(i+1,1:NY1)
  FMSR1z_ij    = FMSR(i,1:NY1)
  FMSR1z_im1j  = FMSR(i-1,1:NY1)

  FMCR1z_ip1j  = FMCR(i+1,1:NY1)
  FMCR1z_ij    = FMCR(i,1:NY1)
  FMCR1z_im1j  = FMCR(i-1,1:NY1)

  CALL MATRIX_UPDATE_Y(rho,k,dt,Lambda,dr,i,dz,NY1,v1z_ip1jkp12,v1z_ijkp12,&
                       v1z_im1jkp12,v1z_ip1jkp1,v1z_ijkp1,v1z_im1jkp1,&
                       FMSR1z_ip1j,FMSR1z_ij,FMSR1z_im1j,&
                       FMCR1z_ip1j,FMCR1z_ij,FMCR1z_im1j,MY1,L1,.FALSE.)

  CALL MATRIX_SOLVER(MY1,L1,v1z_ijkp1_new,NY1-2)
  VELOCITY_kp1_new(i,2:NY1-1) = v1z_ijkp1_new
END DO

! -------------------------- BEGIN OF Y-CONE ------------------------
DO i = NX2_cone+1,NX1_cone+2,-1
  NY2_cone                    = y_cone(i-NX1_cone)
  v1z_c_ip1jkp12(1:NY2_cone)  = VELOCITY_kp12(i+1,1:NY2_cone)
  v1z_c_ijkp12(1:NY2_cone)    = VELOCITY_kp12(i,1:NY2_cone)
  v1z_c_im1jkp12(1:NY2_cone)  = VELOCITY_kp12(i-1,1:NY2_cone)
  v1z_c_ip1jkp1(1:NY2_cone)   = VELOCITY_kp1(i+1,1:NY2_cone)
  v1z_c_ijkp1(1:NY2_cone)     = VELOCITY_kp1(i,1:NY2_cone)
  v1z_c_im1jkp1(1:NY2_cone)   = VELOCITY_kp1(i-1,1:NY2_cone)

  FMSR1z_c_ip1j(1:NY2_cone)   = FMSR(i+1,1:NY2_cone)
  FMSR1z_c_ij(1:NY2_cone)     = FMSR(i,1:NY2_cone)
  FMSR1z_c_im1j(1:NY2_cone)   = FMSR(i-1,1:NY2_cone)

  FMCR1z_c_ip1j(1:NY2_cone)   = FMCR(i+1,1:NY2_cone)
  FMCR1z_c_ij(1:NY2_cone)     = FMCR(i,1:NY2_cone)
  FMCR1z_c_im1j(1:NY2_cone)   = FMCR(i-1,1:NY2_cone)

  CALL MATRIX_UPDATE_Y(rho,k,dt,Lambda,dr,i,dz,NY2_cone,v1z_c_ip1jkp12,v1z_c_ijkp12,&
                       v1z_c_im1jkp12,v1z_c_ip1jkp1,v1z_c_ijkp1,v1z_c_im1jkp1,&
                       FMSR1z_c_ip1j,FMSR1z_c_ij,FMSR1z_c_im1j,&
                       FMCR1z_c_ip1j,FMCR1z_c_ij,FMCR1z_c_im1j,MY1,L1,.FALSE.)

  CALL MATRIX_SOLVER(MY1,L1,v1z_ijkp1_new,NY2_cone-2)
  VELOCITY_kp1_new(i,2:NY2_cone-1) = v1z_ijkp1_new(1:NY2_cone-2)
END DO

DO i = NX1_cone+1,2,-1
  v1z_c_ip1jkp12(1:NY1_cone)  = VELOCITY_kp12(i+1,1:NY1_cone)
  v1z_c_ijkp12(1:NY1_cone)    = VELOCITY_kp12(i,1:NY1_cone)
  v1z_c_im1jkp12(1:NY1_cone)  = VELOCITY_kp12(i-1,1:NY1_cone)
  v1z_c_ip1jkp1(1:NY1_cone)   = VELOCITY_kp1(i+1,1:NY1_cone)
  v1z_c_ijkp1(1:NY1_cone)     = VELOCITY_kp1(i,1:NY1_cone)
  v1z_c_im1jkp1(1:NY1_cone)   = VELOCITY_kp1(i-1,1:NY1_cone)

  FMSR1z_c_ip1j(1:NY1_cone)   = FMSR(i+1,1:NY1_cone)
  FMSR1z_c_ij(1:NY1_cone)     = FMSR(i,1:NY1_cone)
  FMSR1z_c_im1j(1:NY1_cone)   = FMSR(i-1,1:NY1_cone)

  FMCR1z_c_ip1j(1:NY1_cone)   = FMCR(i+1,1:NY1_cone)
  FMCR1z_c_ij(1:NY1_cone)     = FMCR(i,1:NY1_cone)
  FMCR1z_c_im1j(1:NY1_cone)   = FMCR(i-1,1:NY1_cone)

  CALL MATRIX_UPDATE_Y(rho,k,dt,Lambda,dr,i,dz,NY1_cone,v1z_c_ip1jkp12,v1z_c_ijkp12,&
                       v1z_c_im1jkp12,v1z_c_ip1jkp1,v1z_c_ijkp1,v1z_c_im1jkp1,&
                       FMSR1z_c_ip1j,FMSR1z_c_ij,FMSR1z_c_im1j,&
                       FMCR1z_c_ip1j,FMCR1z_c_ij,FMCR1z_c_im1j,MY1,L1,.FALSE.)

  CALL MATRIX_SOLVER(MY1,L1,v1z_ijkp1_new,NY1_cone-2)
  VELOCITY_kp1_new(i,2:NY1_cone-1) = v1z_ijkp1_new(1:NY1_cone-2)
END DO
! --------------------------- END OF Y-CONE -------------------------

!========================= END OF Y-ITERATION (TIME) =========================

CONVERGENCE = .TRUE.

! Settings for testing of convergence (or rather stability):
RMS      = 0.0D0
vel_norm = 1.0D0

S5: DO j = 2,NY1_cone-1   
      DO i = 2,NX2_cone+2
        vel_norm =  (VELOCITY_kp1_new(i,j) + VELOCITY_kp1(i,j))/2.0D0 + EPS
        RMS      = ((VELOCITY_kp1_new(i,j) - VELOCITY_kp1(i,j))/vel_norm)**2.0D0
        IF (RMS > tol) THEN
          CONVERGENCE = .FALSE.
          EXIT S5
        END IF
      END DO
    END DO S5

IF (CONVERGENCE) THEN
S6: DO j = NY1_cone,NY1-1
      NX1_cone_rms = x_cone(j-NY1_cone+1) + 2
      DO i = NX1_cone_rms,NX2_cone+2
        vel_norm =  (VELOCITY_kp1_new(i,j) + VELOCITY_kp1(i,j))/2.0D0 + EPS
        RMS      = ((VELOCITY_kp1_new(i,j) - VELOCITY_kp1(i,j))/vel_norm)**2.0D0
        IF (RMS > tol) THEN
          CONVERGENCE = .FALSE.
          EXIT S6
        END IF
      END DO
    END DO S6
END IF

IF (CONVERGENCE) THEN
S7: DO j = 2,NY1-1
      DO i = NX2_cone+3,NX2-1
        vel_norm =  (VELOCITY_kp1_new(i,j) + VELOCITY_kp1(i,j))/2.0D0 + EPS
        RMS      = ((VELOCITY_kp1_new(i,j) - VELOCITY_kp1(i,j))/vel_norm)**2.0D0
        IF (RMS > tol) THEN
          CONVERGENCE = .FALSE.
          EXIT S7
        END IF
      END DO
    END DO S7
END IF

IF (CONVERGENCE) THEN
S8: DO j = NY1,NY2
      DO i = NX1+1,NX2-1
        vel_norm =  (VELOCITY_kp1_new(i,j) + VELOCITY_kp1(i,j))/2.0D0 + EPS
        RMS      = ((VELOCITY_kp1_new(i,j) - VELOCITY_kp1(i,j))/vel_norm)**2.0D0
        IF (RMS > tol) THEN
          CONVERGENCE = .FALSE.
          EXIT S8
        END IF
      END DO
    END DO S8
END IF

! Updating ..._kp1:
VELOCITY_kp1 = VELOCITY_kp1_new

IF (count == count_max) THEN
  CONVERGENCE       = .TRUE.
  FALSE_CONVERGENCE = .TRUE.
  PRINT *, " WARNING: FALSE CONVERGENCE! TIME STEP k= ", k
  PRINT *, " Maximum amount of successive substitutions is = ", count_max
  PRINT *, " ----------------------------------------------------------- "
  PRINT *, " RECOMMENDATION: Kill this application and reduce the time   "
  PRINT *, " step by an order of magnitude: dt -> dt/10                  "
  PRINT *, " ----------------------------------------------------------- "
END IF

END DO TIME_CONVERGE 
! ================================================================================= !
! ==================== END OF SUCCESSIVE SUBSTITUTION (TIME) ====================== !
! ================================================================================= !
VELOCITY_k = VELOCITY_kp1_new
! ---------------------------------------------------------
CALL ROS_PROFILE(VELOCITY_k,NX1,NX2,NY1,NY2,DUMMY_1,dr,dz,SR,DUMMY_2)
time = (DBLE(k)+1.0D0)*dt
! ---------------------------------------------------------
DO i = 1,NX2
  DO j = 1,NY1
    SR_tmp = SR(i,j)
    CALL FMSR_FMCR(time,SR_tmp,U_o,memory_alpha,memory_beta,H_tmp)
    H(i,j) = H_tmp
  END DO
END DO
DO i = NX1,NX2
  DO j = NY1+1,NY2
    SR_tmp = SR(i,j)
    CALL FMSR_FMCR(time,SR_tmp,U_o,memory_alpha,memory_beta,H_tmp)
    H(i,j) = H_tmp
  END DO
END DO
! ---------------------------------------------------------
ALPHA_I = DEXP(time/memory_alpha)
FMSR    = ALPHA_I*SR*dt_Plastic + FMSR
! The above implicitly consists of two steps:
! 1) FMSR_kp1 = ALPHA_I*SR*dt_Plastic + FMSR_k (hence, time=(k+1)dt)
! 2) FMSR_k   = FMSR_kp1 <- Updating for the next time step, 
!    just as the velocity was updated in the above =>
!    VELOCITY_k = VELOCITY_kp1_new (k <- k+1).
! ---------------------------------------------------------
BETA_I  = DEXP(time/memory_beta)
FMCR    = BETA_I*H*dt_Plastic + FMCR
! The above implicitly consists of two steps:
! 1) FMCR_kp1 = BETA_I*H*dt_Plastic + FMCR_k (hence, time=(k+1)dt)
! 2) FMCR_k   = FMCR_kp1 <- Updating for the next time step,
!    just as the velocity was updated in the above =>
!    VELOCITY_k = VELOCITY_kp1_new (k <- k+1).
! --------------------------------------------------------------------------------- !
! NOTE THAT THE VARIABLE "VELOCITY_k" NOW CONTAINS THE VELOCITY AT THE TIME 
! STEP k+1, C.F. "VELOCITY_k = VELOCITY_kp1_new", JUST ABOVE! 
! ---------------------------------------------------------
! Small output every 0.1 second (vel_corner.dat and vel_upper.dat):
! See Section 8.4.1 about vel_corner.dat and vel_upper.dat.
dt_OUTPUT_torque = 0.1D0 
k_OUTPUT_torque  = IDNINT(dt_OUTPUT_torque/dt)
IF (ABS(MOD((k+1),k_OUTPUT_torque)).LT.small_zero) THEN
  PRINT *, " =================================================================== "
  PRINT *, " WRITING THE COMPACT DATA AT TIME ", DBLE(k+1)*dt
  PRINT *, " SECONDS TO FILE, STAND BY...     "
  CALL WRITE2FILE_torque(VELOCITY_k,FMSR,FMCR,NX1,NX2,NY1,NY2,k+1,dt,Lambda,dr,dz,H3,omega)
  PRINT *, " ...DONE! "
  PRINT *, " =================================================================== "
END IF

! Large output every 1 second (the whole domain of calculation -> vel_k.dat):
dt_OUTPUT = 1.0D0 
k_OUTPUT  = IDNINT(dt_OUTPUT/dt) ! = 100000
IF (ABS(MOD((k+1),k_OUTPUT)).LT.small_zero) THEN
  PRINT *, " =================================================================== "
  PRINT *, " WRITING LARGE DATA AT TIME   ", DBLE(k+1)*dt
  PRINT *, " SECONDS TO FILE, STAND BY... "
  CALL WRITE2FILE_time(VELOCITY_k,NX2,k+1,dt)
  PRINT *, " ...DONE! "
  PRINT *, " =================================================================== "
END IF

END DO TIME_LOOP
! ================================================================================= !
! ======================== End of the time loop (TIME) ============================ !
! ================================================================================= !
IF (WARNING_SIGN) THEN
  PRINT *, " --------------------------------------------------------- "
  PRINT *, " WARNING: For the time independent case, then:             "
  PRINT *, " k = MAX_NUMBER_OF_ITERATIONS; See log.dat                 "
  PRINT *, " --------------------------------------------------------- "
END IF

IF (FALSE_CONVERGENCE) THEN
  PRINT *, " --------------------------------------------------------- "
  PRINT *, " WARNING: For the time independent case, then:             "
  PRINT *, " FALSE CONVERGENCE WAS ACHIEVED. Reduce the time           "
  PRINT *, " step by order of magnitude: dt -> dt/10 and then rerun    "
  PRINT *, " the application.                                          "
  PRINT *, " --------------------------------------------------------- "
END IF

! --------------------------------------------------------------------------------- !
! CLEARING SOME MAJOR VARIABLES FROM THE RANDOM ACCESS MEMORY:
DEALLOCATE(x_cone,y_cone,stat=problem)
IF (problem/=0) THEN 
  PRINT *," MAIN_ROUTINE says: The program could not deallocate space! "
  PRINT *," Error code 6 in main and execution terminated!             "
END IF

DEALLOCATE(MX1,MX2,MY1,MY2,VELOCITY_k,VELOCITY_kp12,VELOCITY_kp12_new,&
           VELOCITY_kp1,VELOCITY_kp1_new,stat=problem)
IF (problem/=0) THEN 
  PRINT *," MAIN_ROUTINE says: The program could not deallocate space! "
  PRINT *," Error code 7 in main!                                      "
END IF

DEALLOCATE(K1,K2,L1,L2,stat=problem)
IF (problem/=0) THEN 
  PRINT *," MAIN_ROUTINE says: The program could not deallocate space! "
  PRINT *," Error code 8 in main!                                      "
END IF

DEALLOCATE(v1r_ijp1k,v1r_ijk,v1r_ijm1k,v1r_ijp1kp12,v1r_ijkp12,&
           v1r_ijm1kp12,v1r_c_ijp1k,v1r_c_ijk,v1r_c_ijm1k,&
           v1r_c_ijp1kp12,v1r_c_ijkp12,v1r_c_ijm1kp12,&
           v2r_ijp1k,v2r_ijk,v2r_ijm1k,v2r_ijp1kp12,&
           v2r_ijkp12,v2r_ijm1kp12,stat=problem)
IF (problem/=0) THEN 
  PRINT *," MAIN_ROUTINE says: The program could not deallocate space! "
  PRINT *," Error code 9 in main!                                      "
END IF

DEALLOCATE(v1r_ijkp12_new,v2r_ijkp12_new,stat=problem)
IF (problem/=0) THEN 
  PRINT *," MAIN_ROUTINE says: The program could not deallocate space! "
  PRINT *," Error code 10 in main!                                     "
END IF

DEALLOCATE(v1z_ip1jkp12,v1z_ijkp12,v1z_im1jkp12,v1z_ip1jkp1,v1z_ijkp1,&
           v1z_im1jkp1,v1z_c_ip1jkp12,v1z_c_ijkp12,v1z_c_im1jkp12,&
           v1z_c_ip1jkp1,v1z_c_ijkp1,v1z_c_im1jkp1,&
           v2z_ip1jkp12,v2z_ijkp12,v2z_im1jkp12,&
           v2z_ip1jkp1,v2z_ijkp1,v2z_im1jkp1,stat=problem)
IF (problem/=0) THEN 
  PRINT *," MAIN_ROUTINE says: The program could not deallocate space! "
  PRINT *," Error code 11 in main!                                     "
END IF

DEALLOCATE(v1z_ijkp1_new,v2z_ijkp1_new,stat=problem)
IF (problem/=0) THEN 
  PRINT *," MAIN_ROUTINE says: The program could not deallocate space! "
  PRINT *," Error code 12 in main!                                     "
END IF

PRINT *, " Number of grid points (not including Dirichlet "
PRINT *, " boundary points) =   ", count_rms
PRINT *, "                      "
PRINT *, " EXECUTION FINISHED!  "
! --------------------------------------------------------------------------------- !
END PROGRAM MAIN_ROUTINE
! --------------------------------------------------------------------------------- !
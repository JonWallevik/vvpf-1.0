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
! The geometry of the viscometer is defined in this part of the software.           !
! --------------------------------------------------------------------------------- !
PROGRAM MAIN_ROUTINE

USE ROTATION
USE MATRIX
USE WRITE_INFORMATION

IMPLICIT NONE

DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: MX1,MX2,MY1,MY2
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)   :: K1,K2,L1,L2

DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: VELOCITY_k,VELOCITY_kp12,&
                                               VELOCITY_kp12_new,&
                                               VELOCITY_kp1_new,VELOCITY_kp1

DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)   :: &
                             v1r_ijp1k,      v1r_ijk,      v1r_ijm1k,&
                             v1r_ijp1kp12,   v1r_ijkp12,   v1r_ijm1kp12,&
                             v1r_c_ijp1k,    v1r_c_ijk,    v1r_c_ijm1k,&
                             v1r_c_ijp1kp12, v1r_c_ijkp12, v1r_c_ijm1kp12,&
                             v2r_ijp1k,      v2r_ijk,      v2r_ijm1k,&
                             v2r_ijp1kp12,   v2r_ijkp12,   v2r_ijm1kp12

DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)   :: v1r_ijkp12_new,v2r_ijkp12_new,&
                                               v1r_c_ijkp12_new

DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)   :: &
                             v1z_ip1jkp12,   v1z_ijkp12,   v1z_im1jkp12,&
                             v1z_ip1jkp1,    v1z_ijkp1,    v1z_im1jkp1,&
                             v1z_c_ip1jkp12, v1z_c_ijkp12, v1z_c_im1jkp12,&
                             v1z_c_ip1jkp1,  v1z_c_ijkp1,  v1z_c_im1jkp1,&
                             v2z_ip1jkp12,   v2z_ijkp12,   v2z_im1jkp12,&
                             v2z_ip1jkp1,    v2z_ijkp1,    v2z_im1jkp1

DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)   :: v1z_ijkp1_new,v2z_ijkp1_new,&
                                               v1z_c_ijkp1_new

INTEGER,ALLOCATABLE,DIMENSION(:)            :: x_corner,y_corner


DOUBLE PRECISION :: dr,dz,dt,dt_Plastic,dt_Newton,ZERO_TIME,rho,omega,R_i,&
                    R_i_corner,R_o,h_gap,h_corner,h_R_i,tol,tol_Newton,tol_Plastic,&
                    tol_RMS,tol_RMS_active,RMS,vel_norm,small_zero,EPS,a,b,Lambda

INTEGER          :: NX,NX_corner,NX1,NX1_corner,NX2,NX2mNX_corner,NY_corner,NY1,&
                    NY1_corner,NY2,NXY_corner,i,j,k,k_OUTPUT_rms,N_Lambda_MAX,&
                    N_Lambda,MAX_NUMBER_OF_ITERATIONS,problem,count,count_max,count_rms

LOGICAL          :: CONVERGENCE,TIME_INDEPENDENCE,WARNING_SIGN,FALSE_CONVERGENCE
                               
CHARACTER        :: IGNORED_INPUT
! --------------------------------------------------------------------------------- !
PRINT *, "      ______________________________________________        "
PRINT *, "        Viscometric-ViscoPlastic-Flow v1.0  (C3P2)          "
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
PRINT '(4X,A49)', "================================================="
PRINT '(4X,A49)', "          Solving for the C3P2-Geometry          "
PRINT '(4X,A49)', "================================================="
PRINT '(4X,A49)', "                                                 "
! --------------------------------------------------------------------------------- !
! ------------------------- BEGIN OF VARIABLE DECLARATION ------------------------- !
! --------------------------------------------------------------------------------- !
rho = 2354D0 ! kg/m^3 ! <- Note that for a time independent calculations, then 
                      !    the value of rho is not physically important (i.e. 
                      !    mass inertia plays no role). 
! ------------
tol_Newton  = 1.0D-3  ! Used as condition for time independence in the Newtonian case. 
                      ! Also used as tolerance for the successive substitution 
                      ! (in this case it acts as a dummy variable since always two
                      ! successive steps are made for the Newtonian case).
tol_Plastic = 1.0D-10 ! For the successive substitution tolerance (Equation 7.73).
tol_RMS     = 1.0D-30 ! Condition for time independence (Equation 7.75).
! ------------
dt_Newton   = 1.0D-1
dt_Plastic  = 1.0D-5 
count_max   = 15      ! Maximum number of successive (substitution) iterations, 
                      ! for each (pseudotransient) time step k.
R_i         = 0.170D0 ! => 17.0 cm  = Inner radius of viscometer.
R_o         = 0.218D0 ! => 21.8 cm  = Outer radius of viscometer.
h_gap       = 12.0D-2 ! The gap between the parallel plates -> 12 cm.
h_R_i       = 12.0D-3 ! Top edge = 12 mm -> 12.0D-3 
h_corner    =  2.0D-2 ! Size of corner in z-direction (for example 2.0D-2 -> 2 cm).
dr          =  2.0D-3 ! => 2.0 mm  = Spacing between grid points in r-direction.
dz          =  2.0D-3 ! => - " -   = Spacing between grid points in z-direction.
! ------------
ZERO_TIME = 4.00D0  
k_OUTPUT_rms = 100    ! -> Information output every dt_OUTPUT_rms times (to console and file).
! ------------
! The term "small_zero" does usually not have to be changed.
small_zero = 0.1D-7   ! -> Used in relation to screen and file output.
EPS        = 1.0D-15  ! -> Used in relation to vel_norm.
! --------------------------------------------------------------------------------- !
! -------------------------- END OF VARIABLE DECLARATION -------------------------- !
! --------------------------------------------------------------------------------- !
NXY_corner  = IDNINT(h_corner/dz)              ! -> 2.0/0.2      = 10  
NX1_corner  = IDNINT(R_i/dr) + 1               ! -> 17.0/0.2 + 1 = 86  
NX1         = NX1_corner - NXY_corner          ! -> 86 - 10      = 76
NX2         = IDNINT(R_o/dr) + 1               ! -> 21.8/0.2 + 1 = 110  
NX          = NX2 - NX1_corner + 1             ! -> 110 - 86 + 1 = 25 
NY1_corner  = IDNINT(h_gap/dz) + 1             ! -> 12.0/0.2 + 1 = 61
NY1         = NY1_corner + NXY_corner          ! -> 61 + 10      = 71
NY2         = NY1 + IDNINT(h_R_i/dz)           ! -> 71 + 4       = 75
R_i         = dr*DBLE(NX1_corner-1)
R_o         = dr*DBLE(NX2-1)
! --------------------------------------------------------------------------------- !
CALL ANGULAR_VELOCITY(0.0D0,dt,omega)
! --------------------------------------------------------------------------------- !
ALLOCATE(x_corner(NY1-NY1_corner),y_corner(NX1_corner-NX1),stat=problem)
IF (problem/=0) THEN 
  PRINT *," MAIN_ROUTINE says: The program could not allocate space! "
  PRINT *," Error code 0 in main and execution terminated!           "
  STOP
END IF
! --------------------------------------------------------------------------------- !
! x_corner = NX1_corner - (/11,10,9,8,7,6,5,4,3,2/)
DO j = NY1_corner,NY1-1 
  k = j - NY1_corner + 1 ! k=1,10
  x_corner(k) = (NX1_corner - NX1 + 2) - k
! x_corner    = (/11,10,9,8,7,6,5,4,3,2/)
END DO
x_corner = NX1_corner - x_corner
DO i = NX1,NX1_corner-1 
  k = i - NX1 + 1 ! k=1,10
  y_corner(k) = (NY1 - NY1_corner + 1) - k
END DO
y_corner = NY1_corner + y_corner
! --------------------------------------------------------------------------------- !
12 FORMAT(8X,"NX1 = ",(I3,1X),"; NX1_corner = ",(I3,1X),"; NX2 = ",(I3,1X))
13 FORMAT(8X,"NY1 = ",(I3,1X),"; NY1_corner = ",(I3,1X),"; NY2 = ",(I3,1X))
14 FORMAT(8X,"R_i = ",F6.4,"m ; R_o = ",F6.4,"m")
15 FORMAT(8X,"NX = ",(I3,1X),"; dr = ",F7.5,"m","; dz = ",F7.5,"m")
16 FORMAT(8X,"h_gap = ",F6.4,"m ; h_R_i = ",F6.4,"m")
17 FORMAT(8X,"h_corner = ",F6.4,"m")
18 FORMAT(8X,"dt_Plastic = ",E9.3,"s; f_o = ",F6.4,"rps")
PRINT '(7X,A26)',"Geometric and time values:"
PRINT 12, NX1,NX1_corner,NX2
PRINT 13, NY1,NY1_corner,NY2
PRINT 14, R_i,R_o
PRINT 15, NX,dr,dz
PRINT 16, dz*(NY1_corner-1),dz*(NY2-NY1)
PRINT 17, dz*(NY1 - NY1_corner)
PRINT 18, dt_Plastic,omega/(2.0D0*ACOS(-1.0D0))
PRINT *, "    "
! --------------------------------------------------------------------------------- !
IF (NX1_corner.GE.NX2) THEN
  PRINT *, " Inner radius 'R_i' is larger or equal to the outer radius 'R_o'! "
  PRINT *, " TERMINAL ERROR!                                                  "
  STOP
END IF

20 FORMAT(2X,"h_R_i = ",F7.5,"m (h_R_i < 3dz)")
21 FORMAT(2X,"Increase h_R_i up to ",F7.5,"m")
IF (h_R_i.LT.6.0D0*dz) THEN
  PRINT 20,  h_R_i
  PRINT 21,  6.0D0*dz
  PRINT *, " TERMINAL ERROR! "
  STOP
END IF

CALL WARNING_FOR_WRITING(NY2)
! --------------------------------------------------------------------------------- !
MAX_NUMBER_OF_ITERATIONS = IDNINT(ZERO_TIME/dt_Plastic)   

PRINT *, "      -------------------------------------------  "
PRINT "( '         MAX_NUMBER_OF_ITERATIONS: ', I10 )        ",&
                   MAX_NUMBER_OF_ITERATIONS
PRINT *, "      -------------------------------------------  "
PRINT *, "                                                   "

PRINT *, "              ___________________________                   "
WRITE(*,"(A)",ADVANCE="NO") "                PRESS 'ENTER' TO CONTINUE"
PRINT *, " "
READ (*,"(A)") IGNORED_INPUT
PRINT *, " "
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
ALLOCATE(MX1(NX2-2,NX2-2),MX2(NX-2,NX-2),MY1(NY1_corner-2,NY1_corner-2),&
         MY2(NY2-2,NY2-2),VELOCITY_k(NX2,NY2),VELOCITY_kp12(NX2,NY2),&
         VELOCITY_kp12_new(NX2,NY2),VELOCITY_kp1(NX2,NY2),&
         VELOCITY_kp1_new(NX2,NY2),stat=problem)
IF (problem/=0) THEN 
  PRINT *," MAIN_ROUTINE says: The program could not allocate space! "
  PRINT *," Error code 1 in main and execution terminated!           "
  STOP
END IF

ALLOCATE(K1(NX2-2),K2(NX-2),L1(NY1_corner-2),L2(NY2-2),stat=problem)
IF (problem/=0) THEN 
  PRINT *," MAIN_ROUTINE says: The program could not allocate space! "
  PRINT *," Error code 2 in main and execution terminated!           "
  STOP
END IF

ALLOCATE(v1r_ijp1k(NX2),           v1r_ijk(NX2),           v1r_ijm1k(NX2),&
         v1r_ijp1kp12(NX2),        v1r_ijkp12(NX2),        v1r_ijm1kp12(NX2),&
         v1r_c_ijp1k(NX2),         v1r_c_ijk(NX2),         v1r_c_ijm1k(NX2),&
         v1r_c_ijp1kp12(NX2),      v1r_c_ijkp12(NX2),      v1r_c_ijm1kp12(NX2),&
         v2r_ijp1k(NX),            v2r_ijk(NX),            v2r_ijm1k(NX),&     
         v2r_ijp1kp12(NX),         v2r_ijkp12(NX),         v2r_ijm1kp12(NX),stat=problem)
IF (problem/=0) THEN 
  PRINT *," MAIN_ROUTINE says: The program could not allocate space! "
  PRINT *," Error code 3 in main and execution terminated!           "
  STOP
END IF

ALLOCATE(v1r_ijkp12_new(NX2-2),v1r_c_ijkp12_new(NX2-2),v2r_ijkp12_new(NX-2),stat=problem)
IF (problem/=0) THEN 
  PRINT *," MAIN_ROUTINE says: The program could not allocate space! "
  PRINT *," Error code 4 in main and execution terminated!           "
  STOP
END IF

ALLOCATE(v1z_ip1jkp12(NY1_corner), v1z_ijkp12(NY1_corner), v1z_im1jkp12(NY1_corner),&
         v1z_ip1jkp1(NY1_corner),  v1z_ijkp1(NY1_corner),  v1z_im1jkp1(NY1_corner),&
         v1z_c_ip1jkp12(NY2),      v1z_c_ijkp12(NY2),      v1z_c_im1jkp12(NY2),&
         v1z_c_ip1jkp1(NY2),       v1z_c_ijkp1(NY2),       v1z_c_im1jkp1(NY2),&
         v2z_ip1jkp12(NY2),        v2z_ijkp12(NY2),        v2z_im1jkp12(NY2),&
         v2z_ip1jkp1(NY2),         v2z_ijkp1(NY2),         v2z_im1jkp1(NY2),stat=problem)
IF (problem/=0) THEN 
  PRINT *," MAIN_ROUTINE says: The program could not allocate space! "
  PRINT *," Error code 5 in main and execution terminated!           "
  STOP
END IF

ALLOCATE(v1z_ijkp1_new(NY1_corner-2),v1z_c_ijkp1_new(NY2-2),&
         v2z_ijkp1_new(NY2-2),stat=problem)
IF (problem/=0) THEN 
  PRINT *," MAIN_ROUTINE says: The program could not allocate space! "
  PRINT *," Error code 6 in main and execution terminated!           "
  STOP
END IF
! --------------------------------------------------------------------------------- !
! Initialization:
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
v1r_c_ijkp12_new  = 0.0D0
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
v1r_c_ijp1k       = 0.0D0  ! v1r_c -> K1(1:x) & MX1(1:x)*v1r_c_new(1:x)=K1(1:x)
v1r_c_ijk         = 0.0D0  ! ..._c -> ..._corner
v1r_c_ijm1k       = 0.0D0
v1r_c_ijp1kp12    = 0.0D0
v1r_c_ijkp12      = 0.0D0
v1r_c_ijm1kp12    = 0.0D0

v2r_ijp1k         = 0.0D0  ! v2r -> K2 & MX2*v2r_new=K2  
v2r_ijk           = 0.0D0
v2r_ijm1k         = 0.0D0
v2r_ijp1kp12      = 0.0D0
v2r_ijkp12        = 0.0D0
v2r_ijm1kp12      = 0.0D0
! ================================================================================= !
MY1               = 0.0D0
MY2               = 0.0D0
L1                = 0.0D0
L2                = 0.0D0

VELOCITY_kp1      = 0.0D0
VELOCITY_kp1_new  = 0.0D0
! --------------------------------------------------------------------------------- !
v1z_ijkp1_new     = 0.0D0
v1z_c_ijkp1_new   = 0.0D0
v2z_ijkp1_new     = 0.0D0

! v1z is used in L1 and MY1 -> L1=L1(v1z) and MY1=MY1(v1z) to solve the system 
! MY1*v1z_new=L1. "v1z" could be called "v1z_old" since it is the velocity 
! from the previous iteration.
v1z_ip1jkp12      = 0.0D0
v1z_ijkp12        = 0.0D0
v1z_im1jkp12      = 0.0D0
v1z_ip1jkp1       = 0.0D0
v1z_ijkp1         = 0.0D0
v1z_im1jkp1       = 0.0D0
! --------------------------------------------------------------------------------- !
v1z_c_ip1jkp12    = 0.0D0  ! v1z_c -> L1(1:y) & MY1(1:y)*v1z_c_new(1:y)=L1(1:y)
v1z_c_ijkp12      = 0.0D0  ! ..._c -> ..._corner
v1z_c_im1jkp12    = 0.0D0
v1z_c_ip1jkp1     = 0.0D0
v1z_c_ijkp1       = 0.0D0
v1z_c_im1jkp1     = 0.0D0

v2z_ip1jkp12      = 0.0D0  ! v2z -> L2 & MY2*v2z_new=L2
v2z_ijkp12        = 0.0D0
v2z_im1jkp12      = 0.0D0
v2z_ip1jkp1       = 0.0D0
v2z_ijkp1         = 0.0D0
v2z_im1jkp1       = 0.0D0
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
VELOCITY_k(NX2,:) = R_o*omega

! ################################################################################# !
! In Section 7.11.1 is a detailed description of the algorithm, which is used in    !
! the following.                                                                    !
! --------------------------------------------------------------------------------- !
! Linear approximation to speed up convergence:
DO i = 2,NX1
  a = VELOCITY_k(i,1)
  b = VELOCITY_k(i,NY1_corner)
  DO j = 2,NY1_corner-1
    VELOCITY_k(i,j) = a - (a - b)*DBLE(j-1)/DBLE(NY1_corner-1)
  END DO
END DO

DO i = NX1+1,NX1_corner
  NY_corner = y_corner(NX1_corner - i + 1)
  a = VELOCITY_k(i,1)
  b = VELOCITY_k(i,NY_corner)
  DO j = 2,NY_corner-1
    VELOCITY_k(i,j) = a - (a - b)*DBLE(j-1)/DBLE(NY_corner-1)
  END DO
END DO 

DO j = 2,NY2-1
  a = VELOCITY_k(NX1_corner,j)
  b = VELOCITY_k(NX2,j)
  DO i = NX1_corner+1,NX2-1
    VELOCITY_k(i,j) = a - (a - b)*DBLE(i-NX1_corner)/DBLE(NX2-NX1_corner)
  END DO
END DO

! Neumann boundary condition:
VELOCITY_k(NX1_corner+1:NX2-1,NY2) = (4.0D0 * &
     VELOCITY_k(NX1_corner+1:NX2-1,NY2-1) - &
     VELOCITY_k(NX1_corner+1:NX2-1,NY2-2))/3.0D0

! CHECK OUT IF VELOCITY_k IS OK:
! CALL WRITE2FILE_k(VELOCITY_k,NX2)
! STOP
! --------------------------------------------------------------------------------- !
VELOCITY_kp12     = VELOCITY_k
VELOCITY_kp12_new = VELOCITY_k
VELOCITY_kp1_new  = VELOCITY_k
VELOCITY_kp1      = VELOCITY_k
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
  tol_RMS_active = tol_Newton
ELSE
  dt  = dt_Plastic
  tol = tol_Plastic
  tol_RMS_active = tol_RMS
END IF

TIME_INDEPENDENCE = .FALSE.
! Initializing time for each CONTINUATION step:
k = 0
! ================================================================================= !
! =========================== Begin of the time loop ============================== !
! ================================================================================= !
ZERO_TIME_LOOP: DO WHILE (.NOT.TIME_INDEPENDENCE)
CONVERGENCE = .FALSE.
! ---------------------------------------------------------
IF (ABS(MOD((k+1),k_OUTPUT_rms)).LT.small_zero) THEN
  PRINT *,"________________________________________________________"
  PRINT *,"                                                        "
  PRINT *,"  PSEUDO-TRANSIENT time step: k+1 = ",k+1
  PRINT *,"--------------------------------------------------------"
END IF
! --------------------------------------------------------------------------------- !
count = 0 
! ================================================================================= !
! ====================== BEGIN OF SUCCESSIVE SUBSTITUTION ========================= !
! ================================================================================= !
! The iteration loop here is because of the non-linearity of the governing
! Equations 7.22 and 7.23. To come around this problem, the successive substitution
! approach is used (see Section 7.8).
SUBSTITUTION: DO WHILE (.NOT.CONVERGENCE)

! If convergence is a problem, then this might help:
! VELOCITY_kp12 = (VELOCITY_k + VELOCITY_kp1_new)/2
count = count + 1
10 FORMAT(4X,"Successive substitution number = ",1(I3,1X))

IF (ABS(MOD((k+1),k_OUTPUT_rms)).LT.small_zero) THEN
  PRINT 10, count
END IF

! ========================== BEGIN OF X-ITERATION ======================= !
! Iteration is made along r-direction (i.e. along the i-direction as 
! in A(i,j)). It starts at the bottom of the viscometer i=(2:NX2-1) at 
! j = 2 and then move upward with increasing j (see Figure 10.29).

! -------------------------- BEGIN OF X-GAP-AREA ------------------------ !
DO j = 2,NY1_corner-1
  v1r_ijp1k    = VELOCITY_k(:,j+1)
  v1r_ijk      = VELOCITY_k(:,j)
  v1r_ijm1k    = VELOCITY_k(:,j-1)
  v1r_ijp1kp12 = VELOCITY_kp12(:,j+1)
  v1r_ijkp12   = VELOCITY_kp12(:,j)
  v1r_ijm1kp12 = VELOCITY_kp12(:,j-1)

  CALL MATRIX_UPDATE_X(rho,k,dt,Lambda,dr,0.0D0,dz,NX2,v1r_ijp1k,v1r_ijk,&
                    v1r_ijm1k,v1r_ijp1kp12,v1r_ijkp12,v1r_ijm1kp12,MX1,K1)
  CALL MATRIX_SOLVER(MX1,K1,v1r_ijkp12_new,NX2-2)

  VELOCITY_kp12_new(2:NX2-1,j) = v1r_ijkp12_new
END DO
! -------------------------- END OF X-GAP-AREA -------------------------- !

! -------------------------- BEGIN OF X-CORNER -------------------------- !
DO j = NY1_corner,NY1-1
  NX_corner                       = x_corner(j - NY1_corner + 1)
  R_i_corner                      = NX_corner*dr
  NX2mNX_corner                   = NX2 - NX_corner
  v1r_c_ijp1k(1:NX2mNX_corner)    = VELOCITY_k(NX_corner+1:NX2,j+1)
  v1r_c_ijk(1:NX2mNX_corner)      = VELOCITY_k(NX_corner+1:NX2,j)
  v1r_c_ijm1k(1:NX2mNX_corner)    = VELOCITY_k(NX_corner+1:NX2,j-1)
  v1r_c_ijp1kp12(1:NX2mNX_corner) = VELOCITY_kp12(NX_corner+1:NX2,j+1)
  v1r_c_ijkp12(1:NX2mNX_corner)   = VELOCITY_kp12(NX_corner+1:NX2,j)
  v1r_c_ijm1kp12(1:NX2mNX_corner) = VELOCITY_kp12(NX_corner+1:NX2,j-1)

  CALL MATRIX_UPDATE_X(rho,k,dt,Lambda,dr,R_i_corner,dz,NX2mNX_corner,v1r_c_ijp1k,&
           v1r_c_ijk,v1r_c_ijm1k,v1r_c_ijp1kp12,v1r_c_ijkp12,v1r_c_ijm1kp12,MX1,K1)
  CALL MATRIX_SOLVER(MX1,K1,v1r_c_ijkp12_new,NX2mNX_corner-2)
  
  VELOCITY_kp12_new(NX_corner+2:NX2-1,j) = v1r_c_ijkp12_new(1:NX2mNX_corner-2)
END DO
! -------------------------- END OF X-CORNER ---------------------------- !

! -------------------------- BEGIN OF X-OPEN-AREA ----------------------- !
DO j = NY1,NY2-1
  v2r_ijp1k    = VELOCITY_k(NX1_corner:NX2,j+1)
  v2r_ijk      = VELOCITY_k(NX1_corner:NX2,j)
  v2r_ijm1k    = VELOCITY_k(NX1_corner:NX2,j-1)
  v2r_ijp1kp12 = VELOCITY_kp12(NX1_corner:NX2,j+1)
  v2r_ijkp12   = VELOCITY_kp12(NX1_corner:NX2,j)
  v2r_ijm1kp12 = VELOCITY_kp12(NX1_corner:NX2,j-1)

  CALL MATRIX_UPDATE_X(rho,k,dt,Lambda,dr,R_i,dz,NX,v2r_ijp1k,v2r_ijk,&
                 v2r_ijm1k,v2r_ijp1kp12,v2r_ijkp12,v2r_ijm1kp12,MX2,K2)
  CALL MATRIX_SOLVER(MX2,K2,v2r_ijkp12_new,NX-2)
  
  VELOCITY_kp12_new(NX1_corner+1:NX2-1,j) = v2r_ijkp12_new
END DO

j = NY2
v2r_ijp1k    = VELOCITY_k(NX1_corner:NX2,j-1)    ! => v(i,j+1) = v(i,j-1)
v2r_ijk      = VELOCITY_k(NX1_corner:NX2,j)
v2r_ijm1k    = VELOCITY_k(NX1_corner:NX2,j-1)    ! => v(i,j+1) = v(i,j-1)
v2r_ijp1kp12 = VELOCITY_kp12(NX1_corner:NX2,j-1) ! => v(i,j+1) = v(i,j-1)
v2r_ijkp12   = VELOCITY_kp12(NX1_corner:NX2,j)
v2r_ijm1kp12 = VELOCITY_kp12(NX1_corner:NX2,j-1) ! => v(i,j+1) = v(i,j-1)  

CALL MATRIX_UPDATE_X(rho,k,dt,Lambda,dr,R_i,dz,NX,v2r_ijp1k,v2r_ijk,&
               v2r_ijm1k,v2r_ijp1kp12,v2r_ijkp12,v2r_ijm1kp12,MX2,K2)
CALL MATRIX_SOLVER(MX2,K2,v2r_ijkp12_new,NX-2)

VELOCITY_kp12_new(NX1_corner+1:NX2-1,j) = v2r_ijkp12_new
! -------------------------- END OF X-OPEN-AREA ------------------------- !

! --------- PAUSE FOR DEBUGGING --------- 
! CALL WRITE2FILE_k(VELOCITY_kp12_new,NX2)
! WRITE(*,"(A)",ADVANCE="NO") " PRESS 'ENTER' TO CONTINUE "
! PRINT *, " "
! READ (*,"(A)") IGNORED_INPUT
! PRINT *, " "

! ========================== END OF X-ITERATION ========================= !

! Updating ..._kp12:
VELOCITY_kp12 = VELOCITY_kp12_new

! ========================== BEGIN OF Y-ITERATION ======================= !
! Iteration is made along z-direction (i.e. along the j-direction as 
! in A(i,j)). It starts at the right side of the viscometer j=(2:NY2-1) at 
! i = NX2-1 and then moves to the left with decreasing i (see Figure 10.29).

! -------------------------- BEGIN OF Y-OPEN-AREA ----------------------- !
DO i = NX2-1,NX1_corner+1,-1
  v2z_ip1jkp12 = VELOCITY_kp12(i+1,:)
  v2z_ijkp12   = VELOCITY_kp12(i,:)
  v2z_im1jkp12 = VELOCITY_kp12(i-1,:)
  v2z_ip1jkp1  = VELOCITY_kp1(i+1,:)
  v2z_ijkp1    = VELOCITY_kp1(i,:)
  v2z_im1jkp1  = VELOCITY_kp1(i-1,:)

  CALL MATRIX_UPDATE_Y(rho,k,dt,Lambda,dr,i,dz,NY2,v2z_ip1jkp12,v2z_ijkp12,&
               v2z_im1jkp12,v2z_ip1jkp1,v2z_ijkp1,v2z_im1jkp1,MY2,L2,.TRUE.)
  CALL MATRIX_SOLVER(MY2,L2,v2z_ijkp1_new,NY2-2)

  VELOCITY_kp1_new(i,2:NY2-1) = v2z_ijkp1_new
END DO

! v(i,j+1) = v(i,j-1) => 
VELOCITY_kp1_new(NX1_corner+1:NX2-1,NY2) = VELOCITY_kp1_new(NX1_corner+1:NX2-1,NY2-2)
! -------------------------- END OF Y-OPEN-AREA ------------------------- !

! -------------------------- BEGIN OF Y-CORNER -------------------------- !
DO i = NX1_corner,NX1+1,-1
  NY_corner                   = y_corner(NX1_corner - i + 1)
  v1z_c_ip1jkp12(1:NY_corner) = VELOCITY_kp12(i+1,1:NY_corner)
  v1z_c_ijkp12(1:NY_corner)   = VELOCITY_kp12(i,1:NY_corner)
  v1z_c_im1jkp12(1:NY_corner) = VELOCITY_kp12(i-1,1:NY_corner)
  v1z_c_ip1jkp1(1:NY_corner)  = VELOCITY_kp1(i+1,1:NY_corner)
  v1z_c_ijkp1(1:NY_corner)    = VELOCITY_kp1(i,1:NY_corner)
  v1z_c_im1jkp1(1:NY_corner)  = VELOCITY_kp1(i-1,1:NY_corner)

  CALL MATRIX_UPDATE_Y(rho,k,dt,Lambda,dr,i,dz,NY_corner,v1z_c_ip1jkp12,v1z_c_ijkp12,&
                v1z_c_im1jkp12,v1z_c_ip1jkp1,v1z_c_ijkp1,v1z_c_im1jkp1,MY2,L2,.FALSE.)
  CALL MATRIX_SOLVER(MY2,L2,v1z_c_ijkp1_new,NY_corner-2)

  VELOCITY_kp1_new(i,2:NY_corner-1) = v1z_c_ijkp1_new(1:NY_corner-2)
END DO
! -------------------------- END OF Y-CORNER ---------------------------- !

! -------------------------- BEGIN OF Y-GAP-AREA ------------------------ !
DO i = NX1,2,-1 
  v1z_ip1jkp12 = VELOCITY_kp12(i+1,1:NY1_corner)
  v1z_ijkp12   = VELOCITY_kp12(i,1:NY1_corner)
  v1z_im1jkp12 = VELOCITY_kp12(i-1,1:NY1_corner)
  v1z_ip1jkp1  = VELOCITY_kp1(i+1,1:NY1_corner)
  v1z_ijkp1    = VELOCITY_kp1(i,1:NY1_corner)
  v1z_im1jkp1  = VELOCITY_kp1(i-1,1:NY1_corner)

  CALL MATRIX_UPDATE_Y(rho,k,dt,Lambda,dr,i,dz,NY1_corner,v1z_ip1jkp12,v1z_ijkp12,&
                     v1z_im1jkp12,v1z_ip1jkp1,v1z_ijkp1,v1z_im1jkp1,MY1,L1,.FALSE.)
  CALL MATRIX_SOLVER(MY1,L1,v1z_ijkp1_new,NY1_corner-2)
  
  VELOCITY_kp1_new(i,2:NY1_corner-1) = v1z_ijkp1_new
END DO
! -------------------------- END OF Y-GAP-AREA -------------------------- !

! --------- PAUSE FOR DEBUGGING --------- 
! CALL WRITE2FILE_k(VELOCITY_kp12_new,NX2)
! WRITE(*,"(A)",ADVANCE="NO") " PRESS 'ENTER' TO CONTINUE "
! PRINT *, " "
! READ (*,"(A)") IGNORED_INPUT
! PRINT *, " "

! ========================== END OF Y-ITERATION ========================= !

CONVERGENCE = .TRUE.

! Settings for testing of convergence:
RMS      = 0.0D0
vel_norm = 1.0D0

S1: DO j = 2,NY1_corner-1   
      DO i = 2,NX1_corner
        vel_norm =  (VELOCITY_kp1_new(i,j) + VELOCITY_kp1(i,j))/2.0D0 + EPS
        RMS      = ((VELOCITY_kp1_new(i,j) - VELOCITY_kp1(i,j))/vel_norm)**2.0D0
        IF (RMS > tol) THEN
          CONVERGENCE = .FALSE.
          EXIT S1
        END IF
      END DO
    END DO S1

IF (CONVERGENCE) THEN
S2: DO j = NY1_corner,NY1-1
      NX_corner = x_corner(j - NY1_corner + 1)
      DO i = NX_corner+2,NX1_corner
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
S3: DO j = NY1_corner,NY2-1
      DO i = NX1_corner+1,NX2-1
        vel_norm =  (VELOCITY_kp1_new(i,j) + VELOCITY_kp1(i,j))/2.0D0 + EPS
        RMS      = ((VELOCITY_kp1_new(i,j) - VELOCITY_kp1(i,j))/vel_norm)**2.0D0
        IF (RMS > tol) THEN
          CONVERGENCE = .FALSE.
          EXIT S3
        END IF
      END DO
    END DO S3
END IF

IF (count == count_max) THEN
  CONVERGENCE       = .TRUE.
  FALSE_CONVERGENCE = .TRUE.
  PRINT *, "WARNING: FALSE CONVERGENCE! TIME STEP k= ", k
  PRINT *, "Maximum amount of successive substitutions is = ", count_max
  PRINT *, " ---------------------------------------------------------- "
  PRINT *, " RECOMMENDATION: Kill this application and reduce the time  "
  PRINT *, " step by an order of magnitude: dt -> dt/10                 "
  PRINT *, " ---------------------------------------------------------- "
END IF

! Updating ..._kp1
VELOCITY_kp1 = VELOCITY_kp1_new

! CALL WRITE2FILE_k(VELOCITY_k,NX2)

END DO SUBSTITUTION 
! ================================================================================= !
! ======================= END OF SUCCESSIVE SUBSTITUTION ========================== !
! ================================================================================= !
count_rms = 0
RMS       = 0.0D0
vel_norm  = 1.0D0

DO j = 2,NY1_corner-1   
  DO i = 2,NX1_corner
    count_rms = count_rms + 1
    vel_norm  =  (VELOCITY_kp1(i,j) + VELOCITY_k(i,j))/2.0D0 + EPS
    RMS       = ((VELOCITY_kp1(i,j) - VELOCITY_k(i,j))/vel_norm)**2.0D0 + RMS
  END DO
END DO

DO j = NY1_corner,NY1-1
  NX_corner = x_corner(j - NY1_corner + 1)
  DO i = NX_corner+2,NX1_corner
    count_rms = count_rms + 1
    vel_norm  =  (VELOCITY_kp1(i,j) + VELOCITY_k(i,j))/2.0D0 + EPS
    RMS       = ((VELOCITY_kp1(i,j) - VELOCITY_k(i,j))/vel_norm)**2.0D0 + RMS
  END DO
END DO

DO j = NY1_corner,NY2-1
  DO i = NX1_corner+1,NX2-1
    count_rms = count_rms + 1
    vel_norm  =  (VELOCITY_kp1(i,j) + VELOCITY_k(i,j))/2.0D0 + EPS
    RMS       = ((VELOCITY_kp1(i,j) - VELOCITY_k(i,j))/vel_norm)**2.0D0 + RMS
  END DO
END DO

RMS = DSQRT(RMS/count_rms)
! ------------
IF (ABS(MOD((k+1),k_OUTPUT_rms)).LT.small_zero) THEN
  PRINT *, "   RMS   =", RMS
  CALL WRITE2FILE_rms(k+1,RMS)
END IF
! ------------
IF (RMS.LT.tol_RMS_active) THEN
  TIME_INDEPENDENCE = .TRUE.
ELSE 
  TIME_INDEPENDENCE = .FALSE.
END IF
! ------------
IF (k == MAX_NUMBER_OF_ITERATIONS) THEN
  TIME_INDEPENDENCE = .TRUE.
  WARNING_SIGN      = .TRUE. 
END IF
! ------------
VELOCITY_k = VELOCITY_kp1_new
k = k + 1
! ------------
! CALL WRITE2FILE_k(VELOCITY_k,NX2)
! STOP
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
  PRINT *, " RECOMMENDATION: Rerun this application with reduced time  "
  PRINT *, " step. Try order of magnitude less: dt -> dt/10            "
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

PRINT *, " Number of grid points:", count_rms
PRINT *, "                       "
PRINT *, " WRITING INFORMATION TO FILE, STAND BY... "
CALL WRITE2FILE_k(VELOCITY_k,NX2)
CALL WRITE2FILE_kp1(VELOCITY_k,NX1,NX1_corner,NX2,NY1,NY1_corner,NY2,dr,dz)
PRINT *, " ...DONE! "

! --------------------------------------------------------------------------------- !
! CLEARING SOME MAJOR VARIABLES FROM THE RANDOM ACCESS MEMORY:
DEALLOCATE(x_corner,y_corner,stat=problem)
IF (problem/=0) THEN 
  PRINT *," MAIN_ROUTINE says: The program could not deallocate space! "
  PRINT *," Error code 7 in main and execution terminated!             "
END IF

DEALLOCATE(MX1,MX2,MY1,MY2,VELOCITY_k,VELOCITY_kp12,VELOCITY_kp12_new,&
           VELOCITY_kp1,VELOCITY_kp1_new,stat=problem)
IF (problem/=0) THEN 
  PRINT *," MAIN_ROUTINE says: The program could not deallocate space! "
  PRINT *," Error code 8 in main!                                      "
END IF

DEALLOCATE(K1,K2,L1,L2,stat=problem)
IF (problem/=0) THEN 
  PRINT *," MAIN_ROUTINE says: The program could not deallocate space! "
  PRINT *," Error code 9 in main!                                      "
END IF

DEALLOCATE(v1r_ijp1k,      v1r_ijk,      v1r_ijm1k,&
           v1r_ijp1kp12,   v1r_ijkp12,   v1r_ijm1kp12,&
           v1r_c_ijp1k,    v1r_c_ijk,    v1r_c_ijm1k,&
           v1r_c_ijp1kp12, v1r_c_ijkp12, v1r_c_ijm1kp12,&
           v2r_ijp1k,      v2r_ijk,      v2r_ijm1k,&
           v2r_ijp1kp12,   v2r_ijkp12,   v2r_ijm1kp12,stat=problem)
IF (problem/=0) THEN 
  PRINT *," MAIN_ROUTINE says: The program could not deallocate space! "
  PRINT *," Error code 10 in main!                                     "
END IF

DEALLOCATE(v1r_ijkp12_new,v1r_c_ijkp12_new,v2r_ijkp12_new,stat=problem)
IF (problem/=0) THEN 
  PRINT *," MAIN_ROUTINE says: The program could not deallocate space! "
  PRINT *," Error code 11 in main!                                     "
END IF

DEALLOCATE(v1z_ip1jkp12,   v1z_ijkp12,   v1z_im1jkp12,&
           v1z_ip1jkp1,    v1z_ijkp1,    v1z_im1jkp1,&
           v1z_c_ip1jkp12, v1z_c_ijkp12, v1z_c_im1jkp12,&
           v1z_c_ip1jkp1,  v1z_c_ijkp1,  v1z_c_im1jkp1,&
           v2z_ip1jkp12,   v2z_ijkp12,   v2z_im1jkp12,&
           v2z_ip1jkp1,    v2z_ijkp1,    v2z_im1jkp1,stat=problem)
IF (problem/=0) THEN 
  PRINT *," MAIN_ROUTINE says: The program could not deallocate space! "
  PRINT *," Error code 12 in main!                                     "
END IF

DEALLOCATE(v1z_ijkp1_new,v1z_c_ijkp1_new,v2z_ijkp1_new,stat=problem)
IF (problem/=0) THEN 
  PRINT *," MAIN_ROUTINE says: The program could not deallocate space! "
  PRINT *," Error code 13 in main!                                     "
END IF

PRINT *, " EXECUTION FINISHED! "
! --------------------------------------------------------------------------------- !
END PROGRAM MAIN_ROUTINE
! --------------------------------------------------------------------------------- !
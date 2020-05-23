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
! File name: update.f90 (MODULE)                                                    !
! This file sets up the system of algebraic Equations 7.28 to 7.31. This file also  !
! contains the Thomas algorithm that is used in solving this system.                !
! --------------------------------------------------------------------------------- !
MODULE MATRIX
  USE RATE_OF_SHEAR
  USE SHEAR_VISCOSITY
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: MATRIX_UPDATE_X,MATRIX_UPDATE_Y,MATRIX_SOLVER
CONTAINS
! ================================================================================= !
SUBROUTINE MATRIX_UPDATE_X(rho,k,dt,Lambda,dr,R_i,dz,NX,vr_ijp1k,vr_ijk,&
                           vr_ijm1k,vr_ijp1kp12,vr_ijkp12,vr_ijm1kp12,&
                           FMSR_ijp1_t,FMSR_ij_t,FMSR_ijm1_t,&
                           FMCR_ijp1_t,FMCR_ij_t,FMCR_ijm1_t,M,K_M)

DOUBLE PRECISION,DIMENSION(:,:),INTENT(OUT) :: M
DOUBLE PRECISION,DIMENSION(:),INTENT(OUT)   :: K_M

DOUBLE PRECISION,DIMENSION(:),INTENT(IN)    :: vr_ijp1k,vr_ijk,&
                                               vr_ijm1k,vr_ijp1kp12,&
                                               vr_ijkp12,vr_ijm1kp12,&
                                               FMSR_ijp1_t,FMSR_ij_t,&
                                               FMSR_ijm1_t,&
                                               FMCR_ijp1_t,FMCR_ij_t,&
                                               FMCR_ijm1_t

DOUBLE PRECISION,INTENT(IN)                 :: dt,dr,R_i,dz,Lambda,rho
INTEGER,INTENT(IN)                          :: k,NX


DOUBLE PRECISION :: BETA,THETA,CHI,rp1,r,rm1

DOUBLE PRECISION :: V_ijk,V_ip1jk,V_im1jk,V_ijp1k,V_ijm1k,&
                    V_ip1jp1k,V_ip1jm1k,V_im1jp1k,V_im1jm1k,&
                    V_ijkp12,V_ip1jkp12,V_im1jkp12,V_ijp1kp12,V_ijm1kp12,&
                    V_ip1jp1kp12,V_ip1jm1kp12,V_im1jp1kp12,V_im1jm1kp12

DOUBLE PRECISION :: SR_ijk,SR_ip12jk,SR_im12jk,&
                    SR_ijp12k,SR_ijm12k,SR_ijkp12,&
                    SR_ip12jkp12,SR_im12jkp12,SR_ijp12kp12,&
                    SR_ijm12kp12

DOUBLE PRECISION :: FMSR_im1j,FMSR_ij,FMSR_ip1j,FMSR_im1jp1,FMSR_ijp1,FMSR_ip1jp1,&
                    FMSR_im1jm1,FMSR_ijm1,FMSR_ip1jm1,FMSR_ip12jp12,FMSR_ip12jm12,&
                    FMSR_im12jp12,FMSR_im12jm12,FMSR_ip12j,FMSR_im12j,&
                    FMSR_ijp12,FMSR_ijm12

DOUBLE PRECISION :: FMCR_im1j,FMCR_ij,FMCR_ip1j,FMCR_im1jp1,FMCR_ijp1,FMCR_ip1jp1,&
                    FMCR_im1jm1,FMCR_ijm1,FMCR_ip1jm1,FMCR_ip12jp12,FMCR_ip12jm12,&
                    FMCR_im12jp12,FMCR_im12jm12,FMCR_ip12j,FMCR_im12j,&
                    FMCR_ijp12,FMCR_ijm12

DOUBLE PRECISION :: ETA_ijk,ETA_ip12jk,ETA_im12jk,ETA_ijp12k,&
                    ETA_ijm12k,ETA_ijkp12,ETA_ip12jkp12,&
                    ETA_im12jkp12,ETA_ijp12kp12,ETA_ijm12kp12

DOUBLE PRECISION :: A_kp12,B_kp12,C_kp12,D_k,E_k,F_k,KK
DOUBLE PRECISION :: TIME_k,TIME_kp12
INTEGER          :: i
! --------------------------------------------------------------------------------- !
BETA      = dt/(2.0D0*dr*rho)
CHI       = dt/(2.0D0*dz*rho)
TIME_k    = DBLE(k)*dt
TIME_kp12 = (DBLE(k) + 0.5D0)*dt

M         =  0.0D0
K_M       =  0.0D0
! --------------------------------------------------------------------------------- !
! i = 1    => i = 2 in main.f90, i.e. near the center of the viscometer.
! i = NX-2 => j = NX2-1 in main.f90, i.e. near the outer cylinder.
DO i = 1,NX-2

  rp1   =  DBLE(i+1)*dr + R_i
  r     =  DBLE(i)*dr + R_i
  rm1   =  DBLE(i-1)*dr + R_i
  THETA =  dt/(r*rho)

  V_im1jk       = vr_ijk(i)
  V_ijk         = vr_ijk(i+1)
  V_ip1jk       = vr_ijk(i+2)
  V_im1jp1k     = vr_ijp1k(i)
  V_ijp1k       = vr_ijp1k(i+1)
  V_ip1jp1k     = vr_ijp1k(i+2)
  V_im1jm1k     = vr_ijm1k(i)
  V_ijm1k       = vr_ijm1k(i+1)
  V_ip1jm1k     = vr_ijm1k(i+2)

  V_im1jkp12    = vr_ijkp12(i)
  V_ijkp12      = vr_ijkp12(i+1)
  V_ip1jkp12    = vr_ijkp12(i+2)
  V_im1jp1kp12  = vr_ijp1kp12(i)
  V_ijp1kp12    = vr_ijp1kp12(i+1)
  V_ip1jp1kp12  = vr_ijp1kp12(i+2)
  V_im1jm1kp12  = vr_ijm1kp12(i)
  V_ijm1kp12    = vr_ijm1kp12(i+1)
  V_ip1jm1kp12  = vr_ijm1kp12(i+2)

  FMSR_im1j     = FMSR_ij_t(i)
  FMSR_ij       = FMSR_ij_t(i+1)
  FMSR_ip1j     = FMSR_ij_t(i+2)
  FMSR_im1jp1   = FMSR_ijp1_t(i)
  FMSR_ijp1     = FMSR_ijp1_t(i+1)
  FMSR_ip1jp1   = FMSR_ijp1_t(i+2)
  FMSR_im1jm1   = FMSR_ijm1_t(i)
  FMSR_ijm1     = FMSR_ijm1_t(i+1)
  FMSR_ip1jm1   = FMSR_ijm1_t(i+2)

  FMSR_ip12jp12 = (FMSR_ij   + FMSR_ip1j + FMSR_ip1jp1   + FMSR_ijp1)/4.0D0
  FMSR_ip12jm12 = (FMSR_ij   + FMSR_ip1j + FMSR_ip1jm1   + FMSR_ijm1)/4.0D0
  FMSR_im12jp12 = (FMSR_ij   + FMSR_im1j + FMSR_im1jp1   + FMSR_ijp1)/4.0D0
  FMSR_im12jm12 = (FMSR_ij   + FMSR_im1j + FMSR_im1jm1   + FMSR_ijm1)/4.0D0

  FMSR_ip12j    = (FMSR_ij   + FMSR_ip1j + FMSR_ip12jp12 + FMSR_ip12jm12)/4.0D0
  FMSR_im12j    = (FMSR_im1j + FMSR_ij   + FMSR_im12jp12 + FMSR_im12jm12)/4.0D0
  FMSR_ijp12    = (FMSR_ij   + FMSR_ijp1 + FMSR_ip12jp12 + FMSR_im12jp12)/4.0D0
  FMSR_ijm12    = (FMSR_ijm1 + FMSR_ij   + FMSR_ip12jm12 + FMSR_im12jm12)/4.0D0

  FMCR_im1j     = FMCR_ij_t(i)
  FMCR_ij       = FMCR_ij_t(i+1)
  FMCR_ip1j     = FMCR_ij_t(i+2)
  FMCR_im1jp1   = FMCR_ijp1_t(i)
  FMCR_ijp1     = FMCR_ijp1_t(i+1)
  FMCR_ip1jp1   = FMCR_ijp1_t(i+2)
  FMCR_im1jm1   = FMCR_ijm1_t(i)
  FMCR_ijm1     = FMCR_ijm1_t(i+1)
  FMCR_ip1jm1   = FMCR_ijm1_t(i+2)

  FMCR_ip12jp12 = (FMCR_ij   + FMCR_ip1j + FMCR_ip1jp1   + FMCR_ijp1)/4.0D0
  FMCR_ip12jm12 = (FMCR_ij   + FMCR_ip1j + FMCR_ip1jm1   + FMCR_ijm1)/4.0D0
  FMCR_im12jp12 = (FMCR_ij   + FMCR_im1j + FMCR_im1jp1   + FMCR_ijp1)/4.0D0
  FMCR_im12jm12 = (FMCR_ij   + FMCR_im1j + FMCR_im1jm1   + FMCR_ijm1)/4.0D0

  FMCR_ip12j    = (FMCR_ij   + FMCR_ip1j + FMCR_ip12jp12 + FMCR_ip12jm12)/4.0D0
  FMCR_im12j    = (FMCR_im1j + FMCR_ij   + FMCR_im12jp12 + FMCR_im12jm12)/4.0D0
  FMCR_ijp12    = (FMCR_ij   + FMCR_ijp1 + FMCR_ip12jp12 + FMCR_im12jp12)/4.0D0
  FMCR_ijm12    = (FMCR_ijm1 + FMCR_ij   + FMCR_ip12jm12 + FMCR_im12jm12)/4.0D0

  CALL ROS(rp1,r,rm1,dr,dz,V_ijk,V_ip1jk,V_im1jk,V_ijp1k,V_ijm1k,&
           V_ip1jp1k,V_ip1jm1k,V_im1jp1k,V_im1jm1k,&
           SR_ijk,SR_ip12jk,SR_im12jk,SR_ijp12k,SR_ijm12k)

  CALL ROS(rp1,r,rm1,dr,dz,V_ijkp12,V_ip1jkp12,V_im1jkp12,V_ijp1kp12,V_ijm1kp12,&
           V_ip1jp1kp12,V_ip1jm1kp12,V_im1jp1kp12,V_im1jm1kp12,&
           SR_ijkp12,SR_ip12jkp12,SR_im12jkp12,SR_ijp12kp12,SR_ijm12kp12)

  ! FMSR is relative to k and since we are calculating ETA_k, the "ALPHA_I*SR*dt"
  ! term must be set equal to zero! -> dt = 0.0D0.
  ! We use dt = 0.0D0 because we are calculating ETA for the time step k. I.e.
  ! we want to have Gamma = ALPHA_II*FMSR and not Gamma = ALPHA_II*(FMSR + ALPHA_I*SR*dt), 
  ! because FMSR is already fully updated relative to the time step k.
  ! [ The same consideration applies for FMCR, BETA_I, BETA_II and H(SR) ]. 
  CALL ETA(0.0D0,TIME_k,Lambda,SR_ijk,SR_ip12jk,SR_im12jk,SR_ijp12k,SR_ijm12k,&
           FMSR_ij,FMSR_ip12j,FMSR_im12j,FMSR_ijp12,FMSR_ijm12,&
           FMCR_ij,FMCR_ip12j,FMCR_im12j,FMCR_ijp12,FMCR_ijm12,&
           ETA_ijk,ETA_ip12jk,ETA_im12jk,ETA_ijp12k,ETA_ijm12k)

  ! FMSR is relative to k and since we are calculating ETA_kp12, the "ALPHA_I*SR*dt"
  ! term must be set equal to "ALPHA_I*SR*(dt/2)"! -> dt = dt/2.
  ! FMSR is relative to k and we are only proceeding from k to k+1/2, hence dt -> dt/2.
  ! [ The same consideration applies for FMCR, BETA_I, BETA_II and H(SR) ]. 
  CALL ETA(dt/2.0D0,TIME_kp12,Lambda,SR_ijkp12,SR_ip12jkp12,SR_im12jkp12,&
           SR_ijp12kp12,SR_ijm12kp12,FMSR_ij,FMSR_ip12j,FMSR_im12j,&
           FMSR_ijp12,FMSR_ijm12,FMCR_ij,FMCR_ip12j,FMCR_im12j,&
           FMCR_ijp12,FMCR_ijm12,ETA_ijkp12,ETA_ip12jkp12,&
           ETA_im12jkp12,ETA_ijp12kp12,ETA_ijm12kp12)

  ! Equations 7.32 to 7.34:
  A_kp12 = (BETA*ETA_ip12jkp12)*(1.0D0/dr - 1.0D0/(rp1+r)) + (THETA*ETA_ijkp12)/(2.0D0*dr)
  B_kp12 = (BETA*ETA_ip12jkp12)*(1.0D0/dr + 1.0D0/(rp1+r)) + &
           (BETA*ETA_im12jkp12)*(1.0D0/dr - 1.0D0/(r+rm1)) + (THETA*ETA_ijkp12)/r 
  C_kp12 = (BETA*ETA_im12jkp12)*(1.0D0/dr + 1.0D0/(r+rm1)) - (THETA*ETA_ijkp12)/(2.0D0*dr)

  ! Equations 7.35 to 7.37:
  D_k    = (CHI*ETA_ijp12k)/dz
  E_k    = (CHI*ETA_ijp12k)/dz + (CHI*ETA_ijm12k)/dz  
  F_k    = (CHI*ETA_ijm12k)/dz

  ! Equation 7.29:
  KK = - D_k*V_ijp1k - (1.0D0 - E_k)*V_ijk - F_k*V_ijm1k
  
  ! Equation 7.28:
  IF (i == 1) THEN
    M(i,i)   = -(1.0D0 + B_kp12)
    M(i,i+1) = A_kp12
    K_M(i)   = KK - C_kp12*V_im1jkp12
  ELSE IF (i == NX-2) THEN
    M(i,i-1) = C_kp12
    M(i,i)   = -(1.0D0 + B_kp12)
    K_M(i)   = KK - A_kp12*V_ip1jkp12
  ELSE 
    M(i,i-1) = C_kp12
    M(i,i)   = -(1.0D0 + B_kp12)
    M(i,i+1) = A_kp12
    K_M(i)   = KK
  END IF

END DO
! --------------------------------------------------------------------------------- !
RETURN
END SUBROUTINE MATRIX_UPDATE_X
! ================================================================================= !
SUBROUTINE MATRIX_UPDATE_Y(rho,k,dt,Lambda,dr,i,dz,NY,vz_ip1jkp12,vz_ijkp12,&
                           vz_im1jkp12,vz_ip1jkp1,vz_ijkp1,vz_im1jkp1,&
                           FMSR_ip1j_t,FMSR_ij_t,FMSR_im1j_t,&
                           FMCR_ip1j_t,FMCR_ij_t,FMCR_im1j_t,M,L,Neumann)

DOUBLE PRECISION,DIMENSION(:,:),INTENT(OUT) :: M
DOUBLE PRECISION,DIMENSION(:),INTENT(OUT)   :: L

DOUBLE PRECISION,DIMENSION(:),INTENT(IN)    :: vz_im1jkp12,vz_ijkp12,&
                                               vz_ip1jkp12,vz_im1jkp1,&
                                               vz_ijkp1,vz_ip1jkp1,&
                                               FMSR_ip1j_t,FMSR_ij_t,&
                                               FMSR_im1j_t,&
                                               FMCR_ip1j_t,FMCR_ij_t,&
                                               FMCR_im1j_t

DOUBLE PRECISION,INTENT(IN)                 :: dt,dr,dz,Lambda,rho
INTEGER,INTENT(IN)                          :: k,NY,i

LOGICAL,INTENT(IN)                          :: Neumann

DOUBLE PRECISION :: BETA,THETA,CHI,rp1,r,rm1

DOUBLE PRECISION :: V_ijkp12,V_ip1jkp12,V_im1jkp12,V_ijp1kp12,V_ijm1kp12,&
                    V_ip1jp1kp12,V_ip1jm1kp12,V_im1jp1kp12,V_im1jm1kp12,&
                    V_ijkp1,V_ip1jkp1,V_im1jkp1,V_ijp1kp1,V_ijm1kp1,&
                    V_ip1jp1kp1,V_ip1jm1kp1,V_im1jp1kp1,V_im1jm1kp1

DOUBLE PRECISION :: SR_ijkp12,SR_ip12jkp12,SR_im12jkp12,&
                    SR_ijp12kp12,SR_ijm12kp12,SR_ijkp1,&
                    SR_ip12jkp1,SR_im12jkp1,SR_ijp12kp1,&
                    SR_ijm12kp1

DOUBLE PRECISION :: FMSR_im1j,FMSR_ij,FMSR_ip1j,FMSR_im1jp1,FMSR_ijp1,FMSR_ip1jp1,&
                    FMSR_im1jm1,FMSR_ijm1,FMSR_ip1jm1,FMSR_ip12jp12,FMSR_ip12jm12,&
                    FMSR_im12jp12,FMSR_im12jm12,FMSR_ip12j,FMSR_im12j,&
                    FMSR_ijp12,FMSR_ijm12

DOUBLE PRECISION :: FMCR_im1j,FMCR_ij,FMCR_ip1j,FMCR_im1jp1,FMCR_ijp1,FMCR_ip1jp1,&
                    FMCR_im1jm1,FMCR_ijm1,FMCR_ip1jm1,FMCR_ip12jp12,FMCR_ip12jm12,&
                    FMCR_im12jp12,FMCR_im12jm12,FMCR_ip12j,FMCR_im12j,&
                    FMCR_ijp12,FMCR_ijm12

DOUBLE PRECISION :: ETA_ijkp12,ETA_ip12jkp12,ETA_im12jkp12,ETA_ijp12kp12,&
                    ETA_ijm12kp12,ETA_ijkp1,ETA_ip12jkp1,&
                    ETA_im12jkp1,ETA_ijp12kp1,ETA_ijm12kp1

DOUBLE PRECISION :: A_kp12,B_kp12,C_kp12,D_kp1,E_kp1,F_kp1,LL
DOUBLE PRECISION :: TIME_kp12,TIME_kp1
INTEGER          :: j
! --------------------------------------------------------------------------------- !
BETA      = dt/(2.0D0*dr*rho)
CHI       = dt/(2.0D0*dz*rho)
TIME_kp12 = (DBLE(k) + 0.5D0)*dt
TIME_kp1  = (DBLE(k) + 1.0D0)*dt

M         =  0.0D0
L         =  0.0D0

rp1       =  DBLE(i)*dr 
r         =  DBLE(i-1)*dr
rm1       =  DBLE(i-2)*dr
THETA     =  dt/(r*rho)
! --------------------------------------------------------------------------------- !
! j = 1    => j = 2 in main.f90, i.e. near the bottom plate of the viscometer.
! j = NY-2 => j = NY2-1 in main.f90, i.e. near the top of the viscometer.
DO j = 1,NY-2

  V_im1jkp12    = vz_im1jkp12(j+1)
  V_ijkp12      = vz_ijkp12(j+1)
  V_ip1jkp12    = vz_ip1jkp12(j+1)
  V_im1jp1kp12  = vz_im1jkp12(j+2)
  V_ijp1kp12    = vz_ijkp12(j+2)
  V_ip1jp1kp12  = vz_ip1jkp12(j+2)
  V_im1jm1kp12  = vz_im1jkp12(j)
  V_ijm1kp12    = vz_ijkp12(j)
  V_ip1jm1kp12  = vz_ip1jkp12(j)

  V_im1jkp1     = vz_im1jkp1(j+1)
  V_ijkp1       = vz_ijkp1(j+1)
  V_ip1jkp1     = vz_ip1jkp1(j+1)
  V_im1jp1kp1   = vz_im1jkp1(j+2)
  V_ijp1kp1     = vz_ijkp1(j+2)
  V_ip1jp1kp1   = vz_ip1jkp1(j+2)
  V_im1jm1kp1   = vz_im1jkp1(j)
  V_ijm1kp1     = vz_ijkp1(j)
  V_ip1jm1kp1   = vz_ip1jkp1(j)

  FMSR_im1j     = FMSR_im1j_t(j+1)
  FMSR_ij       = FMSR_ij_t(j+1)
  FMSR_ip1j     = FMSR_ip1j_t(j+1)
  FMSR_im1jp1   = FMSR_im1j_t(j+2)
  FMSR_ijp1     = FMSR_ij_t(j+2)
  FMSR_ip1jp1   = FMSR_ip1j_t(j+2)
  FMSR_im1jm1   = FMSR_im1j_t(j)
  FMSR_ijm1     = FMSR_ij_t(j)
  FMSR_ip1jm1   = FMSR_ip1j_t(j)

  FMSR_ip12jp12 = (FMSR_ij   + FMSR_ip1j + FMSR_ip1jp1   + FMSR_ijp1)/4.0D0
  FMSR_ip12jm12 = (FMSR_ij   + FMSR_ip1j + FMSR_ip1jm1   + FMSR_ijm1)/4.0D0
  FMSR_im12jp12 = (FMSR_ij   + FMSR_im1j + FMSR_im1jp1   + FMSR_ijp1)/4.0D0
  FMSR_im12jm12 = (FMSR_ij   + FMSR_im1j + FMSR_im1jm1   + FMSR_ijm1)/4.0D0

  FMSR_ip12j    = (FMSR_ij   + FMSR_ip1j + FMSR_ip12jp12 + FMSR_ip12jm12)/4.0D0
  FMSR_im12j    = (FMSR_im1j + FMSR_ij   + FMSR_im12jp12 + FMSR_im12jm12)/4.0D0
  FMSR_ijp12    = (FMSR_ij   + FMSR_ijp1 + FMSR_ip12jp12 + FMSR_im12jp12)/4.0D0
  FMSR_ijm12    = (FMSR_ijm1 + FMSR_ij   + FMSR_ip12jm12 + FMSR_im12jm12)/4.0D0

  FMCR_im1j     = FMCR_im1j_t(j+1)
  FMCR_ij       = FMCR_ij_t(j+1)
  FMCR_ip1j     = FMCR_ip1j_t(j+1)
  FMCR_im1jp1   = FMCR_im1j_t(j+2)
  FMCR_ijp1     = FMCR_ij_t(j+2)
  FMCR_ip1jp1   = FMCR_ip1j_t(j+2)
  FMCR_im1jm1   = FMCR_im1j_t(j)
  FMCR_ijm1     = FMCR_ij_t(j)
  FMCR_ip1jm1   = FMCR_ip1j_t(j)

  FMCR_ip12jp12 = (FMCR_ij   + FMCR_ip1j + FMCR_ip1jp1   + FMCR_ijp1)/4.0D0
  FMCR_ip12jm12 = (FMCR_ij   + FMCR_ip1j + FMCR_ip1jm1   + FMCR_ijm1)/4.0D0
  FMCR_im12jp12 = (FMCR_ij   + FMCR_im1j + FMCR_im1jp1   + FMCR_ijp1)/4.0D0
  FMCR_im12jm12 = (FMCR_ij   + FMCR_im1j + FMCR_im1jm1   + FMCR_ijm1)/4.0D0

  FMCR_ip12j    = (FMCR_ij   + FMCR_ip1j + FMCR_ip12jp12 + FMCR_ip12jm12)/4.0D0
  FMCR_im12j    = (FMCR_im1j + FMCR_ij   + FMCR_im12jp12 + FMCR_im12jm12)/4.0D0
  FMCR_ijp12    = (FMCR_ij   + FMCR_ijp1 + FMCR_ip12jp12 + FMCR_im12jp12)/4.0D0
  FMCR_ijm12    = (FMCR_ijm1 + FMCR_ij   + FMCR_ip12jm12 + FMCR_im12jm12)/4.0D0

  CALL ROS(rp1,r,rm1,dr,dz,V_ijkp12,V_ip1jkp12,V_im1jkp12,V_ijp1kp12,V_ijm1kp12,&
           V_ip1jp1kp12,V_ip1jm1kp12,V_im1jp1kp12,V_im1jm1kp12,&
           SR_ijkp12,SR_ip12jkp12,SR_im12jkp12,SR_ijp12kp12,SR_ijm12kp12)

  CALL ROS(rp1,r,rm1,dr,dz,V_ijkp1,V_ip1jkp1,V_im1jkp1,V_ijp1kp1,V_ijm1kp1,&
           V_ip1jp1kp1,V_ip1jm1kp1,V_im1jp1kp1,V_im1jm1kp1,&
           SR_ijkp1,SR_ip12jkp1,SR_im12jkp1,SR_ijp12kp1,SR_ijm12kp1)

  ! FMSR is relative to k and since we are calculating ETA_kp12, the "ALPHA_I*SR*dt"
  ! term must be set equal to "ALPHA_I*SR*(dt/2)"! -> dt = dt/2.
  ! FMSR is relative to k and we are only proceeding from k to k+1/2, hence dt -> dt/2.
  ! [ The same consideration applies for FMCR, BETA_I, BETA_II and H(SR) ]. 
  CALL ETA(dt/2.0D0,TIME_kp12,Lambda,SR_ijkp12,SR_ip12jkp12,SR_im12jkp12,&
           SR_ijp12kp12,SR_ijm12kp12,FMSR_ij,FMSR_ip12j,FMSR_im12j,&
           FMSR_ijp12,FMSR_ijm12,FMCR_ij,FMCR_ip12j,FMCR_im12j,&
           FMCR_ijp12,FMCR_ijm12,ETA_ijkp12,ETA_ip12jkp12,&
           ETA_im12jkp12,ETA_ijp12kp12,ETA_ijm12kp12)

  ! FMSR is relative to k and since we are calculating ETA_kp1, the "ALPHA_I*SR*dt"
  ! term must be set equal to "ALPHA_I*SR*dt"! -> dt = dt.
  ! FMSR is relative to k and we are now proceeding from k to k+1, hence dt -> dt.
  ! [ The same consideration applies for FMCR, BETA_I, BETA_II and H(SR) ]. 
  CALL ETA(dt,TIME_kp1,Lambda,SR_ijkp1,SR_ip12jkp1,SR_im12jkp1,SR_ijp12kp1,SR_ijm12kp1,&
           FMSR_ij,FMSR_ip12j,FMSR_im12j,FMSR_ijp12,FMSR_ijm12,&
           FMCR_ij,FMCR_ip12j,FMCR_im12j,FMCR_ijp12,FMCR_ijm12,&
           ETA_ijkp1,ETA_ip12jkp1,ETA_im12jkp1,ETA_ijp12kp1,ETA_ijm12kp1)

  ! Equations 7.32 to 7.34:
  A_kp12 = (BETA*ETA_ip12jkp12)*(1.0D0/dr - 1.0D0/(rp1+r)) + (THETA*ETA_ijkp12)/(2.0D0*dr)
  B_kp12 = (BETA*ETA_ip12jkp12)*(1.0D0/dr + 1.0D0/(rp1+r)) + &
           (BETA*ETA_im12jkp12)*(1.0D0/dr - 1.0D0/(r+rm1)) + (THETA*ETA_ijkp12)/r 
  C_kp12 = (BETA*ETA_im12jkp12)*(1.0D0/dr + 1.0D0/(r+rm1)) - (THETA*ETA_ijkp12)/(2.0D0*dr)

  ! Equations 7.35 to 7.37:
  D_kp1  = (CHI*ETA_ijp12kp1)/dz
  E_kp1  = (CHI*ETA_ijp12kp1)/dz + (CHI*ETA_ijm12kp1)/dz  
  F_kp1  = (CHI*ETA_ijm12kp1)/dz
  
  ! Equation 7.31:
  LL = - A_kp12*V_ip1jkp12 - (1.0D0 - B_kp12)*V_ijkp12 - C_kp12*V_im1jkp12
  
  ! Equation 7.30:
  IF (j == 1) THEN
    M(j,j)   = -(1.0D0 + E_kp1)
    M(j,j+1) = D_kp1
    L(j)     = LL - F_kp1*V_ijm1kp1
  ELSE IF ((j.EQ.NY-2).AND.Neumann) THEN ! (Equation 7.46)
    M(j,j-1) = 2.0D0*F_kp1
    M(j,j)   = -(1.0D0 + E_kp1)
    L(j)     = LL
  ELSE IF (j == NY-2) THEN
    M(j,j-1) = F_kp1
    M(j,j)   = -(1.0D0 + E_kp1)
    L(j)     = LL - D_kp1*V_ijp1kp1
  ELSE 
    M(j,j-1) = F_kp1
    M(j,j)   = -(1.0D0 + E_kp1)
    M(j,j+1) = D_kp1
    L(j)     = LL
  END IF

END DO
! --------------------------------------------------------------------------------- !
RETURN
END SUBROUTINE MATRIX_UPDATE_Y
! ================================================================================= !
SUBROUTINE MATRIX_SOLVER(M,D,v,dim)
! Subroutine that solves the trigonal system "M * v = D" with the 
! Thomas algorithm (also known as the "Crout reduction for tridiagonal linear
! systems" - algorithm).
! M ->  Left side of the linear system (a tridiagonal array).
! D ->  Right side of the linear system (a vector).
! v ->  The variable to be solved: v = (M)^(-1) * D.

DOUBLE PRECISION,DIMENSION(:,:),INTENT(INOUT) :: M
DOUBLE PRECISION,DIMENSION(:),INTENT(INOUT)   :: D
DOUBLE PRECISION,DIMENSION(:),INTENT(OUT)     :: v
INTEGER,INTENT(IN)                            :: dim

DOUBLE PRECISION                              :: coef
INTEGER                                       :: i
! --------------------------------------------------------------------------------- !
v = -1.0D0
! --------------------------------------------------------------------------------- !
! Forward and then back substitution:
DO i = 1,dim-1
  coef       = M(i+1,(i-1)+1)/M(i,i)
  M(i+1,i+1) = M(i+1,i+1) - coef*M(i,i+1)
  D(i+1)     = D(i+1) - coef*D(i)
END DO
v(dim) = D(dim)/M(dim,dim)
DO i = dim-1,1,-1
  v(i) = (D(i) - M(i,i+1)*v(i+1))/M(i,i)
END DO
! --------------------------------------------------------------------------------- !
RETURN
END SUBROUTINE MATRIX_SOLVER
! ================================================================================= !
END MODULE MATRIX
! --------------------------------------------------------------------------------- !
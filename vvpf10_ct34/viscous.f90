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
! File name: viscous.f90 (MODULE)                                                   !
! In this file, the shear viscosity function ETA = ETA(SR,t,...) is defined and     !
! calculated. This information is requested by update.f90.                          !
! --------------------------------------------------------------------------------- !
MODULE SHEAR_VISCOSITY
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ETA,VISCOSITY,FMSR_FMCR
CONTAINS
! ================================================================================= !
SUBROUTINE ETA(dt,time,Lambda,SR_ij,SR_ip12j,SR_im12j,SR_ijp12,SR_ijm12,&
               FMSR_ij,FMSR_ip12j,FMSR_im12j,FMSR_ijp12,FMSR_ijm12,&
               FMCR_ij,FMCR_ip12j,FMCR_im12j,FMCR_ijp12,FMCR_ijm12,&
               ETA_ij,ETA_ip12j,ETA_im12j,ETA_ijp12,ETA_ijm12)
! --------------------------------------------------------------------------------- !
DOUBLE PRECISION,INTENT(IN)    :: dt,time,Lambda,&
                                  SR_ij,SR_ip12j,&
                                  SR_im12j,SR_ijp12,SR_ijm12,&
                                  FMSR_ij,FMSR_ip12j,&
                                  FMSR_im12j,FMSR_ijp12,FMSR_ijm12,&
                                  FMCR_ij,FMCR_ip12j,&
                                  FMCR_im12j,FMCR_ijp12,FMCR_ijm12

DOUBLE PRECISION,INTENT(OUT)   :: ETA_ij,ETA_ip12j,ETA_im12j,ETA_ijp12,&
                                  ETA_ijm12
! --------------------------------------------------------------------------------- !
CALL VISCOSITY(dt,time,Lambda,SR_ij,FMSR_ij,FMCR_ij,ETA_ij)
CALL VISCOSITY(dt,time,Lambda,SR_ip12j,FMSR_ip12j,FMCR_ip12j,ETA_ip12j)
CALL VISCOSITY(dt,time,Lambda,SR_im12j,FMSR_im12j,FMCR_im12j,ETA_im12j)
CALL VISCOSITY(dt,time,Lambda,SR_ijp12,FMSR_ijp12,FMCR_ijp12,ETA_ijp12)
CALL VISCOSITY(dt,time,Lambda,SR_ijm12,FMSR_ijm12,FMCR_ijm12,ETA_ijm12)
! --------------------------------------------------------------------------------- !
RETURN
END SUBROUTINE ETA
! ================================================================================= !
SUBROUTINE VISCOSITY(dt,time,Lambda_tmp,SR,FMSR,FMCR,ETA)

DOUBLE PRECISION,INTENT(IN)    :: dt,time,Lambda_tmp,SR,FMSR,FMCR
DOUBLE PRECISION,INTENT(OUT)   :: ETA
DOUBLE PRECISION               :: mu,tau,a_1,a_2,mu_tmp,tau_tmp,delta,Lambda,&
                                  B3_n3_23,U_o,U_3,memory_alpha,memory_beta,H,&
                                  Gamma,Theta,BETA_I,BETA_II,ALPHA_I,ALPHA_II                                  
! --------------------------------------------------------------------------------- !
CALL FMSR_FMCR(time,SR,U_o,memory_alpha,memory_beta,H)
! --------------------------------------------------------------------------------- !
! "Lambda >= 0" (Lambda_tmp.GT.-0.5D0) means that time independent calculations is
! present, with a constant shear viscosity ETA = constant. In this case, "Lambda" is
! related to the Continuation Method (see Section 7.8).
! "Lambda = -1" (Lambda_tmp.LT.-0.5D0) means that time dependent calculations 
! (basically thixotropic) have begun with time dependent shear viscosity ETA = ETA(r,z,t).
! --------------------------------------------------------------------------------- !
IF (Lambda_tmp.GT.-0.5D0) THEN        ! -> Time independent calculations.
  Lambda  = Lambda_tmp                !    Lambda = continuation parameter.
  Gamma   = 0.0D0
  Theta   = 0.0D0
ELSE IF (Lambda_tmp.LT.-0.5D0) THEN   ! -> Time dependent calculations.
  Lambda  = 1.0D0                     !    Lambda = dummy variable.
  ! ----
  ALPHA_I  = DEXP(time/memory_alpha)
  ALPHA_II = DEXP(-time/memory_alpha)
  Gamma    = ALPHA_II*(FMSR + ALPHA_I*SR*dt) ! Equations 7.69 to 7.70 (Equation 9.3)
  ! ----
  BETA_I   = DEXP(time/memory_beta)
  BETA_II  = DEXP(-time/memory_beta)
  Theta    = BETA_II*(FMCR + BETA_I*H*dt)    ! Equations 7.71 to 7.72 (Equation 9.4) 
  Theta    = 0.0D0  ! <- This condition is only used for the VHMW Na-case. With 
                    ! "Theta = 0.0D0", the re-coagulation is set equal to zero, as is
                    ! done in Sections 9.4 and 9.5. In Sections 9.6 to 9.8, the 
                    ! re-coagulation is non-zero and hence, Theta is calculated 
                    ! according to Equation 9.4 as "BETA_II*(FMCR + BETA_I*H*dt)".
END IF
! --------------------------------------------------------------------------------- !
delta    =  0.005D0 ! <- The regularization parameter (see Section 7.9).
! ---------------------------------------------------------
! Values from Table 9.1 are shown here (t = 72 min and t = 102 min):
mu       =  0.650D0 
a_1      =  1.100D0
B3_n3_23 = 30.000D0
tau      =  0.000D0 
a_2      =  0.800D0
! ---------------------------------------------------------
! Equation 9.5:
U_3 = (U_o*(Theta*Gamma + 1.0D0) + Theta)/((Theta + 1.0D0)*(Gamma + 1.0D0))
! ---------------------------------------------------------
mu_tmp  = mu  + a_1*(B3_n3_23*(U_3**(2.0D0/3.0D0)))   ! Equation 9.8
tau_tmp = tau + a_2*(B3_n3_23*(U_3**(2.0D0/3.0D0)))   ! Equation 9.9
ETA     = mu_tmp + (tau_tmp*Lambda)/(SR + delta)      ! Equation 9.7
! --------------------------------------------------------------------------------- !
RETURN 
END SUBROUTINE VISCOSITY
! ================================================================================= !
SUBROUTINE FMSR_FMCR(time,SR,U_o,memory_alpha,memory_beta,H)

DOUBLE PRECISION,INTENT(IN)    :: time,SR
DOUBLE PRECISION,INTENT(OUT)   :: U_o,memory_alpha,memory_beta,H
DOUBLE PRECISION               :: k_1,k_2,k_3
! --------------------------------------------------------------------------------- !
! As Theta = 0.0D0, then memory_beta is not used (Sections 9.4 and 9.5, only).
U_o          =  1.0D0
memory_alpha = 30.0D0 
memory_beta  = 18.0D0
! ---------------------------------------------------------
! As Theta = 0.0D0, then k_1, k_2 and k_3 are not used (Sections 9.4 and 9.5, only).
k_1 = 0.005D0
k_2 = 0.100D0
k_3 = 0.005D0
! ---------------------------------------------------------
! Equations 9.10 and 9.11:
IF (time.EQ.0.0D0) THEN 
  H = k_1*(1 - U_o)/4.0D0
ELSE IF ((time.GT. 0.0D0).AND.(time.LT.25.0D0)) THEN 
  H = k_1/(SR**2 + 1.0D0)
ELSE IF ((time.GE.25.0D0).AND.(time.LT.45.0D0)) THEN 
  H = k_2/(SR**2 + 1.0D0)
ELSE IF ((time.GE.45.0D0).AND.(time.LE.50.0D0)) THEN 
  H = k_3/(SR**2 + 1.0D0)
END IF
! --------------------------------------------------------------------------------- !
RETURN 
END SUBROUTINE FMSR_FMCR
! ================================================================================= !
END MODULE SHEAR_VISCOSITY
! --------------------------------------------------------------------------------- !
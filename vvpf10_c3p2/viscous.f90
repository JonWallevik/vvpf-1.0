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
! In this file, the shear viscosity function ETA = ETA(SR) is defined and           !
! calculated. This information is requested by update.f90.                          !
! --------------------------------------------------------------------------------- !
MODULE SHEAR_VISCOSITY
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ETA,VISCOSITY
CONTAINS
! ================================================================================= !
SUBROUTINE ETA(TIME,Lambda,ROS_ij,ROS_ip12j,ROS_im12j,ROS_ijp12,&
               ROS_ijm12,ETA_ij,ETA_ip12j,ETA_im12j,ETA_ijp12,ETA_ijm12)

DOUBLE PRECISION,INTENT(IN)    :: TIME,Lambda,ROS_ij,ROS_ip12j,ROS_im12j,&
                                  ROS_ijp12,ROS_ijm12

DOUBLE PRECISION,INTENT(OUT)   :: ETA_ij,ETA_ip12j,ETA_im12j,ETA_ijp12,&
                                  ETA_ijm12
! --------------------------------------------------------------------------------- !
CALL VISCOSITY(TIME,Lambda,ROS_ij,ETA_ij)
CALL VISCOSITY(TIME,Lambda,ROS_ip12j,ETA_ip12j)
CALL VISCOSITY(TIME,Lambda,ROS_im12j,ETA_im12j)
CALL VISCOSITY(TIME,Lambda,ROS_ijp12,ETA_ijp12)
CALL VISCOSITY(TIME,Lambda,ROS_ijm12,ETA_ijm12)
! --------------------------------------------------------------------------------- !
RETURN
END SUBROUTINE ETA
! ================================================================================= !
SUBROUTINE VISCOSITY(TIME_tmp,Lambda_tmp,ROS,ETA)

DOUBLE PRECISION,INTENT(IN)    :: TIME_tmp,Lambda_tmp,ROS
DOUBLE PRECISION,INTENT(OUT)   :: ETA
DOUBLE PRECISION               :: mu,tau,delta,Lambda
! --------------------------------------------------------------------------------- !
! Lambda => Continuation Method (see Section 7.8).
Lambda = Lambda_tmp
! ---------------------------------------------------------
mu     =  20.0D0  ! Plastic viscosity [Pa.s].
tau    = 150.0D0  ! Yield value [Pa].
delta  =  2.0D-3  ! <- The regularization parameter (see Section 7.9).
! ---------------------------------------------------------
ETA = mu + (tau*Lambda)/(ROS + delta)
! --------------------------------------------------------------------------------- !
RETURN 
END SUBROUTINE VISCOSITY
! ================================================================================= !
END MODULE SHEAR_VISCOSITY
! --------------------------------------------------------------------------------- !
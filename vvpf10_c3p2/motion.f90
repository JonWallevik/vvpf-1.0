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
! The information about the (constant) angular velocity "omega" is requested by the !
! routine main.f90.                                                                 !
! --------------------------------------------------------------------------------- !
MODULE ROTATION
  IMPLICIT NONE
  PRIVATE
  PUBLIC ::  ANGULAR_VELOCITY
CONTAINS
! ================================================================================= !
SUBROUTINE ANGULAR_VELOCITY(double_prec_k,dt,omega)

DOUBLE PRECISION,INTENT(IN)     :: double_prec_k,dt
DOUBLE PRECISION,INTENT(OUT)    :: omega
DOUBLE PRECISION                :: PI,f
! --------------------------------------------------------------------------------- !
! The terms dt and double_prec_k are a legacy from Appendix A.2.
PI    = DACOS(-1.0D0)
f     = 0.5D0
omega = 2*PI*f ! rad/s (for example, omega = 3.0D0 rad/s) 
! --------------------------------------------------------------------------------- !
END SUBROUTINE ANGULAR_VELOCITY
! ================================================================================= !
END MODULE ROTATION
! --------------------------------------------------------------------------------- !
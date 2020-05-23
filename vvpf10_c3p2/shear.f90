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
! File name: shear.f90 (MODULE) [is the same to what is shown in Appendix A.2.6]    !
! This routine calculates the shear rate SR(i,j) from the computed velocity profile !
! VELOCITY_k(i,j). It is the program update.f90 that makes the request.             !
! See Section 7.5 about the formulas for the shear rate (SR). Note that ROS and SR  !
! means the same thing: ROS = rate of shear = SR = shear rate.                      !
! --------------------------------------------------------------------------------- !
MODULE SHEAR_RATE
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: SR
CONTAINS
! ================================================================================= !
SUBROUTINE SR(rp1,r,rm1,dr,dz,V_ij,V_ip1j,V_im1j,V_ijp1,V_ijm1,V_ip1jp1,V_ip1jm1,&
              V_im1jp1,V_im1jm1,SR_ij,SR_ip12j,SR_im12j,SR_ijp12,SR_ijm12)

DOUBLE PRECISION,INTENT(IN)  ::  rp1,r,rm1,dr,dz,V_ij,&
                                 V_ip1j,    V_im1j,    V_ijp1,    V_ijm1,&
                                 V_ip1jp1,  V_ip1jm1,  V_im1jp1,  V_im1jm1
DOUBLE PRECISION,INTENT(OUT) ::  SR_ij,     SR_ip12j,  SR_im12j,  SR_ijp12, SR_ijm12
DOUBLE PRECISION             ::  SR1_ij,    SR2_ij,    SR1_ip12j, SR2_ip12j,&
                                 SR1_im12j, SR2_im12j, SR1_ijp12,&
                                 SR2_ijp12, SR1_ijm12, SR2_ijm12
! --------------------------------------------------------------------------------- !
SR1_ij    = (V_ip1j - V_im1j)/(2.0D0*dr) - V_ij/r
SR2_ij    = (V_ijp1 - V_ijm1)/(2.0D0*dz)
SR_ij     = DSQRT(SR1_ij**2.0D0 + SR2_ij**2.0D0)

SR1_ip12j = (V_ip1j - V_ij)/dr - (V_ip1j + V_ij)/(rp1 + r)
SR2_ip12j = (V_ip1jp1 + V_ijp1 - V_ip1jm1 - V_ijm1)/(4.0D0*dz)
SR_ip12j  = DSQRT(SR1_ip12j**2.0D0 + SR2_ip12j**2.0D0)

SR1_im12j = (V_ij - V_im1j)/dr - (V_ij + V_im1j)/(r + rm1)
SR2_im12j = (V_ijp1 + V_im1jp1 - V_ijm1 - V_im1jm1)/(4.0D0*dz)
SR_im12j  = DSQRT(SR1_im12j**2.0D0 + SR2_im12j**2.0D0)

SR1_ijp12 = (V_ip1jp1 + V_ip1j - V_im1jp1 - &
              V_im1j)/(4.0D0*dr) - (V_ijp1 + V_ij)/(2.0D0*r)
SR2_ijp12 = (V_ijp1 - V_ij)/dz 
SR_ijp12  = DSQRT(SR1_ijp12**2.0D0 + SR2_ijp12**2.0D0)

SR1_ijm12 = (V_ip1j + V_ip1jm1 - V_im1j - & 
              V_im1jm1)/(4.0D0*dr) - (V_ij + V_ijm1)/(2.0D0*r)
SR2_ijm12 = (V_ij - V_ijm1)/dz 
SR_ijm12  = DSQRT(SR1_ijm12**2.0D0 + SR2_ijm12**2.0D0)
! --------------------------------------------------------------------------------- !
RETURN
END SUBROUTINE SR
! ================================================================================= !
END MODULE SHEAR_RATE
! --------------------------------------------------------------------------------- !
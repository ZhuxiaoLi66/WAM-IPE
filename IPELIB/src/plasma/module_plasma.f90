!dbg20120501: add v// to perp transport
!20110911: note: openmp was tried on jet but did not work: only thread 0 was used not the other thread...although other threads did exist...needs more investigation...
! DATE: 08 September, 2011
!********************************************
!***      Copyright 2011 NAOMI MARUYAMA   ***
!***      ALL RIGHTS RESERVED             ***
!********************************************
! LICENSE AGREEMENT Ionosphere Plasmasphere Electrodynamics (IPE) model
! DEVELOPER: Dr. Naomi Maruyama
! CONTACT INFORMATION:
! E-MAIL : Naomi.Maruyama@noaa.gov
! PHONE  : 303-497-4857
! ADDRESS: 325 Broadway, Boulder, CO 80305
!--------------------------------------------  
      MODULE module_PLASMA

      USE module_precision
      USE module_IPE_dimension,ONLY: IPDIM,ISTOT
      IMPLICIT NONE
      include "gptl.inc"
      REAL(KIND=real_prec8),DIMENSION(ISTOT,IPDIM),PUBLIC :: plasma_1d
      INTEGER (KIND=int_prec),PUBLIC:: utime_save

      END MODULE module_PLASMA

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
! --- PRIVATE ---
!
! --- PUBLIC ---
!solver REALL needs these parameters from the previous time step
      !.. TE_TI_k(3,J) = Te, TE_TI_k(2,J) = Ti = TE_TI_k(2,J) [kelvin]
!dbg20120501      TYPE :: plasma_data_1d
!dbg20120501        REAL(KIND=real_prec), DIMENSION(      IPDIM) :: Te_k  !TE_TI(3)
!dbg20120501        REAL(KIND=real_prec), DIMENSION(ISPET,IPDIM) :: Ti_k  !TE_TI(1:2)
!dbg20120501        REAL(KIND=real_prec), DIMENSION(ISPEC,IPDIM) :: N_m3
!dbg20110927        REAL(KIND=real_prec), DIMENSION(ISPEV,IPDIM) :: V_ms1
!dbg20110927        REAL(KIND=real_prec), DIMENSION(      IPDIM) :: heating_rate_e_cgs ![eV cm-3 s-1] !EHT(3)
!dbg20110927        REAL(KIND=real_prec), DIMENSION(ISPET,IPDIM) :: heating_rate_i_cgs  !EHT(1:2)
!dbg20110923        REAL(KIND=real_prec), DIMENSION(      IPDIM) :: NO_m3 
!dbg20120501      END TYPE  plasma_data_1d
!      TYPE(plasma_data_1d),ALLOCATABLE,TARGET,PUBLIC :: plasma_3d(:,:)

!neutral needs these parameters from the ionosphere in addition to N&T
!dbg20120501      TYPE :: plasma_data_1d4n
!dbg20120501        REAL(KIND=real_prec), DIMENSION(ISPEV) :: V_ms1
!???      REAL(KIND=real_prec) :: NHEAT
!dbg20120501      END TYPE  plasma_data_1d4n
!dbg20120501      TYPE(plasma_data_1d4n),ALLOCATABLE,TARGET,PUBLIC :: plasma_3d4n(:,:) !(NPTS2D, NMP)

!dbg20120501      TYPE(plasma_data_1d), PUBLIC :: n0_1d !N&T after perpendicular transport
!dbg20120501
!     REAL(KIND=real_prec),ALLOCATABLE,TARGET,PUBLIC :: plasma_3d(:,:,:,:)!ISTOT,NPTS2D,NLP,NMP
      REAL(KIND=real_prec),DIMENSION(ISTOT,IPDIM),PUBLIC :: plasma_1d
!1:9:den:o+,h+,he+,n+
!10 :te
!11:12:ti
!13:16:vi:o+,h+,he+,n+


!only for debug, o+
!d      REAL(KIND=real_prec), DIMENSION(NPTS2D), PUBLIC :: n0_2dbg 
!      REAL(KIND=real_prec), DIMENSION(3,NPTS2D,NMP), PUBLIC :: TE_TI_k
!      REAL(KIND=real_prec), DIMENSION(ISPEC,NPTS2D,NMP), PUBLIC :: XIONN_m3
!      REAL(KIND=real_prec), DIMENSION(ISPEC,NPTS2D,NMP), PUBLIC :: XIONV_ms1

! save ut so that other subroutines can refer to it
      INTEGER (KIND=int_prec),PUBLIC:: utime_save

!nm20121003:subroutine plasma is separated into module_sub_plasma.f90
      END MODULE module_PLASMA

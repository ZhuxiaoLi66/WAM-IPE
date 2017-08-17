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
!NOTE: dbg20120228: only temporary used to estimate the error created by the perp transport for debug purpose
      subroutine dbg_estimate_trans_error ( utime )
        USE module_precision
        USE module_IPE_dimension,ONLY:NMP,NLP,IPDIM
        USE module_FIELD_LINE_GRID_MKS,ONLY: plasma_3d
        USE module_input_parameters,ONLY : parallelBuild
        IMPLICIT NONE
        INTEGER (KIND=int_prec),  INTENT(IN)  :: utime    !universal time [sec]

!---local variables---
        INTEGER (KIND=int_prec) :: lp,mp
!!SMS$DISTRIBUTE(dh,NLP,NMP) BEGIN
        REAL (KIND=real_prec),dimension(NMP,NLP) :: opmaxpc,opminpc
!!SMS$DISTRIBUTE END
        REAL (KIND=real_prec),PARAMETER :: fixedval=100.0

if(parallelBuild) then
  print*,'sw_dbg_perp_trans=true does not work for parallel runs'
  print*,'Stopping in dbg_estimate_trans_error'
  stop
endif

do lp=1,NLP
do mp=1,NMP
opmaxpc(mp,lp)=( &
     MAXVAL( plasma_3d(mp,lp)%N_m3( 1,1:IPDIM) )-fixedval &
              & )/fixedval
opminpc(mp,lp)=( &
     MINVAL( plasma_3d(mp,lp)%N_m3( 1,1:IPDIM) )-fixedval &
              & )/fixedval
enddo
enddo

write(5000,fmt=*) utime,MAXVAL(opmaxpc),MINVAL(opminpc)


END       subroutine dbg_estimate_trans_error 

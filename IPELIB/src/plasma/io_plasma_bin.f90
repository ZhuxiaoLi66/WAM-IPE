!CAUTION!!!!!!!the plasma i-o orders/file names have been changed !!!

!CAUTION!!!!!!!the plasma i-o orders/file names have been changed !!!
!date:Thu Sep 29 18:58:05 GMT 2011
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
SUBROUTINE io_plasma_bin ( switch, utime, timestamp_for_IPE )
USE module_precision
USE module_IO,ONLY: LUN_PLASMA1,LUN_PLASMA2,lun_min1,lun_min2,lun_ut,record_number_plasma,lun_max1

USE module_FIELD_LINE_GRID_MKS,ONLY: JMIN_IN,JMAX_IS,plasma_3d,JMIN_ING,JMAX_ISG,VEXBup, &
                                     vn_ms1_4output,tn_k,on_m3,n2n_m3,o2n_m3, maxFluxTube

USE module_IPE_dimension,ONLY: NMP,NLP,NPTS2D,ISPEC,ISPEV,IPDIM,ISPET,ISTOT
USE module_input_parameters,ONLY:record_number_plasma_start,mype, &
                                 sw_record_number,stop_time,start_time, &
                                 duration,mpstop, sw_output_wind, sw_use_wam_fields_for_restart

USE module_physical_constants,ONLY:zero
! ghgm - now need the open_file module....
USE module_open_file,ONLY: open_file

!USE module_myIPE_Init, only: model_start_time, model_stop_time

IMPLICIT NONE

INTEGER (KIND=int_prec ),INTENT(IN) :: switch !2:read; 1:write
INTEGER (KIND=int_prec ),INTENT(IN) :: utime !universal time [sec]
CHARACTER(len=12), INTENT(IN)       :: timestamp_for_IPE
INTEGER (KIND=int_prec ),PARAMETER  :: n_max=10000
INTEGER (KIND=int_prec )            :: stat_alloc
INTEGER (KIND=int_prec )            :: jth,mp,lp,npts
INTEGER (KIND=int_prec )            :: lun,in,is
INTEGER (KIND=int_prec )            :: n_read,n_read_min, utime_dum,record_number_plasma_dum
INTEGER (KIND=int_prec )            :: n_count
INTEGER (KIND=int_prec )            :: ipts !dbg20120501
INTEGER                             :: iErr
REAL    (KIND=real_prec)            :: dumm(NPTS2D,NMP)
CHARACTER(len=200)                  :: restart_directory




IF ( switch==1 ) THEN !1:Output the 16 plasma* files


!SMS$SERIAL(<plasma_3d, tn_k, vn_ms1_4output, on_m3, n2n_m3, o2n_m3, IN>:default=ignore) BEGIN

  record_number_plasma = record_number_plasma+1

  OPEN( UNIT = 5999, &
        FILE = 'ipe_grid_plasma_params.'//timestamp_for_IPE, &
        FORM = 'UNFORMATTED', &
        STATUS = 'REPLACE', &
        IOSTAT = iErr )
  IF( iErr /= 0 )THEN
    PRINT*, 'sub-io_plasma_bin : Error Opening file ipe_grid_plasma_params.'//timestamp_for_IPE//' for writing.'
    STOP
  ENDIF
  

  j_loop1: DO jth=1,ISTOT !=(ISPEC+3+ISPEV)
  
    mp_loop1:DO mp=1,mpstop
      lp_loop1:DO lp=1,nlp
      
        IN=JMIN_IN(lp)
        IS=JMAX_IS(lp)
        npts = IS-IN+1 
        dumm(JMIN_ING(lp):JMAX_ISG(lp),mp) = plasma_3d(IN:IS,lp,mp,jth)
        
      ENDDO lp_loop1!lp
    ENDDO mp_loop1!mp

    ! ghgm - also write the 16 plasma files to a single unit.....
    write (unit=5999) (dumm(:,mp),mp=1,mpstop)

  ENDDO j_loop1
  
  CLOSE(unit=5999)


    OPEN( UNIT = 5998, &
        FILE = 'ipe_grid_neutral_params.'//timestamp_for_IPE, &
        FORM = 'UNFORMATTED', &
        STATUS = 'REPLACE', &
        IOSTAT = iErr )
    IF( iErr /= 0 )THEN
      PRINT*, 'sub-io_plasma_bin : Error Opening file ipe_grid_neutral_params.'//timestamp_for_IPE
      STOP
    ENDIF

    WRITE (unit=5998) tn_k(1:MaxFluxTube,1:NLP,1:NMP)
    WRITE (unit=5998) vn_ms1_4output(1:MaxFluxTube,1:NLP,1:NMP,1:3)
    WRITE (unit=5998) on_m3(1:MaxFluxTube,1:NLP,1:NMP) 
    WRITE (unit=5998) n2n_m3(1:MaxFluxTube,1:NLP,1:NMP)
    WRITE (unit=5998) o2n_m3(1:MaxFluxTube,1:NLP,1:NMP)
    CLOSE (unit=5998)


!SMS$SERIAL END


ELSE IF ( switch==2 ) THEN !2:RESTART: 

   CALL getenv("RESDIR", restart_directory)
   print *, 'GHGM IO_PLASMA restart_directory=',TRIM(restart_directory)

!SMS$SERIAL (<plasma_3d, tn_k, vn_ms1_4output, on_m3, n2n_m3, o2n_m3, OUT> : DEFAULT=IGNORE) BEGIN

  OPEN( UNIT = 5997, &
        FILE = trim(restart_directory)//'ipe_grid_plasma_params', &
        FORM = 'UNFORMATTED', &
        STATUS = 'OLD', &
        IOSTAT = iErr )
  IF( iErr /= 0 )THEN
    PRINT*, 'sub-io_plasma_bin : Error Opening file '//trim(restart_directory)//'ipe_grid_plasma_params'
    STOP
  ENDIF
    
  j_loop3: DO jth=1,(ISPEC+3)
    read (unit=5997) dumm
    mp_loop3:do mp=1,NMP
      lp_loop3:do lp=1,NLP
        IN   = JMIN_IN(lp)
        IS   = JMAX_IS(lp)
        npts = IS-IN+1 
        plasma_3d(IN:IS,lp,mp,jth) = dumm(JMIN_ING(lp):JMAX_ISG(lp),mp) 
      end do lp_loop3
    end do mp_loop3
  END DO j_loop3

  CLOSE (5997)

  OPEN( UNIT = 5996, &
        FILE = trim(restart_directory)//'ipe_grid_neutral_params', &
        FORM = 'UNFORMATTED', &
        STATUS = 'OLD', &
        IOSTAT = iErr )
  IF( iErr /= 0 )THEN
    PRINT*, 'sub-io_plasma_bin : Error Opening file='//trim(restart_directory)//'ipe_grid_neutral_params'
    STOP
  ENDIF


print*, mype,' shape=',shape(tn_k)

  READ (unit=5996) tn_k(1:MaxFluxTube,1:NLP,1:NMP)

print*,mype,'reading tn_k finished'
  READ (unit=5996) vn_ms1_4output(1:MaxFluxTube,1:NLP,1:NMP,1:3)

print*,mype,'reading vn_ms1_4 finished'

  READ (unit=5996) on_m3(1:MaxFluxTube,1:NLP,1:NMP) 
print*,mype,'reading on_m3 finished'

  READ (unit=5996) n2n_m3(1:MaxFluxTube,1:NLP,1:NMP) 

print*,mype,'reading n2n_m3 finished'

  READ (unit=5996) o2n_m3(1:MaxFluxTube,1:NLP,1:NMP) 
print*,mype,'reading o2n_m3 finished'

  CLOSE(5996)

!SMS$SERIAL END


END IF !( switch==1 ) THEN

END SUBROUTINE io_plasma_bin

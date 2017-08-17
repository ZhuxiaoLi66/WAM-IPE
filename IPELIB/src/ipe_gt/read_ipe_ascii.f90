program read_ipe_ascii
implicit none
INTEGER, PARAMETER :: NPTS2D =44438
INTEGER, PARAMETER :: NMP =80  !# of magnetic longitude sectors
INTEGER, PARAMETER :: NLP =170 !# of field line grids in one longitude sector
INTEGER, PARAMETER :: real_prec  = 4
REAL (KIND=real_prec),DIMENSION(NPTS2D,NMP) :: dumm
integer :: n_read,lun,lun1
integer,PARAMETER :: rec_start=5 !utime=14400
integer,PARAMETER :: rec_stop=7  !utime=16200
integer:: NMP_all,NLP_all,NPTS2D_dum
integer,PARAMETER:: NMP0=1
integer,PARAMETER:: NMP1=80
integer,dimension(NLP)::IN,IS
REAL (KIND=real_prec),DIMENSION(NPTS2D) :: Z,GL
REAL (KIND=real_prec),DIMENSION(NPTS2D,NMP) :: GCOLAT,GLON
REAL (KIND=real_prec),DIMENSION(1:NMP+1) :: mlon_rad
REAL(KIND=real_prec), PARAMETER :: pi = 3.1415926536
CHARACTER(LEN=*),PARAMETER:: which_format='ascii' !'bin'
CHARACTER(LEN=*),PARAMETER:: flnm0='plasma_grid'
CHARACTER(LEN=*),PARAMETER:: flnm1='stup'
CHARACTER(LEN=100):: flnm2
CHARACTER (LEN=100) :: string_tmp
INTEGER,PARAMETER::sw_read_grid=.true.
INTEGER,PARAMETER::sw_read_ipe=.true.
LOGICAL,PARAMETER::sw_output2asc=.false.
INTEGER,PARAMETER::lun_min=0
INTEGER,PARAMETER::lun_max=12
CHARACTER(LEN=*),PARAMETER:: FMT_asc="(20E12.4)"

! --- read grid
IF ( sw_read_grid ) THEN
!grid
lun=100
lun1=101
IF ( which_format=='bin') THEN
  OPEN(UNIT=lun,FILE=flnm0,STATUS='old',FORM='unformatted')
ELSE IF ( which_format=='ascii') THEN
OPEN(UNIT=lun,FILE=flnm0//'.'//which_format,STATUS='old',FORM='formatted')
END IF

IF ( which_format=='bin') THEN
    READ(UNIT=lun) NMP_all
    READ(UNIT=lun) NLP_all
    READ(UNIT=lun) NPTS2D_dum
    READ(UNIT=lun) IN(1:NLP_all)
    READ(UNIT=lun) IS(1:NLP_all)
    READ(UNIT=lun) mlon_rad( 1: NMP_all+1 ) !rad
    READ(UNIT=lun) Z(     1:NPTS2D_dum) !meter
    READ(UNIT=lun) GL(    1:NPTS2D_dum) !rad
    READ(UNIT=lun) GCOLAT(1:NPTS2D_dum, NMP0:NMP1)  !rad
    READ(UNIT=lun) GLON(  1:NPTS2D_dum, NMP0:NMP1) !rad
ELSE IF ( which_format=='ascii') THEN

    READ(UNIT=lun,FMT="(i10)") NMP_all
    READ(UNIT=lun,FMT="(i10)") NLP_all
    READ(UNIT=lun,FMT="(i10)") NPTS2D_dum
    READ(UNIT=lun,FMT="(20i10)") IN(1:NLP_all)
    READ(UNIT=lun,FMT="(20i10)") IS(1:NLP_all)
    READ(UNIT=lun,FMT="(20E12.4)") mlon_rad( 1: NMP_all+1 ) !rad
    READ(UNIT=lun,FMT="(20E12.4)") Z(     1:NPTS2D_dum) !meter
    READ(UNIT=lun,FMT="(20E12.4)") GL(    1:NPTS2D_dum) !rad
    READ(UNIT=lun,FMT="(20E12.4)") GCOLAT(1:NPTS2D_dum, NMP0:NMP1)  !rad
    READ(UNIT=lun,FMT="(20E12.4)") GLON(  1:NPTS2D_dum, NMP0:NMP1) !rad

END IF
    CLOSE(UNIT=lun)
    print *,'read plasma_grid finished'

!dbg20111003:
print *,'NMP_all',NMP_all
print *,'NLP_all',NLP_all
print *,'NPTS2D_dum',NPTS2D_dum

print *,'IN',IN(2)
print *,'IS',IS(1)

print *,'mlon[deg]',mlon_rad(2)*180./pi
 print *,'Z[meter]',z(1118)
!print *,'Z[meter]',z(1:1200)
print *,'GL[deg]',GL(1118)*180./pi
print *,'GLAT[deg]',90.-GCOLAT(1118,2)*180./pi
print *,'GLON[deg]',GLON(1118,2)*180./pi
END IF !( sw_read_grid ) THEN

! --- read ipe output
IF ( sw_read_ipe ) THEN

lun_loop: DO lun=lun_min,lun_max
!note:
! 0: O+ density [m-3]
! 1: H+ density
! 2: He+ density
! 3: N+ density
! 4: NO+ density
! 5: O2+ density
! 6: N2+ density
! 7: O+(2D) density
! 8: O+(2P) density

! Ne(electron density)= sum of all[0-8] the above ion densities

! 9: Te: electron temperature [k]
!10: Ti: ion temperature [k]
!11: (not used...)
!12: Vo+: O+ field aligned velocity [m/s]


           IF ( lun < 10 ) THEN
              WRITE( string_tmp, FMT="('0',i1)" ) lun
           ELSE IF ( lun < 100 ) THEN
              WRITE( string_tmp, FMT="(i2)" ) lun
           END IF
           flnm2 =flnm1//TRIM(string_tmp)
           print *,lun,'filename',flnm2




IF ( which_format=='bin') THEN
  OPEN(UNIT=lun,FILE=flnm2,STATUS='old',FORM='unformatted')
ELSE IF ( which_format=='ascii') THEN
  OPEN(UNIT=lun,FILE=TRIM(flnm2)//'.'//which_format,STATUS='old',FORM='formatted')
END IF
lun1=lun+10
IF ( sw_output2asc )  OPEN(UNIT=lun1,FILE=TRIM(flnm2)//'.ascii',STATUS='unknown',FORM='formatted')

read_loop: DO n_read=rec_start, rec_stop

IF ( which_format=='bin') THEN
  READ (UNIT=lun ) dumm
ELSE IF ( which_format=='ascii') THEN
  READ (UNIT=lun,FMT=FMT_asc ) dumm
END IF
print *,'n_read=',n_read,'!dbg! read dumm finished'
print *,lun,'th species ',maxval(dumm),minval(dumm)


IF ( sw_output2asc )  WRITE (UNIT=lun1,FMT=FMT_asc ) dumm

end do read_loop
CLOSE(lun)

IF ( sw_output2asc )  CLOSE(lun1)

end do lun_loop !: DO i=lun_min,lun_max
END IF !( sw_read_ipe ) THEN


print *,'end program read_ipe_ascii'
end program read_ipe_ascii

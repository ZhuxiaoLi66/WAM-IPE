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
      MODULE module_open_output_files
      USE module_precision
      USE module_IPE_dimension,ONLY: ISPEC,ISPEV
      IMPLICIT NONE

!nm20121003:module parameters are separated into module_io.f90!

      PRIVATE
      PUBLIC :: open_output_files

      CONTAINS
!---------------------------
        SUBROUTINE open_output_files ( )
        USE module_input_parameters,ONLY: NYEAR,NDAY,HPEQ_flip,sw_output_plasma_grid &
           ,record_number_plasma_start,sw_output_fort167,sw_output_wind,mype,sw_use_wam_fields_for_restart
        USE module_IO,ONLY: &
  filename,FORM_dum,STATUS_dum &
, LUN_pgrid,LUN_LOG &
, LUN_FLIP1,LUN_FLIP2,LUN_FLIP3,LUN_FLIP4 &
, LUN_PLASMA0, LUN_PLASMA1,LUN_PLASMA2, LUN_UT &
, lun_min1,lun_max1,lun_min2,lun_max2 &
, record_number_plasma,luntmp1,luntmp2,luntmp3

        USE module_open_file,ONLY: open_file

        IMPLICIT NONE
        CHARACTER (LEN=100) :: string_tmp
        INTEGER (KIND=int_prec)::i
!
        IF ( NDAY < 100 ) THEN
          WRITE ( string_tmp, FMT="(i4,A2,i2)" ) NYEAR,'_0',NDAY
        ELSE         ! NDAY>=100
          WRITE ( string_tmp, FMT="(i4,A1,i3)" ) NYEAR,'_' ,NDAY
        END IF
        WRITE( UNIT=LUN_LOG, FMT=*) string_tmp


!--- unit=9
!nm20120303        filename ='logfile'//TRIM(string_tmp)//'.log'
        filename ='input_par'
        FORM_dum ='formatted  ' 
        STATUS_dum ='unknown'
        CALL open_file ( filename, LUN_LOG, FORM_dum, STATUS_dum )  

        IF ( sw_output_fort167 ) THEN
!SMS$IGNORE begin
           !--- unit=167
           LUN_FLIP1=167
           WRITE ( filename, FMT="('fort.',i3)" ) LUN_FLIP1  !'fort.167'
           if(mype==0)print *,'fort.167?', filename,'unit_number=',LUN_FLIP1
           FORM_dum ='formatted  ' 
           STATUS_dum ='unknown'
CALL open_file ( filename, LUN_FLIP1, FORM_dum, STATUS_dum ) 

           !--- unit=168
           LUN_FLIP2=168
           WRITE ( filename, FMT="('fort.',i3)" ) LUN_FLIP2  !'fort.168'
           if(mype==0)print *,'fort.168?', filename,'unit_number=',LUN_FLIP2
CALL open_file ( filename, LUN_FLIP2, FORM_dum, STATUS_dum )

           !--- unit=170
           LUN_FLIP3=170
           WRITE ( filename, FMT="('fort.',i3)" ) LUN_FLIP3  !'fort.170'
           if(mype==0)print *,'fort.170?', filename,'unit_number=',LUN_FLIP3
CALL open_file ( filename, LUN_FLIP3, FORM_dum, STATUS_dum )

           !--- unit=171
           LUN_FLIP4=171
           WRITE ( filename, FMT="('fort.',i3)" ) LUN_FLIP4  !'fort.171'
           if(mype==0)print *,'fort.171?', filename,'unit_number=',LUN_FLIP4
CALL open_file ( filename, LUN_FLIP4, FORM_dum, STATUS_dum )
print*,mype,'check unit#',LUN_FLIP1,LUN_FLIP3,LUN_FLIP2,LUN_FLIP4
!SMS$IGNORE end
        END IF !( sw_output_fort167


IF ( sw_output_plasma_grid ) THEN
        LUN_PLASMA0=98
        FORM_dum = 'unformatted' 
        filename = 'plasma_grid'
        CALL open_file ( filename, LUN_PLASMA0, FORM_dum, STATUS_dum )
END IF !( sw_output_plasma_grid ) THEN

!nm20110923        LUN_PLASMA2=99
!nm20110923        filename = 'plasma_startup0'
!nm20110923        CALL open_file ( filename, LUN_PLASMA2, FORM_dum, STATUS_dum )

        LUN_UT=lun_min1-1 !=99
!nm20120303        filename ='ut_rec.log'
!        filename ='ut_rec'
!        FORM_dum ='formatted  ' 
!        STATUS_dum ='unknown'
!        CALL open_file ( filename, LUN_UT, FORM_dum, STATUS_dum ) 



        END SUBROUTINE open_output_files
!---------------------------
END MODULE module_open_output_files


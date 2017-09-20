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
        USE module_input_parameters,ONLY: NYEAR,NDAY,HPEQ_flip,sw_debug,sw_output_plasma_grid &
&, record_number_plasma_start,sw_output_fort167,sw_output_wind,mype,peFort167 &
&, sw_use_wam_fields_for_restart
        USE module_IO,ONLY: &
  filename,FORM_dum,STATUS_dum &
, LUN_pgrid,LUN_LOG &
, LUN_FLIP1,LUN_FLIP2,LUN_FLIP3,LUN_FLIP4 &
, LUN_PLASMA0, LUN_PLASMA1,LUN_PLASMA2, LUN_UT, LUN_UT2 &
, lun_min1,lun_max1,lun_min2,lun_max2 &
, record_number_plasma,luntmp1,luntmp2,luntmp3 &
, lun_ipe_grid_neut_params_ut &
, lun_ipe_grid_neut_wind &
, lun_ipe_grid_neut_temp &
, lun_ipe_grid_neut_O_den &
, lun_ipe_grid_neut_N2_den &
, lun_ipe_grid_neut_O2_den &
, LUN_WAM_RESTART0,LUN_WAM_RESTART1,lun_wam_tn &
, LUN_WAM_RESTART3,LUN_WAM_RESTART4,LUN_WAM_RESTART5
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
        filename ='ut_rec'
        FORM_dum ='formatted  ' 
        STATUS_dum ='unknown'
        CALL open_file ( filename, LUN_UT, FORM_dum, STATUS_dum ) 


!--- plasma output: unit=100~~115
        FORM_dum = 'unformatted' 
        DO i=lun_min1,lun_max1
           if(sw_debug) print *,'plasma: unit=',i        
           LUN_PLASMA1(i)=i
           IF ( (i-lun_min1) < 10 ) THEN
              WRITE( string_tmp, FMT="('0',i1)" )(i-lun_min1)
           ELSE IF ( (i-lun_min1) < 100 ) THEN
              WRITE( string_tmp, FMT="(i2)" )(i-lun_min1)
           END IF
           filename ='plasma'//TRIM(string_tmp)
           if(sw_debug) print *,(i-lun_min1),'filename',filename
           CALL open_file ( filename, LUN_PLASMA1(i), FORM_dum, STATUS_dum )
        END DO
        record_number_plasma = record_number_plasma_start - 1

        IF ( HPEQ_flip==0.0 ) THEN
!--- unit=1181 : input history file
!          LUN_PLASMA12=1181
!          STATUS_dum ='old'
!          filename = 'plasma_startup1'
!          CALL open_file ( filename, LUN_PLASMA12, FORM_dum, STATUS_dum )

          LUN_UT2=lun_min2-1 !=199
!nm20120303          filename ='startup_ut_rec.log'
          filename ='stup_ut_rec'
          FORM_dum ='formatted  ' 
          STATUS_dum ='old'
          CALL open_file ( filename, LUN_UT2, FORM_dum, STATUS_dum )

!--- unit=200~~215
          FORM_dum = 'unformatted' 
!note: reading Vi will be needed for neutral coupling
          DO i=lun_min2,lun_max2   !t +ISPEV)
             !if(sw_debug) 
             print *,'startup: unit=',i        
             LUN_PLASMA2(i)=i
             IF ( (i-lun_min2) < 10 ) THEN
                WRITE( string_tmp, FMT="('0',i1)" )(i-lun_min2)
             ELSE IF ( (i-lun_min2) < 100 ) THEN
                WRITE( string_tmp, FMT="(i2)" )(i-lun_min2)
             END IF
!nm20120303             filename ='startup'//TRIM(string_tmp)
             filename ='stup'//TRIM(string_tmp)
             !if(sw_debug) 
             print *,(i-lun_min2),'filename',filename
             CALL open_file ( filename, LUN_PLASMA2(i), FORM_dum, STATUS_dum )
          END DO
        END IF  ! ( HPEQ_flip==0.0 ) THEN


!nm20141001 wind output
        IF ( sw_output_wind ) THEN
!--- unit=6000 ut for wind
           lun_ipe_grid_neut_params_ut = 6000
           filename='ipe_grid_neut_params_ut'
           print *,'fort.6000? ', filename, lun_ipe_grid_neut_params_ut
           FORM_dum ='formatted  ' 
           STATUS_dum ='unknown'
           CALL open_file ( filename, lun_ipe_grid_neut_params_ut, FORM_dum, STATUS_dum )

!--- unit=6001 wind
           lun_ipe_grid_neut_wind = 6001
           filename='ipe_grid_neut_wind'
           print *,'fort.6001? ', filename, lun_ipe_grid_neut_wind
           FORM_dum ='unformatted' 
           STATUS_dum ='unknown'

!--- unit=6002 tn
           lun_ipe_grid_neut_temp = 6002
           filename='ipe_grid_neut_temp'
           print *,'fort.6002? ', filename, lun_ipe_grid_neut_temp
           FORM_dum ='unformatted' 
           STATUS_dum ='unknown'
           CALL open_file ( filename, lun_ipe_grid_neut_temp, FORM_dum, STATUS_dum ) 

!--- unit=6003 on: neutral atomic oxygen density
           lun_ipe_grid_neut_O_den = 6003
           filename='ipe_grid_neut_O_den'
           print *,'fort.6003? ', filename, lun_ipe_grid_neut_O_den
           FORM_dum ='unformatted' 
           STATUS_dum ='unknown'
           CALL open_file ( filename, lun_ipe_grid_neut_O_den, FORM_dum, STATUS_dum ) 

!--- unit=6004 on: neutral molecular nitrogen density
           lun_ipe_grid_neut_N2_den = 6004
           filename='ipe_grid_neut_N2_den'
           print *,'fort.6004? ', filename, lun_ipe_grid_neut_N2_den
           FORM_dum ='unformatted' 
           STATUS_dum ='unknown'
           CALL open_file ( filename, lun_ipe_grid_neut_N2_den, FORM_dum, STATUS_dum ) 

!--- unit=6005 on: neutral molecular oxygen density
           lun_ipe_grid_neut_O2_den = 6005
           filename='ipe_grid_neut_O2_den'
           print *,'fort.6004? ', filename, lun_ipe_grid_neut_O2_den
           FORM_dum ='unformatted' 
           STATUS_dum ='unknown'
           CALL open_file ( filename, lun_ipe_grid_neut_O2_den, FORM_dum, STATUS_dum ) 

        END IF !( sw_output_wind

        IF ( sw_use_wam_fields_for_restart ) THEN
!--- unit=5000 ut for wind
!          LUN_WAM_RESTART0=5000
!          filename='ut_out4wind'
!          print *,'fort.5000? ', filename, LUN_WAM_RESTART0
!          FORM_dum ='formatted  ' 
!          STATUS_dum ='old'
!          CALL open_file ( filename, LUN_WAM_RESTART0, FORM_dum, STATUS_dum )

!--- unit=5001 wind
           LUN_WAM_RESTART1=15001
           filename='/scratch3/NCEPDEV/swpc/noscrub/George.Millward/wam-ipe_new_svn/wam_data/ipe_grid_neut_wind'
           print *,'fort.15001? ', filename, LUN_WAM_RESTART1
           FORM_dum ='unformatted' 
           STATUS_dum ='old'
           CALL open_file ( filename, LUN_WAM_RESTART1, FORM_dum, STATUS_dum ) 

!--- unit=5002 tn
           lun_wam_tn=15002
           filename='/scratch3/NCEPDEV/swpc/noscrub/George.Millward/wam-ipe_new_svn/wam_data/ipe_grid_neut_temp'
           print *,'fort.15002? ', filename, lun_wam_tn
           FORM_dum ='unformatted' 
           STATUS_dum ='old'
           CALL open_file ( filename, lun_wam_tn, FORM_dum, STATUS_dum ) 

!--- unit=5003 on: neutral atomic oxygen density
           LUN_WAM_RESTART3=15003
           filename='/scratch3/NCEPDEV/swpc/noscrub/George.Millward/wam-ipe_new_svn/wam_data/ipe_grid_neut_O_den'
           print *,'fort.15003? ', filename, LUN_WAM_RESTART3
           FORM_dum ='unformatted' 
           STATUS_dum ='old'
           CALL open_file ( filename, LUN_WAM_RESTART3, FORM_dum, STATUS_dum ) 

!--- unit=5004 on: neutral molecular nitrogen density
           LUN_WAM_RESTART4=15004
           filename='/scratch3/NCEPDEV/swpc/noscrub/George.Millward/wam-ipe_new_svn/wam_data/ipe_grid_neut_N2_den'
           print *,'fort.15004? ', filename, LUN_WAM_RESTART4
           FORM_dum ='unformatted' 
           STATUS_dum ='old'
           CALL open_file ( filename, LUN_WAM_RESTART4, FORM_dum, STATUS_dum ) 

!--- unit=5005 on: neutral molecular oxygen density
           LUN_WAM_RESTART5=15005
           filename='/scratch3/NCEPDEV/swpc/noscrub/George.Millward/wam-ipe_new_svn/wam_data/ipe_grid_neut_O2_den'
           print *,'fort.15005? ', filename, LUN_WAM_RESTART5
           FORM_dum ='unformatted' 
           STATUS_dum ='old'
           CALL open_file ( filename, LUN_WAM_RESTART5, FORM_dum, STATUS_dum ) 

        END IF !( sw_use_wam_fields_for_restart

        END SUBROUTINE open_output_files
!---------------------------
END MODULE module_open_output_files


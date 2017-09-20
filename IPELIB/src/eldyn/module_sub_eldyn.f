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
!Aug2011: the original code was provided from Fei Wu from WAM version
!nm20110906: modified to implement to IPE
!copied from /.../testall.f
!ylonm(1:nmlon=180)
!ylatm(1:nmlat=90)
!      program ts_efield
      MODULE module_sub_eldyn
      USE module_precision
!----------------------
!c idea
!      subroutine idea_geteb(im,ix,dayno,utsec,f107,kp,maglat,maglon,
!     &essa,ee1,ee2)
      USE efield !,ONLY:iday,imo,iday_m,iyear,ut,kp,by,bz,f107d
      USE module_get_efield,ONLY:get_efield
!c     use date_def
!c     use physcons, pi => con_pi
      IMPLICIT NONE

!nm20121003:module parameters are separated into module_eldyn.f

      PRIVATE
      PUBLIC :: eldyn
      CONTAINS
!---
      SUBROUTINE eldyn ( utime )
      USE module_precision
      USE module_cal_monthday
      USE module_input_parameters,ONLY:NYEAR,NDAY,start_time,mype       &
     &,ip_freq_output,sw_debug,F107D_ipe => F107D                       &
     &,lpi,bz_eld,by_eld,kp_eld,sw_bnd_wei,bnd_wei_eld,lat_sft_eld      &
     &,sw_ctip_input,utime0LPI
      USE module_physical_constants,ONLY:rtd
!nm20121003:
      USE module_eldyn,ONLY:theta90_rad,j0,Ed1_90,Ed2_90,j1,coslam_m
      use efield, only:nmlat,ylatm,bnd_wei,lat_sft,ilat_sft,ef_max
      USE module_IPE_dimension,ONLY: NLP
      IMPLICIT NONE
      INTEGER (KIND=int_prec),INTENT(IN)   :: utime !universal time [sec]
!---local
      real :: kp ! 
      integer :: j
!c initiate
!c calculate efield only if diff time step
!      if(utsec.ne.utsec_last) then
!        utsec_last=utsec
      iday = NDAY !254 ! dayno                   ! day of year
!c       imo=idate(2)
!c       iday_m=idate(3) 
!!!need to create a routine to calculate month/day from NDAY!!!
!      imo=3                     !month
      iyear = NYEAR 
!nm20121127: calculate month/day from iyear and iday
!dbg20170916
      print*,'sub-eld:iyear',iyear,iday,utime
      call cal_monthday ( iyear,iday, imo,iday_m )
!nm20130402: temporarily hard-code the iday to get b4bconfirmed.
!      iday_m=15                 !day of month 

!!! F107D is global both in module efield & ipe input
      f107d = F107D_ipe         !f107
      ut = REAL(utime,real_prec)/3600.0
      if(ut>=24.) ut=MOD(ut,24.)
      kp = kp_eld  !=1.                   !???
!      bz = .433726 - kp*(.0849999*kp + .0810363)                        &
!     &        + f107d*(.00793738 - .00219316*kp)
 
!dbg20170916
      print*,'sub-eld: utime=', utime,' start_time=', start_time
      if ( utime==start_time ) then
        print *,'iday',iday, 'imo',imo,' iday_m',iday_m,' iyear',iyear
        print *,' kp',kp
        print *,'IPE By=',by,' Bz=',bz,' F107d=',f107d
      end if

!nm20151104
      if ( sw_ctip_input ) then
        LPI = INT( ( utime - utime0LPI ) / 60. ) +1
!t        if(sw_debug)
        print*,'sub-eld: LPI=',lpi
!t        if(sw_debug)
        print*,'sub-eld: utime',utime,'dt_m=',((utime-utime0LPI)/60.)
      else
        LPI=1
      end if
      bz = bz_eld(LPI)
      by = by_eld(LPI)

      if(sw_debug)print *,'Bz=',bz,' By=',by

!nm20151105
!-------------------------------------------------------------------
! find latitudinal shift
!-------------------------------------------------------------------
      if ( sw_bnd_wei==1 ) then 
         bnd_wei = bnd_wei_eld
         lat_sft = lat_sft_eld
         do j = 0,nmlat
            ilat_sft = j
            if( lat_sft <= ylatm(j) ) then
               exit
            end if
         end do 
         print *,'bnd_wei=',bnd_wei,' lat_sft=',lat_sft                    &
     &        ,' ilat_sft=',ilat_sft,' ef_max=',ef_max
      end if                    !( sw_bnd_wei==1 ) then 
!------------------------------------------------------------------

      call get_efield
!dbg        print*,'www'
!        ed11=ed1
!        ed22=ed2
      IF ( utime==start_time ) THEN
          write(unit=2003,FMT='(20f10.4)')ylatm
          write(unit=2004,FMT='(20f10.4)')ylonm
      END IF
!      endif !if(utsec.ne.utsec_last) then


! get ED1/2(nmp=80 X nlp=170) at 90km from potent(181x91)at 130km
      IF ( utime==start_time ) then
        j0=-999 !missing_value in module_find_nei...
      endif
      CALL GET_EFIELD90km ( utime )
      if ( sw_debug )  print *,'GET_EFIELD90km finished'
      IF ( utime==start_time ) THEN 
        write(unit=2007,FMT='(20f10.4)') (90.-theta90_rad*rtd)    
      ENDIF

      write(unit=1002,FMT='(i8,3F8.2)')utime,(maxval(potent)*1.0E-03)   & !vmax [kV]
     &,(minval(potent)*1.0E-03)                                         & !vmin
     &,(( maxval(potent)-minval(potent) )*1.0E-03)  !cpcp
      IF ( MOD( (utime-start_time),ip_freq_output)==0 ) THEN
!SMS$SERIAL(<ed1_90,ed2_90,IN>:default=ignore) BEGIN
        write(unit=2000,FMT='(20E12.4)')potent !V
        write(unit=2001,FMT='(20E12.4)')ed1 *1.0E+03 !V/m-->mV/m
        write(unit=2002,FMT='(20E12.4)')ed2 *1.0E+03 !V/m-->mV/m
        write(unit=2008,FMT='(20E12.4)')ed1_90 *1.0E+03 !V/m-->mV/m
        write(unit=2009,FMT='(20E12.4)')ed2_90 *1.0E+03 !V/m-->mV/m
        write(unit=2010,FMT='(I12)')utime !sec
!SMS$SERIAL END
      END IF ! ( mod( (utime
!c
      return
!      end
!      end program ts_efield
!
      END SUBROUTINE eldyn
      END MODULE module_sub_eldyn

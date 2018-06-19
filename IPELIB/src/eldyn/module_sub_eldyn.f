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
      USE efield_ipe !,ONLY:iday,imo,iday_m,iyear,ut,kp,by,bz,f107d
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
!SMS$IGNORE BEGIN
#ifdef HAVE_MPI
      USE module_input_parameters,ONLY:NYEAR,NDAY,start_time,mype       &
     &,ip_freq_output,lpi,kp_eld,sw_bnd_wei,bnd_wei_eld        &
     &,lat_sft_eld,sw_ctip_input,utime0LPI,f107_new,f107d_new           &
     &,input_params_begin,input_params_interval
#else
      USE module_input_parameters,ONLY:NYEAR,NDAY,start_time,mype,      &
     & ip_freq_output,lpi,kp_eld,sw_bnd_wei,bnd_wei_eld,       &
     & lat_sft_eld,sw_ctip_input,utime0LPI,f107_new,f107d_new,          &
     & input_params_begin,input_params_interval
#endif
!SMS$IGNORE END
      USE module_physical_constants,ONLY:rtd
!nm20121003:
      USE module_eldyn,ONLY:theta90_rad,j0,Ed1_90,Ed2_90,j1,coslam_m
      USE efield_ipe, only:nmlat,ylatm,bnd_wei,lat_sft,ilat_sft,ef_max
      USE module_IPE_dimension,ONLY: NLP
      IMPLICIT NONE
      INTEGER (KIND=int_prec),INTENT(IN)   :: utime !universal time [sec]
!---local
      real :: kp ! 
      integer :: j
      iday = NDAY !254 ! dayno                   ! day of year
      iyear = NYEAR 
      call cal_monthday ( iyear,iday, imo,iday_m )

      if ( sw_ctip_input ) then
        LPI = INT( ( utime - utime0LPI ) / real(input_params_interval) ) &
     &+ 1 + input_params_begin
      else
        LPI=1
      end if

!!! F107D is global both in module efield & ipe input
      f107d = F107_new(LPI)         !f107
      ut = REAL(utime,real_prec)/3600.0
      if(ut>=24.) ut=MOD(ut,24.)
      kp = kp_eld(LPI)  !=1.                   !???
!      bz = .433726 - kp*(.0849999*kp + .0810363)                        &
!     &        + f107d*(.00793738 - .00219316*kp)
 
!-------------------------------------------------------------------
! find latitudinal shift
!-------------------------------------------------------------------
      if ( sw_bnd_wei ) then 
         bnd_wei = bnd_wei_eld
         lat_sft = lat_sft_eld
         do j = 0,nmlat
            ilat_sft = j
            if( lat_sft <= ylatm(j) ) then
               exit
            end if
         end do 
      end if                    !( sw_bnd_wei==1 ) then 
!------------------------------------------------------------------

      call get_efield


! get ED1/2(nmp=80 X nlp=170) at 90km from potent(181x91)at 130km
      IF ( utime==start_time ) then
        j0=-999 !missing_value in module_find_nei...
      endif
      CALL GET_EFIELD90km ( utime )

      return
!
      END SUBROUTINE eldyn
      END MODULE module_sub_eldyn

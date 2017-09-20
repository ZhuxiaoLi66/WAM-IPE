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
!20120829: copied from magfield.f and splitted

!
      module module_sunloc
!
! Sub sunloc is called once per timestep from advance to determine
!   sun's longitudes for current ut model time.
! sunlons: sun's longitude in dipole coordinates (see sub sunloc)
! (this was dlons in earlier versions)
!
      PUBLIC:: sunloc
!
      contains
!-----------------------------------------------------------------------
      subroutine sunloc(iyr,iday,secs)
      USE module_precision
      USE module_physical_constants,ONLY:dtr      
      USE module_magfield,ONLY:sunlons
!nm20111003:      use cons_module,only: dtr  ! degrees to radians (pi/180)
! am 10/04
! tiegcm uses the sun's longitude in dipole coordinates
! we changed the approximation of the sun's location in geographic
! coodinates from
!      glats=asin(.398749*sin(2.*PI*(iday-80)/365.))
!      glons=pi*(1.-2.*secs/86400.)
! to use the apex routines based on formulas in Astronomical Almanac
! difference is around 6/7 min 
! This is called every timestep from advance.
!
! Args:
      integer (KIND=int_prec),intent(in) :: iyr,  ! year
     |   iday ! day of year
      real (KIND=real_prec),intent(in) :: secs    ! ut in seconds
!
! Local:
      integer :: ihr,imn
      real :: sec,date,vp,xmlon ! apex magnetic longitude
     |  sbsllat,    ! geographic latitude of subsolar point (degrees)
     |  sbsllon,    ! geographic longitude of subsolar point (degrees)
     |  colat,      ! Geocentric colatitude of geomagnetic dipole north pole (deg)
     |  elon        ! East longitude of geomagnetic dipole north pole (deg)
      
      ihr = int(secs/3600.)
      imn = int((secs - float(ihr)*3600.)/60.)
      sec = secs - float(ihr)*3600. - float(imn)*60.
      
!  calculate subsol point: given universal time
!          input: iyr,iday,ihr,imn,sec
!          output: sbsllat,sbsllon 
!                  
      call subsol(iyr,iday,ihr,imn,sec ,sbsllat,sbsllon)
      
      date = float(iyr) + float(iday)/365. + float(ihr)/24./365. +
     |  float(imn)/60./24./365.+ sec/60./60./24./365.
      call cofrm(date)
      call dypol(colat,elon,vp)
      
! calculate geomagn. diploe longitude
!        input: aloni,sbsllat,sbsllon,colat,elon
!        output: xmlon  
      call solgmlon(sbsllat,sbsllon,colat,elon,xmlon) 
!
      sunlons(1) = xmlon*dtr
!      sunlons(1) = 1.247029293721108 ! for cpmapring with TIEGCM code
!     
      end subroutine sunloc
!-----------------------------------------------------------------------
      end module module_sunloc

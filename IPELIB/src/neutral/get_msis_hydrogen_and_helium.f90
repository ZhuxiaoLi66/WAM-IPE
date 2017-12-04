      subroutine get_msis_hydrogen (npts, glon_deg, glat_deg, alt_km, &                                
                               iyear, iday, ut_hour, f107D_dum, f107A_dum, AP_dum, &
                               h_density_m3, he_density_m3)

      USE module_precision
      USE module_unit_conversion,ONLY: M3_TO_CM3
      implicit none

      integer, intent(in) :: npts, iyear, iday
      real(8), intent(in) :: ut_hour, f107D_dum, f107A_dum, AP_dum(7)
      real(8), dimension(npts), intent(in) :: glon_deg, glat_deg, alt_km

      REAL(KIND=real_prec), dimension(npts), intent(out) :: &
                        h_density_m3, he_density_m3

      real(4) :: sec, alt, glat, glon, stl, f107a_msis, f107d_msis
      real(4) :: ap_msis(7), d(9), t(2)
      real(4) :: ap_hwm(2)

      integer :: iyd, mass

      integer :: i

      ap_msis(:) = AP_dum(1) ! (daily) magnetic index
      ap_hwm (:) = AP_dum(1) ! (daily) magnetic index
      mass       = 48 ! mass number is calculated for all
      f107a_msis = f107A_dum  ! 81 day average of f10.7 flux (centered on day ddd) 
      f107d_msis = f107D_dum  ! daily f10.7 flux for previous day 

      iyd = 99000 + iday
      sec = ut_hour * 3600. 

      do i = 1, npts 

      glon = glon_deg(i)
      glat = glat_deg(i)
      alt  = alt_km  (i)

      stl  = sec/3600. + glon/15.
!dbg      if ( stl>=24.) stl=stl-24.

      call gtd7(iyd,sec,alt,glat,glon,stl,f107a_msis,f107d_msis,ap_msis,mass,d,t)

! convert from cm-3 to m-3
      h_density_m3 (i) = d(1)/M3_TO_CM3 !*1.e6
      he_density_m3  (i) = d(7)/M3_TO_CM3

      enddo

      end subroutine get_msis_hydrogen

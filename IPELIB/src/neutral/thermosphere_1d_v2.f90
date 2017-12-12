      subroutine get_thermosphere (npts, &
                               iyear, iday, ut_hour, f107D_dum, f107A_dum, AP_dum, &
                               glon_deg, glat_deg, alt_km, &
                               he_density_m3, o_density_m3, &
                               o2_density_m3, n2_density_m3, &
                               h_density_m3, n4s_density_m3, tn_k_dum, tinf_k_dum, &
                               vn_ms1_dum)

      USE module_precision
      USE module_unit_conversion,ONLY: M3_TO_CM3
      implicit none

      integer, intent(in) :: npts

      integer, intent(in) :: iyear, iday

      real(8), intent(in) :: ut_hour, f107D_dum, f107A_dum, AP_dum(7)

      real(8), dimension(npts), intent(in) ::  &
                         glon_deg, glat_deg, alt_km

      REAL(KIND=real_prec), dimension(npts), intent(out) ::  &
                         he_density_m3, &         
                         o_density_m3, &         
                         o2_density_m3, &       
                         n2_density_m3, &      
                         h_density_m3, &
                         n4s_density_m3, &    
                         tn_k_dum, &                
                         tinf_k_dum
      REAL(KIND=real_prec), dimension(3,npts), intent(out) ::  vn_ms1_dum

      real(4) :: sec, alt, glat, glon, stl, f107a_msis, f107d_msis
      real(4) :: ap_msis(7), d(9), t(2)
      real(4) :: ap_hwm(2), w(2)

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
!dbg20110923
!dbg      if ( stl>=24.) stl=stl-24.

      ! horizontal & vertical wind

      w(:)=0.0000000
!dbg20110729
!dbg      print *,'call gws5',iyd,sec,alt,glat,glon,stl,f107a_msis,f107d_msis,ap_hwm(1),w
      call gws5(iyd,sec,alt,glat,glon,stl,f107a_msis,f107d_msis,ap_hwm(1),w)
!dbg      print *,'(2)w=',w

      vn_ms1_dum(1,i) =   w(2)  !eastward
      vn_ms1_dum(2,i) =   w(1)  !northward
      vn_ms1_dum(3,i) =   0.0   !upward

      ! composition & temperature

!dbg20110923
!dbg print *,'call gtd7',iyd,sec,alt,glat,glon,stl,f107a_msis,f107d_msis,ap_msis,mass
      call gtd7(iyd,sec,alt,glat,glon,stl,f107a_msis,f107d_msis,ap_msis,mass,d,t)
!dbg print *,'msis output',d,t

! convert from cm-3 to m-3
      he_density_m3 (i) = d(1)/M3_TO_CM3 !*1.e6
      o_density_m3  (i) = d(2)/M3_TO_CM3
      n2_density_m3 (i) = d(3)/M3_TO_CM3
      o2_density_m3 (i) = d(4)/M3_TO_CM3
      h_density_m3  (i) = d(7)/M3_TO_CM3
      n4s_density_m3(i) = d(8)/M3_TO_CM3

      tinf_k_dum(i)         = t(1) !exospheric temperature
      tn_k_dum(i)           = t(2) !temperature at alt 

      enddo

      end subroutine get_thermosphere

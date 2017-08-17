MODULE modSizeFixedGridThermo

IMPLICIT NONE

PRIVATE  ! Set everything to private access, except those
         ! marked as public

PUBLIC :: nFixedGridThermoHeights , &
          !nFixedGridThermoLats, &   not used - gt lats used instead
          !nFixedGridThermoLons, &   not used - gt lons used instead
          fixedThermoHeights_km

SAVE

!----------------------------------
! sizes of Thermospheric fixed grid 
!----------------------------------
integer, parameter :: nFixedGridThermoHeights = 31

! lrm This was fixed_heights_km(interface_hts)
REAL(kind=8) :: fixedThermoHeights_km(nFixedGridThermoHeights)

DATA fixedThermoHeights_km /90., 95., 100., 105., 110., 115., 120., 125., &
         150., 175., 200., 225., 250., 275., 300., 325., 350., 375., 400.,  &
         450., 500., 550., 600., 700., 800., 900., 1000., &
         2000., 4000., 6370., 9000./

END MODULE modSizeFixedGridThermo

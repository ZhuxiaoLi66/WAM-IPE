MODULE modSizeFixedGridIono

IMPLICIT NONE

PRIVATE  ! Set everything to private access, except those
         ! marked as public

PUBLIC :: nFixedGridIonoHeights, &
          nFixedGridIonoLats, &
          nFixedGridIonoLons

SAVE

!----------------------------------
! sizes of Ionospheric fixed grid 
!----------------------------------
integer, parameter :: nFixedGridIonoHeights = 183
integer, parameter :: nFixedGridIonoLats = 91
integer, parameter :: nFixedGridIonoLons = 90

END MODULE modSizeFixedGridIono

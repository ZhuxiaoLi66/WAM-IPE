MODULE modSizeThermo

! module for sizes of the thermospheric pressure grid

IMPLICIT NONE

PRIVATE  ! Set everything to private access, except those
         ! marked as public

PUBLIC :: GT_ht_dim, &
          GT_lat_dim, &
          GT_lon_dim


SAVE

!--------------------------------------
! sizes of Thermospheric pressure grid 
!--------------------------------------
integer, parameter :: GT_ht_dim = 15
integer, parameter :: GT_lat_dim = 91
integer, parameter :: GT_lon_dim = 20

 

END MODULE modSizeThermo

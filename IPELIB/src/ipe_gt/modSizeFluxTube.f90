MODULE modSizeFluxTube

IMPLICIT NONE

PRIVATE  ! Set everything to private access, except those
         ! marked as public

PUBLIC :: NPTS, NMP, NLP

SAVE

!----------------------------------
! sizes of flux tube grid
!----------------------------------
! changed to fit new plasma grid 20120717lrm
INTEGER, PARAMETER :: NPTS = 44514 ! Total number of gridpoints along flux tubes
INTEGER, PARAMETER :: NMP  = 80  ! number of longitude sectors
INTEGER, PARAMETER :: NLP  = 170  ! number of tubes (N. pole to equator)

END MODULE modSizeFluxTube

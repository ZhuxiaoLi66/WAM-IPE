program driver_grid

USE module_FIELD_LINE_GRID_MKS,ONLY: init_plasma_grid
IMPLICIT NONE

! set up plasma grids by reading file
CALL init_plasma_grid ( )

end program driver_grid

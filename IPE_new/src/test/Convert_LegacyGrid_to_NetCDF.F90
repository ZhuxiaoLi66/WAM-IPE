PROGRAM Convert_LegacyGrid_to_NetCDF

USE IPE_Grid_Class
USE IPE_Model_Parameters_Class


IMPLICIT NONE

  TYPE( IPE_Grid )             :: grid
  TYPE( IPE_Model_Parameters ) :: params 
  
  CHARACTER(500) :: gridFile
  LOGICAL :: initSuccess, readSuccess


    CALL InitializeFromCommandLine( gridFile, initSuccess )

    IF( initSuccess )THEN

      CALL params % Build( readSuccess )

      IF( readSuccess )THEN

        

      ELSE

        PRINT*, '  Re-run Legacy2NetCDF with the generated IPE.inp file' 

      ENDIF
 
    ENDIF

CONTAINS

 SUBROUTINE InitializeFromCommandLine( gridFile, initSuccess )
   IMPLICIT NONE
   CHARACTER(*), INTENT(out) :: gridFile
   LOGICAL, INTENT(out)      :: initSuccess
   ! Local
   INTEGER        :: nArg, argID
   CHARACTER(500) :: argname
   LOGICAL        :: fileExists
   LOGICAL        :: gridGiven, gridFileGiven, helpRequested

     gridFileGiven = .FALSE.
     gridGiven     = .FALSE.
     helpRequested = .FALSE.
     initSuccess   = .FALSE.

     ! Default grid file
     gridFile      = './ipe_grid'

     nArg = command_argument_count( )

     IF( nArg > 0 )THEN
 
       DO argID = 1, nArg
  
         CALL get_command_argument( argID, argName )

         SELECT CASE( TRIM(argName) )

           CASE("--legacy-grid")

             gridGiven = .TRUE.

           CASE("--help")

             helpRequested = .TRUE.

           CASE DEFAULT

             IF( gridGiven )THEN

               gridFileGiven = .TRUE.
               gridFile  = TRIM(argName)

             ENDIF

         END SELECT
       ENDDO

     ENDIF

     IF( helpRequested )THEN


       PRINT*, '  legacy2netcdf ' 
       PRINT*, ' '
       PRINT*, '  Usage : legacy2netcdf [OPTIONS] '
       PRINT*, ' '
       PRINT*, '    Options '
       PRINT*, ' '
       PRINT*, '      --help '
       PRINT*, '          Displays this message '
       PRINT*, ' '
       PRINT*, '      --legacy-grid <grid-file>'
       PRINT*, '          Specifies to use <grid-file> as the legacy grid file '
       PRINT*, '          to be converted to NetCDF. If this option is not  '
       PRINT*, '          provide, ./ipe_grid is assumed.  '
       PRINT*, ' '

     ENDIF

     IF( .NOT. helpRequested )THEN 

       INQUIRE( FILE = TRIM(gridFile), &
                EXIST = fileExists )

       IF( .NOT.(fileExists) )THEN
         PRINT*, '  Grid file :'//TRIM(gridFile)//': not found.'
         PRINT*, '  Setup Failed.'
         initSuccess = .FALSE.
       ELSE
         initSuccess = .TRUE.
       ENDIF

     ENDIF

 END SUBROUTINE InitializeFromCommandLine

END PROGRAM Convert_LegacyGrid_to_NetCDF

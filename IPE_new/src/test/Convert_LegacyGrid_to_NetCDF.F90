PROGRAM Convert_LegacyGrid_to_NetCDF

USE IPE_Grid_Class
USE IPE_Model_Parameters_Class


IMPLICIT NONE

  TYPE( IPE_Grid )             :: grid
  TYPE( IPE_Model_Parameters ) :: params 
  
  CHARACTER(500) :: gridFile, ncFile
  LOGICAL :: initSuccess, readSuccess


    CALL InitializeFromCommandLine( gridFile, ncFile, initSuccess )

    IF( initSuccess )THEN

      CALL params % Build( readSuccess )

      IF( readSuccess )THEN

        CALL grid % Build( params % nFluxTube, params % NLP, params % NMP, params % NPTS2D )
       
        CALL grid % Read_IPE_Grid( gridFile )        

        CALL grid % Write_IPE_Grid_NetCDF( ncFile )

        CALL grid % Trash( )

      ELSE

        PRINT*, ' '
        PRINT*, '  Re-run Legacy2NetCDF with the generated IPE.inp file' 
        PRINT*, ' '

      ENDIF
 
    ENDIF

CONTAINS

 SUBROUTINE InitializeFromCommandLine( gridFile, outFile, initSuccess )
   IMPLICIT NONE
   CHARACTER(*), INTENT(out) :: gridFile, outFile
   LOGICAL, INTENT(out)      :: initSuccess
   ! Local
   INTEGER        :: nArg, argID
   CHARACTER(500) :: argname
   LOGICAL        :: fileExists
   LOGICAL        :: gridGiven, gridFileGiven, helpRequested
   LOGICAL        :: outGiven, outFileGiven

     gridFileGiven = .FALSE.
     gridGiven     = .FALSE.
     outGiven      = .FALSE.
     outFileGiven  = .FALSE.
     helpRequested = .FALSE.
     initSuccess   = .FALSE.

     ! Default grid file
     gridFile      = './ipe_grid'
     outFile       = './IPE_Grid.nc'

     nArg = command_argument_count( )

     IF( nArg > 0 )THEN
 
       DO argID = 1, nArg
  
         CALL get_command_argument( argID, argName )

         SELECT CASE( TRIM(argName) )

           CASE("--legacy-grid")

             gridGiven = .TRUE.

           CASE("--nc-file")

             outGiven = .TRUE.

           CASE("--help")

             helpRequested = .TRUE.

           CASE DEFAULT

             IF( gridGiven )THEN

               gridFileGiven = .TRUE.
               gridFile  = TRIM(argName)
               gridGiven = .FALSE.

             ENDIF

             IF( outGiven )THEN

               outFileGiven = .TRUE.
               outFile      = TRIM(argName)
               outGiven     = .FALSE.

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
       PRINT*, '          provided, ./ipe_grid is assumed.  '
       PRINT*, ' '
       PRINT*, '      --nc-file <output-nc-file>'
       PRINT*, '          Specifies the name of the output netcdf file to use.'
       PRINT*, '          If this option is not provided, IPE_Grid.nc is '
       PRINT*, '          assumed. '
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
         PRINT*, ' '
         PRINT*, ' Using Legacy File :', TRIM(gridFile)
         PRINT*, ' '
         initSuccess = .TRUE.
       ENDIF

     ENDIF

 END SUBROUTINE InitializeFromCommandLine

END PROGRAM Convert_LegacyGrid_to_NetCDF

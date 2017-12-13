! IPEtoHeightGrid.f90
!
!  Authors : George Millward  (george.millward@noaa.gov)
!            Joe Schoonover   (joseph.schoonover@noaa.gov)
!
!  Cooperative Institute for Research in Environmental Sciences (CIRES)
!  National Oceanographic and Atmospheric Administration (NOAA)
!  Space Weather Prediction Center (SWPC)
!
!


PROGRAM IPEtoHeightGrid

USE module_precision
USE module_IPE_dimension,ONLY: NMP, NLP, NPTS2D, ISTOT
USE module_input_parameters
USE module_FIELD_LINE_GRID_MKS,ONLY: MaxFluxTube, JMIN_IN, JMAX_IS, JMIN_ING, JMAX_ISG
USE module_init_plasma_grid


  IMPLICIT NONE




  REAL(kind=real_prec8)                                 :: dtotinv, factor, d, fac
  REAL(kind=real_prec8),DIMENSION(:,:,:,:), ALLOCATABLE :: facfac_interface,dd_interface
  REAL(kind=real_prec),DIMENSION(:,:), ALLOCATABLE      :: oplus_ft
  REAL(kind=real_prec),DIMENSION(:,:), ALLOCATABLE      :: hplus_ft
  REAL(kind=real_prec),DIMENSION(:,:), ALLOCATABLE      :: heplus_ft
  REAL(kind=real_prec),DIMENSION(:,:), ALLOCATABLE      :: nplus_ft
  REAL(kind=real_prec),DIMENSION(:,:), ALLOCATABLE      :: noplus_ft
  REAL(kind=real_prec),DIMENSION(:,:), ALLOCATABLE      :: o2plus_ft
  REAL(kind=real_prec),DIMENSION(:,:), ALLOCATABLE      :: n2plus_ft
  REAL(kind=real_prec),DIMENSION(:,:), ALLOCATABLE      :: oplus2d_ft
  REAL(kind=real_prec),DIMENSION(:,:), ALLOCATABLE      :: oplus2p_ft
  REAL(kind=real_prec),DIMENSION(:,:), ALLOCATABLE      :: tn_ft
  REAL(kind=real_prec),DIMENSION(:,:,:), ALLOCATABLE    :: vn_ft
  REAL(kind=real_prec),DIMENSION(:,:), ALLOCATABLE      :: on_ft
  REAL(kind=real_prec),DIMENSION(:,:), ALLOCATABLE      :: n2n_ft
  REAL(kind=real_prec),DIMENSION(:,:), ALLOCATABLE      :: o2n_ft
  REAL(kind=real_prec8),DIMENSION(:,:,:), ALLOCATABLE   :: oplus_fixed

  CHARACTER(100) :: filename
  CHARACTER(500) :: plasmaFile
  CHARACTER(500) :: neutralFile1
  CHARACTER(500) :: neutralFile2
  CHARACTER(500) :: gridFile
  CHARACTER(500) :: outputDir
  

  LOGICAL :: plasmaFileGiven, neutralFileGiven, outputDirGiven, gridFileGiven


    CALL InitializeFromCommandLine( )

    CALL read_input_parameters ( )

    PRINT*, "AllocateArrays"
    CALL AllocateArrays( )

    PRINT*, "ReadGridAndPlasmaFiles"
    CALL ReadGridAndPlasmaFiles( )
 
!   PRINT*, "InterpolateOntoFixedHeightGrid"
!   CALL InterpolateOntoFixedHeightGrid( )

!   PRINT*,"WriteToNetCDF"
!   CALL WriteToNetCDF( )

    PRINT*,"CleanupArrays"
    CALL CleanupArrays( )



CONTAINS
 SUBROUTINE InitializeFromCommandLine( )
   IMPLICIT NONE

   INTEGER        :: nArg, argID
   CHARACTER(500) :: argname
   LOGICAL        :: fileExists
   LOGICAL        :: plasmaGiven, neutralsGiven,  gridGiven, outputGiven

   CHARACTER(500) :: argv(3)
   

!    plasmaFileGiven   = .FALSE. 
!    gridFileGiven     = .FALSE.
!    neutralFileGiven = .FALSE.
!    outputDirGiven    = .FALSE.

!    plasmaGiven    = .FALSE.
!    neutralsGiven  = .FALSE.
!    gridGiven      = .FALSE.
!    outputGiven    = .FALSE.

!    nArg = command_argument_count( )

     IF( .true. ) THEN
!    IF( nArg > 0 )THEN
!
!       DO argID = 1, nArg
! 
!         CALL get_command_argument( argID, argName )

!         SELECT CASE( TRIM(argName) )

!            CASE("--plasma-file")
!               plasmaGiven = .TRUE.
!               plasmaFileGiven = .TRUE.
!            CASE("--neutral-file1")
!               neutralFileGiven = .TRUE.
!               neutralFile1  = TRIM(argName) ! Capture the directory name
!               print *,'argname = ' // trim(argName)
!            CASE("--neutral-file2")
!               neutralFileGiven = .TRUE.
!               neutralFile2  = TRIM(argName) ! Capture the directory name
!               print *,'argname = ' // trim(argName)
!            CASE("--grid-file")
!               gridGiven = .TRUE.
!            CASE("--output-dir")
!               outputGiven = .TRUE.
!               outputDirGiven = .TRUE.
!            CASE DEFAULT

!              IF( plasmaGiven )THEN

!                 plasmaFileGiven = .TRUE.
!                 plasmaFile  = TRIM(argName) ! Capture the directory name
!                 plasmaGiven = .FALSE.

!              ELSEIF( gridGiven )THEN

!                 gridFileGiven = .TRUE.
!                 gridFile  = TRIM(argName)
!                 gridGiven = .FALSE.

!              ELSEIF( outputGiven )THEN

!                 outputDirGiven = .TRUE.
!                 outputDir   = TRIM(argName)
!                 outputGiven = .FALSE.
!              ENDIF

!         END SELECT 
!       ENDDO

!       IF( .NOT. plasmaFileGiven .OR. &
!       IF( .NOT. gridFileGiven )THEN
!          PRINT*, '  i2hg : A grid file is required'
!          PRINT*, '  i2hg : Setup Failed.'
!          STOP
!       ENDIF

        call getarg(1, neutralFile1)
        print *,'neutralFile1 = '//trim(neutralFile1)
        call getarg(2, neutralFile2)
        print *,'neutralFile2 = '//trim(neutralFile2)

       
        INQUIRE( FILE = TRIM(neutralFile1), &
                 EXIST = fileExists )
        IF( .NOT.(fileExists) )THEN
          PRINT*, ' Input file :'//TRIM(neutralFile1)//': not found.'
          PRINT*, '  i2hg : Setup Failed.'
        ENDIF
        INQUIRE( FILE = TRIM(neutralFile2), &
                 EXIST = fileExists )
        IF( .NOT.(fileExists) )THEN
          PRINT*, ' Input file :'//TRIM(neutralFile2)//': not found.'
          PRINT*, '  i2hg : Setup Failed.'
        ENDIF
        
        INQUIRE( FILE = TRIM(gridFile), &
                 EXIST = fileExists )
        IF( .NOT.(fileExists) )THEN
          PRINT*, ' Grid file :'//TRIM(gridFile)//': not found.'
          PRINT*, '  i2hg : Setup Failed.'
        ENDIF

        IF( .NOT. outputDirGiven )THEN
           outputDir = './'
        ENDIF
        
     ELSE

        ! List possible options
        PRINT*, '  i2hg : IPE to Height Grid'
        PRINT*, '    A tool for interpolating IPE plasma output onto a fixed height grid.'
        PRINT*, '--------------------------------------------------------------'
        PRINT*, '  Usage : i2hg <inputs>                                      '
        PRINT*, ''
        PRINT*, ' Required inputs'
        PRINT*, ''
        PRINT*, ' --plasma-file /full/path/to/plasma/file '
        PRINT*, ''
        PRINT*, ' --grid-file /full/path/to/grid/file'
        PRINT*, ''
        PRINT*, ' Optional inputs'
        PRINT*, ''
        PRINT*, ' --neutral-file /full/path/to/neutrals/file '
        PRINT*, ''
        PRINT*, ' If the neutrals file is not provided, these parameters will  '
        PRINT*, ' not be included in the netcdf output.'
        PRINT*, ''
        PRINT*, ' --output-directory /full/path/to/output directory '
        PRINT*, ''
        PRINT*, ' If the output directory is not specified, your current       '
        PRINT*, ' directory will be used.  '
        PRINT*, ' '
        PRINT*, '--------------------------------------------------------------'
        PRINT*, ' This executable takes READs in the plasma file that was      '
        PRINT*, ' output by IPE, and interpolates the ions to a fixed height   '
        PRINT*, ' grid. NmF2, HmF2, and total electron content are also        '
        PRINT*, ' calculated and all outputs are written to a NetCDF file.     '
        PRINT*, '--------------------------------------------------------------'
        STOP
        
     ENDIF   


 END SUBROUTINE InitializeFromCommandLine
!
 SUBROUTINE AllocateArrays( )
   IMPLICIT NONE

      ALLOCATE ( oplus_ft(NPTS2D,80), &
                 hplus_ft(NPTS2D,80), &  
                 heplus_ft(NPTS2D,80), &  
                 nplus_ft(NPTS2D,80), &  
                 noplus_ft(NPTS2D,80), &  
                 o2plus_ft(NPTS2D,80), &  
                 n2plus_ft(NPTS2D,80), &  
                 oplus2d_ft(NPTS2D,80), &  
                 oplus2p_ft(NPTS2D,80), &  
                 tn_ft(NPTS2D,80), &  
                 vn_ft(NPTS2D,80,1:3), &  
                 on_ft(NPTS2D,80), &  
                 n2n_ft(NPTS2D,80), &  
                 o2n_ft(NPTS2D,80) ) 
            

 END SUBROUTINE AllocateArrays
!
 SUBROUTINE CleanupArrays( )
   IMPLICIT NONE

      DEALLOCATE ( oplus_ft, &
                   hplus_ft, &  
                   heplus_ft, &  
                   nplus_ft, &  
                   noplus_ft, &  
                   o2plus_ft, &  
                   n2plus_ft, &  
                   oplus2d_ft, &  
                   oplus2p_ft, &  
                   tn_ft, &  
                   vn_ft, &  
                   on_ft, &  
                   n2n_ft, &  
                   o2n_ft ) 
              

 END SUBROUTINE CleanupArrays
!
 SUBROUTINE ReadGridandPlasmaFiles()
   IMPLICIT NONE

   REAL(real_prec), DIMENSION(:,:,:), ALLOCATABLE :: field1, field2
   REAL(real_prec), DIMENSION(:,:,:,:), ALLOCATABLE :: v1, v2
   REAL(real_prec), ALLOCATABLE :: tn_k(:,:,:)
   REAL(real_prec), ALLOCATABLE :: vn_ms1(:,:,:,:)
   REAL(real_prec), ALLOCATABLE :: on_ms1(:,:,:)
   REAL(real_prec), ALLOCATABLE :: n2n_ms1(:,:,:)
   REAL(real_prec), ALLOCATABLE :: o2n_m3(:,:,:)

   CHARACTER(LEN=*), PARAMETER :: vcomp(3) = (/ 'u', 'v', 'w' /)
   CHARACTER(LEN=*), PARAMETER :: gases(3) = (/ 'O', 'N2', 'O2' /)
   INTEGER :: i 
 
    CALL init_plasma_grid( )

    ALLOCATE( field1(MaxFluxTube,NLP,NMP), &
              field2(MaxFluxTube,NLP,NMP), &
              v1(MaxFluxTube,NLP,NMP,1:3), &
              v2(MaxFluxTube,NLP,NMP,1:3) )

     OPEN( UNIT = 20, &
           FILE = TRIM(neutralFile1), &
           FORM = 'UNFORMATTED', &
           STATUS = 'OLD' )

     OPEN( UNIT = 40, &
           FILE = TRIM(neutralFile2), &
           FORM = 'UNFORMATTED', &
           STATUS = 'OLD' )

     REWIND(20)
     REWIND(40)

    
     READ(20) field1
     READ(40) field2
     CALL ComputeDiff(field1, field2, "temp_neutral")

     READ(20) v1
     READ(40) v2
     do i = 1, 3
       field1 = v1(:,:,:,i)
       field2 = v2(:,:,:,i)
       CALL ComputeDiff(field1, field2, "wind_neutral: "//trim(vcomp(i)))
     end do

     do i = 1, 3
       READ(20) field1
       READ(40) field2
       CALL ComputeDiff(field1, field2, trim(gases(i)) // "_Density")
     end do

     CLOSE(20)
     CLOSE(40)

    DEALLOCATE( field1, field2, v1, v2 )

 END SUBROUTINE ReadGridandPlasmaFiles

 SUBROUTINE ComputeDiff(field1, field2, label)
   REAL(real_prec), DIMENSION(:,:,:), INTENT(IN) :: field1, field2
   CHARACTER(LEN=*), INTENT(IN) :: label

   
   REAL(real_prec) :: dmin, dmax

   dmin = MINVAL(ABS(field1-field2))
   dmax = MAXVAL(ABS(field1-field2))

   WRITE(*,*) 'DIFF REP : '//TRIM(label), dmin, dmax, MAXVAL(field1), MAXVAL(field2)

 END SUBROUTINE ComputeDiff

END PROGRAM IPEtoHeightGrid

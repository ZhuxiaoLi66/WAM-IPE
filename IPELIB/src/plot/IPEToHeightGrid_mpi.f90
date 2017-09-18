! IPEtoHeightGrid.f90
!
!  Authors : George Millward  (george.millward@noaa.gov)
!            Joe Schoonover   (joseph.schoonover@noaa.gov)
!            Adam Kubaryk     (adam.kubaryk@noaa.gov)
!
!  Cooperative Institute for Research in Environmental Sciences (CIRES)
!  National Oceanographic and Atmospheric Administration (NOAA)
!  Space Weather Prediction Center (SWPC)
!
!


PROGRAM IPEtoHeightGrid

USE module_precision
USE netcdf

  IMPLICIT NONE
  INTEGER, PARAMETER :: nheights=183
  INTEGER, PARAMETER :: nlon=90 
  INTEGER, PARAMETER :: nlat=91
  REAL(real_prec8), PARAMETER:: fillValue = -999999.9999_real_prec8



  INTEGER,DIMENSION(3,nheights,nlat,nlon)  :: ii1, ii2, ii3, ii4
  INTEGER,DIMENSION(3,nheights,nlat,nlon)  :: ii1_interface, ii2_interface,&
                                              ii3_interface, ii4_interface
  REAL(kind=real_prec8),DIMENSION(nheights)  :: e_profile
  REAL(kind=real_prec8),DIMENSION(nheights)  :: z
  REAL(kind=real_prec8),DIMENSION(nlon) :: x
  REAL(kind=real_prec8),DIMENSION(nlat) :: y 
  INTEGER :: height_km
  INTEGER :: i, iheight, iup, ido, ih, ip
  INTEGER :: l, m, mp, lp, in1, in2
  INTEGER :: ifile, the_index
  INTEGER :: my_pe, num_pe, iprovided, iret
  INTEGER :: NMP, NLP, NPTS2D

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
  REAL(kind=real_prec8),DIMENSION(:,:,:), ALLOCATABLE   :: hplus_fixed
  REAL(kind=real_prec8),DIMENSION(:,:,:), ALLOCATABLE   :: heplus_fixed
  REAL(kind=real_prec8),DIMENSION(:,:,:), ALLOCATABLE   :: nplus_fixed
  REAL(kind=real_prec8),DIMENSION(:,:,:), ALLOCATABLE   :: noplus_fixed
  REAL(kind=real_prec8),DIMENSION(:,:,:), ALLOCATABLE   :: o2plus_fixed
  REAL(kind=real_prec8),DIMENSION(:,:,:), ALLOCATABLE   :: n2plus_fixed
  REAL(kind=real_prec8),DIMENSION(:,:,:), ALLOCATABLE   :: oplus2d_fixed
  REAL(kind=real_prec8),DIMENSION(:,:,:), ALLOCATABLE   :: oplus2p_fixed
  REAL(kind=real_prec8),DIMENSION(:,:,:), ALLOCATABLE   :: tn_fixed
  REAL(kind=real_prec8),DIMENSION(:,:,:,:), ALLOCATABLE :: vn_fixed
  REAL(kind=real_prec8),DIMENSION(:,:,:), ALLOCATABLE   :: on_fixed
  REAL(kind=real_prec8),DIMENSION(:,:,:), ALLOCATABLE   :: n2n_fixed
  REAL(kind=real_prec8),DIMENSION(:,:,:), ALLOCATABLE   :: o2n_fixed
  REAL(kind=real_prec8),DIMENSION(:,:,:), ALLOCATABLE   :: electron_density_fixed
  REAL(kind=real_prec8),DIMENSION(:,:), ALLOCATABLE     :: electron_density_fixed_300km
  REAL(kind=real_prec8),DIMENSION(:,:), ALLOCATABLE     :: total_electron_content
  REAL(kind=real_prec8),DIMENSION(:,:), ALLOCATABLE     :: hmf2
  REAL(kind=real_prec8),DIMENSION(:,:), ALLOCATABLE     :: nmf2
  REAL(kind=real_prec8),DIMENSION(3)                    :: oplus_interpolated
  REAL(kind=real_prec8),DIMENSION(3)                    :: hplus_interpolated
  REAL(kind=real_prec8),DIMENSION(3)                    :: heplus_interpolated
  REAL(kind=real_prec8),DIMENSION(3)                    :: nplus_interpolated
  REAL(kind=real_prec8),DIMENSION(3)                    :: noplus_interpolated
  REAL(kind=real_prec8),DIMENSION(3)                    :: o2plus_interpolated
  REAL(kind=real_prec8),DIMENSION(3)                    :: n2plus_interpolated
  REAL(kind=real_prec8),DIMENSION(3)                    :: oplus2d_interpolated
  REAL(kind=real_prec8),DIMENSION(3)                    :: oplus2p_interpolated
  REAL(kind=real_prec8),DIMENSION(3)                    :: tn_interpolated
  REAL(kind=real_prec8),DIMENSION(3,3)                  :: vn_interpolated
  REAL(kind=real_prec8),DIMENSION(3)                    :: on_interpolated
  REAL(kind=real_prec8),DIMENSION(3)                    :: n2n_interpolated
  REAL(kind=real_prec8),DIMENSION(3)                    :: o2n_interpolated

  CHARACTER(100) :: filename
  CHARACTER(500) :: plasmaFile
  CHARACTER(500) :: neutralFile
  CHARACTER(500) :: gridFile
  CHARACTER(500) :: outputDir
  

  LOGICAL :: plasmaFileGiven, neutralFileGiven, outputDirGiven, gridFileGiven

    include 'mpif.h'

    CALL MPI_Init_thREAD(MPI_THREAD_SERIALIZED,iprovided,iret)
    CALL MPI_Comm_rank(MPI_COMM_WORLD,my_pe,iret)
    CALL MPI_Comm_size(MPI_COMM_WORLD,num_pe,iret)

    CALL InitializeFromCommandLine( )

    CALL ReadInputParameters( )

    CALL AllocateArrays( )

    PRINT*, "ReadGridAndPlasmaFiles"
    CALL ReadGridAndPlasmaFiles( )
 
    PRINT*, "InterpolateOntoFixedHeightGrid"
    CALL InterpolateOntoFixedHeightGrid( )

    PRINT*,"WriteToNetCDF"
    CALL WriteToNetCDF( )

    CALL CleanupArrays( )

    CALL MPI_Finalize(iret)

CONTAINS
 SUBROUTINE InitializeFromCommandLine( )
   IMPLICIT NONE

   INTEGER        :: nArg, argID, io
   CHARACTER(500) :: argname, neutralFileFile, plasmaFileFile
   LOGICAL        :: fileExists
   LOGICAL        :: plasmaGiven, neutralsGiven,  gridGiven, outputGiven

     plasmaFileGiven   = .FALSE. 
     gridFileGiven     = .FALSE.
     neutralFileGiven = .FALSE.
     outputDirGiven    = .FALSE.

     plasmaGiven    = .FALSE.
     neutralsGiven  = .FALSE.
     gridGiven      = .FALSE.
     outputGiven    = .FALSE.

     nArg = command_argument_count( )

     IF( nArg > 0 )THEN
 
        DO argID = 1, nArg
  
          CALL get_command_argument( argID, argName )

          SELECT CASE( TRIM(argName) )

             CASE("--plasma-file")
                plasmaGiven = .TRUE.
                plasmaFileGiven = .TRUE.
             CASE("--neutral-file")
                neutralsGiven = .TRUE.
                neutralFileGiven = .TRUE.
             CASE("--grid-file")
                gridGiven = .TRUE.
             CASE("--output-dir")
                outputGiven = .TRUE.
                outputDirGiven = .TRUE.
             CASE DEFAULT

               IF( plasmaGiven )THEN

                  plasmaFileGiven = .TRUE.
                  plasmaFileFile  = TRIM(argName) ! capture list of plasma files
                  OPEN(1, file=plasmaFileFile)
                  DO i=0,my_pe
                     READ(1,'(a)',advance='NO',iostat=io) plasmaFile
                  END DO
                  CLOSE(1)
                  plasmaGiven = .FALSE.

               ELSEIF( neutralsGiven )THEN

                  neutralFileGiven = .TRUE.
                  neutralFileFile  = TRIM(argName) ! capture list of neutral files
                  OPEN(2, file=neutralFileFile)
                  DO i=0,my_pe
                     READ(2,'(a)',advance='NO',iostat=io) neutralFile
                  END DO
                  CLOSE(2)
                  neutralsGiven = .FALSE.

               ELSEIF( gridGiven )THEN

                  gridFileGiven = .TRUE.
                  gridFile  = TRIM(argName)
                  gridGiven = .FALSE.

               ELSEIF( outputGiven )THEN

                  outputDirGiven = .TRUE.
                  outputDir   = TRIM(argName)
                  outputGiven = .FALSE.
               ENDIF

          END SELECT 
        ENDDO

        IF( .NOT. plasmaFileGiven .OR. &
            .NOT. gridFileGiven )THEN
           PRINT*, '  i2hg : A plasma file and grid file are required'
           PRINT*, '  i2hg : Setup Failed.'
           STOP
        ENDIF
       
        
        INQUIRE( FILE = TRIM(plasmaFile), &
                 EXIST = fileExists )
        IF( .NOT.(fileExists) )THEN
          PRINT*, ' Input file :'//TRIM(plasmaFile)//': not found.'
          PRINT*, '  i2hg : Setup Failed.'
        ENDIF
        
        INQUIRE( FILE = TRIM(gridFile), &
                 EXIST = fileExists )
        IF( .NOT.(fileExists) )THEN
          PRINT*, ' Grid file :'//TRIM(gridFile)//': not found.'
          PRINT*, '  i2hg : Setup Failed.'
        ENDIF


        IF( neutralFileGiven )THEN
        
           INQUIRE( FILE = TRIM(neutralFile), &
                    EXIST = fileExists )
           IF( .NOT.(fileExists) )THEN
             PRINT*, ' Input file :'//TRIM(neutralFile)//': not found.'
             PRINT*, '  i2hg : Setup Failed.'
           ENDIF

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
 SUBROUTINE ReadInputParameters( )
   IMPLICIT NONE
   CHARACTER(LEN=*),PARAMETER :: INPTNMLT='IPE.inp'
   INTEGER ierr
   NAMELIST/IPEDIMS/NLP,NMP,NPTS2D
   OPEN(1,FILE=INPTNMLT,STATUS='OLD',IOSTAT=ierr)
   READ(1,NML=IPEDIMS,IOSTAT=ierr)
   CLOSE(1)
 END SUBROUTINE ReadInputParameters
!
 SUBROUTINE AllocateArrays( )
   IMPLICIT NONE

      ALLOCATE ( facfac_interface(1:3,nheights,nlat,nlon), &
                 dd_interface(1:3,nheights,nlat,nlon), &
                 oplus_ft(NPTS2D,80), &
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
                 o2n_ft(NPTS2D,80), &  
                 oplus_fixed(nlon,nlat,nheights), &
                 hplus_fixed(nlon,nlat,nheights), &
                 heplus_fixed(nlon,nlat,nheights), &
                 nplus_fixed(nlon,nlat,nheights), &
                 noplus_fixed(nlon,nlat,nheights), &
                 o2plus_fixed(nlon,nlat,nheights), &
                 n2plus_fixed(nlon,nlat,nheights), &
                 oplus2d_fixed(nlon,nlat,nheights), &
                 oplus2p_fixed(nlon,nlat,nheights), &
                 tn_fixed(nlon,nlat,nheights), &
                 vn_fixed(nlon,nlat,nheights,1:3), &
                 on_fixed(nlon,nlat,nheights), &
                 n2n_fixed(nlon,nlat,nheights), &
                 o2n_fixed(nlon,nlat,nheights), &
                 electron_density_fixed(nlon,nlat,nheights), &
                 electron_density_fixed_300km(nlon,nlat), &
                 total_electron_content(nlon,nlat), &
                 hmf2(nlon,nlat), nmf2(nlon,nlat) )
            

 END SUBROUTINE AllocateArrays
!
 SUBROUTINE CleanupArrays( )
   IMPLICIT NONE

      DEALLOCATE ( facfac_interface, &
                   dd_interface, &
                   oplus_ft, &
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
                   o2n_ft, &  
                   oplus_fixed, &
                   hplus_fixed, &
                   heplus_fixed, &
                   nplus_fixed, &
                   noplus_fixed, &
                   o2plus_fixed, &
                   n2plus_fixed, &
                   oplus2d_fixed, &
                   oplus2p_fixed, &
                   tn_fixed, &
                   vn_fixed, &
                   on_fixed, &
                   n2n_fixed, &
                   o2n_fixed, &
                   electron_density_fixed, &
                   electron_density_fixed_300km, &
                   total_electron_content, &
                   hmf2, nmf2 )
              

 END SUBROUTINE CleanupArrays
!
 SUBROUTINE ReadGridandPlasmaFiles( )
   IMPLICIT NONE
   REAL(real_prec), ALLOCATABLE :: tn_k(:,:,:)
   REAL(real_prec), ALLOCATABLE :: vn_ms1(:,:,:,:)
   REAL(real_prec), ALLOCATABLE :: on_ms1(:,:,:)
   REAL(real_prec), ALLOCATABLE :: n2n_ms1(:,:,:)
   REAL(real_prec), ALLOCATABLE :: o2n_m3(:,:,:)
   INTEGER,DIMENSION(:,:), ALLOCATABLE :: JMIN_IN_ALL, JMAX_IS_ALL
   INTEGER,DIMENSION(:  ), ALLOCATABLE :: JMIN_IN, JMAX_IS, JMIN_ING, JMAX_ISG
   INTEGER ierr, MaxFluxTube
 
    allocate(JMIN_IN_ALL(NMP,NLP),JMAX_IS_ALL(NMP,NLP),JMIN_IN(NLP),JMAX_IS(NLP),JMIN_ING(NLP),JMAX_ISG(NLP))
    if ( my_pe .eq. 0 ) then
    OPEN(UNIT=20,file=TRIM(gridFile),form='unformatted',status='old')
       READ(20) facfac_interface
       READ(20) dd_interface
       READ(20) ii1_interface
       READ(20) ii2_interface
       READ(20) ii3_interface
       READ(20) ii4_interface
    CLOSE(20)

    OPEN(UNIT=20, file='ipe_grid', form='formatted', status='old')
       READ (20, FMT=*) JMIN_IN_all, JMAX_IS_all
       JMIN_ING = JMIN_IN_all(1,:)
       JMAX_ISG = JMAX_IS_all(1,:)
       JMIN_IN = 1
       JMAX_IS = JMAX_ISG - JMIN_ING + 1
    CLOSE(20)
    end if
    call MPI_Bcast(facfac_interface,3*nheights*nlat*nlon,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(dd_interface    ,3*nheights*nlat*nlon,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(ii1_interface   ,3*nheights*nlat*nlon,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(ii2_interface   ,3*nheights*nlat*nlon,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(ii3_interface   ,3*nheights*nlat*nlon,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(ii4_interface   ,3*nheights*nlat*nlon,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(JMIN_IN ,NLP,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(JMIN_ING,NLP,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(JMAX_IS ,NLP,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(JMAX_ISG,NLP,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    MaxFluxTube = maxval(JMAX_ISG-JMIN_ING+1)
    
    ALLOCATE( tn_k(MaxFluxTube,NLP,NMP), &
              vn_ms1(MaxFluxTube,NLP,NMP,1:3), &
              on_ms1(MaxFluxTube,NLP,NMP), &
              n2n_ms1(MaxFluxTube,NLP,NMP), &
              o2n_m3(MaxFluxTube,NLP,NMP) )


     OPEN(UNIT   = 20, &
          FILE   = TRIM(plasmaFile),&
          FORM   = 'UNFORMATTED', &
          STATUS = 'OLD')

        READ(20) oplus_ft
        READ(20) hplus_ft
        READ(20) heplus_ft
        READ(20) nplus_ft
        READ(20) noplus_ft
        READ(20) o2plus_ft
        READ(20) n2plus_ft
        READ(20) oplus2d_ft
        READ(20) oplus2p_ft

     CLOSE(20)

     OPEN( UNIT = 20, &
           FILE = TRIM(neutralFile), &
           FORM = 'UNFORMATTED', &
           STATUS = 'OLD' )


     READ(20) tn_k
     READ(20) vn_ms1
     READ(20) on_ms1
     READ(20) n2n_ms1
     READ(20) o2n_m3 

     DO mp=1,NMP
       DO lp=1,NLP

        tn_ft(JMIN_ING(lp):JMAX_ISG(lp),mp)     = tn_k(JMIN_IN(lp):JMAX_IS(lp),lp,mp)
        vn_ft(JMIN_ING(lp):JMAX_ISG(lp),mp,1:3) = vn_ms1(JMIN_IN(lp):JMAX_IS(lp),lp,mp,1:3)
        on_ft(JMIN_ING(lp):JMAX_ISG(lp),mp)     = on_ms1(JMIN_IN(lp):JMAX_IS(lp),lp,mp)
        n2n_ft(JMIN_ING(lp):JMAX_ISG(lp),mp)    = n2n_ms1(JMIN_IN(lp):JMAX_IS(lp),lp,mp)
        o2n_ft(JMIN_ING(lp):JMAX_ISG(lp),mp)    = o2n_m3(JMIN_IN(lp):JMAX_IS(lp),lp,mp)

      ENDDO 
     ENDDO



#ifdef DEBUG
  PRINT*, MAXVAL(oplus_ft), MINVAL(oplus_ft)
  PRINT*, MAXVAL(hplus_ft), MINVAL(hplus_ft)
  PRINT*, MAXVAL(heplus_ft), MINVAL(heplus_ft)
  PRINT*, MAXVAL(nplus_ft), MINVAL(nplus_ft)
  PRINT*, MAXVAL(noplus_ft), MINVAL(noplus_ft)
  PRINT*, MAXVAL(o2plus_ft), MINVAL(o2plus_ft)
  PRINT*, MAXVAL(n2plus_ft), MINVAL(n2plus_ft)
  PRINT*, 'Neutrals'
  PRINT*, MAXVAL(vn_ms1(:,:,:,1)), MINVAL(vn_ms1(:,:,:,1))
  PRINT*, MAXVAL(vn_ms1(:,:,:,2)), MINVAL(vn_ms1(:,:,:,2))
  PRINT*, MAXVAL(vn_ms1(:,:,:,3)), MINVAL(vn_ms1(:,:,:,3))
#endif


    DEALLOCATE( tn_k, &
                vn_ms1, &
                on_ms1, &
                n2n_ms1, &
                o2n_m3  )


 END SUBROUTINE ReadGridandPlasmaFiles

 SUBROUTINE InterpolateOntoFixedHeightGrid( )
  IMPLICIT NONE
  REAL(real_prec8) :: dlat, dlon
  REAL(real_prec8) :: a, b, c, x1, x2, x3, y1, y2, y3
  INTEGER         :: i1, i2, i3, i_max_location

    dlat = 180.0_real_prec/REAL( nlat-1,real_prec8 ) 
    dlon = 360.0_real_prec/REAL( nlon-1,real_prec8 ) 

    do l=1,nlat
      y(l) = -90.0_real_prec + REAL(l-1,real_prec8)*dlat
    enddo
    do m=1,nlon
      x(m) = REAL(m-1,real_prec8)*dlon
    enddo


    do iheight=1,nheights
      do l=1,nlat
          do m=1,nlon

              height_km = ((iheight - 1) * 5) + nlat
              z(iHeight) = REAL( height_km, real_prec8 ) 
              oplus_fixed(m,l,iheight) = 0.0
              hplus_fixed(m,l,iheight) = 0.0
              heplus_fixed(m,l,iheight) = 0.0
              nplus_fixed(m,l,iheight) = 0.0
              noplus_fixed(m,l,iheight) = 0.0
              o2plus_fixed(m,l,iheight) = 0.0
              n2plus_fixed(m,l,iheight) = 0.0
              oplus2d_fixed(m,l,iheight) = 0.0
              oplus2p_fixed(m,l,iheight) = 0.0
              dtotinv = 0.0

              do i=1,3
                  factor=facfac_interface(i,iheight,l,m)
                  mp=ii1_interface(i,iheight,l,m)
                  lp=ii2_interface(i,iheight,l,m)
                  in2=ii3_interface(i,iheight,l,m)
                  in1=ii4_interface(i,iheight,l,m)

                  ! Ions
                  oplus_interpolated(i) = ((oplus_ft(in2,mp)-oplus_ft(in1,mp))*factor) + oplus_ft(in1,mp)
                  hplus_interpolated(i) = ((hplus_ft(in2,mp)-hplus_ft(in1,mp))*factor) + hplus_ft(in1,mp)
                  heplus_interpolated(i) = ((heplus_ft(in2,mp)-heplus_ft(in1,mp))*factor) + heplus_ft(in1,mp)
                  nplus_interpolated(i) = ((nplus_ft(in2,mp)-nplus_ft(in1,mp))*factor) + nplus_ft(in1,mp)
                  noplus_interpolated(i) = ((noplus_ft(in2,mp)-noplus_ft(in1,mp))*factor) + noplus_ft(in1,mp)
                  o2plus_interpolated(i) = ((o2plus_ft(in2,mp)-o2plus_ft(in1,mp))*factor) + o2plus_ft(in1,mp)
                  n2plus_interpolated(i) = ((n2plus_ft(in2,mp)-n2plus_ft(in1,mp))*factor) + n2plus_ft(in1,mp)
                  oplus2d_interpolated(i) = ((oplus2d_ft(in2,mp)-oplus2d_ft(in1,mp))*factor) + oplus2d_ft(in1,mp)
                  oplus2p_interpolated(i) = ((oplus2p_ft(in2,mp)-oplus2p_ft(in1,mp))*factor) + oplus2p_ft(in1,mp)

                  ! Neutrals
                  tn_interpolated(i)   = ((tn_ft(in2,mp)-tn_ft(in1,mp))*factor) + tn_ft(in1,mp)
                  vn_interpolated(i,1) = ((vn_ft(in2,mp,1)-vn_ft(in1,mp,1))*factor) + vn_ft(in1,mp,1)
                  vn_interpolated(i,2) = ((vn_ft(in2,mp,2)-vn_ft(in1,mp,2))*factor) + vn_ft(in1,mp,2)
                  vn_interpolated(i,3) = ((vn_ft(in2,mp,3)-vn_ft(in1,mp,3))*factor) + vn_ft(in1,mp,3)
                  on_interpolated(i)   = ((on_ft(in2,mp)-on_ft(in1,mp))*factor) + on_ft(in1,mp)
                  n2n_interpolated(i)  = ((n2n_ft(in2,mp)-n2n_ft(in1,mp))*factor) + n2n_ft(in1,mp)
                  o2n_interpolated(i)  = ((o2n_ft(in2,mp)-o2n_ft(in1,mp))*factor) + o2n_ft(in1,mp)



                  d = dd_interface(i,iheight,l,m)
                  dtotinv = (1./d) + dtotinv

                  ! Ions
                  oplus_fixed(m,l,iheight) = (oplus_interpolated(i)/d) + oplus_fixed(m,l,iheight)
                  hplus_fixed(m,l,iheight) = (hplus_interpolated(i)/d) + hplus_fixed(m,l,iheight)
                  heplus_fixed(m,l,iheight) = (heplus_interpolated(i)/d) + heplus_fixed(m,l,iheight)
                  nplus_fixed(m,l,iheight) = (nplus_interpolated(i)/d) + nplus_fixed(m,l,iheight)
                  noplus_fixed(m,l,iheight) = (noplus_interpolated(i)/d) + noplus_fixed(m,l,iheight)
                  o2plus_fixed(m,l,iheight) = (o2plus_interpolated(i)/d) + o2plus_fixed(m,l,iheight)
                  n2plus_fixed(m,l,iheight) = (n2plus_interpolated(i)/d) + n2plus_fixed(m,l,iheight)
                  oplus2d_fixed(m,l,iheight) = (oplus2d_interpolated(i)/d) + oplus2d_fixed(m,l,iheight)
                  oplus2p_fixed(m,l,iheight) = (oplus2p_interpolated(i)/d) + oplus2p_fixed(m,l,iheight)

                  ! Neutrals
                  tn_fixed(m,l,iheight)     = (tn_interpolated(i)/d) + tn_fixed(m,l,iheight)
                  vn_fixed(m,l,iheight,1:3) = (vn_interpolated(i,1:3)/d) + vn_fixed(m,l,iheight,1:3)
                  on_fixed(m,l,iheight)     = (on_interpolated(i)/d) + on_fixed(m,l,iheight)
                  n2n_fixed(m,l,iheight)    = (n2n_interpolated(i)/d) + n2n_fixed(m,l,iheight)
                  o2n_fixed(m,l,iheight)    = (o2n_interpolated(i)/d) + o2n_fixed(m,l,iheight)

              ENDDO ! do 400 i=1,3

              ! Ions
              oplus_fixed(m,l,iheight)   = oplus_fixed(m,l,iheight)/dtotinv
              hplus_fixed(m,l,iheight)   = hplus_fixed(m,l,iheight)/dtotinv
              heplus_fixed(m,l,iheight)  = heplus_fixed(m,l,iheight)/dtotinv
              nplus_fixed(m,l,iheight)   = nplus_fixed(m,l,iheight)/dtotinv
              noplus_fixed(m,l,iheight)  = noplus_fixed(m,l,iheight)/dtotinv
              o2plus_fixed(m,l,iheight)  = o2plus_fixed(m,l,iheight)/dtotinv
              n2plus_fixed(m,l,iheight)  = n2plus_fixed(m,l,iheight)/dtotinv
              oplus2d_fixed(m,l,iheight) = oplus2d_fixed(m,l,iheight)/dtotinv
              oplus2p_fixed(m,l,iheight) = oplus2p_fixed(m,l,iheight)/dtotinv

              ! Neutrals
              tn_fixed(m,l,iheight)     = tn_fixed(m,l,iheight)/dtotinv
              vn_fixed(m,l,iheight,1:3) = vn_fixed(m,l,iheight,1:3)/dtotinv
              on_fixed(m,l,iheight)     = on_fixed(m,l,iheight)/dtotinv
              n2n_fixed(m,l,iheight)    = n2n_fixed(m,l,iheight)/dtotinv
              o2n_fixed(m,l,iheight)    = o2n_fixed(m,l,iheight)/dtotinv

          ENDDO
      ENDDO
  ENDDO

#ifdef DEBUG
  PRINT*, MAXVAL(oplus_fixed), MINVAL(oplus_fixed)
  PRINT*, MAXVAL(hplus_fixed), MINVAL(hplus_fixed)
  PRINT*, MAXVAL(heplus_fixed), MINVAL(heplus_fixed)
  PRINT*, MAXVAL(nplus_fixed), MINVAL(nplus_fixed)
  PRINT*, MAXVAL(noplus_fixed), MINVAL(noplus_fixed)
  PRINT*, MAXVAL(o2plus_fixed), MINVAL(o2plus_fixed)
  PRINT*, MAXVAL(n2plus_fixed), MINVAL(n2plus_fixed)
#endif
  electron_density_fixed = oplus_fixed + hplus_fixed + heplus_fixed + nplus_fixed + noplus_fixed & 
                           + o2plus_fixed + n2plus_fixed + oplus2d_fixed + oplus2p_fixed

  total_electron_content = 0.0
  DO iheight=1,nheights
    total_electron_content(1:nlon,1:nlat) = total_electron_content(1:nlon,1:nlat) + &
                                       (electron_density_fixed(1:nlon,1:nlat,iheight) * 5000.0)                                           
  ENDDO
  total_electron_content = total_electron_content*10.0_real_prec**(-16)

! NmF2 and hmF2....
  DO l = 1 , nlat
    DO  m = 1 , nlon

      e_profile = electron_density_fixed(m,l,1:183)
      i_max_location = maxval(maxloc(e_profile(1:183)))
    
      y1 = e_profile(i_max_location - 1)
      y2 = e_profile(i_max_location)
      y3 = e_profile(i_max_location + 1)
      x1 = (float(i_max_location - 2) * 5.0) + 90.0
      x2 = (float(i_max_location - 1) * 5.0) + 90.0
      x3 = (float(i_max_location) * 5.0) + 90.0
    
      c = (x3*y1 - x3*y2 + x1*y2 + x2*y3 - x2*y1 - x1*y3) / &
          (x3*x1*x1 - x3*x2*x2 + x1*x2*x2 - x2*x1*x1 + x2*x3*x3 - x1*x3*x3)
      b = (y2 - (c*x2*x2) + (c*x1*x1) - y1) / (x2 - x1)
      a = y1 - (b*x1) - (c*x1*x1)
      hmf2(m,l) = (0.0 - b) / (2*c)
      nmf2(m,l) = a + (b*hmf2(m,l)) + (c*hmf2(m,l)*hmf2(m,l))

    
    ENDDO
  ENDDO

 END SUBROUTINE InterpolateOntoFixedHeightGrid
!
 SUBROUTINE WriteToNetCDF(  )
   IMPLICIT NONE
   ! Local
   INTEGER :: NF90_PREC
   INTEGER :: ncid
   INTEGER :: x_dimid, y_dimid, z_dimid, time_dimid
   INTEGER :: x_varid, y_varid, z_varid, time_varid
   INTEGER :: oplus_varid, hplus_varid, heplus_varid
   INTEGER :: nplus_varid, noplus_varid, o2plus_varid
   INTEGER :: n2plus_varid, tec_varid, nmf2_varid, hmf2_varid
   INTEGER :: edens_varid
   INTEGER :: tn_varid, u_varid, v_varid, w_varid, on_varid, n2n_varid, o2n_varid

   CHARACTER( LEN(TRIM(plasmaFile)) ) :: shortenedFile
   CHARACTER( 12 )                    :: timeStamp


      NF90_PREC = NF90_DOUBLE

      shortenedFile = TRIM(plasmaFile)
      timestamp     = shortenedFile( LEN(shortenedFile)-11:LEN(shortenedFile) )

      CALL Check( nf90_create( TRIM(outputDir)//'IPE_State.'//TRIM(timestamp)//'.nc',&
                               NF90_NETCDF4,ncid))

      CALL Check( nf90_def_dim( ncid, "z", nheights, z_dimid ) ) 
      CALL Check( nf90_def_dim( ncid, "longitude", nlon, x_dimid ) ) 
      CALL Check( nf90_def_dim( ncid, "latitude", nlat, y_dimid ) ) 
      CALL Check( nf90_def_dim( ncid, "time", NF90_UNLIMITED, time_dimid ) )

      ! Create variables -- here we need to create arrays for the dimensions
      CALL Check( nf90_def_var( ncid, "z", NF90_PREC, z_dimid, z_varid ) )
      CALL Check( nf90_put_att( ncid, z_varid, "long_name", "Height above sea level" ) )
      CALL Check( nf90_put_att( ncid, z_varid, "units", "km" ) )
      CALL Check( nf90_put_att( ncid, z_varid, "_FillValue", fillValue) )
      CALL Check( nf90_put_att( ncid, z_varid, "missing_value", fillValue) )


      CALL Check( nf90_def_var( ncid, "latitude", NF90_PREC, y_dimid, y_varid ) )
      CALL Check( nf90_put_att( ncid, y_varid, "long_name", "Geographic Latitude" ) )
      CALL Check( nf90_put_att( ncid, y_varid, "units", "degrees_north" ) )
      CALL Check( nf90_put_att( ncid, y_varid, "_FillValue", fillValue) )
      CALL Check( nf90_put_att( ncid, y_varid, "missing_value", fillValue) )


      CALL Check( nf90_def_var( ncid, "longitude", NF90_PREC, x_dimid, x_varid ) )
      CALL Check( nf90_put_att( ncid, x_varid, "long_name", "Geographic Longitude" ) )
      CALL Check( nf90_put_att( ncid, x_varid, "units", "degrees_east" ) )
      CALL Check( nf90_put_att( ncid, x_varid, "_FillValue", fillValue) )
      CALL Check( nf90_put_att( ncid, x_varid, "missing_value", fillValue) )

      CALL Check( nf90_def_var( ncid, "time", NF90_PREC, time_dimid, time_varid ) )
      CALL Check( nf90_put_att( ncid, time_varid, "long_name", "time" ) )
      CALL Check( nf90_put_att( ncid, time_varid, "units", "hours since "//timestamp(1:4)//"-"//timestamp(5:6)//"-"//timestamp(7:8)//" "//timestamp(9:10)//":"//timestamp(11:12)//":00" ) )
      CALL Check( nf90_put_att( ncid, time_varid, "_FillValue", fillValue) )
      CALL Check( nf90_put_att( ncid, time_varid, "missing_value", fillValue) )

      CALL Check( nf90_def_var( ncid, "O+", NF90_PREC,&
                               (/ x_dimid, y_dimid, z_dimid, time_dimid /),&
                               oplus_varid ) )
      CALL Check( nf90_put_att( ncid, oplus_varid, "long_name","Density of positively chared Oxygen ions" ) )
      CALL Check( nf90_put_att( ncid, oplus_varid, "units","[unknown]" ) )
      CALL Check( nf90_put_att( ncid, oplus_varid, "_FillValue",fillValue) )
      CALL Check( nf90_put_att( ncid, oplus_varid,"coordinates", "latitude longitude z" ) )


      CALL Check( nf90_def_var( ncid, "H+", NF90_PREC,&
                               (/ x_dimid, y_dimid, z_dimid, time_dimid /),&
                               hplus_varid ) )
      CALL Check( nf90_put_att( ncid, hplus_varid, "long_name","Density of positively chared Hydrogen ions" ) )
      CALL Check( nf90_put_att( ncid, hplus_varid, "units","[unknown]" ) )
      CALL Check( nf90_put_att( ncid, hplus_varid, "_FillValue",fillValue) )
      CALL Check( nf90_put_att( ncid, hplus_varid,"coordinates", "latitude longitude z" ) )


      CALL Check( nf90_def_var( ncid, "He+", NF90_PREC,&
                               (/ x_dimid, y_dimid, z_dimid, time_dimid /),&
                               heplus_varid ) )
      CALL Check( nf90_put_att( ncid, heplus_varid, "long_name","Density of positively chared Helium ions" ) )
      CALL Check( nf90_put_att( ncid, heplus_varid, "units","[unknown]" ) )
      CALL Check( nf90_put_att( ncid, heplus_varid, "_FillValue",fillValue) )
      CALL Check( nf90_put_att( ncid, heplus_varid,"coordinates", "latitude longitude z" ) )


      CALL Check( nf90_def_var( ncid, "N+", NF90_PREC,&
                               (/ x_dimid, y_dimid, z_dimid, time_dimid /),&
                               nplus_varid ) )
      CALL Check( nf90_put_att( ncid, nplus_varid, "long_name","Density of positively charged Nitrogen ions" ) )
      CALL Check( nf90_put_att( ncid, nplus_varid, "units","[unknown]" ) )
      CALL Check( nf90_put_att( ncid, nplus_varid, "_FillValue",fillValue) )
      CALL Check( nf90_put_att( ncid, nplus_varid,"coordinates", "latitude longitude z" ) )


      CALL Check( nf90_def_var( ncid, "NO+", NF90_PREC,&
                               (/ x_dimid, y_dimid, z_dimid, time_dimid /),&
                               noplus_varid ) )
      CALL Check( nf90_put_att( ncid, noplus_varid, "long_name","Density of Nitrilooxnium." ) )
      CALL Check( nf90_put_att( ncid, noplus_varid, "units","[unknown]" ) )
      CALL Check( nf90_put_att( ncid, noplus_varid, "_FillValue",fillValue) )
      CALL Check( nf90_put_att( ncid, noplus_varid,"coordinates", "latitude longitude z" ) )


      CALL Check( nf90_def_var( ncid, "O2+", NF90_PREC,&
                               (/ x_dimid, y_dimid, z_dimid, time_dimid /),&
                               o2plus_varid ) )
      CALL Check( nf90_put_att( ncid, o2plus_varid, "long_name","Density of Dioxygenyl." ) )
      CALL Check( nf90_put_att( ncid, o2plus_varid, "units","[unknown]" ) )
      CALL Check( nf90_put_att( ncid, o2plus_varid, "_FillValue",fillValue) )
      CALL Check( nf90_put_att( ncid, o2plus_varid,"coordinates", "latitude longitude z" ) )


      CALL Check( nf90_def_var( ncid, "N2+", NF90_PREC,&
                               (/ x_dimid, y_dimid, z_dimid, time_dimid /),&
                               n2plus_varid ) )
      CALL Check( nf90_put_att( ncid, n2plus_varid, "long_name","Density of positively charged Nitrogen ions" ) )
      CALL Check( nf90_put_att( ncid, n2plus_varid, "units","[unknown]" ) )
      CALL Check( nf90_put_att( ncid, n2plus_varid, "_FillValue",fillValue) )
      CALL Check( nf90_put_att( ncid, n2plus_varid,"coordinates", "latitude longitude z" ) )


      CALL Check( nf90_def_var( ncid, "e", NF90_PREC,&
                               (/ x_dimid, y_dimid, z_dimid, time_dimid /),&
                               edens_varid ) )
      CALL Check( nf90_put_att( ncid, edens_varid, "long_name","Total electron density" ) )
      CALL Check( nf90_put_att( ncid, edens_varid, "units","[unknown]" ) )
      CALL Check( nf90_put_att( ncid, edens_varid, "_FillValue",fillValue) )
      CALL Check( nf90_put_att( ncid, edens_varid,"coordinates", "latitude longitude z" ) )


      CALL Check( nf90_def_var( ncid, "TEC", NF90_PREC,&
                               (/ x_dimid, y_dimid, time_dimid /),&
                               tec_varid ) )
      CALL Check( nf90_put_att( ncid, tec_varid, "long_name","Total electron content" ) )
      CALL Check( nf90_put_att( ncid, tec_varid, "units","TECu" ) )
      CALL Check( nf90_put_att( ncid, tec_varid, "_FillValue",fillValue) )
      CALL Check( nf90_put_att( ncid, tec_varid,"coordinates", "latitude longitude" ) )


      CALL Check( nf90_def_var( ncid, "nmf2", NF90_PREC,&
                               (/ x_dimid, y_dimid, time_dimid /),&
                               nmf2_varid ) )
      CALL Check( nf90_put_att( ncid, nmf2_varid, "long_name","Maximum electron density in the F layer." ) )
      CALL Check( nf90_put_att( ncid, nmf2_varid, "units","[unknown]" ) )
      CALL Check( nf90_put_att( ncid, nmf2_varid, "_FillValue",fillValue) )
      CALL Check( nf90_put_att( ncid, nmf2_varid,"coordinates", "latitude longitude" ) )


      CALL Check( nf90_def_var( ncid, "hmf2", NF90_PREC,&
                               (/ x_dimid, y_dimid, time_dimid /),&
                               hmf2_varid ) )
      CALL Check( nf90_put_att( ncid, hmf2_varid, &
                                "long_name","Height of the maximum electron density in the F layer." ) )
      CALL Check( nf90_put_att( ncid, hmf2_varid, "units","[unknown]" ) )
      CALL Check( nf90_put_att( ncid, hmf2_varid, "_FillValue",fillValue) )
      CALL Check( nf90_put_att( ncid, hmf2_varid,"coordinates", "latitude longitude" ) )
         

      CALL Check( nf90_def_var( ncid, "tn", NF90_PREC,&
                               (/ x_dimid, y_dimid, z_dimid, time_dimid /),&
                               tn_varid ) )
      CALL Check( nf90_put_att( ncid, tn_varid, "long_name","Thermosphere temperature" ) )
      CALL Check( nf90_put_att( ncid, tn_varid, "units","[K]" ) )
      CALL Check( nf90_put_att( ncid, tn_varid, "_FillValue",fillValue) )
      CALL Check( nf90_put_att( ncid, tn_varid,"coordinates", "latitude longitude z" ) )

      CALL Check( nf90_def_var( ncid, "vn_zonal", NF90_PREC,&
                               (/ x_dimid, y_dimid, z_dimid, time_dimid /),&
                               u_varid ) )
      CALL Check( nf90_put_att( ncid, u_varid, "long_name","Geographic Zonal wind (positive eastward)" ) )
      CALL Check( nf90_put_att( ncid, u_varid, "units","[m/s]" ) )
      CALL Check( nf90_put_att( ncid, u_varid, "_FillValue",fillValue) )
      CALL Check( nf90_put_att( ncid, u_varid,"coordinates", "latitude longitude z" ) )

      CALL Check( nf90_def_var( ncid, "vn_meridional", NF90_PREC,&
                               (/ x_dimid, y_dimid, z_dimid, time_dimid /),&
                               v_varid ) )
      CALL Check( nf90_put_att( ncid, v_varid, "long_name","Geographic Meridional wind (positive northward)" ) )
      CALL Check( nf90_put_att( ncid, v_varid, "units","[m/s]" ) )
      CALL Check( nf90_put_att( ncid, v_varid, "_FillValue",fillValue) )
      CALL Check( nf90_put_att( ncid, v_varid,"coordinates", "latitude longitude z" ) )

      CALL Check( nf90_def_var( ncid, "vn_vertical", NF90_PREC,&
                               (/ x_dimid, y_dimid, z_dimid, time_dimid /),&
                               w_varid ) )
      CALL Check( nf90_put_att( ncid, w_varid, "long_name","Geographic Vertical wind (positive upward)" ) )
      CALL Check( nf90_put_att( ncid, w_varid, "units","[m/s]" ) )
      CALL Check( nf90_put_att( ncid, w_varid, "_FillValue",fillValue) )
      CALL Check( nf90_put_att( ncid, w_varid,"coordinates", "latitude longitude z" ) )

      CALL Check( nf90_def_var( ncid, "O", NF90_PREC,&
                               (/ x_dimid, y_dimid, z_dimid, time_dimid /),&
                               on_varid ) )
      CALL Check( nf90_put_att( ncid, on_varid, "long_name","Oxygen density" ) )
      CALL Check( nf90_put_att( ncid, on_varid, "units","[#/m^3]" ) )
      CALL Check( nf90_put_att( ncid, on_varid, "_FillValue",fillValue) )
      CALL Check( nf90_put_att( ncid, on_varid,"coordinates", "latitude longitude z" ) )

      CALL Check( nf90_def_var( ncid, "N2", NF90_PREC,&
                               (/ x_dimid, y_dimid, z_dimid, time_dimid /),&
                               n2n_varid ) )
      CALL Check( nf90_put_att( ncid, n2n_varid, "long_name","Molecular Nitrogen density" ) )
      CALL Check( nf90_put_att( ncid, n2n_varid, "units","[#/m^3]" ) )
      CALL Check( nf90_put_att( ncid, n2n_varid, "_FillValue",fillValue) )
      CALL Check( nf90_put_att( ncid, n2n_varid,"coordinates", "latitude longitude z" ) )

      CALL Check( nf90_def_var( ncid, "O2", NF90_PREC,&
                               (/ x_dimid, y_dimid, z_dimid, time_dimid /),&
                               o2n_varid ) )
      CALL Check( nf90_put_att( ncid, o2n_varid, "long_name","Molecular Oxygen density" ) )
      CALL Check( nf90_put_att( ncid, o2n_varid, "units","[#/m^3]" ) )
      CALL Check( nf90_put_att( ncid, o2n_varid, "_FillValue",fillValue) )
      CALL Check( nf90_put_att( ncid, o2n_varid,"coordinates", "latitude longitude z" ) )
      ! End the Define Mode
      CALL Check( nf90_enddef(ncid) )
      ! 
      CALL Check( nf90_put_var( ncid, z_varid, z ) )
      CALL Check( nf90_put_var( ncid, y_varid, y ) )
      CALL Check( nf90_put_var( ncid, x_varid, x ) )
      CALL Check( nf90_put_var( ncid, time_varid, 0.0 ) )

      CALL Check( nf90_put_var( ncid, oplus_varid, oplus_fixed ) )
      CALL Check( nf90_put_var( ncid, hplus_varid, hplus_fixed ) )
      CALL Check( nf90_put_var( ncid, heplus_varid, heplus_fixed ) )
      CALL Check( nf90_put_var( ncid, nplus_varid, nplus_fixed ) )
      CALL Check( nf90_put_var( ncid, noplus_varid, noplus_fixed ) )
      CALL Check( nf90_put_var( ncid, o2plus_varid, o2plus_fixed ) )
      CALL Check( nf90_put_var( ncid, n2plus_varid, n2plus_fixed ) )
      CALL Check( nf90_put_var( ncid, edens_varid, electron_density_fixed ) )
      CALL Check( nf90_put_var( ncid, tec_varid, total_electron_content ) )
      CALL Check( nf90_put_var( ncid, nmf2_varid, nmf2 ) )
      CALL Check( nf90_put_var( ncid, hmf2_varid, hmf2 ) )
      CALL Check( nf90_put_var( ncid, tn_varid, tn_fixed ) )
      CALL Check( nf90_put_var( ncid, u_varid, vn_fixed(:,:,:,1) ) )
      CALL Check( nf90_put_var( ncid, v_varid, vn_fixed(:,:,:,2) ) )
      CALL Check( nf90_put_var( ncid, w_varid, vn_fixed(:,:,:,3) ) )
      CALL Check( nf90_put_var( ncid, on_varid, on_fixed ) )
      CALL Check( nf90_put_var( ncid, n2n_varid, n2n_fixed ) )
      CALL Check( nf90_put_var( ncid, o2n_varid, o2n_fixed ) )

      CALL Check( nf90_close( ncid ) )

 END SUBROUTINE WriteToNetCDF
!
 SUBROUTINE Check(status)
   IMPLICIT NONE
   INTEGER, INTENT (in) :: status
    
    IF(status /= nf90_noerr) THEN 
      PRINT *, trim(nf90_strerror(status))
      STOP "NetCDF Error, Stopped"
    ENDIF
END SUBROUTINE Check 
END PROGRAM IPEtoHeightGrid

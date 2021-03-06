! ModelDataInstances_Class.f90
! 
!  Copyright (2017), Joseph Schoonover, Cooperative Institute for Research in
!  Environmental Sciences, NOAA, (joseph.schoonover@noaa.gov)
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

MODULE ModelDataInstances_Class

!src/common/
USE Module_Precision

IMPLICIT NONE

   INTEGER, PARAMETER      :: strLen = 30
   INTEGER, PARAMETER      :: ioChunkSize=100000 ! The number of array elements to write in a single binary file record

   CHARACTER(6), PARAMETER :: strFMT = "(A30)"

   TYPE ModelDataInstance
      CHARACTER(strLen)               :: moduleName
      CHARACTER(strLen)               :: subroutineName
      CHARACTER(strLen)               :: statusCheckName
      INTEGER                         :: lineNumber
      INTEGER                         :: instanceID
      INTEGER                         :: nObs
      INTEGER                         :: arraySize
      REAL(real_prec), ALLOCATABLE    :: array(:) 
      TYPE(ModelDataInstance),POINTER :: next

   END TYPE ModelDataInstance

   TYPE ModelDataInstances
      TYPE(ModelDataInstance), POINTER :: head, current, tail

      CONTAINS

      PROCEDURE :: Build       => Build_ModelDataInstances
      PROCEDURE :: Trash       => Trash_ModelDataInstances 
      PROCEDURE :: AddInstance => AddInstance_ModelDataInstances 

      PROCEDURE :: Update              => Update_ModelDataInstances
      PROCEDURE :: SetNames            => SetNames_ModelDataInstances
      PROCEDURE :: GetNames            => GetNames_ModelDataInstances
      PROCEDURE :: PointToInstance     => PointToInstance_ModelDataInstances
      PROCEDURE :: ThereAreNoInstances
      PROCEDURE :: CompareWith         => CompareWith_ModelDataInstances
      PROCEDURE :: CalculateStorageCost => CalculateStorageCost_ModelDataInstances

        
      PROCEDURE :: Write_ModelDataInstances
      PROCEDURE :: Read_ModelDataInstances

    END TYPE ModelDataInstances

 INTEGER, PRIVATE, PARAMETER :: defaultNameLength = 40
 CONTAINS
!
!
!==================================================================================================!
!---------------------------- CONSTRUCTOR/DESTRUCTOR ROUTINES -------------------------------------!
!==================================================================================================!
!
 SUBROUTINE Build_ModelDataInstances( theInstances  )

   IMPLICIT NONE
   CLASS( ModelDataInstances ) :: theInstances
   
     theInstances % head => NULL( )
     ! Point the tail to null
     theInstances % tail => NULL()
     ! Set the current position to Null
     theInstances % current => NULL( )


 END SUBROUTINE Build_ModelDataInstances
!
 SUBROUTINE Trash_ModelDataInstances( theInstances )

   IMPLICIT NONE
   CLASS( ModelDataInstances ) :: theInstances
   ! LOCAL 
   TYPE( ModelDataInstance ), POINTER :: pNext
   
     ! Set the current position of the list to the head
     theInstances % current => theInstances % head
     
     ! Scroll through the list until the current position is nullified
     DO WHILE ( ASSOCIATED( theInstances % current ) )

        ! temporarily point to the next in the list
        pNext => theInstances % current % next 

        ! Deallocate memory pointed to by current position
        DEALLOCATE( theInstances % current % array )
        DEALLOCATE( theInstances % current ) 

        ! Update current position
        theInstances % current => pNext 

     ENDDO
      

 END SUBROUTINE Trash_ModelDataInstances
!
!
!==================================================================================================!
!---------------------------------- ACCESSOR ROUTINES ---------------------------------------------!
!==================================================================================================!
!
! 
 SUBROUTINE SetNames_ModelDataInstances( theInstances, moduleName, subroutineName, statusCheckName, instanceID )

   IMPLICIT NONE
   CLASS( ModelDataInstances ), INTENT(inout) :: theInstances
   CHARACTER(*), INTENT(IN)                   :: moduleName, subroutineName, statusCheckName
   INTEGER, OPTIONAL, INTENT(IN)              :: instanceID
      
      IF( PRESENT( instanceID ) )THEN
         CALL theInstances % PointToInstance( instanceID )  
      ENDIF
      
      theInstances % current % moduleName      = TRIM(moduleName)
      theInstances % current % subroutineName  = TRIM(subroutineName)
      theInstances % current % statusCheckName = TRIM(statusCheckName)

 END SUBROUTINE SetNames_ModelDataInstances
!
 SUBROUTINE GetNames_ModelDataInstances( theInstances, moduleName, subroutineName, statusCheckName, instanceID )

   IMPLICIT NONE
   CLASS( ModelDataInstances ), INTENT(in) :: theInstances
   CHARACTER(*), INTENT(out)               :: moduleName, subroutineName, statusCheckName
   INTEGER, OPTIONAL, INTENT(IN)           :: instanceID
      
      IF( PRESENT( instanceID ) )THEN
         CALL theInstances % PointToInstance( instanceID )
      ENDIF
      
      moduleName      = theInstances % current % moduleName
      subroutineName  = theInstances % current % subroutineName
      statusCheckName = theInstances % current % statusCheckName

 END SUBROUTINE GetNames_ModelDataInstances
!
 SUBROUTINE Update_ModelDataInstances( theInstances, &
                                       moduleName, &
                                       subroutineName, &
                                       statusCheckName, &
                                       lineNumber, &
                                       arraySize, &
                                       array )
   IMPLICIT NONE
   CLASS( ModelDataInstances ), INTENT(inout) :: theInstances
   CHARACTER(*), INTENT(in)                   :: moduleName
   CHARACTER(*), INTENT(in)                   :: subroutineName
   CHARACTER(*), INTENT(in)                   :: statusCheckName
   INTEGER, INTENT(in)                        :: lineNumber, arraySize
   REAL(real_prec), INTENT(in)                :: array(1:arraySize)
   ! Local
   LOGICAL :: success

      CALL theInstances % PointToInstance( CharToIntHashFunction(statusCheckName), success )
      IF( success ) THEN
         theInstances % current % array = array
         theInstances % current % nObs  = theInstances % current % nObs + 1  
      ELSE
         CALL theInstances % AddInstance( moduleName, &
                                          subroutineName, &
                                          statusCheckName, &
                                          lineNumber, &
                                          arraySize, &
                                          array )
        ! PRINT*, 'ModelDataInstances_Class.f90 : Update_ModelDataInstances '
        ! PRINT*, 'Instance "'//TRIM(statusCheckName)//'" not found. STOPPING!'
        ! STOP
      ENDIF
 END SUBROUTINE Update_ModelDataInstances
!
!==================================================================================================!
!-------------------------------- Linked-List Type Operations -------------------------------------!
!==================================================================================================!
!
! 
 FUNCTION ThereAreNoInstances( theInstances ) RESULT( TorF )
  IMPLICIT NONE
  CLASS( ModelDataInstances ) :: theInstances
  LOGICAL              :: TorF

     TorF = .NOT.( ASSOCIATED( theInstances % head  ) )
     
 END FUNCTION ThereAreNoInstances
!
 SUBROUTINE AddInstance_ModelDataInstances( theInstances, &
                                            moduleName, &
                                            subroutineName, &
                                            statusCheckName, &
                                            lineNumber,& 
                                            arraySize, &
                                            array )
 
   IMPLICIT NONE
   CLASS( ModelDataInstances ), INTENT(inout) :: theInstances
   CHARACTER(*)                               :: moduleName
   CHARACTER(*)                               :: subroutineName
   CHARACTER(*)                               :: statusCheckName
   INTEGER                                    :: lineNumber, arraySize
   REAL(real_prec), OPTIONAL                  :: array(1:arraySize)
   ! LOCAL
   INTEGER :: allocationStatus
   TYPE( ModelDataInstance ), POINTER :: pNext

     ! Check to see if this list is empty
     IF( theInstances % ThereAreNoInstances() )THEN
     
        ALLOCATE( pNext, STAT = allocationStatus )
        IF( allocationStatus /=0 )THEN
           PRINT*, 'MODULE ModelDataInstances_Class.f90 : S/R AddInstance : Memory not allocated for next entry in list.'
        ENDIF      
      
        theInstances % head => pNext
        theInstances % tail => pNext
        ! Point the current position to the head
        theInstances % current => pNext
        ! Set the data
        CALL theInstances % SetNames( moduleName, subroutineName, statusCheckName  )
        
        theInstances % current % instanceID = CharToIntHashFunction( theInstances % current % statusCheckName )
        theInstances % current % nObs       = 1
        theInstances % current % arraySize  = arraySize
        ALLOCATE( theInstances % current % array(1:arraySize) )
        IF( PRESENT( array ) )THEN
          theInstances % current % array      = array
        ELSE
          theInstances % current % array      = 0.0_real_prec
        ENDIF        

        ! Point the next to null and the tail to current
        theInstances % current % next => NULL( )

     ELSE ! the list is not empty
    
        ! Then we allocate space for the next item in the list    
        ALLOCATE( pNext, STAT = allocationStatus )
        IF( allocationStatus /=0 )THEN
           PRINT*, 'MODULE ModelDataInstances_Class.f90 : S/R AddInstance : Memory not allocated for next entry in list.'
        ENDIF      

        ! Temporarily point to the tail
        theInstances % current => theInstances % tail
        ! Set the "next" pointer 
        theInstances % current % next => pNext

        ! Reassign the tail
        theInstances % tail    => pNext
        theInstances % current => pNext
        theInstances % current % next => NULL( )
  
        ! Fill in the data
        CALL theInstances % SetNames( moduleName, subroutineName, statusCheckName  )
        
        ! Fill in the key information
        theInstances % current % instanceID = CharToIntHashFunction( theInstances % current % statusCheckName )
        theInstances % current % nObs       = 1
        theInstances % current % arraySize  = arraySize
        ALLOCATE( theInstances % current % array(1:arraySize) )
        IF( PRESENT( array ) )THEN
          theInstances % current % array      = array
        ELSE
          theInstances % current % array      = 0.0_real_prec
        ENDIF        
        
        
     ENDIF

 END SUBROUTINE AddInstance_ModelDataInstances
!
 SUBROUTINE PointToInstance_ModelDataInstances( theInstances, instanceID, instanceFound )
 
   IMPLICIT NONE
   CLASS( ModelDataInstances )    :: theInstances
   INTEGER, INTENT(IN)            :: instanceID
   LOGICAL, OPTIONAL, INTENT(OUT) :: instanceFound
   
   
      theInstances % current => theInstances % head ! Point to the head of the list
      IF( PRESENT( instanceFound ) )THEN
        instanceFound = .FALSE.
      ENDIF
      
      DO WHILE(ASSOCIATED(theInstances % current))
      
         IF( theInstances % current % instanceID == instanceID )THEN
            IF( PRESENT( instanceFound ) )THEN
              instanceFound = .TRUE.
            ENDIF
            EXIT
         ENDIF
         theInstances % current => theInstances % current % next
       
      ENDDO
      
      
 END SUBROUTINE PointToInstance_ModelDataInstances
!
 SUBROUTINE CompareWith_ModelDataInstances( theInstances, otherInstances )
   IMPLICIT NONE
   CLASS( ModelDataInstances ), INTENT(inout) :: theInstances
   TYPE( ModelDataInstances ), INTENT(inout)  :: otherInstances
   ! Local
   LOGICAL         :: instanceFound   
   REAL(real_prec) :: relDiff, thisMag, scaleFac
   INTEGER         :: i

      ! Rewind the main list
      theInstances % current => theInstances % head 

      PRINT*, ' Comparing Instances '
      PRINT*, '-----------------------------------------------------'
      DO WHILE( ASSOCIATED( theInstances % current ) )

         PRINT*, '  Check Name : '//TRIM(theInstances % current % statusCheckName) 
         PRINT '(A,1x,I20)', '   Check ID   : ',theInstances % current % instanceID 
         PRINT*, '  Module     : '//TRIM(theInstances % current % moduleName) 
         PRINT*, '  Subroutine : '//TRIM(theInstances % current % subroutineName) 
         ! Point the "otherInstances" to the record with the same hash ID
         CALL otherInstances % PointToInstance( theInstances % current % instanceID, &
                                                instanceFound )
         IF( instanceFound )THEN

            ! When an instance is found, we need to compare the two arrays that
            ! are stored within
            IF( theInstances % current % arraySize == otherInstances % current % arraySize ) THEN

               relDiff = 0.0_real_prec
               thisMag = 0.0_real_prec
               scaleFac = 1.0_real_prec
               DO i = 1, theInstances % current % arraySize 
                  relDiff = relDiff + scaleFac*( theInstances % current % array(i) - &
                                                 otherInstances % current % array(i) )**2
 
                  thisMag = thisMag + scaleFac*( 0.5_real_prec*( theInstances % current % array(i) + &
                                                         otherInstances % current % array(i) ) )**2
                  IF( thisMag > 10.0_real_prec**6 )THEN
                     scaleFac = scaleFac*10.0_real_prec**(-6)
                     relDiff  = relDiff*10.0_real_prec**(-6)
                     thisMag  = thisMag*10.0_real_prec**(-6)
                  ENDIF
               ENDDO

               IF( sqrt(thisMag) <= 10.0_real_prec**(-10) )THEN
                  PRINT '(A)', '  >>> Solution magnitude is small.'
               ENDIF
               IF( scaleFac /= 1.0_real_prec )THEN
                  PRINT '(A)', '  >>> Solution magnitude is large '
                  PRINT '(A,2x,E11.4)', '  >>> Applied scale factor :', scaleFac 
               ENDIF 

               PRINT '(A,2x,E11.4)', '  >>> Solution magnitude : ', sqrt( thisMag ) 
               PRINT '(A,2x,E11.4)', '  >>> Absolute Difference : ', sqrt( relDiff ) 
               PRINT '(A,2x,E11.4," %")', '  >>> Relative Difference : ', sqrt( relDiff )/sqrt( thisMag )*100.0_real_prec 

            ELSE
               PRINT*, '  >>> *Instance found, but array sizes do not match.'
               PRINT*, '  >>> *Check your instrumentation for this instance.'
            ENDIF

         ELSE

            PRINT*, '  >>> *No instance found for check name :'//TRIM(theInstances % current % statusCheckName)

         ENDIF

         PRINT*, '-----------------------------------------------------'
         theInstances % current => theInstances % current % next
      ENDDO
    

 END SUBROUTINE CompareWith_ModelDataInstances
!
 SUBROUTINE CalculateStorageCost_ModelDataInstances( theInstances )
    IMPLICIT NONE
    CLASS( ModelDataInstances ), INTENT(inout) :: theInstances
    ! Local
    INTEGER :: storageCost
 
       theInstances % current => theInstances % head
       storageCost = 0

       DO WHILE( ASSOCIATED( theInstances % current ) )
          storageCost = storageCost + theInstances % current %arraySize*real_prec
          theInstances % current => theInstances % current % next
       ENDDO 

       PRINT '(A,1x,E11.4,1x,A)', "   MDI Storage Cost : ", &
                                 REAL(storageCost,real_prec)/10.0_real_prec**9,"GB"

 END SUBROUTINE CalculateStorageCost_ModelDataInstances
!
 SUBROUTINE Write_ModelDataInstances( theInstances, baseFileName )

   IMPLICIT NONE
   CLASS( ModelDataInstances ), INTENT(INOUT) :: theInstances 
   CHARACTER(*), INTENT(IN)                   :: baseFileName
   ! LOCAL
   INTEGER :: k, fUnit, fUnit2, recID, i, chunkSizeUsed, rStart, nChunks
   CHARACTER(3) :: countChar 
   REAL(real_prec) :: bufferArray(1:ioChunkSize)

      theInstances % current => theInstances % head
      WRITE( countChar, '(I3.3)' ) theInstances % current % nObs
      OPEN( UNIT = NewUnit(fUnit), &
            FILE = TRIM(baseFileName)//'.mdi.hdr', &
            FORM = 'FORMATTED', &
            ACCESS = 'SEQUENTIAL', &
            STATUS = 'REPLACE', &
            ACTION = 'WRITE' )
            
      OPEN( UNIT = NewUnit(fUnit2), &
            FILE = TRIM(baseFileName)//'.'//countChar//'.mdi', &
            FORM = 'BINARY', &
            ACCESS = 'STREAM', &
            ACTION = 'WRITE', &
            STATUS = 'REPLACE', &
            CONVERT = 'BIG_ENDIAN')

      k     = 0
      recID = 0
      DO WHILE( ASSOCIATED(theInstances % current) )
         k = k+1


         WRITE(fUnit,*) TRIM( theInstances % current % moduleName )
         WRITE(fUnit,*) TRIM( theInstances % current % subroutineName )
         WRITE(fUnit,*) TRIM( theInstances % current % statusCheckName )
         WRITE(fUnit,*) theInstances % current % lineNumber
         WRITE(fUnit,*) theInstances % current % arraySize
         WRITE(fUnit,*) theInstances % current % instanceID
         WRITE(fUnit,*) '------------------------------------------------------------'

         rStart = 1
         nChunks = theInstances % current % arraySize/ioChunkSize
         IF( nChunks*ioChunkSize < theInstances % current % arraySize )THEN
           nChunks = nChunks + 1
         ENDIF
         DO i = 1, nChunks
         
            chunkSizeUsed = MIN( ioChunkSize, theInstances % current % arraySize - rStart )
          
            bufferArray(1:ioChunkSize)   = 0.0_real_prec
            bufferArray(1:chunkSizeUsed) = theInstances % current % array(rStart:rStart+chunkSizeUsed-1) 

            recID = recID + 1

            WRITE( fUnit2 ) bufferArray(1:ioChunkSize)

            rStart = rStart + chunkSizeUsed

         ENDDO
         theInstances % current => theInstances % current % next

      ENDDO
      CLOSE(fUnit)
      CLOSE(fUnit2)


 END SUBROUTINE Write_ModelDataInstances
!
 SUBROUTINE Read_ModelDataInstances( theInstances, baseFileName, obsCount, fileExists )

   IMPLICIT NONE
   CLASS( ModelDataInstances ), INTENT(INOUT) :: theInstances 
   CHARACTER(*), INTENT(IN)                   :: baseFileName
   INTEGER, INTENT(in)                        :: obsCount
   LOGICAL, INTENT(out)                       :: fileExists
   ! LOCAL
   INTEGER           :: k, fUnit, fUnit2, recID, i, ioErr, rStart, chunkSizeUsed, nChunks
   CHARACTER(3)      :: countChar 
   CHARACTER(strLen) :: moduleName, subroutineName, statusCheckName, dummyChar
   INTEGER           :: lineNumber, arraySize, instanceID
   REAL(real_prec)   :: bufferArray(1:ioChunkSize)

      WRITE( countChar, '(I3.3)' ) obsCount

      INQUIRE( FILE = TRIM(baseFileName)//'.'//countChar//'.mdi', &
               EXIST= fileExists )
      IF( .NOT.( fileExists ) )THEN
         PRINT*, 'ModelDataInstances_Class.f90 : Read_ModelDataInstances '
         PRINT*, 'File '//TRIM(baseFileName)//'.'//countChar//'.mdi not found.'
         RETURN
      ENDIF

      OPEN( UNIT = NewUnit(fUnit), &
            FILE = TRIM(baseFileName)//'.mdi.hdr', &
            FORM = 'FORMATTED', &
            ACCESS = 'SEQUENTIAL', &
            STATUS = 'OLD', &
            ACTION = 'READ' )
            
      OPEN( UNIT = NewUnit(fUnit2), &
            FILE = TRIM(baseFileName)//'.'//countChar//'.mdi', &
            FORM = 'BINARY', &
            ACCESS = 'STREAM', &
            ACTION = 'READ', &
            STATUS = 'OLD', &
            CONVERT = 'BIG_ENDIAN' )

      k     = 0
      recID = 0
      ioErr = 0
      DO WHILE( ioErr == 0 )
         k = k+1


         READ(fUnit,strFMT,IOSTAT=ioErr) moduleName 
         IF( ioErr < 0 ) EXIT
         READ(fUnit,strFMT) subroutineName
         READ(fUnit,strFMT) statusCheckName
         READ(fUnit,*) lineNumber
         READ(fUnit,*) arraySize
         READ(fUnit,*) instanceID
         READ(fUnit,strFMT) dummyChar

         IF( obsCount == 1 )THEN
           CALL theInstances % AddInstance( moduleName, &
                                            subroutineName, &
                                            statusCheckName, &
                                            lineNumber, &
                                            arraySize )
         ELSE
            CALL theInstances % PointToInstance( CharToIntHashFunction(statusCheckName) )
         ENDIF             

       

         rStart = 1
         nChunks = theInstances % current % arraySize/ioChunkSize
         IF( nChunks*ioChunkSize < theInstances % current % arraySize )THEN
           nChunks = nChunks + 1
         ENDIF
         DO i = 1 , nChunks 
         
            chunkSizeUsed = MIN( ioChunkSize, theInstances % current % arraySize - rStart )

            recID = recID + 1

            READ( fUnit2 ) bufferArray(1:ioChunkSize)

            theInstances % current % array(rStart:rStart+chunkSizeUsed-1) =bufferArray(1:chunkSizeUsed)
            rStart = rStart + chunkSizeUsed
         ENDDO

      ENDDO
      CLOSE(fUnit)
      CLOSE(fUnit2)


 END SUBROUTINE Read_ModelDataInstances
!
 FUNCTION CharToIntHashFunction( inputChar ) RESULT( hash )
    IMPLICIT NONE
    CHARACTER(*) :: inputChar
    INTEGER      :: hash
    ! Local
   ! CHARACTER(strLen) :: localChar
    INTEGER           :: i

      ! localChar = UpperCase( inputChar )
      
       hash = 5381
       DO i = 1, LEN_TRIM(inputChar)
          hash = ( ISHFT(hash,5)+hash ) + ICHAR( inputChar(i:i) )
       ENDDO

 END FUNCTION CharToIntHashFunction
!
 FUNCTION UpperCase( str ) RESULT( upper )

    IMPLICIT NONE
    CHARACTER(strLen), INTENT(IN) :: str
    CHARACTER(strLen)           :: Upper

    INTEGER :: ic, i

    CHARACTER(27), PARAMETER :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ '
    CHARACTER(27), PARAMETER :: low = 'abcdefghijklmnopqrstuvwxyz '

    DO i = 1, LEN_TRIM(str)
        ic = INDEX(low, str(i:i))
        IF (ic > 0)THEN
         Upper(i:i) = cap(ic:ic)
        ELSE
         Upper(i:i) = str(i:i)
        ENDIF
    ENDDO

 END FUNCTION UpperCase
!
 INTEGER FUNCTION NewUnit(thisunit)
 
   IMPLICIT NONE
   INTEGER, INTENT(out), OPTIONAL :: thisunit
   ! Local
   INTEGER, PARAMETER :: unitMin=100, unitMax=1000
   LOGICAL :: isopened
   INTEGER :: iUnit

     newunit=-1

     DO iUnit=unitMin, unitMax
        ! Check to see IF this UNIT is opened
        INQUIRE(UNIT=iUnit,opened=isopened)
        IF( .not. isopened )THEN
           newunit=iUnit
           EXIT
        ENDIF
     ENDDO

     IF (PRESENT(thisunit)) thisunit = newunit
 
END FUNCTION NewUnit
END MODULE ModelDataInstances_Class

! ModelDataInstances_Driver.f90
! 
!  Copyright (2017), Joseph Schoonover, Cooperative Institute for Research in
!  Environmental Sciences, NOAA, (joseph.schoonover@noaa.gov)
!
!  Author : Joseph Schoonover
!  Creation Date : July 25, 2017
!
! ------------------------------------------------------------------------------------------------ !
!
!  This example driver tests the functionality of the
!  ModelDataInstances_Class.f90 module and provides examples the shows how to
!  make use of the ModelDataInstances data structure and type-bound procedures
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

PROGRAM ModelDataInstances_Driver

 USE Module_Precision
 USE ModelDataInstances_Class


   TYPE(ModelDataInstances) :: mdi
   REAL(real_prec)          :: sampleArray1(1:100)       
   REAL(real_prec)          :: sampleArray2(1:20,1:20)       
   REAL(real_prec)          :: sampleArray3(1:10,1:10,1:10)       
   INTEGER                  :: lineNumber

     ! Initialize the ModelDataInstances linked list, by pointing all of the
     ! pointers to Null
     CALL mdi % Build( )

     ! Here, the "sampleArray's" are set to initialized values
     sampleArray1 = 0.0_real_prec
     sampleArray2 = 0.0_real_prec
     sampleArray3 = 0.0_real_prec

     ! Grab the line number (if desired)
     lineNumber = 39
     ! Here, we add instances for logging the initial condition
     ! All strings that are passed must be 30 characters or less.
     !
     ! This section of code shows how to use the ModelDataInstances data
     ! structure to add new instances for storing the state of single and
     ! multidimensional arrays.
     ! When "AddInstance" is called, more space is allocated to store the array
     ! and for keeping track of this unique instance. The array is stored within
     ! the data structure until the user calls the "Trash" routine.



     ! In practice, you will likely do some work that changes the values stored
     ! in each of the arrays. Here, the call to "RANDOM_NUMBER" is a proxy for
     ! doing such work.

     CALL RANDOM_NUMBER( sampleArray1 ) ! Fill SampleArray1 with random data
     CALL mdi % AddInstance( "ModelDataInstances_Driver.f90", & ! Module name -here we use the main program name 
                             "Main", &                          ! Subroutine name - here we use "Main" indicating we are not in a subroutine
                             "Randomize Array1", &              ! Unique name of the status check. Pick something descriptive of what it is.
                             lineNumber, &                      ! Line number near where this instance was added
                             SIZE(sampleArray1), &              ! The total number of elements in the array
                             sampleArray1 )                     ! The array that you want to capture

     CALL RANDOM_NUMBER( sampleArray2 ) ! Fill SampleArray2 with random data
     CALL mdi % AddInstance( "ModelDataInstances_Driver.f90", & ! Module name -here we use the main program name 
                             "Main", &                          ! Subroutine name - here we use "Main" indicating we are not in a subroutine
                             "Randomize Array2", &              ! Unique name of the status check. Pick something descriptive of what it is.
                             lineNumber, &                      ! Line number near where this instance was added
                             SIZE(sampleArray2), &              ! The total number of elements in the array
                             PACK(sampleArray2,.TRUE.) )        ! The array that you want to capture

     CALL RANDOM_NUMBER( sampleArray3 ) ! Fill SampleArray3 with random data
     CALL mdi % AddInstance( "ModelDataInstances_Driver.f90", & ! Module name -here we use the main program name 
                             "Main", &                          ! Subroutine name - here we use "Main" indicating we are not in a subroutine
                             "Randomize Array3", &              ! Unique name of the status check. Pick something descriptive of what it is.
                             lineNumber, &                      ! Line number near where this instance was added
                             SIZE(sampleArray3), &              ! The total number of elements in the array
                             PACK(sampleArray3,.TRUE.) )        ! The array that you want to capture

     ! When writing to file, a file "base name" needs to be passed.
     ! This call generates, here, test.mdi.hdr and test.0001.mdi
     CALL mdi % Write_ModelDataInstances( 'test' )  

     ! It is possible that a particular instance falls
     ! within a do-loop. The type-bound procedure "Update" keeps track of the
     ! count. Keep in mind that, unless written to file between updates, stale
     ! data can be lost.
     !
     ! In the current implementation, it is much safer to only track data that
     ! are updated the same number of times
     DO i = 1, 3
       CALL RANDOM_NUMBER( sampleArray1 ) ! Fill SampleArray1 with random data

       CALL mdi % Update( "ModelDataInstances_Driver.f90", & ! Module name -here we use the main program name 
                          "Main", &                          ! Subroutine name - here we use "Main" indicating we are not in a subroutine
                          "Randomize Array1", &              ! Unique name of the status check. Pick something descriptive of what it is.
                          lineNumber, &                      ! Line number near where this instance was added
                          SIZE(sampleArray1), &              ! The total number of elements in the array
                          sampleArray1 )        ! The array that you want to capture

       CALL RANDOM_NUMBER( sampleArray2 ) ! Fill SampleArray2 with random data
       CALL mdi % Update( "ModelDataInstances_Driver.f90", & ! Module name -here we use the main program name 
                          "Main", &                          ! Subroutine name - here we use "Main" indicating we are not in a subroutine
                          "Randomize Array2", &              ! Unique name of the status check. Pick something descriptive of what it is.
                          lineNumber, &                      ! Line number near where this instance was added
                          SIZE(sampleArray2), &              ! The total number of elements in the array
                          PACK(sampleArray2,.TRUE.) )        ! The array that you want to capture

       CALL RANDOM_NUMBER( sampleArray3 ) ! Fill SampleArray3 with random data
       CALL mdi % Update( "ModelDataInstances_Driver.f90", & ! Module name -here we use the main program name 
                          "Main", &                          ! Subroutine name - here we use "Main" indicating we are not in a subroutine
                          "Randomize Array3", &              ! Unique name of the status check. Pick something descriptive of what it is.
                          lineNumber, &                      ! Line number near where this instance was added
                          SIZE(sampleArray3), &              ! The total number of elements in the array
                          PACK(sampleArray3,.TRUE.) )        ! The array that you want to capture

       ! To avoid losing track of data on each loop, we write the data structure
       ! to file each time around.
       CALL mdi % Write_ModelDataInstances( 'test' )  
     ENDDO

     CALL mdi % CalculateStorageCost( )


     ! Here, we clear all of the memory held by the data structure.
     ! If this driver is run underneath valgrind, there should be no memory
     ! leaks and no lost or unreachable data.
     CALL mdi % Trash( )

END PROGRAM ModelDataInstances_Driver

MODULE moduleDriverDebug

IMPLICIT NONE
!-----------------------------
! debug Interpolation files
!-----------------------------

CHARACTER(LEN=*), PARAMETER :: debugThermoInterpFileName = 'interpOut.dat'
INTEGER, parameter :: unitCheckThermoInterpBefore = 7600  ! for writing out thermosphere values
INTEGER, parameter :: unitCheckThermoInterpAfter = 7500  ! for writing out thermosphere values
INTEGER, parameter :: unitCheckThermoInterpFluxT = 7400  ! for writing out thermosphere values
INTEGER, parameter :: unitCheckPressureHeightThermo = 8800
INTEGER, parameter :: unitPressureGridHeight = 7810
!-----------------------------------------------------------------------------
! File unit number for checking Ionospheric values to the pressure grid
! and file unit number for getting the average height for each pressure level
!-----------------------------------------------------------------------------
INTEGER, parameter :: unitCheckPressureInterp = 8900
INTEGER, parameter :: unitCheckPressureHeight = 7820


PRIVATE

PUBLIC :: checkGridIPE
PUBLIC :: check2GridIPE
PUBLIC :: checkThermo
PUBLIC :: checkThermoArrays
PUBLIC :: checkInterp
PUBLIC :: checkFixedGeo
PUBLIC :: openInterpFiles
PUBLIC :: openPressureInterpFiles
PUBLIC :: writeInterpThermo
PUBLIC :: writeThermo
PUBLIC :: writeThermoFixed
PUBLIC :: writeThermoToIono

CONTAINS

!==========================================================================

Function openInterpFiles(debugDir) result(j)

   IMPLICIT NONE

   character(*), intent(in) :: debugDir
   INTEGER :: j ! output

   ! Begin code ------------------------

   !--------------------------------------------------------------
   ! Open files for writing out interpreted variables b/f & after
   !--------------------------------------------------------------
       OPEN (unitCheckThermoInterpBefore, FILE=TRIM(debugDir)//TRIM('SouthWindBefore.txt'), STATUS='REPLACE', ERR=100)
       OPEN (unitCheckThermoInterpAfter, FILE=TRIM(debugDir)//TRIM('SouthWindAfter.txt'), STATUS='REPLACE', ERR=100)
       OPEN (unitCheckThermoInterpFluxT, FILE=TRIM(debugDir)//TRIM('SouthWindFluxTube.txt'), STATUS='REPLACE', ERR=100)

       OPEN (unitCheckThermoInterpBefore+1, FILE=TRIM(debugDir)//TRIM('EastWindBefore.txt'), STATUS='REPLACE', ERR=100)
       OPEN (unitCheckThermoInterpAfter+1, FILE=TRIM(debugDir)//TRIM('EastWindAfter.txt'), STATUS='REPLACE', ERR=100)
       OPEN (unitCheckThermoInterpFluxT+1, FILE=TRIM(debugDir)//TRIM('EastWindFluxTube.txt'), STATUS='REPLACE', ERR=100)

       OPEN (unitCheckThermoInterpBefore+2, FILE=TRIM(debugDir)//TRIM('WVZBefore.txt'), STATUS='REPLACE')
       OPEN (unitCheckThermoInterpAfter+2, FILE=TRIM(debugDir)//TRIM('WVZAfter.txt'), STATUS='REPLACE')
       OPEN (unitCheckThermoInterpFluxT+2, FILE=TRIM(debugDir)//TRIM('WVZFluxTube.txt'), STATUS='REPLACE')


       OPEN (unitCheckThermoInterpBefore+3, FILE=TRIM(debugDir)//TRIM('RMTBefore.txt'), STATUS='REPLACE')
       OPEN (unitCheckThermoInterpAfter+3, FILE=TRIM(debugDir)//TRIM('RMTAfter.txt'), STATUS='REPLACE')
       OPEN (unitCheckThermoInterpFluxT+3, FILE=TRIM(debugDir)//TRIM('RMTFluxTube.txt'), STATUS='REPLACE')


       OPEN (unitCheckThermoInterpBefore+4, FILE=TRIM(debugDir)//TRIM('TemperatureBefore.txt'), STATUS='REPLACE')
       OPEN (unitCheckThermoInterpAfter+4, FILE=TRIM(debugDir)//TRIM('TemperatureAfter.txt'), STATUS='REPLACE')
       OPEN (unitCheckThermoInterpFluxT+4, FILE=TRIM(debugDir)//TRIM('TemperatureFluxTube.txt'), STATUS='REPLACE')


       OPEN (unitCheckThermoInterpBefore+5, FILE=TRIM(debugDir)//TRIM('OBefore.txt'), STATUS='REPLACE')
       OPEN (unitCheckThermoInterpAfter+5, FILE=TRIM(debugDir)//TRIM('OAfter.txt'), STATUS='REPLACE')
       OPEN (unitCheckThermoInterpFluxT+5, FILE=TRIM(debugDir)//TRIM('OFluxTube.txt'), STATUS='REPLACE')


       OPEN (unitCheckThermoInterpBefore+6, FILE=TRIM(debugDir)//TRIM('O2Before.txt'), STATUS='REPLACE')
       OPEN (unitCheckThermoInterpAfter+6, FILE=TRIM(debugDir)//TRIM('O2After.txt'), STATUS='REPLACE')
       OPEN (unitCheckThermoInterpFluxT+6, FILE=TRIM(debugDir)//TRIM('O2FluxTube.txt'), STATUS='REPLACE')


       OPEN (unitCheckThermoInterpBefore+7, FILE=TRIM(debugDir)//TRIM('N2Before.txt'), STATUS='REPLACE')
       OPEN (unitCheckThermoInterpAfter+7, FILE=TRIM(debugDir)//TRIM('N2After.txt'), STATUS='REPLACE')
       OPEN (unitCheckThermoInterpFluxT+7, FILE=TRIM(debugDir)//TRIM('N2FluxTube.txt'), STATUS='REPLACE')


       OPEN (unitCheckThermoInterpBefore+8, FILE=TRIM(debugDir)//TRIM('QionBefore.txt'), STATUS='REPLACE')
       OPEN (unitCheckThermoInterpAfter+8, FILE=TRIM(debugDir)//TRIM('QionAfter.txt'), STATUS='REPLACE')
       OPEN (unitCheckThermoInterpFluxT+8, FILE=TRIM(debugDir)//TRIM('QionFluxTube.txt'), STATUS='REPLACE')

       OPEN (unitCheckPressureHeightThermo, &
             FILE=TRIM(debugDir)//TRIM('averagePressureHeightGridfromThermo.txt'), STATUS='REPLACE')

       ! Write all heights out to a file 
       OPEN (unitPressureGridHeight, FILE=TRIM(debugDir)//TRIM('PressureHeightGrid.txt'), STATUS='REPLACE')

   j = 0
   RETURN

100 j = -1
    RETURN

END FUNCTION openInterpFiles

!============================================================================================================


Function openPressureInterpFiles(debugDir) result(j)

   IMPLICIT NONE

   character(*), intent(in) :: debugDir
   INTEGER :: j ! output

   ! Begin code -----------------------------------------

   OPEN (unitCheckPressureInterp, FILE=TRIM(debugDir)//TRIM('OplusPressureGrid.txt'), STATUS='REPLACE')
   OPEN (unitCheckPressureInterp+1, FILE=TRIM(debugDir)//TRIM('HplusPressureGrid.txt'), STATUS='REPLACE')
   OPEN (unitCheckPressureInterp+2, FILE=TRIM(debugDir)//TRIM('NplusPressureGrid.txt'), STATUS='REPLACE')
   OPEN (unitCheckPressureInterp+3, FILE=TRIM(debugDir)//TRIM('NOplusPressureGrid.txt'), STATUS='REPLACE')
   OPEN (unitCheckPressureInterp+4, FILE=TRIM(debugDir)//TRIM('O2plusPressureGrid.txt'), STATUS='REPLACE')
   OPEN (unitCheckPressureInterp+5, FILE=TRIM(debugDir)//TRIM('N2plusPressureGrid.txt'), STATUS='REPLACE')
   OPEN (unitCheckPressureInterp+6, FILE=TRIM(debugDir)//TRIM('NePressureGrid.txt'), STATUS='REPLACE')
   OPEN (unitCheckPressureInterp+7, FILE=TRIM(debugDir)//TRIM('TePressureGrid.txt'), STATUS='REPLACE')
   OPEN (unitCheckPressureInterp+8, FILE=TRIM(debugDir)//TRIM('TiPressureGrid.txt'), STATUS='REPLACE')
   OPEN (unitCheckPressureInterp+9, FILE=TRIM(debugDir)//TRIM('nHratePressureGrid.txt'), STATUS='REPLACE')

   OPEN (unitCheckPressureHeight, FILE=TRIM(debugDir)//TRIM('averagePressureHeightGrid.txt'), STATUS='REPLACE')

   j = 0
   RETURN

101 j = -1
    RETURN


END FUNCTION openPressureInterpFiles

! ==================================================================================
FUNCTION writeInterpThermo(time, Oplus, Hplus, Nplus, NOplus, O2plus,&
                           N2plus, Ne, Te, Ti_Oplus, neutralHeatingRates, height)   result(j)

   USE module_precision
   USE modSizeThermo

   IMPLICIT NONE

   !----------------------------------------------
   ! Ionospheric Variables in thermospheric grid
   !----------------------------------------------
   INTEGER (KIND=int_prec), intent(in) :: Time
   REAL*8, intent(in) :: Oplus(GT_ht_dim, GT_lat_dim, GT_lon_dim)
   REAL*8, intent(in) :: Hplus(GT_ht_dim, GT_lat_dim, GT_lon_dim)  ! NOT USED YET lrm20120508
   REAL*8, intent(in) :: N2plus(GT_ht_dim, GT_lat_dim, GT_lon_dim) ! NOT USED YET lrm20120508
   REAL*8, intent(in) :: Te(GT_ht_dim, GT_lat_dim, GT_lon_dim)             ! NOT USED YET lrm20120508
   REAL*8, intent(in) :: Nplus(GT_ht_dim, GT_lat_dim, GT_lon_dim)  ! NOT USED YET lrm20120508
   REAL*8, intent(in) :: Ne(GT_ht_dim, GT_lat_dim, GT_lon_dim)
   REAL*8, intent(in) :: NOplus(GT_ht_dim, GT_lat_dim, GT_lon_dim)
   REAL*8, intent(in) :: O2plus(GT_ht_dim, GT_lat_dim, GT_lon_dim)
   REAL*8, intent(in) :: Ti_Oplus(GT_ht_dim, GT_lat_dim, GT_lon_dim)
   !INTEGER, intent(in) :: numHrate
   REAL(kind=8), intent(in) :: neutralHeatingRates(GT_ht_dim, GT_lat_dim, GT_lon_dim)
   REAL(kind=8), intent(in) :: height(GT_ht_dim, GT_lat_dim, GT_lon_dim)

   INTEGER :: ii
   INTEGER :: j ! output

   ! Begin code -------------------------------------------------------------------

   write(unitCheckPressureInterp,*)     Oplus
   write(unitCheckPressureInterp + 1,*) Hplus
   write(unitCheckPressureInterp + 2,*) Nplus
   write(unitCheckPressureInterp + 3,*) NOplus
   write(unitCheckPressureInterp + 4,*) O2plus
   write(unitCheckPressureInterp + 5,*) N2plus
   write(unitCheckPressureInterp + 6,*) Ne
   write(unitCheckPressureInterp + 7,*) Te
   write(unitCheckPressureInterp + 8,*) Ti_Oplus
   write(unitCheckPressureInterp + 9,*) neutralHeatingRates

   ! write out the pressure grid heights
   write(unitCheckPressureHeight,*) height

   !--------------------------------------------------------------------
   ! Calculate and write out the average height for each pressure level
   !--------------------------------------------------------------------
   write(unitCheckPressureHeight, FMT="(I12)" ) Time
   do ii = 1, GT_ht_dim

      !avgHeightofPressure = Height(ii,:,:)/(GT_lat_dim*GT_lon_dim)

      write(unitCheckPressureHeight, FMT="(E12.4)") SUM(Height(ii,:,:))/(GT_lat_dim*GT_lon_dim)
      !write(unitCheckPressureHeight, FMT="(A1)") " "

   enddo
   j = 0
END FUNCTION writeInterpThermo

!================================================================================================================


! ==================================================================================

FUNCTION writeThermo(time, wind_southwards, wind_eastwards, &
                     wvz, rmt, Temperature, &
                     O_density, O2_density, N2_density, &
                     qion3d, height)   result(j)

   USE module_precision
   USE modSizeThermo, ONLY : GT_ht_dim, GT_lat_dim, GT_lon_dim

   IMPLICIT NONE

   INTEGER (KIND=int_prec), intent(in) :: Time
   REAL(kind=8), intent(in) :: wind_southwards(GT_ht_dim, GT_lat_dim, GT_lon_dim)
   REAL(kind=8), intent(in) :: wind_eastwards(GT_ht_dim, GT_lat_dim, GT_lon_dim)
   REAL(kind=8), intent(in) :: wvz(GT_ht_dim, GT_lat_dim, GT_lon_dim)

   REAL(kind=8), intent(in) :: rmt(GT_ht_dim, GT_lat_dim, GT_lon_dim)
   REAL(kind=8), intent(in) :: Temperature(GT_ht_dim, GT_lat_dim, GT_lon_dim)
   REAL(kind=8), intent(in) :: height(GT_ht_dim, GT_lat_dim, GT_lon_dim)
   REAL(kind=8), intent(in) :: O_density(GT_ht_dim, GT_lat_dim, GT_lon_dim)
   REAL(kind=8), intent(in) :: O2_density(GT_ht_dim, GT_lat_dim, GT_lon_dim)
   REAL(kind=8), intent(in) :: N2_density(GT_ht_dim, GT_lat_dim, GT_lon_dim)
   REAL*8, intent(in) :: qion3d(GT_ht_dim,GT_lat_dim,GT_lon_dim)


   INTEGER :: ii
   INTEGER :: j ! output


   ! Begin code -------------------------------------------

   ! write out arrays to file for plotting
   write(unitCheckThermoInterpBefore,*) wind_southwards
   write(unitCheckThermoInterpBefore + 1,*) wind_eastwards
   write(unitCheckThermoInterpBefore + 2,*) wvz
   write(unitCheckThermoInterpBefore + 3,*) rmt
   write(unitCheckThermoInterpBefore + 4,*) Temperature
   write(unitCheckThermoInterpBefore + 5,*) O_density
   write(unitCheckThermoInterpBefore + 6,*) O2_density
   write(unitCheckThermoInterpBefore + 7,*) N2_density
   write(unitCheckThermoInterpBefore + 8,*) qion3d

   !--------------------------------------------------------------------
   ! Calculate and write out the average height for each pressure level
   !--------------------------------------------------------------------
   write(unitCheckPressureHeightThermo, FMT="(I12)" ) Time
   do ii = 1, GT_ht_dim

        write(unitCheckPressureHeightThermo, FMT="(E12.4)") SUM(Height(ii,:,:))/(GT_lat_dim*GT_lon_dim)

   enddo

   j = 0
END FUNCTION writeThermo

!=======================================================================================

FUNCTION writeThermoFixed(V_South_FixedHeight, V_East_FixedHeight, &
                          V_Upward_FixedHeight, tts_fixed_ht, &
                          O_density_fixed_ht, O2_density_fixed_ht, &
                          N2_density_fixed_ht, qion3d_fixed_ht)  result(j)

USE modSizeFixedGridThermo, ONLY : nFixedGridThermoHeights
USE modSizeThermo, ONLY : GT_lat_dim, GT_lon_dim
IMPLICIT NONE

!-------------------------------- 
! Fixed grid neutral atmosphere
!--------------------------------                                              
REAL(kind=8) :: O_density_fixed_ht(nFixedGridThermoHeights, GT_lat_dim, GT_lon_dim)
REAL(kind=8) :: O2_density_fixed_ht(nFixedGridThermoHeights, GT_lat_dim, GT_lon_dim)
REAL(kind=8) :: N2_density_fixed_ht(nFixedGridThermoHeights, GT_lat_dim, GT_lon_dim)
REAL(kind=8) :: V_East_FixedHeight(nFixedGridThermoHeights, GT_lat_dim, GT_lon_dim)
REAL(kind=8) :: V_South_FixedHeight(nFixedGridThermoHeights, GT_lat_dim, GT_lon_dim)
REAL(kind=8) :: V_Upward_FixedHeight(nFixedGridThermoHeights, GT_lat_dim, GT_lon_dim)
REAL(kind=8) :: tts_fixed_ht(nFixedGridThermoHeights, GT_lat_dim, GT_lon_dim)
REAL(kind=8) :: qion3d_fixed_ht(nFixedGridThermoHeights, GT_lat_dim, GT_lon_dim)
REAL(kind=8) :: elx_fixed_ht(nFixedGridThermoHeights, GT_lat_dim, GT_lon_dim)
REAL(kind=8) :: ely_fixed_ht(nFixedGridThermoHeights, GT_lat_dim, GT_lon_dim)

   INTEGER :: j ! output
! Begin code ----------------------------------------------------------------

      write(unitCheckThermoInterpAfter,*) V_South_FixedHeight
      write(unitCheckThermoInterpAfter+1,*) V_East_FixedHeight
      write(unitCheckThermoInterpAfter+2,*) V_Upward_FixedHeight
      write(unitCheckThermoInterpAfter+3,*) ! should be rmt , but the composition is not being used ***
      write(unitCheckThermoInterpAfter+4,*) tts_fixed_ht
      write(unitCheckThermoInterpAfter+5,*) O_density_fixed_ht
      write(unitCheckThermoInterpAfter+6,*) O2_density_fixed_ht
      write(unitCheckThermoInterpAfter+7,*) N2_density_fixed_ht
      write(unitCheckThermoInterpAfter+8,*) qion3d_fixed_ht
   j = 0
END FUNCTION writeThermoFixed

!=================================================================================

FUNCTION writeThermoToIono(V_South_plasma, V_East_plasma, V_Upward_plasma, &
                           TN_plasma_input_3d, O_plasma_input_3d, O2_plasma_input_3d, &
                           N2_plasma_input_3d) result(j)

USE modSizeFluxTube, ONLY : NPTS, NMP, NLP
IMPLICIT NONE

REAL(kind=8), intent(in) :: V_east_plasma(npts, nmp), &
                V_south_plasma(npts, nmp), &
                V_upward_plasma(npts, nmp)
REAL(kind=8), intent(in) :: TN_plasma_input_3d(npts, nmp), O_plasma_input_3d(npts, nmp), &
                O2_plasma_input_3d(npts, nmp), N2_plasma_input_3d(npts, nmp)

INTEGER :: j


! Begin code --------------------------------------------------

        write(unitCheckThermoInterpFluxT,*) V_South_plasma
        write(unitCheckThermoInterpFluxT+1,*) V_East_plasma
        write(unitCheckThermoInterpFluxT+2,*) V_Upward_plasma
        write(unitCheckThermoInterpFluxT+3,*) ! should be rmt , but the composition is not being used ***
        write(unitCheckThermoInterpFluxT+4,*) TN_plasma_input_3d
        write(unitCheckThermoInterpFluxT+5,*) O_plasma_input_3d
        write(unitCheckThermoInterpFluxT+6,*) O2_plasma_input_3d
        write(unitCheckThermoInterpFluxT+7,*) N2_plasma_input_3d
        write(unitCheckThermoInterpFluxT+8,*) !qion3d_

   j = 0

END FUNCTION writeThermoToIono

!===============================================================================

SUBROUTINE checkGridIPE(directory, fileName, unitNumber, NMP, NLP, &
                        inGIP1d, isGIP1d, mlon_rad)

IMPLICIT NONE

character(*), intent(in) :: directory
character(*), intent(in) :: fileName
INTEGER :: unitNumber

INTEGER, intent(in) :: NLP, NMP
INTEGER, intent(in) :: inGIP1d(NLP)
INTEGER, intent(in) :: isGIP1d(NLP)
REAL (KIND=8), intent(in) :: mlon_rad(1:NMP+1)


!--------------------
! local variables 
!--------------------
!------------------------------------------------------
! File unit number for checking IPE grid values :
!------------------------------------------------------
INTEGER ::  FileOpenStatus  ! for checking if we can open the file
INTEGER :: ii

CHARACTER(100) :: fullFileName

! BEGIN CODE ====================================================================

fullFileName = TRIM(directory)//TRIM(fileName)

!================================================
! Open for checking values from the IPE grid file
!================================================
OPEN (UNIT=unitNumber, FILE=fullFileName, &
                 ACTION="WRITE", IOSTAT=FileOpenStatus)
IF (FileOpenStatus > 0) then
    print *,"checkGridIPE : Cannot open "//fullFileName//" ***********"
    STOP 
ENDIF

print *,'checkGridIPE : writing to ',fullFileName


WRITE (unitNumber, '(A)' ) 'inGIP1d : '
do ii = 1, NLP, 10
   WRITE (unitNumber, '(10I8)' ) inGIP1d(ii:ii+9)
enddo 

WRITE (unitNumber, '(A)' ) ' '
WRITE (unitNumber, '(A)' ) 'isGIP1d : '
do ii = 1, NLP, 10
   WRITE (unitNumber, '(10I8)' ) isGIP1d(ii:ii+9)
enddo 

WRITE (unitNumber, '(A)' ) ' '
WRITE (unitNumber, '(A)' ) 'mlon_rad : '
do ii = 1, NMP, 10
   WRITE (unitNumber, '(10ES15.5)' ) mlon_rad(ii:ii+9)
enddo 
WRITE (unitNumber, '(ES15.5)' ) mlon_rad(NMP+1)



WRITE (unitNumber, '(A,F15.5)') 'mlon_rad(2) = ', mlon_rad(2)
WRITE (unitNumber, '(A)') ' '

END SUBROUTINE checkGridIPE



! =======================================================================================


SUBROUTINE check2GridIPE(unitNumber, NPTS, NMP, &
                         glat_plasma_3d, glond_plasma_3d, pz_plasma_3d)

IMPLICIT NONE

INTEGER :: unitNumber

INTEGER, intent(in) :: NPTS, NMP


REAL(kind=8), intent(in) :: glat_plasma_3d(npts,nmp), glond_plasma_3d(npts,nmp)
REAL(kind=8), intent(in) :: pz_plasma_3d(npts,nmp)


!--------------------
! local variables 
!--------------------

INTEGER :: ii, jj 

! BEGIN CODE ====================================================================

ii = 1118
jj = 1

WRITE (unitNumber, '(A,i4,a1,i4,a4,F15.5)') 'Pz_plasma_3d(', ii, ',' , jj, ') = ', Pz_plasma_3d(ii,jj)

jj = 2
WRITE (unitNumber, '(A,i4,a1,i4,a4,F15.5)') 'GLAT_plasma_3d(' , ii, ',' , jj, ') = ', GLAT_plasma_3d(ii,jj)
WRITE (unitNumber, '(A,i4,a1,i4,a4,F15.5)') 'GLONd_plasma_3d(', ii, ',' , jj, ' ) = ', GLONd_plasma_3d(ii,jj)

WRITE (unitNumber, '(A)') ' '

ii = 700
jj =  20
!WRITE (unitNumber, '(A,F15.5)') 'Pz_plasma_3d(700,20) = ', Pz_plasma_3d(ii,jj)
!WRITE (unitNumber, '(A,F15.5)') 'GLAT_plasma_3d(700,20) = ', GLAT_plasma_3d(ii,jj)
!WRITE (unitNumber, '(A,F15.5)') 'GLONd_plasma_3d(700,20) = ', GLONd_plasma_3d(ii,jj)

WRITE (unitNumber, '(A,i4,a1,i4,a4,F15.5)') 'Pz_plasma_3d(' , ii, ',' , jj, ') = ', Pz_plasma_3d(ii,jj)
WRITE (unitNumber, '(A,i4,a1,i4,a4,F15.5)') 'GLAT_plasma_3d(' , ii, ',' , jj, ') = ', GLAT_plasma_3d(ii,jj)
WRITE (unitNumber, '(A,i4,a1,i4,a4,F15.5)') 'GLONd_plasma_3d(', ii, ',' , jj, ' ) = ', GLONd_plasma_3d(ii,jj)


END SUBROUTINE check2GridIPE


!=============================================================================================================

SUBROUTINE checkThermo(directory, fileName, unitNumber, &
       GT_input_dataset, GT_output_dataset, nday, &
       Universal_Time_seconds, solar_declination_angle_radians,  &
       houghSize, hough11, hough22, hough23, hough24, hough25, &
       amplitude, tidalPhase, GT_ht_dim, GT_lat_dim, GT_lon_dim, &
       wind_southwards_ms1_FROM_GT, &
       wind_eastwards_ms1_FROM_GT, wvz_FROM_GT, rmt_FROM_GT, &
       Temperature_K_FROM_GT, ht_FROM_GT, Ne_density_FOR_GT)


USE moduleAmplitude  ! has amplitude type
USE moduleTidalPhase ! has tidalPhase type

IMPLICIT NONE

character(*), intent(in) :: directory
character(*), intent(in) :: fileName
INTEGER, intent(in) :: unitNumber
CHARACTER(*), intent(in) :: GT_input_dataset , GT_output_dataset
INTEGER, intent(in) :: nday
REAL*8, intent(in) :: universal_time_seconds 
REAL*8, intent(in) :: solar_declination_angle_radians
INTEGER, intent(in) :: houghSize
REAL*8, intent(in) :: hough11(houghSize)
REAL*8, intent(in) :: hough22(houghSize)
REAL*8, intent(in) :: hough23(houghSize)
REAL*8, intent(in) :: hough24(houghSize)
REAL*8, intent(in) :: hough25(houghSize)
TYPE(amplitudeType), intent(in) :: amplitude
TYPE(tidalPhaseType), intent(in) :: tidalPhase

INTEGER, intent(in) :: GT_ht_dim, GT_lat_dim, GT_lon_dim
REAL(kind=8), intent(in) :: wind_southwards_ms1_FROM_GT(GT_ht_dim, GT_lat_dim, GT_lon_dim)
REAL(kind=8), intent(in) :: wind_eastwards_ms1_FROM_GT(GT_ht_dim, GT_lat_dim, GT_lon_dim)
REAL(kind=8), intent(in) :: wvz_FROM_GT(GT_ht_dim, GT_lat_dim, GT_lon_dim)
REAL(kind=8), intent(in) :: rmt_FROM_GT(GT_ht_dim, GT_lat_dim, GT_lon_dim)
REAL(kind=8), intent(in) :: Temperature_K_FROM_GT(GT_ht_dim, GT_lat_dim, GT_lon_dim)
REAL(kind=8), intent(in) :: ht_FROM_GT(GT_ht_dim, GT_lat_dim, GT_lon_dim)

REAL*8, intent(in) :: Ne_density_FOR_GT(GT_ht_dim, GT_lat_dim, GT_lon_dim)


! Local Variables ----------
CHARACTER(100) :: fullFileName
INTEGER ::  FileOpenStatus  ! for checking if we can open the file

! BEGIN CODE ====================================================================

fullFileName = TRIM(directory)//TRIM(fileName)

!================================================
! Open for checking values from GT 
!================================================
OPEN (UNIT=unitNumber, FILE=fullFileName, &
                 ACTION="WRITE", IOSTAT=FileOpenStatus)
IF (FileOpenStatus > 0) then
    print *,"checkThermo : Cannot open "//fullFileName//" ***********"
    STOP 
ENDIF

print *,'checkThermo : writing to ',fullFileName

WRITE (unitNumber, '(A)') GT_input_dataset
WRITE (unitNumber, '(A)') GT_output_dataset
WRITE (unitNumber, '(A,I5)') 'nday = ', nday
WRITE (unitNumber, '(A, F20.7)' ) 'Universal_Time_seconds = ',  Universal_Time_seconds
WRITE (unitNumber, '(A, 3x, F10.5)' ) 'solar_declination_angle_radians = ', &
                    solar_declination_angle_radians
WRITE (unitNumber, '(A, 181F10.5)' ) 'hough11 = ',hough11
WRITE (unitNumber, '(A)' ) ' '
WRITE (unitNumber, '(A, 181F10.5)' ) 'hough22 = ',hough22
WRITE (unitNumber, '(A)' ) ' '
WRITE (unitNumber, '(A, 181F10.5)' ) 'hough23 = ',hough23
WRITE (unitNumber, '(A)' ) ' '
WRITE (unitNumber, '(A, 181F10.5)' ) 'hough24 = ',hough24
WRITE (unitNumber, '(A)' ) ' '
WRITE (unitNumber, '(A, 181F10.5)' ) 'hough25 = ',hough25

WRITE (unitNumber, '(A)' ) ' '
WRITE (unitNumber, '(A, 5F10.5)' ) 'amplitude = ', amplitude

WRITE (unitNumber, '(A)' ) ' '
WRITE (unitNumber, '(A, 5F10.5)' ) 'tidalPhase = ', tidalPhase


CALL checkThermoArrays(unitNumber, &
       GT_ht_dim, GT_lat_dim, GT_lon_dim, &
       wind_southwards_ms1_FROM_GT, &
       wind_eastwards_ms1_FROM_GT, wvz_FROM_GT, rmt_FROM_GT, &
       Temperature_K_FROM_GT, ht_FROM_GT)


WRITE (unitNumber, '(A)' ) ' '
WRITE (unitNumber, '(A, 20F25.5)' ) 'Ne_density_FOR_GT(1,1,:) = ', Ne_density_FOR_GT(1,1,:)
WRITE (unitNumber, '(A, 20F25.5)' ) 'Ne_density_FOR_GT(11,51,:) = ', Ne_density_FOR_GT(11,51,:)
WRITE (unitNumber, '(A, 20F25.5)' ) 'Ne_density_FOR_GT(15,91,:) = ', Ne_density_FOR_GT(15,91,:)


END SUBROUTINE checkThermo

!=================================================================================================



SUBROUTINE checkThermoArrays(unitNumber, &
       GT_ht_dim, GT_lat_dim, GT_lon_dim, &
       wind_southwards_ms1_FROM_GT, &
       wind_eastwards_ms1_FROM_GT, wvz_FROM_GT, rmt_FROM_GT, &
       Temperature_K_FROM_GT, ht_FROM_GT)


IMPLICIT NONE

INTEGER, intent(in) :: unitNumber

INTEGER, intent(in) :: GT_ht_dim, GT_lat_dim, GT_lon_dim
REAL(kind=8), intent(in) :: wind_southwards_ms1_FROM_GT(GT_ht_dim, GT_lat_dim, GT_lon_dim)
REAL(kind=8), intent(in) :: wind_eastwards_ms1_FROM_GT(GT_ht_dim, GT_lat_dim, GT_lon_dim)
REAL(kind=8), intent(in) :: wvz_FROM_GT(GT_ht_dim, GT_lat_dim, GT_lon_dim)
REAL(kind=8), intent(in) :: rmt_FROM_GT(GT_ht_dim, GT_lat_dim, GT_lon_dim)
REAL(kind=8), intent(in) :: Temperature_K_FROM_GT(GT_ht_dim, GT_lat_dim, GT_lon_dim)
REAL(kind=8), intent(in) :: ht_FROM_GT(GT_ht_dim, GT_lat_dim, GT_lon_dim)


! Local Variables ----------


! BEGIN CODE ====================================================================


WRITE (unitNumber, '(A)' ) ' '
WRITE (unitNumber, '(A, 20F10.5)' ) 'wind_southwards_ms1_FROM_GT(1,1,:) = ', &
                      wind_southwards_ms1_FROM_GT(1,1,:)
WRITE (unitNumber, '(A, 20F10.5)' ) 'wind_southwards_ms1_FROM_GT(11,51,:) = ', &
                      wind_southwards_ms1_FROM_GT(11,51,:)
WRITE (unitNumber, '(A, 20F10.5)' ) 'wind_southwards_ms1_FROM_GT(15,91,:) = ', &
                      wind_southwards_ms1_FROM_GT(15,91,:)
WRITE (unitNumber, '(A)' ) ' '
WRITE (unitNumber, '(A, 20F10.5)' ) 'wind_eastwards_ms1_FROM_GT(1,1,:) = ', &
                      wind_eastwards_ms1_FROM_GT(1,1,:)
WRITE (unitNumber, '(A, 20F10.5)' ) 'wind_eastwards_ms1_FROM_GT(11,51,:) = ', &
                      wind_eastwards_ms1_FROM_GT(11,51,:)
WRITE (unitNumber, '(A, 20F10.5)' ) 'wind_eastwards_ms1_FROM_GT(15,91,:) = ', &
                      wind_eastwards_ms1_FROM_GT(15,91,:)

WRITE (unitNumber, '(A)' ) ' '
WRITE (unitNumber, '(A, 20F10.5)' ) 'wvz_FROM_GT(1,1,:) = ', wvz_FROM_GT(1,1,:)
WRITE (unitNumber, '(A, 20F10.5)' ) 'wvz_FROM_GT(11,51,:) = ', wvz_FROM_GT(11,51,:)
WRITE (unitNumber, '(A, 20F10.5)' ) 'wvz_FROM_GT(15,91,:) = ', wvz_FROM_GT(15,91,:)

WRITE (unitNumber, '(A)' ) ' '
WRITE (unitNumber, '(A, 20F10.5)' ) 'rmt_FROM_GT(1,1,:) = ',rmt_FROM_GT(1,1,:)
WRITE (unitNumber, '(A, 20F10.5)' ) 'rmt_FROM_GT(11,51,:) = ',rmt_FROM_GT(11,51,:)
WRITE (unitNumber, '(A, 20F10.5)' ) 'rmt_FROM_GT(15,91,:) = ',rmt_FROM_GT(15,91,:)

WRITE (unitNumber, '(A)' ) ' '
WRITE (unitNumber, '(A, 20F15.5)' ) 'Temperature_K_FROM_GT(1,1,:) = ',Temperature_K_FROM_GT(1,1,:)
WRITE (unitNumber, '(A, 20F15.5)' ) 'Temperature_K_FROM_GT(11,51,:) = ',Temperature_K_FROM_GT(11,51,:)
WRITE (unitNumber, '(A, 20F15.5)' ) 'Temperature_K_FROM_GT(15,91,:) = ',Temperature_K_FROM_GT(15,91,:)

WRITE (unitNumber, '(A)' ) ' '
WRITE (unitNumber, '(A, 20F15.5)' ) 'ht_FROM_GT(1,1,:) = ',ht_FROM_GT(1,1,:)
WRITE (unitNumber, '(A, 20F15.5)' ) 'ht_FROM_GT(11,51,:) = ',ht_FROM_GT(11,51,:)
WRITE (unitNumber, '(A, 20F15.5)' ) 'ht_FROM_GT(15,91,:) = ',ht_FROM_GT(15,91,:)




END SUBROUTINE checkThermoArrays


!=================================================================================================


SUBROUTINE checkInterp(directory, unitNumber, &
                       interface_hts, GT_lat_dim, GT_lon_dim, &
                       O_density_fixed_ht, O2_density_fixed_ht, N2_density_fixed_ht, &
                       Vx_fixed_ht, Vy_fixed_ht, wvz_fixed_ht, &
                       tts_fixed_ht, qion3d_fixed_ht, elx_fixed_ht, &
                       ely_fixed_ht)


IMPLICIT NONE

character(*), intent(in) :: directory
!character(*), intent(in) :: fileName
INTEGER, intent(in) :: unitNumber

INTEGER, intent(in) ::  interface_hts, GT_lat_dim, GT_lon_dim


REAL(kind=8), intent(in) :: O_density_fixed_ht(interface_hts, GT_lat_dim, GT_lon_dim)
REAL(kind=8), intent(in) :: O2_density_fixed_ht(interface_hts, GT_lat_dim, GT_lon_dim)
REAL(kind=8), intent(in) :: N2_density_fixed_ht(interface_hts, GT_lat_dim, GT_lon_dim)
REAL(kind=8), intent(in) :: Vx_fixed_ht(interface_hts, GT_lat_dim, GT_lon_dim)
REAL(kind=8), intent(in) :: Vy_fixed_ht(interface_hts, GT_lat_dim, GT_lon_dim)
REAL(kind=8), intent(in) :: Wvz_fixed_ht(interface_hts, GT_lat_dim, GT_lon_dim)
REAL(kind=8), intent(in) :: tts_fixed_ht(interface_hts, GT_lat_dim, GT_lon_dim)
REAL(kind=8), intent(in) :: qion3d_fixed_ht(interface_hts, GT_lat_dim, GT_lon_dim)
REAL(kind=8), intent(in) :: elx_fixed_ht(interface_hts, GT_lat_dim, GT_lon_dim)
REAL(kind=8), intent(in) :: ely_fixed_ht(interface_hts, GT_lat_dim, GT_lon_dim)



! Local Variables ---------------
CHARACTER(100) :: fullFileName
INTEGER ::  FileOpenStatus  ! for checking if we can open the file
INTEGER ::  ii, jj

! BEGIN CODE ====================================================================

!fullFileName = TRIM(directory)//TRIM(fileName)
fullFileName = TRIM(directory)//TRIM(debugThermoInterpFileName)

!======================================================
! Open file for checking values from the interpolation
!======================================================
OPEN (UNIT=unitNumber, FILE=fullFileName, &
                 ACTION="WRITE", IOSTAT=FileOpenStatus)
IF (FileOpenStatus > 0) then
    print *,"checkInterp : Cannot open "//fullFileName//" ***********"
    STOP 
ENDIF

!========================================
! Open for checking values from the interpolation
!========================================
!OPEN (UNIT=unitNumber, FILE='/mtb/save/wx20lrm/DATA/IONOSPHERE/gtgipIPE/OUTPUT/interpOut.dat', &
!      ACTION="WRITE", IOSTAT=FileOpenStatus)
!IF (FileOpenStatus > 0) STOP "Cannot open  interpOut.dat file ********"


do ii = 1, interface_hts
   do jj = 1, GT_lat_dim

       WRITE (unitNumber, '(A,I2,A,I2,A)' ) 'O_density_fixed_ht(',  ii, ',' ,jj, ',:) = '
       WRITE (unitNumber, '(20ES20.10)' ) O_density_fixed_ht(ii,jj,:)

    end do  ! jj
end do ! ii


do ii = 1, interface_hts
   do jj = 1, GT_lat_dim

      WRITE (unitNumber, '(A,I2,A,I2,A)' ) 'O2_density_fixed_ht(',  ii, ',' ,jj, ',:) = '
      WRITE (unitNumber, '(20ES20.10)' ) O2_density_fixed_ht(ii,jj,:)

   end do  ! jj
end do ! ii


do ii = 1, interface_hts
   do jj = 1, GT_lat_dim


      WRITE (unitNumber, '(A,I2,A,I2,A)' ) 'N2_density_fixed_ht(',  ii, ',' ,jj, ',:) = '
      WRITE (unitNumber, '(20ES20.10)' ) N2_density_fixed_ht(ii,jj,:)

   end do  ! jj
end do ! ii

do ii = 1, interface_hts
   do jj = 1, GT_lat_dim


      WRITE (unitNumber, '(A,I2,A,I2,A)' ) 'Vx_fixed_ht(',  ii, ',' ,jj, ',:) = '
      WRITE (unitNumber, '(20ES20.10)' ) Vx_fixed_ht(ii,jj,:)

   end do  ! jj
end do ! ii


do ii = 1, interface_hts
   do jj = 1, GT_lat_dim

      WRITE (unitNumber, '(A,I2,A,I2,A)' ) 'Vy_fixed_ht(',  ii, ',' ,jj, ',:) = '
      WRITE (unitNumber, '(20ES20.10)' ) Vy_fixed_ht(ii,jj,:)

   end do  ! jj
end do ! ii

do ii = 1, interface_hts
   do jj = 1, GT_lat_dim

      WRITE (unitNumber, '(A,I2,A,I2,A)' ) 'wvz_fixed_ht(',  ii, ',' ,jj, ',:) = '
      WRITE (unitNumber, '(20ES20.10)' ) wvz_fixed_ht(ii,jj,:)

   end do  ! jj
end do ! ii

do ii = 1, interface_hts
   do jj = 1, GT_lat_dim


      WRITE (unitNumber, '(A,I2,A,I2,A)' ) 'tts_fixed_ht(',  ii, ',' ,jj, ',:) = '
      WRITE (unitNumber, '(20ES20.10)' ) tts_fixed_ht(ii,jj,:)

   end do  ! jj
end do ! ii


do ii = 1, interface_hts
   do jj = 1, GT_lat_dim

      WRITE (unitNumber, '(A,I2,A,I2,A)' ) 'qion3d_fixed_ht(',  ii, ',' ,jj, ',:) = '
      WRITE (unitNumber, '(20ES20.10)' ) qion3d_fixed_ht(ii,jj,:)

    end do  ! jj
end do ! ii

do ii = 1, interface_hts
   do jj = 1, GT_lat_dim

      WRITE (unitNumber, '(A,I2,A,I2,A)' ) 'elx_fixed_ht(',  ii, ',' ,jj, ',:) = '
      WRITE (unitNumber, '(20ES20.10)' ) elx_fixed_ht(ii,jj,:)

   end do  ! jj
end do ! ii

do ii = 1, interface_hts
   do jj = 1, GT_lat_dim

      WRITE (unitNumber, '(A,I2,A,I2,A)' ) 'ely_fixed_ht(',  ii, ',' ,jj, ',:) = '
      WRITE (unitNumber, '(20ES20.10)' ) ely_fixed_ht(ii,jj,:)

   end do  ! jj
end do ! ii


END SUBROUTINE checkInterp

!====================================================================================================


SUBROUTINE checkFixedGeo(directory, fileName, unitNumber, &
                         npts, nmp, &
                         TN_plasma_input_3d, O_plasma_input_3d, O2_plasma_input_3d, &
                         N2_plasma_input_3d, GLAt_plasma_3d, GLOnd_plasma_3d, &
                         pz_plasma_3d, V_east_plasma, V_south_plasma, &
                         V_upward_plasma, ilon1_3d_fixed_ht, ilon2_3d_fixed_ht, &
                         ilat1_3d_fixed_ht, ilat2_3d_fixed_ht, ispecial_3d_fixed_ht, &
                         ihl_3d_fixed_ht, ihu_3d_fixed_ht, isFirstCallFixedHeight)


IMPLICIT NONE

character(*), intent(in) :: directory
character(*), intent(in) :: fileName
INTEGER, intent(in) :: unitNumber

INTEGER, intent(in) :: npts, nmp

REAL(kind=8), intent(in) :: TN_plasma_input_3d(npts, nmp), O_plasma_input_3d(npts, nmp), &
                            O2_plasma_input_3d(npts, nmp), N2_plasma_input_3d(npts, nmp)

REAL(kind=8), intent(in)  :: glat_plasma_3d(npts,nmp), glond_plasma_3d(npts,nmp), pz_plasma_3d(npts,nmp)

REAL(kind=8), intent(in) :: V_east_plasma(npts, nmp), &
                            V_south_plasma(npts, nmp), &
                            V_upward_plasma(npts, nmp)

INTEGER, intent(in) :: ilon1_3d_fixed_ht(npts,nmp), ilon2_3d_fixed_ht(npts,nmp)
INTEGER, intent(in) :: ilat1_3d_fixed_ht(npts,nmp), ilat2_3d_fixed_ht(npts,nmp)
INTEGER, intent(in) :: ispecial_3d_fixed_ht(npts,nmp)
INTEGER, intent(in) :: ihl_3d_fixed_ht(npts,nmp), ihu_3d_fixed_ht(npts,nmp)
LOGICAL, intent(in)  :: isFirstCallFixedHeight

! Local Variables ---------------
CHARACTER(100) :: fullFileName
INTEGER ::  FileOpenStatus  ! for checking if we can open the file


! BEGIN CODE ====================================================================

fullFileName = TRIM(directory)//TRIM(fileName)

!======================================================
! Open file for checking values from GT
!======================================================
OPEN (UNIT=unitNumber, FILE=fullFileName, &
                 ACTION="WRITE", IOSTAT=FileOpenStatus)
IF (FileOpenStatus > 0) then
    print *,"checkInterp : Cannot open "//fullFileName//" ***********"
    STOP 
ENDIF


    WRITE (unitNumber, '(A)' ) ' '
    WRITE (unitNumber, '(A, 80ES20.10)' ) 'TN_plasma_input_3d(1,:) = ', TN_plasma_input_3d(1,:)
    WRITE (unitNumber, '(A, 80ES20.10)' ) 'TN_plasma_input_3d(700,:) = ', TN_plasma_input_3d(700,:)
    WRITE (unitNumber, '(A, 80ES20.10)' ) 'TN_plasma_input_3d(NPTS,:) = ', TN_plasma_input_3d(NPTS,:)

    WRITE (unitNumber, '(A)' ) ' '
    WRITE (unitNumber, '(A, 80ES20.10)' ) 'O_plasma_input_3d(1,:) = ', O_plasma_input_3d(1,:)
    WRITE (unitNumber, '(A, 80ES20.10)' ) 'O_plasma_input_3d(700,:) = ', O_plasma_input_3d(700,:)
    WRITE (unitNumber, '(A, 80ES20.10)' ) 'O_plasma_input_3d(NPTS,:) = ', O_plasma_input_3d(NPTS,:)

    WRITE (unitNumber, '(A)' ) ' '
    WRITE (unitNumber, '(A, 80ES20.10)' ) 'O2_plasma_input_3d(1,:) = ', O2_plasma_input_3d(1,:)
    WRITE (unitNumber, '(A, 80ES20.10)' ) 'O2_plasma_input_3d(700,:) = ', O2_plasma_input_3d(700,:)
    WRITE (unitNumber, '(A, 80ES20.10)' ) 'O2_plasma_input_3d(NPTS,:) = ', O2_plasma_input_3d(NPTS,:)

    WRITE (unitNumber, '(A)' ) ' '
    WRITE (unitNumber, '(A, 80ES20.10)' ) 'N2_plasma_input_3d(1,:) = ', N2_plasma_input_3d(1,:)
    WRITE (unitNumber, '(A, 80ES20.10)' ) 'N2_plasma_input_3d(700,:) = ', N2_plasma_input_3d(700,:)
    WRITE (unitNumber, '(A, 80ES20.10)' ) 'N2_plasma_input_3d(NPTS,:) = ', N2_plasma_input_3d(NPTS,:)

    WRITE (unitNumber, '(A)' ) ' '
    WRITE (unitNumber, '(A, 80ES20.10)' ) 'GLAt_plasma_3d(1,:) = ', GLAt_plasma_3d(1,:)
    WRITE (unitNumber, '(A, 80ES20.10)' ) 'GLAt_plasma_3d(700,:) = ', GLAt_plasma_3d(700,:)
    WRITE (unitNumber, '(A, 80ES20.10)' ) 'GLAt_plasma_3d(NPTS,:) = ', GLAt_plasma_3d(NPTS,:)

    WRITE (unitNumber, '(A)' ) ' '
    WRITE (unitNumber, '(A, 80ES20.10)' ) 'GLOnd_plasma_3d(1,:) = ', GLOnd_plasma_3d(1,:)
    WRITE (unitNumber, '(A, 80ES20.10)' ) 'GLOnd_plasma_3d(700,:) = ', GLOnd_plasma_3d(700,:)
    WRITE (unitNumber, '(A, 80ES20.10)' ) 'GLOnd_plasma_3d(NPTS,:) = ', GLOnd_plasma_3d(NPTS,:)


    WRITE (unitNumber, '(A)' ) ' '
    WRITE (unitNumber, '(A, 80ES20.10)' ) 'pz_plasma_3d(1,:) = ', pz_plasma_3d(1,:)
    WRITE (unitNumber, '(A, 80ES20.10)' ) 'pz_plasma_3d(700,:) = ', pz_plasma_3d(700,:)
    WRITE (unitNumber, '(A, 80ES20.10)' ) 'pz_plasma_3d(NPTS,:) = ', pz_plasma_3d(NPTS,:)

    WRITE (unitNumber, '(A)' ) ' '
    WRITE (unitNumber, '(A, 80ES20.10)' ) 'V_east_plasma(1,:) = ', V_east_plasma(1,:)
    WRITE (unitNumber, '(A, 80ES20.10)' ) 'V_east_plasma(700,:) = ', V_east_plasma(700,:)
    WRITE (unitNumber, '(A, 80ES20.10)' ) 'V_east_plasma(NPTS,:) = ', V_east_plasma(NPTS,:)


    WRITE (unitNumber, '(A)' ) ' '
    WRITE (unitNumber, '(A, 80ES20.10)' ) 'V_south_plasma(1,:) = ', V_south_plasma(1,:)
    WRITE (unitNumber, '(A, 80ES20.10)' ) 'V_south_plasma(700,:) = ', V_south_plasma(700,:)
    WRITE (unitNumber, '(A, 80ES20.10)' ) 'V_south_plasma(NPTS,:) = ', V_south_plasma(NPTS,:)


    WRITE (unitNumber, '(A)' ) ' '
    WRITE (unitNumber, '(A, 80ES20.10)' ) 'V_upward_plasma(1,:) = ', V_upward_plasma(1,:)
    WRITE (unitNumber, '(A, 80ES20.10)' ) 'V_upward_plasma(700,:) = ', V_upward_plasma(700,:)
    WRITE (unitNumber, '(A, 80ES20.10)' ) 'V_upward_plasma(NPTS,:) = ', V_upward_plasma(NPTS,:)

    WRITE (unitNumber, '(A)' ) ' '
    WRITE (unitNumber, '(A, 80I7)' ) 'ilon1_3d_fixed_ht(1,:) = ', ilon1_3d_fixed_ht(1,:)
    WRITE (unitNumber, '(A, 80I7)' ) 'ilon1_3d_fixed_ht(700,:) = ', ilon1_3d_fixed_ht(700,:)
    WRITE (unitNumber, '(A, 80I7)' ) 'ilon1_3d_fixed_ht(NPTS,:) = ', ilon1_3d_fixed_ht(NPTS,:)

    WRITE (unitNumber, '(A)' ) ' '
    WRITE (unitNumber, '(A, 80I7)' ) 'ilon2_3d_fixed_ht(1,:) = ', ilon2_3d_fixed_ht(1,:)
    WRITE (unitNumber, '(A, 80I7)' ) 'ilon2_3d_fixed_ht(700,:) = ', ilon2_3d_fixed_ht(700,:)
    WRITE (unitNumber, '(A, 80I7)' ) 'ilon2_3d_fixed_ht(NPTS,:) = ', ilon2_3d_fixed_ht(NPTS,:)

    WRITE (unitNumber, '(A)' ) ' '
    WRITE (unitNumber, '(A, 80I7)' ) 'ilat1_3d_fixed_ht(1,:) = ', ilat1_3d_fixed_ht(1,:)
    WRITE (unitNumber, '(A, 80I7)' ) 'ilat1_3d_fixed_ht(700,:) = ', ilat1_3d_fixed_ht(700,:)
    WRITE (unitNumber, '(A, 80I7)' ) 'ilat1_3d_fixed_ht(NPTS,:) = ', ilat1_3d_fixed_ht(NPTS,:)

    WRITE (unitNumber, '(A)' ) ' '
    WRITE (unitNumber, '(A, 80I7)' ) 'ilat2_3d_fixed_ht(1,:) = ', ilat2_3d_fixed_ht(1,:)
    WRITE (unitNumber, '(A, 80I7)' ) 'ilat2_3d_fixed_ht(700,:) = ', ilat2_3d_fixed_ht(700,:)
    WRITE (unitNumber, '(A, 80I7)' ) 'ilat2_3d_fixed_ht(NPTS,:) = ', ilat2_3d_fixed_ht(NPTS,:)

    WRITE (unitNumber, '(A)' ) ' '
    WRITE (unitNumber, '(A, 80I7)' ) 'ispecial_3d_fixed_ht(1,:) = ', ispecial_3d_fixed_ht(1,:)
    WRITE (unitNumber, '(A, 80I7)' ) 'ispecial_3d_fixed_ht(700,:) = ', ispecial_3d_fixed_ht(700,:)
    WRITE (unitNumber, '(A, 80I7)' ) 'ispecial_3d_fixed_ht(NPTS,:) = ', ispecial_3d_fixed_ht(NPTS,:)

    WRITE (unitNumber, '(A)' ) ' '
    WRITE (unitNumber, '(A, 80I7)' ) 'ihl_3d_fixed_ht(1,:) = ', ihl_3d_fixed_ht(1,:)
    WRITE (unitNumber, '(A, 80I7)' ) 'ihl_3d_fixed_ht(700,:) = ', ihl_3d_fixed_ht(700,:)
    WRITE (unitNumber, '(A, 80I7)' ) 'ihl_3d_fixed_ht(NPTS,:) = ', ihl_3d_fixed_ht(NPTS,:)

    WRITE (unitNumber, '(A)' ) ' '
    WRITE (unitNumber, '(A, 80I7)' ) 'ihu_3d_fixed_ht(1,:) = ', ihu_3d_fixed_ht(1,:)
    WRITE (unitNumber, '(A, 80I7)' ) 'ihu_3d_fixed_ht(700,:) = ', ihu_3d_fixed_ht(700,:)
    WRITE (unitNumber, '(A, 80I7)' ) 'ihu_3d_fixed_ht(NPTS,:) = ', ihu_3d_fixed_ht(NPTS,:)

    WRITE (unitNumber, '(A)' ) ' '
    WRITE (unitNumber, '(A, I3)' ) 'isFirstCallFixedHeight = ', isFirstCallFixedHeight


END SUBROUTINE checkFixedGeo


!====================================================================================================
END MODULE moduleDriverDebug

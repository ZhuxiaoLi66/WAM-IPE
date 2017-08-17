program interface_to_fixed_height_single_file
  IMPLICIT NONE
  INTEGER, PARAMETER :: nheights=183
  INTEGER :: height_km
  REAL(kind=8),Dimension(3,nheights,91,90) ::  facfac_interface,dd_interface
  INTEGER,DIMENSION(3,nheights,91,90) :: ii1_interface, ii2_interface,&
&         ii3_interface, ii4_interface

  INTEGER, PARAMETER :: NPTS2D=44514
  INTEGER, PARAMETER :: NMP=80
  INTEGER, PARAMETER :: NLP=170
  REAL(kind=8) :: dtotinv, factor, d, fac
  INTEGER :: i , iheight , iup, ido, ih, ip, lun1, lun2, lun3, lun5
  INTEGER :: l , m ,  mp , lp , in1, in2
  INTEGER :: ifile , the_index
  INTEGER :: ilat, ilon, i_max_location
  INTEGER,DIMENSION(3,nheights,91,90) :: ii1, ii2, ii3, ii4
  REAL(kind=4),dimension(44514,80) :: oplus_ft
  REAL(kind=4),dimension(44514,80) :: hplus_ft
  REAL(kind=4),dimension(44514,80) :: heplus_ft
  REAL(kind=4),dimension(44514,80) :: nplus_ft
  REAL(kind=4),dimension(44514,80) :: noplus_ft
  REAL(kind=4),dimension(44514,80) :: o2plus_ft
  REAL(kind=4),dimension(44514,80) :: n2plus_ft
  REAL(kind=4),dimension(44514,80) :: oplus2d_ft
  REAL(kind=4),dimension(44514,80) :: oplus2p_ft
  REAL(kind=8),Dimension(nheights,91,90) :: oplus_fixed
  REAL(kind=8),Dimension(nheights,91,90) :: hplus_fixed
  REAL(kind=8),Dimension(nheights,91,90) :: heplus_fixed
  REAL(kind=8),Dimension(nheights,91,90) :: nplus_fixed
  REAL(kind=8),Dimension(nheights,91,90) :: noplus_fixed
  REAL(kind=8),Dimension(nheights,91,90) :: o2plus_fixed
  REAL(kind=8),Dimension(nheights,91,90) :: n2plus_fixed
  REAL(kind=8),Dimension(nheights,91,90) :: oplus2d_fixed
  REAL(kind=8),Dimension(nheights,91,90) :: oplus2p_fixed
  REAL(kind=8),Dimension(nheights,91,90) :: electron_density_fixed
  REAL(kind=8),Dimension(91,90) :: electron_density_fixed_300km
  REAL(kind=8),Dimension(91,90) :: total_electron_content
  REAL(kind=8),Dimension(nheights) :: electron_density_profile
  REAL(kind=8),Dimension(91,90) :: nmf2
  REAL(kind=8),Dimension(91,90) :: hmf2
  REAL(kind=8),Dimension(3) :: oplus_interpolated
  REAL(kind=8),Dimension(3) :: hplus_interpolated
  REAL(kind=8),Dimension(3) :: heplus_interpolated
  REAL(kind=8),Dimension(3) :: nplus_interpolated
  REAL(kind=8),Dimension(3) :: noplus_interpolated
  REAL(kind=8),Dimension(3) :: o2plus_interpolated
  REAL(kind=8),Dimension(3) :: n2plus_interpolated
  REAL(kind=8),Dimension(3) :: oplus2d_interpolated
  REAL(kind=8),Dimension(3) :: oplus2p_interpolated
  REAL(kind=8) :: x1,x2,x3,y1,y2,y3,a,b,c
  INTEGER :: i1,i2,i3
  CHARACTER (len=2) :: number
  CHARACTER (len=20) :: filename
  character(len=8) :: fmt 
  CHARACTER(len=255), dimension(880) :: list_of_plasma_files
  CHARACTER(len=255) :: a_plasma_file
  CHARACTER(len=255), dimension(880) :: list_of_timestamps
  CHARACTER(len=255) :: results_directory
  CHARACTER(len=255) :: input_plasma_file
  CHARACTER(len=255) :: output_plasma_file
  fmt = '(I2.2)' 

print *,''
print *,'************* IPE_POSTPROCESSOR *************'
print *,''

CALL getenv("RESULTS_DIRECTORY", results_directory)
CALL getenv("INPUT_PLASMA_FILE", input_plasma_file)
CALL getenv("OUTPUT_PLASMA_FILE", output_plasma_file)

lun1=101
print *, 'Reading GIP grid interpolation file'
!open(UNIT=lun1,file='GIP_Fixed_GEO_grid_lowres_corrected',form='formatted',status='old')
!    read(lun1,*) facfac_interface
!    read(lun1,*) dd_interface
!    read(lun1,*) ii1_interface
!    read(lun1,*) ii2_interface
!    read(lun1,*) ii3_interface
!    read(lun1,*) ii4_interface
open(UNIT=lun1,file='GIP_Fixed_GEO_grid_lowres_corrected.bin',form='unformatted',status='old')
    read(lun1) facfac_interface
    read(lun1) dd_interface
    read(lun1) ii1_interface
    read(lun1) ii2_interface
    read(lun1) ii3_interface
    read(lun1) ii4_interface
close(lun1)

lun5=105
open(UNIT=lun5,file='GIP_Fixed_GEO_grid_lowres_corrected.bin',form='unformatted',status='unknown')
    write(lun5) facfac_interface
    write(lun5) dd_interface
    write(lun5) ii1_interface
    write(lun5) ii2_interface
    write(lun5) ii3_interface
    write(lun5) ii4_interface
close(lun5)

lun2=102
print *, 'Data directory : ',trim(results_directory)
print *, 'Reading input file   : ', trim(input_plasma_file)
open(UNIT=lun2,file=trim(results_directory)//trim(input_plasma_file),form='unformatted',status='old')
     read(lun2) oplus_ft
     read(lun2) hplus_ft
     read(lun2) heplus_ft
     read(lun2) nplus_ft
     read(lun2) noplus_ft
     read(lun2) o2plus_ft
     read(lun2) n2plus_ft
     read(lun2) oplus2d_ft
     read(lun2) oplus2p_ft
close(lun2)

!g
!g  loop over all heights and the lats/longs......
!g
  do 100 iheight=1,nheights
      height_km = ((iheight - 1) * 5) + 90
!      print *, iheight, height_km
      do 200 l=1,90
          do 300 m=1,91
!g
!g  Initialise our output parameters....
!g
              oplus_fixed(iheight,m,l) = 0.0
              hplus_fixed(iheight,m,l) = 0.0
              heplus_fixed(iheight,m,l) = 0.0
              nplus_fixed(iheight,m,l) = 0.0
              noplus_fixed(iheight,m,l) = 0.0
              o2plus_fixed(iheight,m,l) = 0.0
              n2plus_fixed(iheight,m,l) = 0.0
              oplus2d_fixed(iheight,m,l) = 0.0
              oplus2p_fixed(iheight,m,l) = 0.0
              dtotinv = 0.0
!g
!g  The number of contributing flux-tube points
!g  is always 3....
!g
              do 400 i=1,3
!g
!g  The first interpolation uses facfac
!g
                  factor=facfac_interface(i,iheight,m,l)
                  mp=ii1_interface(i,iheight,m,l)
                  lp=ii2_interface(i,iheight,m,l)
                  in2=ii3_interface(i,iheight,m,l)
                  in1=ii4_interface(i,iheight,m,l)

                  oplus_interpolated(i) = ((oplus_ft(in2,mp)-oplus_ft(in1,mp))*factor) + oplus_ft(in1,mp)
                  hplus_interpolated(i) = ((hplus_ft(in2,mp)-hplus_ft(in1,mp))*factor) + hplus_ft(in1,mp)
                  heplus_interpolated(i) = ((heplus_ft(in2,mp)-heplus_ft(in1,mp))*factor) + heplus_ft(in1,mp)
                  nplus_interpolated(i) = ((nplus_ft(in2,mp)-nplus_ft(in1,mp))*factor) + nplus_ft(in1,mp)
                  noplus_interpolated(i) = ((noplus_ft(in2,mp)-noplus_ft(in1,mp))*factor) + noplus_ft(in1,mp)
                  o2plus_interpolated(i) = ((o2plus_ft(in2,mp)-o2plus_ft(in1,mp))*factor) + o2plus_ft(in1,mp)
                  n2plus_interpolated(i) = ((n2plus_ft(in2,mp)-n2plus_ft(in1,mp))*factor) + n2plus_ft(in1,mp)
                  oplus2d_interpolated(i) = ((oplus2d_ft(in2,mp)-oplus2d_ft(in1,mp))*factor) + oplus2d_ft(in1,mp)
                  oplus2p_interpolated(i) = ((oplus2p_ft(in2,mp)-oplus2p_ft(in1,mp))*factor) + oplus2p_ft(in1,mp)

!g
!g Now we calculate the parameters at each point.....
!g
                  d = dd_interface(i,iheight,m,l)
                  dtotinv = (1./d) + dtotinv

                  oplus_fixed(iheight,m,l) = (oplus_interpolated(i)/d) + oplus_fixed(iheight,m,l)
                  hplus_fixed(iheight,m,l) = (hplus_interpolated(i)/d) + hplus_fixed(iheight,m,l)
                  heplus_fixed(iheight,m,l) = (heplus_interpolated(i)/d) + heplus_fixed(iheight,m,l)
                  nplus_fixed(iheight,m,l) = (nplus_interpolated(i)/d) + nplus_fixed(iheight,m,l)
                  noplus_fixed(iheight,m,l) = (noplus_interpolated(i)/d) + noplus_fixed(iheight,m,l)
                  o2plus_fixed(iheight,m,l) = (o2plus_interpolated(i)/d) + o2plus_fixed(iheight,m,l)
                  n2plus_fixed(iheight,m,l) = (n2plus_interpolated(i)/d) + n2plus_fixed(iheight,m,l)
                  oplus2d_fixed(iheight,m,l) = (oplus2d_interpolated(i)/d) + oplus2d_fixed(iheight,m,l)
                  oplus2p_fixed(iheight,m,l) = (oplus2p_interpolated(i)/d) + oplus2p_fixed(iheight,m,l)

              400 ENDDO ! do 400 i=1,3

              oplus_fixed(iheight,m,l) = oplus_fixed(iheight,m,l)/dtotinv
              hplus_fixed(iheight,m,l) = hplus_fixed(iheight,m,l)/dtotinv
              heplus_fixed(iheight,m,l) = heplus_fixed(iheight,m,l)/dtotinv
              nplus_fixed(iheight,m,l) = nplus_fixed(iheight,m,l)/dtotinv
              noplus_fixed(iheight,m,l) = noplus_fixed(iheight,m,l)/dtotinv
              o2plus_fixed(iheight,m,l) = o2plus_fixed(iheight,m,l)/dtotinv
              n2plus_fixed(iheight,m,l) = n2plus_fixed(iheight,m,l)/dtotinv
              oplus2d_fixed(iheight,m,l) = oplus2d_fixed(iheight,m,l)/dtotinv
              oplus2p_fixed(iheight,m,l) = oplus2p_fixed(iheight,m,l)/dtotinv

          300 ENDDO
      200 ENDDO
  100 ENDDO

electron_density_fixed = oplus_fixed + hplus_fixed + heplus_fixed + nplus_fixed + noplus_fixed & 
                       + o2plus_fixed + n2plus_fixed + oplus2d_fixed + oplus2p_fixed

! Electron Density at a fixed height of 300km....
electron_density_fixed_300km(1:91,1:90) = electron_density_fixed(43,1:91,1:90)

! Total Electron Content....
total_electron_content = 0.0
do 500 iheight=1,nheights
  total_electron_content(1:91,1:90) = total_electron_content(1:91,1:90) + (electron_density_fixed(iheight,1:91,1:90) * 5000.0) 
500 enddo

! NmF2 and hmF2....
do 660 ilon = 1 , 90
do 650 ilat = 1 , 91
  electron_density_profile = electron_density_fixed(1:183,ilat,ilon)
  i_max_location = maxval(maxloc(electron_density_profile))

  i1 = i_max_location - 1
  i2 = i_max_location
  i3 = i_max_location + 1
  x1 = (float(i1 - 1) * 5.0) + 90.0
  x2 = (float(i2 - 1) * 5.0) + 90.0
  x3 = (float(i3 - 1) * 5.0) + 90.0
  y1 = electron_density_profile(i1)
  y2 = electron_density_profile(i2)
  y3 = electron_density_profile(i3)

  c = (x3*y1 - x3*y2 + x1*y2 + x2*y3 - x2*y1 - x1*y3) / (x3*x1*x1 - x3*x2*x2 + x1*x2*x2 - x2*x1*x1 + x2*x3*x3 - x1*x3*x3)         
  b = (y2 - (c*x2*x2) + (c*x1*x1) - y1) / (x2 - x1)
  a = y1 - (b*x1) - (c*x1*x1)
  hmf2(ilat,ilon) = (0.0 - b) / (2*c)
  nmf2(ilat,ilon) = a + (b*hmf2(ilat,ilon)) + (c*hmf2(ilat,ilon)*hmf2(ilat,ilon))                      

!  nmf2(ilat,ilon) = electron_density_profile(i_max_location)
!  hmf2(ilat,ilon) = (float(i_max_location - 1) * 5.0) + 90.0
!  print *, ilat, ilon, i_max_location, nmf2(ilat,ilon), hmf2(ilat,ilon)
650 enddo
660 enddo


lun3=103
print *, 'Writing output file  : ', trim(output_plasma_file)
open(UNIT=lun3,file=trim(results_directory)//trim(output_plasma_file),form='formatted',status='unknown')
write(lun3,"(20e12.4)") total_electron_content
write(lun3,"(20e12.4)") nmf2
write(lun3,"(20e12.4)") hmf2
write(lun3,"(20e12.4)") electron_density_fixed_300km
close(lun3)
print *,''
print *, '******************* DONE ********************'
print *,''

end program interface_to_fixed_height_single_file

!Jan2011:original code was provided from Astrid from WACCM.
!Aug2011:this code was provided from Fei Wu from the WAM version.
!--------------------------------------------  
      module module_read_acoef
!--------------------------------------------------------------------- 
!   L.Scherliess ludger@gaim.cass.usu.edu
! - Weimer96 model
! Author: A. Maute Dec 2003  am 12/30/03 
!------------------------------------------------------------------------------ 

      implicit none

!nm20121003:module parameters are separated into efield.f!
      public :: read_acoef

      contains
!-----------------------------------------------------------------------
      subroutine read_acoef (efield_lflux_file, efield_hflux_file)
!----------------------------------------------------------------     
! Purpose:  
!    1. read in coefficients A_klmn^lf for solar cycle minimum and
!      A_klmn^hf for maximum 
!    2. adjust S_a (f107d) such that if S_a<80 or S_a > 220 it has reasonable numbers
!      S_aM = [atan{(S_a-65)^2/90^2}-a90]/[a180-a90]
!      a90  = atan [(90-65)/90]^2
!      a180 = atan [(180-65)/90]^2
!    3. inter/extrapolation of the coefficient to the actual flux which is
!      given by the user
!      A_klmn = S_aM [A_klmn^hf-A_klmn^lf]/90. + 2*A_klmn^lf-A_klmn^hf
!
! Method:
!
! Author: A. Maute Nov 2003  am 11/19/03
!---------------------------------------------------------------

!     use ioFileMod,     only : getfil
!     use units,         only : getunit, freeunit
!nm20121003
      USE efield !,ONLY:

      character(len=*), intent(in) :: efield_lflux_file
      character(len=*), intent(in) :: efield_hflux_file

      integer  :: i,ios,unit,istat
      character (len=256):: locfn

!SMS$SERIAL(<a_lf,OUT> : default=ignore) BEGIN
!------------------------------------------------------------------    
!  get coefficients file for solar minimum: 
!-----------------------------------------------------------------                                                                   
!     unit     = getunit()
      unit     = 11
!     call getfil( efield_lflux_file, locfn, 0 )
      locfn=efield_lflux_file

!------------------------------------------------------------------    
! open datafile with coefficients A_klnm
!------------------------------------------------------------------     
!     write(iulog,*) 'read_acoef: open file ',trim(locfn),
!    &' unit ',unit
!dbg20110906
      print *,'unit',unit
      open(unit=unit,file=trim(locfn),status = 'old',iostat = ios)
!     if(ios.gt.0) then
!     write(iulog,*) 
!    &'read_acoef: error in opening coeff_lf file',
!    &' unit ',unit
!       call endrun
!     end if

!----------------------------------------------------------------------------                                                                   
! read datafile with coefficients A_klnm
!--------------------------------------------------------------------   
!     write(iulog,*) 'read_acoef: read file ',trim(locfn),' unit ',unit
      read(unit,*,iostat = ios) a_lf
!     if(ios.gt.0) then
!     write(iulog,*) 
!    &'read_acoef: error in reading coeff_lf file',' unit ',unit
!       call endrun
!     end if

!--------------------------------------------------------------------  
! close & free unit      
!--------------------------------------------------------------------  
      close(unit)
!     call freeunit(unit)
!     write(iulog,*) 'read_acoef: free unit ',unit
!SMS$SERIAL END

!SMS$SERIAL(<a_hf,OUT> : default=ignore) BEGIN

!--------------------------------------------------------------------  
!  get coefficients file for solar maximum: 
!--------------------------------------------------------------------
!     unit     = getunit()
      unit     = 10
!     call getfil( efield_hflux_file, locfn, 0 )
      locfn= efield_hflux_file

!-------------------------------------------------------------------
! open datafile with coefficients A_klnm
!------------------------------------------------------------------
!     write(iulog,*) 'read_acoef: open file ',trim(locfn),' unit ',unit
      open(unit=unit,file=trim(locfn),status = 'old',iostat = ios)
!     if(ios.gt.0) then
!      write(iulog,*) 
!    &'read_acoef: error in opening coeff_hf file',' unit ',unit
!       call endrun
!     end if

!-----------------------------------------------------------------
! read datafile with coefficients A_klnm
!----------------------------------------------------------------
!     write(iulog,*) 'read_acoef: read file ',trim(locfn)
      read(unit,*,iostat = ios) a_hf
!     if(ios.gt.0) then
!      write(iulog,*) 
!    &'read_acoef: error in reading coeff_hf file',' unit ',unit
!       call endrun
!     end if

!---------------------------------------------------------------
! close & free unit      
!-------------------------------------------------------------- 
      close(unit)
!     call freeunit(unit)
!     write(iulog,*) 'read_acoef: free unit ',unit
!SMS$SERIAL END

      end subroutine read_acoef
!-----------------------------------------------------------------------
      end module module_read_acoef

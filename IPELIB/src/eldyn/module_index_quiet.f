!Jan2011:original code was provided from Astrid from WACCM.
!Aug2011:this code was provided from Fei Wu from the WAM version.
!--------------------------------------------  
      module module_index_quiet
!--------------------------------------------------------------------- 
! description: calculates the electric potential for a given year,
!      day of year,UT, F10.7, B_z(K_p)
! - low/midlatitudes electric potential is from an empirical model from
!   L.Scherliess ludger@gaim.cass.usu.edu
! - high latitude electric potential is from Weimer96 model
! - the transition zone is smoothed
! - output is the horizontal global electric field in magnetic coordinates direction
!  at every magnetic local time grid point expressed in degrees (0 deg-0MLT; 360deg 24 MLT)
!
! input 
!      integer :: iday,     ! day number of year
!                 iyear     ! year
!      real:: ut,       ! universal time 
!                 F10.7,    ! solar flux       (see ionosphere module)
!                 bz        ! component of IMF (see ionosphere module)
! output
!      real ::               &
!       ed1(0:nmlon,0:nmlat),    &  ! zonal electric field Ed1  [V/m] 
!       ed2(0:nmlon,0:nmlat)        ! meridional electric field Ed2/sin I_m  [V/m]  
!
! notes:
!
! - !to be done (commented out): input S_a F10.7/ Kp from WACCM and calculate B_z 
!    from these inputs
! - assume regular geomagnetic grid 
! - uses average year 365.24 days/year 30.6001 day/mo s. Weimer
! - get_tilt works only for iyear >= 1900
! - Weimer model 1996, Dan Weimer (not with the updates from B.Emery)
! - fixed parameters: B_z, B_y units nT  CHANGE THIS
!                     F10.7
! - we assume that the reference height is 300km for the emperical potential model
! - as a first approximation the electric field is constant in height
!   WATCH what is the upper boundary condition in WACCM
! - for all the calculation done here we set the reference height to the same 
!   value as in tiegcm (hr=130km)
! - 12/15/03 input value iseasav : replaced by day -> month and day of month
! - 12/15/03 S_aM calculated according to Scherliess draft paper and added
!   S_aM(corrected) = 90*(S_aM+1) to get variation in fig 1 Scherliess draft
!
! Author: A. Maute Dec 2003  am 12/30/03 
!------------------------------------------------------------------------------ 

c     use shr_kind_mod,  only: r8 => shr_kind_r8
c     use physconst,     only: pi
c     use abortutils,    only: endrun
c     use cam_logfile,   only: iulog
   
      implicit none

!nm20121003:module parameters are separated into efield.f!
      public :: index_quiet

      contains

                                                                      



                                                                     

      subroutine index_quiet
!-----------------------------------------------------------------
! Purpose: set up index for factors f_m(mlt),f_l(UT),f_-k(d) to
!    describe the electric potential Phi for the empirical model   
!
! Method:
!    Phi = sum_k sum_l sum_m sum_n [ A_klmn * P_n^m *f_m(mlt)*f_l(UT)*f_-k(d)]
!    - since the electric potential is symmetric about the equator
!      n+m odd terms are set zero resp. not used
!    - in the summation for calculation Phi the index have the following
!      range n=1,12 and m=-n,n, k=0,2 l=-2,2
!
! Author: A. Maute Nov 2003  am 11/18/03
!----------------------------------------------------------------       
!nm20121003
      USE efield !,ONLY:
      implicit none

!----------------------------------------------------------------      
!	... local variables
!----------------------------------------------------------------                                                                   
      integer :: i, j, k, l, n, m

      i = 0 	! initialize
      j = 1 
      do k = 2,0,-1
        do l = -2,2
          if( k == 2 .and. abs(l) == 2 ) then
             cycle
          end if
          do n = 1,12
            do m = -18,18 
              if( abs(m) <= n ) then		    !  |m| < n
                if( (((n-m)/2)*2) == (n-m) ) then   ! only n+m even
             	  if( n-abs(m) <= 9 ) then	    ! n-|m| <= 9 why?
             	    kf(i) = 2-k
             	    lf(i) = l
             	    nf(i) = n
             	    mf(i) = m
             	    jf(i) = j
             	    i	  = i + 1	 ! counter
                  end if
                end if
              end if
            end do ! m
          end do ! n
        end do ! l
      end do ! k

      imax = i - 1  
      if(imax /= ni ) then    ! check if imax == ni 
c       write(iulog,'(a19,i5,a18,i5)') 'index_quiet: imax= ',imax,  
c    &     ' not equal to ni =',ni 
        stop
      end if							
c     if(debug) write(iulog,*) 'imax=',imax

      end subroutine index_quiet                                                           



      end module module_index_quiet

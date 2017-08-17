!Jan2011:original code was provided from Astrid from WACCM.
!Aug2011:this code was provided from Fei Wu from the WAM version.
!sep2012: efield.f was separated into each routin for SMS compatibility.
!--------------------------------------------  
      module module_svbksb
!--------------------------------------------------------------------- 
!   L.Scherliess ludger@gaim.cass.usu.edu
! - Weimer96
! Author: A. Maute Dec 2003  am 12/30/03 
!------------------------------------------------------------------------------ 

      implicit none

!nm20121003:module parameters are separated into efield.f!
      public :: svbksb

      contains
!-----------------------------------------------------------------------
! purpose: solves a*x = b
!
! method:     
! solves a*x = b for a vector x, where a is specified by the arrays
! u,w,v as returned by svdcmp. m and n
! are the logical dimensions of a, and will be equal for square matrices.
! mp and np are the physical dimensions of a. b(1:m) is the input right-hand 
! side. x(1:n) is the output solution vector. no input quantities are 
! destroyed, so the routine may be called sequentially with different b
!
! author:  a. maute dec 2002   
! (* copyright (c) 1985 numerical recipes software -- svbksb *!
! from numerical recipes 1986 pp. 57 or can be find on web-sites
!-------------------------------------------------------------------------      

      subroutine svbksb( u, w, v, m, n, mp, np, b, x )
!------------------------------------------------------------------------- 
!	... dummy arguments
!------------------------------------------------------------------------- 
      implicit none
      integer, parameter :: nmax = 1600
      integer, intent(in)   :: m
      integer, intent(in)   :: n
      integer, intent(in)   :: mp
      integer, intent(in)   :: np
      real , intent(in)  :: u(mp,np)
      real , intent(in)  :: w(np)
      real , intent(in)  :: v(np,np)
      real , intent(in)  :: b(mp)
      real , intent(out) :: x(np)

!------------------------------------------------------------------------- 
!	... local variables
!------------------------------------------------------------------------- 
      integer  :: i, j, jj
      real :: s
      real :: tmp(nmax)

      do j = 1,n
        s = 0. 
        if( w(j) /= 0. ) then
          do i = 1,m
            s = s + u(i,j)*b(i)
          end do
          s = s/w(j)
        endif
        tmp(j) = s
      end do

      do j = 1,n
        s = 0. 
        do jj = 1,n
          s = s + v(j,jj)*tmp(jj)
        end do
        x(j) = s
      end do

      end subroutine svbksb

!-----------------------------------------------------------------------
      end module module_svbksb

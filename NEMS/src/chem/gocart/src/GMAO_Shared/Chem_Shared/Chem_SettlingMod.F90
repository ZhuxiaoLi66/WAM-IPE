!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  Chem_SettlingMod --- Gravitional Sedimentation & Settling Speed
!
! !INTERFACE:
!

   module  Chem_SettlingMod

! !USES:

   use Chem_Mod
   use Chem_ConstMod, only: grav        ! Constants !
   use Chem_UtilMod

   use m_mpout

   implicit none

! !PUBLIC TYPES:
!
   PRIVATE

!
! !PUBLIC MEMBER FUNCTIONS:
!

   PUBLIC  Chem_Settling
   PUBLIC  Chem_CalcVsettle

!
! !DESCRIPTION:
!
!  This module implements the (pre-ESMF) DU Grid Component. 
!
! !REVISION HISTORY:
!
!  16Sep2003 da Silva  First crack.
!
!EOP
!-------------------------------------------------------------------------
CONTAINS
!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Chem_Settling - 
!
! !INTERFACE:
!

   subroutine Chem_Settling ( i1, i2, j1, j2, km, nbeg, nend, nbins, flag, &
                              radiusInp, rhopInp, cdt, w_c, tmpu, rhoa, &
                              hsurf, hghte, fluxout, rc, &
                              vsettleOut, correctionMaring )

! !USES:

  implicit NONE

! !INPUT PARAMETERS:

   integer, intent(in) :: i1, i2, j1, j2, km, nbeg, nend, nbins
   integer, intent(in) :: flag     ! flag to control particle swelling (see note)
   real, intent(in)    :: cdt
   real, pointer, dimension(:)     :: radiusInp, rhopInp
   real, pointer, dimension(:,:)   :: hsurf
   real, pointer, dimension(:,:,:) :: tmpu, rhoa, hghte

! !OUTPUT PARAMETERS:

   type(Chem_Bundle), intent(inout) :: w_c        ! Chemical tracer fields
   type(Chem_Array), pointer        :: fluxout(:) ! Mass lost by settling
                                                  ! to surface, kg/m2/s
   integer, intent(out)             :: rc         ! Error return code:
                                                  !  0 - all is well
                                                  !  1 - 
!  Optionally output the settling velocity calculated
   type(Chem_Array), pointer, optional, dimension(:)  :: vsettleOut

!  Optionally correct the settling velocity following Maring et al, 2003
   logical, optional, intent(in)    :: correctionMaring

   character(len=*), parameter :: myname = 'Settling'

! !DESCRIPTION: Gravitational settling of aerosol between vertical
!               layers.  Assumes input radius in [m] and density (rhop) 
!               in [kg m-3]. If flag is set, use the Fitzgerald 1975 (flag = 1)
!               or Gerber 1985 (flag = 2) parameterization to update the 
!               particle radius for the calculation (local variables radius
!               and rhop).
!
! !REVISION HISTORY:
!
!  17Sep2004  Colarco   Strengthen sedimentation flux out at surface
!                       by setting removal to be valid from middle of
!                       surface layer
!  06Nov2003  Colarco   Based on Ginoux
!  23Jan2003  da Silva  Standardization
!
!EOP
!-------------------------------------------------------------------------

! !Local Variables
   integer  ::  i, j, k, iit, n
   real, parameter ::  rhow = 1000.  ! Density of water [kg m-3]
   real :: pm(i1:i2,j1:j2,km)                    ! midpoint air pressure [pa]
   real :: dz(i1:i2,j1:j2,km)                    ! layer thickness [m]
   real :: dc(i1:i2,j1:j2)                       ! change in mr due to sfc loss
!   real :: pdog(i1:i2,j1:j2)                     ! air mass factor dp/g [kg m-2]
!   real :: pdog_m1(i1:i2,j1:j2)                  ! air mass factor dp/g [kg m-2]
   real :: vsettle(i1:i2,j1:j2,km)               ! fall speed [m s-1]
   real :: q_save(i1:i2,j1:j2)                   ! save the dust mmr [kg kg-1]
!   real :: q_before(i1:i2,j1:j2)                 ! save the dust mmr [kg kg-1]
   real :: cmass_before(i1:i2,j1:j2)
   real :: cmass_after(i1:i2,j1:j2)
   real :: diff_coef                 ! Brownian diffusion coefficient [m2 s-1]
   real :: pdog                      ! ratio of air mass factors
   real :: q_before                  ! save the dust mmr [kg kg-1]

!  The following parameters relate to the swelling of seasalt like particles
!  following Fitzgerald, Journal of Applied Meteorology, 1975.
   real, parameter :: epsilon = 1.   ! soluble fraction of deliqeuscing particle
   real, parameter :: alphaNaCl = 1.35
   real :: alpha, alpha1, alpharat, beta, theta, f1, f2

!  parameter from Gerber 1985 (units require radius in cm, see rcm)
   real :: rcm
   real, parameter :: c1=0.7674, c2=3.079, c3=2.573e-11, c4=-1.424
!  parameters for ammonium sulfate
   real, parameter :: SU_c1=0.4809, SU_c2=3.082, SU_c3=3.110e-11, SU_c4=-1.428


!  parameters from Maring et al, 2003
   real, parameter :: v_upwardMaring = 0.33e-2   ! upward velocity, [m s-1]
   real, parameter :: diameterMaring = 7.30e-6   ! particle diameter, [m]

!
   real :: sum
   real :: sat, rrat
   real :: radius, rhop   ! particle radius and density passed to
                          ! fall velocity calculation
   real    :: dt_settle, minTime, qmin, qmax, mydata
   integer :: nSubSteps, dk, ijl

   rc = 0

   ijl  = ( i2 - i1 + 1 ) * ( j2 - j1 + 1 )

!  Loop over vertical coordinate to get pressures, etc. (only need to do once)
   pm(i1:i2,j1:j2,1) = w_c%grid%ptop + 0.5*w_c%delp(i1:i2,j1:j2,1)
   do k = 2, km
    pm(i1:i2,j1:j2,k) = pm(i1:i2,j1:j2,k-1) &
         +0.5*(w_c%delp(i1:i2,j1:j2,k)+w_c%delp(i1:i2,j1:j2,k-1))
   end do

!  Handle the fact that hghte may be in the range [1,km+1] or [0,km]
!  -----------------------------------------------------------------
   dk = lbound(hghte,3) - 1  ! This is either 0 or 1
  
!  Layer thickness from hydrostatic equation
   k = km
   dz(:,:,k) = hghte(:,:,k+dk)-hsurf(:,:)
   do k = km-1, 1, -1
    dz(:,:,k) = hghte(:,:,k+dk) - hghte(:,:,k+dk+1)
   enddo

!  Loop over the number of dust bins
   do n = 1, nbins

    radius = radiusInp(n)
    rhop = rhopInp(n)

!   Reset a (large) minimum time to cross a grid cell in settling
    minTime = cdt

    if( associated(fluxout(n)%data2d) ) fluxout(n)%data2d(i1:i2,j1:j2) = 0.0
    cmass_before(:,:) = 0.0
    cmass_after(:,:) = 0.0

!   If radius le 0 then get out of loop
    if(radius .le. 0.) cycle

    do k = 1, km
     do j = j1, j2
      do i = i1, i2

!      Find the column dry mass before sedimentation
       cmass_before(i,j) = cmass_before(i,j) &
        + w_c%qa(nbeg+n-1)%data3d(i,j,k)*w_c%delp(i,j,k)/grav

!      Adjust the particle size for relative humidity effects
       sat = max(w_c%rh(i,j,k),tiny(1.0)) ! to avoid zero FPE

!      Fitzgerald
       if(flag .eq. 1 .and. sat .ge. 0.80) then
!       parameterization blows up for RH > 0.995, so set that as max
!       rh needs to be scaled 0 - 1
        sat = min(0.995,sat)
!       Calculate the alpha and beta parameters for the wet particle
!       relative to amonium sulfate
        beta = exp( (0.00077*sat) / (1.009-sat) )
        if(sat .le. 0.97) then
         theta = 1.058
        else
         theta = 1.058 - (0.0155*(sat-0.97)) /(1.02-sat**1.4)
        endif
        alpha1 = 1.2*exp( (0.066*sat) / (theta-sat) )
        f1 = 10.2 - 23.7*sat + 14.5*sat**2.
        f2 = -6.7 + 15.5*sat - 9.2*sat**2.
        alpharat = 1. - f1*(1.-epsilon) - f2*(1.-epsilon**2.)
        alpha = alphaNaCl * (alpha1*alpharat)
!       radius is the radius of the wet particle
        radius = alpha * radiusInp(n)**beta
        rrat = (radiusInp(n)/radius)**3.
        rhop = rrat*rhopInp(n) + (1.-rrat)*rhow
       elseif(flag .eq. 2) then   ! Gerber
        sat = min(0.995,sat)
        rcm = radiusInp(n)*100.
        radius = 0.01 * (   c1*rcm**c2 / (c3*rcm**c4-alog10(sat)) &
                          + rcm**3.)**(1./3.)
        rrat = (radiusInp(n)/radius)**3.
        rhop = rrat*rhopInp(n) + (1.-rrat)*rhow
       elseif(flag .eq. 3) then   
!       Gerber parameterization for Ammonium Sulfate
        sat = min(0.995,sat)
        rcm = radiusInp(n)*100.
        radius = 0.01 * (   SU_c1*rcm**SU_c2 / (SU_c3*rcm**SU_c4-alog10(sat)) &
                      + rcm**3.)**(1./3.)
        rrat = (radiusInp(n)/radius)**3.
        rhop = rrat*rhopInp(n) + (1.-rrat)*rhow
       elseif(flag .eq. 4) then
!       Petters and Kreidenweis (ACP2007) parameterization
        sat = min(0.99,sat)
        radius = (radiusInp(n)**3 * (1+1.19*sat/(1-sat)))**(1./3.)
        rrat = (radiusInp(n)/radius)**3
        rhop = rrat*rhopInp(n) + (1.-rrat)*rhow
       endif

!      Calculate the settling velocity
       call Chem_CalcVsettle(radius, rhop, rhoa(i,j,k), &
                        tmpu(i,j,k), diff_coef, vsettle(i,j,k))
      end do
     end do
    end do

    if(present(correctionMaring)) then
     if ((correctionMaring) .and. (radiusInp(n) .le. (0.5*diameterMaring))) then
       vsettle = max(1.0e-9, vsettle - v_upwardMaring)
     endif
    endif

    if(present(vsettleOut)) then
     vsettleOut(n)%data3d = vsettle
    endif

!   Determine global max/min time to cross grid cell
    call pmaxmin ( 'Chem_Settling: dt', dz(i1:i2,j1:j2,1:km)/vsettle(i1:i2,j1:j2,1:km), &
                                        qmin, qmax, ijl, km, 0. )
    minTime = min(minTime,qmin)


!   Now, how many iterations do we need to do?
    if ( minTime < 0 ) then
         nSubSteps = 0
         call mpout_log(myname,'no Settling because minTime = ', minTime )
    else if(minTime .ge. cdt) then
     nSubSteps = 1
     dt_settle = cdt
    else
     nSubSteps = cdt/minTime+1
     dt_settle = cdt/nSubSteps
    endif

!   Loop over sub-timestep
    do iit = 1, nSubSteps

!    Update dust mixing ratio (backward Euler scheme)
!    Note that you can not actually transport mixing ratio, so instead
!    we scale by the ratio of air mass factors between 2 levels
!    This approach should be "pretty" mass conservative: the flux
!    into a layer from above is calculated from the saved q
!    before the q is adjusted, so it is consistant.

!     q_save = w_c%qa(nbeg+n-1)%data3d(i1:i2,j1:j2,1)
!
!     w_c%qa(nbeg+n-1)%data3d(i1:i2,j1:j2,1) = w_c%qa(nbeg+n-1)%data3d(i1:i2,j1:j2,1) &
!            / (1.+dt_settle*vsettle(i1:i2,j1:j2,1)/dz(i1:i2,j1:j2,1))
     do j=j1,j2
       do i=i1,i2
         q_save(i,j) = w_c%qa(nbeg+n-1)%data3d(i,j,1) 
         w_c%qa(nbeg+n-1)%data3d(i,j,1) = q_save(i,j)/(1.+dt_settle*vsettle(i,j,1)/dz(i,j,1) )
       enddo
     enddo

     do k = 2, km


!     Air mass factors of layers k, k-1
!      pdog = w_c%delp(i1:i2,j1:j2,k)/grav
!      pdog_m1 = w_c%delp(i1:i2,j1:j2,k-1)/grav
!      q_before = w_c%qa(nbeg+n-1)%data3d(i1:i2,j1:j2,k)
!      w_c%qa(nbeg+n-1)%data3d(i1:i2,j1:j2,k) = &
!         1./(1.+dt_settle*vsettle(i1:i2,j1:j2,k)/dz(i1:i2,j1:j2,k)) &
!       * (       w_c%qa(nbeg+n-1)%data3d(i1:i2,j1:j2,k) &
!           + ( dt_settle*vsettle(i1:i2,j1:j2,k-1)/dz(i1:i2,j1:j2,k-1) &
!              * q_save*pdog_m1(i1:i2,j1:j2)/pdog(i1:i2,j1:j2) &
!             ) &
!         )

      do j=j1,j2
        do i=i1,j2
           pdog = w_c%delp(i,j,k-1)/w_c%delp(i,j,k)
           q_before = w_c%qa(nbeg+n-1)%data3d(i,j,k)
           w_c%qa(nbeg+n-1)%data3d(i,j,k) =                       &
             1./(1.+dt_settle*vsettle(i,j,k)/dz(i,j,k))           &
             *( q_before +( dt_settle*vsettle(i,j,k-1)/dz(i,j,k-1)  &
                     *q_save(i,j)*pdog ) )     
           q_save(i,j) = q_before
        enddo
      enddo

!      q_save = q_before

     end do ! k

    end do  ! iit

!   Find the column dry mass after sedimentation and thus the loss flux
    do k = 1, km
     do j = j1, j2
      do i = i1, i2
       cmass_after(i,j) = cmass_after(i,j) &
        + w_c%qa(nbeg+n-1)%data3d(i,j,k)*w_c%delp(i,j,k)/grav
      enddo
     enddo
    enddo

    if( associated(fluxout(n)%data2d) ) then
     do j = j1, j2
      do i = i1, i2
       fluxout(n)%data2d(i,j) &
        = (cmass_before(i,j) - cmass_after(i,j))/cdt
      enddo
     enddo
    endif

   end do   ! n

 end  subroutine Chem_Settling


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Chem_CalcVsettle - Calculate the aerosol settling velocity
!
! !INTERFACE:
!

   subroutine Chem_CalcVsettle ( radius, rhop, rhoa, tmpu, &
                                 diff_coef, vsettle )

! !USES:

  implicit NONE

! !INPUT PARAMETERS:

   real, intent(in)    :: radius              ! Particle radius [m]
   real, intent(in)    :: rhop                ! Particle density [kg m-3]
   real, intent(in)    :: rhoa                ! Layer air density [kg m-3]
   real, intent(in)    :: tmpu                ! Layer temperature [K]

! !OUTPUT PARAMETERS:

   real, intent(out)   :: diff_coef               ! Brownian diffusion 
                                                  ! coefficient [m2 s-1]
   real, intent(out)   :: vsettle                 ! Layer fall speed [m s-1]

   character(len=*), parameter :: myname = 'Vsettle'

! !DESCRIPTION: Calculates the aerosol settling velocity and Brownian diffusion
!               coefficient
!               Follows discussions in Seinfeld and Pandis, Pruppacher and
!               Klett, and the coding in CARMA (Toon et al., 1988)
!               Should work satisfactorily for al reasonable sized aerosols
!               (up to Reynolds number 300)
!
! !REVISION HISTORY:
!
!  06Nov2003  Colarco   Initial version.
!  23Jan2003  da Silva  Standardization
!
!EOP
!-------------------------------------------------------------------------

! !Local Variables
   real*8 rmu                       ! Dynamic viscosity [kg m-1 s-1]
   real*8 vt                        ! Thermal velocity of air molecule [m s-1]
   real*8 rmfp                      ! Air molecule mean free path [m]
   real*8 bpm                       ! Cunningham slip correction factor
   real*8 rkn                       ! Knudsen number
   real*8 re, x, y                  ! reynolds number and parameters
   real, parameter :: kb = 1.3807e-23 ! Boltzmann constant [kg m2 s-1 K-1 mol-1]
   real, parameter :: m_air = 4.8096e-26 ! Mass of <avg> air molecule [kg]
   real, parameter :: pi = 3.141529265

!  Dynamic viscosity from corrected Sutherland Equation
   rmu = 1.8325e-5*(416.16/(tmpu+120.))*(tmpu/296.16)**1.5

!  Thermal velocity of air molecule
   vt = sqrt(8.*kb*tmpu/pi/m_air)

!  Air molecule mean free path
   rmfp = 2.*rmu/rhoa/vt

!  Knudsen number
   rkn = rmfp/radius

!  Cunningham slip correction factor
   bpm = 1. + 1.246*rkn + 0.42*rkn*exp(-0.87/rkn)

!  Brownian diffusion coefficient
   diff_coef = kb*tmpu*bpm/3./pi/rmu/(2.*radius)

!  Fall speed (assumes Reynolds # < 0.01)
   vsettle = 2./9.*rhop*radius**2.*grav*bpm/rmu

!  Check the Reynold number to see if we need a drag correction
!  First guess at Reynold's number using Stoke's calculation
   re = 2.*rhoa*radius*vsettle/rmu

!  If Re > 0.01 then apply drag correction following Pruppacher and
!  Klett regime 2 (eq. 10-142).  Assuming reasonable aerosols we
!  do not consider that particle Re may exceed 300.
   if(re .gt. 0.01) then
    x = log(24.*re/bpm)
    y = -3.18657 + 0.992696   *x     - .00153193   *x**2. &
                 - 0.000987059*x**3. - .000578878  *x**4. &
                 + 8.55176E-05*x**5. -  3.27815E-06*x**6.
    re = exp(y)*bpm
    vsettle = rmu*re/2./rhoa/radius
   endif

#ifdef NEMS
!     Set a minimum value of settling velocity; and scale with a factor
      vsettle = 1.35 * max(vsettle,1.e-4)
#endif


   end subroutine Chem_CalcVsettle

   end module Chem_SettlingMod

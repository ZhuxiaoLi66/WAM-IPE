!dbg2010923: to reduce memory1: apexE etc
!dbg20110927: to reduce memory2:
! DATE: 08 September, 2011
!********************************************
!***      Copyright 2011 NAOMI MARUYAMA   ***
!***      ALL RIGHTS RESERVED             ***
!********************************************
! LICENSE AGREEMENT Ionosphere Plasmasphere Electrodynamics (IPE) model
! DEVELOPER: Dr. Naomi Maruyama
! CONTACT INFORMATION:
! E-MAIL : Naomi.Maruyama@noaa.gov
! PHONE  : 303-497-4857
! ADDRESS: 325 Broadway, Boulder, CO 80305
!--------------------------------------------  
      MODULE module_FIELD_LINE_GRID_MKS
      USE module_precision
      USE module_IPE_dimension,ONLY: NMP,NLP,ISTOT
      IMPLICIT NONE

! --- PRIVATE ---
!
! --- PUBLIC ---
!nm20121115      INTEGER (KIND=int_prec ),PUBLIC :: mp_save,lp_save,MaxFluxTube
      INTEGER (KIND=int_prec ),PUBLIC :: MaxFluxTube
      REAL    (KIND=real_prec8 ),PUBLIC :: minAltitude,maxAltitude
      REAL    (KIND=real_prec),PUBLIC :: minTheta   ,maxTheta
      REAL(KIND=real_prec),PARAMETER,PUBLIC :: ht90  = 90.0E+03  !reference height in meter
!... read in parameters
!SMS$DISTRIBUTE(dh,2,1) BEGIN
      INTEGER(KIND=int_prec),ALLOCATABLE :: JMIN_IN_all(:,:),JMAX_IS_all(:,:)!.. first and last indices on field line grid
!SMS$DISTRIBUTE END
      INTEGER(KIND=int_prec),ALLOCATABLE,PUBLIC :: JMIN_IN (:),JMAX_IS (:)   !.. first and last indices on field line grid
      INTEGER(KIND=int_prec),ALLOCATABLE,PUBLIC :: JMIN_ING(:),JMAX_ISG(:)   !.. first and last indices on field line grid
      INTEGER(KIND=int_prec),ALLOCATABLE,PUBLIC :: midpnt  (:)

      TYPE :: plasma_grid
         REAL(KIND=real_prec) :: SL !.. distance of point from northern hemisphere foot point [meter]
         REAL(KIND=real_prec) :: BM !.. magnetic field strength [T]
         REAL(KIND=real_prec) :: GR !.. Gravity [m2 s-1]
         REAL(KIND=real_prec) :: Q  
         REAL(KIND=real_prec) :: GCOLAT !.. geographic co-latitude [rad]
         REAL(KIND=real_prec) :: GLON   !.. geographic longitude [rad]
      END TYPE plasma_grid

      INTEGER (KIND=int_prec) :: ISL=1,IBM=2,IGR=3,IQ=4,IGCOLAT=5,IGLON=6
      INTEGER (KIND=int_prec) :: sendCount
!SMS$DISTRIBUTE(dh,,2) BEGIN
      REAL(KIND=real_prec),ALLOCATABLE,PUBLIC :: plasma_mp (:,:,:  )!MaxFluxTube,MPendMax,ISTOT
!SMS$DISTRIBUTE END

      REAL(KIND=real_prec),ALLOCATABLE,PUBLIC :: plasma_mpG(:,:,:,:)!MaxFluxTube,MPendMax,ISTOT,Ntot G means Global
!SMS$DISTRIBUTE(dh,2,3) BEGIN
      REAL(KIND=real_prec),ALLOCATABLE,PUBLIC,TARGET :: plasma_grid_3d(:,:,:,:)!MaxFluxTube,NLP,NMP,6
      REAL(KIND=real_prec),ALLOCATABLE,PUBLIC,TARGET :: plasma_3d     (:,:,:,:)!MaxFluxTube,NLP,NMP,ISTOT
      REAL(KIND=real_prec),ALLOCATABLE,PUBLIC,TARGET :: plasma_3d_old (:,:,:,:)!MaxFluxTube,NLP,NMP,ISTOT
      REAL(KIND=real_prec),ALLOCATABLE,PUBLIC,TARGET :: WamField      (:,:,:,:)!MaxFluxTube,NLP,NMP,NSTOT=7
      REAL(KIND=real_prec),ALLOCATABLE,PUBLIC        :: vn_ms1_4output(:,:,:,:)!MaxFluxTube,NLP,NMP,3
      REAL(KIND=real_prec),ALLOCATABLE,PUBLIC        :: Un_ms1        (:,:,:,:)!MaxFluxTube,NLP,NMP,3:3 Ue1 Eq5.6 in magnetic frame last dim = apexD1-3
      REAL(KIND=real_prec),ALLOCATABLE,PUBLIC        :: apexD       (:,:,:,:,:)!MaxFluxTube,NLP,NMP,3,3.. Eq(3.8-10 ) Richmond 1995
      REAL(KIND=real_prec),ALLOCATABLE,PUBLIC        :: apexE       (:,:,:,:,:)!MaxFluxTube,NLP,NMP,3,2.. Eq(3.11-12) Richmond 1995
      REAL(KIND=real_prec),ALLOCATABLE,PUBLIC        :: apexDscalar (:,:,:)    !MaxFluxTube,NLP,NMP.. Eq(3.15) Richmond 1995
      REAL(KIND=real_prec),ALLOCATABLE,PUBLIC        :: l_mag       (:,:,:,:,:)!MaxFluxTube,NLP,NMP,3,2.. Eq(3.11-12) Richmond 1995
      REAL(KIND=real_prec),ALLOCATABLE,PUBLIC        :: hrate_mks3d   (:,:,:,:)!MaxFluxTube,NLP,NMP,7 each component of Neutral heating rate(eV/kg/s)
      REAL(KIND=real_prec),ALLOCATABLE,PUBLIC,dimension(:,:,:) :: ON_m3,HN_m3,N2N_m3,O2N_m3,HE_m3,N4S_m3,TN_k,TINF_k

      ! ------------- Coupling fields to be passed to WAM ------------- !
      REAL(KIND=real_prec),ALLOCATABLE,PUBLIC        :: eastward_momentum_tendency(:,:,:) ! MaxFluxTube,NLP,NMP
      REAL(KIND=real_prec),ALLOCATABLE,PUBLIC        :: northward_momentum_tendency(:,:,:) ! MaxFluxTube,NLP,NMP
      REAL(KIND=real_prec),ALLOCATABLE,PUBLIC        :: solar_heating_rate(:,:,:) ! MaxFluxTube,NLP,NMP
      REAL(KIND=real_prec),ALLOCATABLE,PUBLIC        :: joule_heating_rate(:,:,:) ! MaxFluxTube,NLP,NMP
      REAL(KIND=real_prec),ALLOCATABLE,PUBLIC        :: auroral_heating_rate(:,:,:) ! MaxFluxTube,NLP,NMP
      REAL(KIND=real_prec),ALLOCATABLE,PUBLIC        :: o2_dissociation_rate(:,:,:) ! MaxFluxTube,NLP,NMP
      ! ---------------------------------------------------------------
!
!SMS$DISTRIBUTE END


!V_ExB m/s at the apex height
!SMS$DISTRIBUTE(dh,1,2) BEGIN
      REAL(KIND=real_prec),ALLOCATABLE,PUBLIC,DIMENSION(:,:) :: Be3    ! .. [T] Eq(4.13) Richmond 1995: "Ed1, Ed2, and Be3 are constant along magnetic field lines" !nm20130830
      REAL(KIND=real_prec),ALLOCATABLE,PUBLIC,DIMENSION(:,:) :: VEXBup !DIMENSION(NLP,NMP)
      REAL(KIND=real_prec),ALLOCATABLE,PUBLIC,DIMENSION(:,:) :: VEXBe  !DIMENSION(NLP,NMP)
      REAL(KIND=real_prec),ALLOCATABLE,PUBLIC,DIMENSION(:,:) :: VEXBth !DIMENSION(NLP,NMP)
!SMS$DISTRIBUTE END

      REAL(KIND=real_prec),ALLOCATABLE,PUBLIC,TARGET :: plasma_grid_Z (:,:)  !.. altitude [meter] (MaxFluxTube,NLP)
      REAL(KIND=real_prec),ALLOCATABLE,PUBLIC,TARGET :: plasma_grid_mag_colat(:,:)  !.. magnetic co-latitude Eq(6.1) [rad]
      REAL(KIND=real_prec),ALLOCATABLE               :: r_meter2D     (:,:)  !.. distance from the center of the Earth[meter]
      REAL(KIND=real_prec),ALLOCATABLE,PUBLIC,TARGET :: mlon_rad      (:)    !mag longitude in [rad]

      REAL(KIND=real_prec),PARAMETER  ,PUBLIC        :: dlonm90km = 4.5 ![deg]
!------------
      REAL(KIND=real_prec), ALLOCATABLE,     PUBLIC ::  Pvalue(:)  !.. p coordinate (L-shell)


!-------------
      REAL    (KIND=real_prec),ALLOCATABLE, PUBLIC :: poleVal(:,:) ! MaxFluxTube,ISTOT
      INTEGER (KIND=int_prec) ,ALLOCATABLE, PUBLIC :: DISPLS(:),MPends(:),recvCounts(:) ! nprocs
!
! components (east, north, up) of base vectors
      TYPE :: geographic_coords
         REAL(KIND=real_prec) :: east
         REAL(KIND=real_prec) :: north
         REAL(KIND=real_prec) :: up
      END TYPE geographic_coords

      INTEGER (KIND=int_prec) :: east=1,north=2,up=3
!

!---------------------------
      END MODULE module_FIELD_LINE_GRID_MKS

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
!
      SUBROUTINE allocate_arrays ( switch )
      USE module_precision
      USE module_IPE_dimension,ONLY: NMP,NLP,ISTOT
      USE module_FIELD_LINE_GRID_MKS,ONLY: &
     & plasma_grid_3d,plasma_3d,plasma_mp,plasma_mpG,r_meter2D,ON_m3,HN_m3,N2N_m3,O2N_m3&
     &,apexD,apexE,VEXBup,VEXBe,VEXBth,MaxFluxTube,HE_m3,N4S_m3,TN_k,TINF_K,Un_ms1 &
     &,Be3, Pvalue, JMIN_IN, JMAX_IS,hrate_mks3d,midpnt &
     &,mlon_rad, plasma_grid_Z, plasma_grid_GL, plasma_3d_old &
     &,apexDscalar, l_mag, WamField, poleVal,DISPLS,MPends,recvCounts &
     &,vn_ms1_4output
  
      USE module_input_parameters,ONLY: sw_neutral_heating_flip,mpHaloSize,nprocs &
!nm20170424 wind output corrected
     &, sw_neutral

      IMPLICIT NONE
      INTEGER (KIND=int_prec),INTENT(IN) :: switch
      INTEGER (KIND=int_prec) :: stat_alloc

! (0) ALLOCATE arrays
      IF ( switch==0 ) THEN
        print *,'ALLOCATing ARRAYS',ISTOT,MaxFluxTube,NLP,NMP
        allocate( plasma_grid_3d(MaxFluxTube,NLP,NMP,6    ) &
     &,           plasma_grid_Z (MaxFluxTube,NLP          ) &
     &,           plasma_grid_GL(MaxFluxTube,NLP          ) &
     &,           r_meter2D     (MaxFluxTube,NLP          ) &
     &,           plasma_3d     (MaxFluxTube,NLP,NMP,ISTOT) &
     &,           plasma_3d_old (MaxFluxTube,NLP,NMP,ISTOT) &
     &,           poleVal       (MaxFluxTube        ,ISTOT) &
     &,           apexD         (MaxFluxTube,NLP,NMP,3,1:3) &
     &,           apexE         (MaxFluxTube,NLP,NMP,3,2  ) &
     &,           apexDscalar   (MaxFluxTube,NLP,NMP      ) &
     &,           l_mag         (MaxFluxTube,NLP,NMP,3,2  ) )

!---neutral

        allocate( ON_m3 (MaxFluxTube,NLP,NMP)     &
     &,           HN_m3 (MaxFluxTube,NLP,NMP)     &
     &,           N2N_m3(MaxFluxTube,NLP,NMP)     &
     &,           O2N_m3(MaxFluxTube,NLP,NMP)     &
     &,           HE_m3 (MaxFluxTube,NLP,NMP)     &
     &,           N4S_m3(MaxFluxTube,NLP,NMP)     &
     &,           TN_k  (MaxFluxTube,NLP,NMP)     &
     &,           TINF_K(MaxFluxTube,NLP,NMP)     &
     &,           Un_ms1(MaxFluxTube,NLP,NMP,3:3) &
     &,           vn_ms1_4output(MaxFluxTube,NLP,NMP,3) )

!nm20170424 wind output corrected
!t        if ( sw_neutral==0.or.sw_neutral==1 ) then
           allocate( WamField(MaxFluxTube,NLP,NMP,7) )
!t        end if

        IF ( sw_neutral_heating_flip==1 ) THEN
          ALLOCATE(hrate_mks3d(MaxFluxTube,NLP,NMP,7),STAT=stat_alloc)
          IF ( stat_alloc==0 ) THEN
            print *,' hrate_mks3d ALLOCATION SUCCESSFUL!!!'
          ELSE !stat_alloc/=0
            print *,"!STOP hrate_mks3d ALLOCATION FAILED!:NHEAT",stat_alloc
            STOP
          END IF
        END IF !( sw_neutral_heating_flip==1 )

        ALLOCATE ( Be3     (  NLP,NMP  ) &
     &,            VEXBup  (  NLP,NMP  ) &
     &,            VEXBe   (  NLP,NMP  ) &
     &,            VEXBth  (  NLP,NMP  ) &
     &,            Pvalue  (  NLP      ) &
     &,            JMIN_IN (  NLP      ) &
     &,            JMAX_IS (  NLP      ) &
     &,            midpnt  (  NLP      ) &
!     &,            mlon_rad(      NMP+1) &
     &,            mlon_rad(1-mpHaloSize:NMP+mpHaloSize) &
     &,            STAT=stat_alloc       )
        IF ( stat_alloc==0 ) THEN
           print *,'ALLOCATION Be3 etc. SUCCESSFUL!!!'
        ELSE !stat_alloc/=0
           print *,switch,"!STOP! ALLOCATION Be3 etc. FAILED!:",stat_alloc
           STOP
        END IF

       allocate ( DISPLS(nprocs),MPends(nprocs),recvCounts(nprocs) &
     &,            STAT=stat_alloc ) 
       IF ( stat_alloc==0 ) THEN
          print *,'ALLOCATION using nprocs SUCCESSFUL!!!'
       ELSE !stat_alloc/=0
          print *,switch,"!STOP! ALLOCATION using nprocs FAILD!:",stat_alloc
          STOP
       END IF

!SMS$IGNORE BEGIN
       VEXBup = 0.0
       VEXBe  = 0.0
!SMS$IGNORE END

! (1) DEALLOCATE arrays
    ELSE IF ( switch==1 ) THEN
       print *,'DE-ALLOCATing ARRAYS'
! field line grid
       DEALLOCATE ( &
     &    plasma_grid_3d &
!---
     &,        apexD     &
!dbg20110923     &,        apexE &
!---
     &,          Be3     &
     &,       Pvalue     &
     &,      JMIN_IN     &
     &,      JMAX_IS     &
!---
     &,      mlon_rad    &
!---plasma
     &,  plasma_3d       &
     &,  plasma_3d_old   &
     &,  poleVal         &
     &,  VEXBup          &
     &,  VEXBe           )
 
       deallocate(           &
     &            ON_m3      &
     &,           HN_m3      &
     &,           N2N_m3     &
     &,           O2N_m3     &
     &,           HE_m3      &
     &,           N4S_m3     &
     &,           TN_k       &
     &,           TINF_K     &
     &,           Un_ms1     &
     &,           vn_ms1_4output )

!nm20170424 wind output corrected
!t       if ( sw_neutral==0.or.sw_neutral==1 ) then 
          DEallocate( WamField )
!t       end if

!---neutral heating
       IF ( sw_neutral_heating_flip==1 ) THEN
          DEALLOCATE ( hrate_mks3d )
       END IF !( sw_neutral_heating_flip==1 ) THEN

       deallocate( DISPLS,MPends,recvCounts )

    END IF !( switch==1 ) THEN

  END SUBROUTINE allocate_arrays

!nm20111118: heating rate output available!!!
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
MODULE module_heating_rate
      USE module_precision
      USE module_IPE_dimension,ONLY: NLP,NMP
      IMPLICIT NONE
!
      !.. EHT_cgs(3,J) = e heating rate, EHT_cgs(1,J) = ion heating rate, EHT_cgs(2,J) unused (eV cm-3 s-1)
!20110516: temporary moved to module plasma for the new structure for mpi test
!      REAL(KIND=real_prec), DIMENSION(3,NPTS2D,NMP), PUBLIC :: EHT_cgs
      !.. TE_TI_k(3,J) = Te, TE_TI_k(2,J) = Ti = TE_TI_k(2,J) [kelvin]
!dbg20110927      REAL(KIND=real_prec), DIMENSION(NPTS2D, NMP), PUBLIC :: NHEAT_cgs !.. Neutral heating rate (eV/cm^3/s) 
!tmp20110404: neutral heating rate in eV kg-1 s-1
!t      REAL(KIND=real_prec), DIMENSION(NPTS2D,NMP), PUBLIC :: NHEAT_mks !.. Neutral heating rate (eV/kg/s) 
!nm20111118: moved to module_FIELD_LINE_GRID_MKS
!nm20111118:      INTEGER,PUBLIC :: mp_save,lp_save !mp,lp values to be referred to from outside of subroutine plasma


      PRIVATE
      PUBLIC :: get_neutral_heating_rate

      CONTAINS
!-----------------------------
! convert the neutral heating rate from eV/cm3/s --> eV/kg/s -->J/kg/s
      SUBROUTINE get_neutral_heating_rate ( hrate_cgs, lp,mp )
      USE module_physical_constants,ONLY:mass_kg,AMU
      USE module_unit_conversion,ONLY:M3_to_CM3,eV2J
      USE module_FIELD_LINE_GRID_MKS,ONLY: plasma_grid_Z,MaxFluxTube,hrate_mks3d,ON_m3,HN_m3,N2N_m3,O2N_m3,HE_m3,N4S_m3,JMIN_IN,JMAX_IS
      IMPLICIT NONE
!------------------------
      DOUBLE PRECISION, DIMENSION(22,MaxFluxTube), INTENT(IN) :: hrate_cgs !.. each component of the Neutral heating rate (eV/kg/s) 
      INTEGER(KIND=int_prec),INTENT(IN) :: mp
      INTEGER(KIND=int_prec),INTENT(IN) :: lp
      REAL(KIND=real_prec) :: total_rho
      INTEGER(KIND=int_prec) :: i,jth

!nm20121020!SMS$PARALLEL(dh, lp) BEGIN
!nm20121020      do lp=1,NLP
        i_loop: DO i=JMIN_IN(lp),JMAX_IS(lp)
          IF(plasma_grid_Z(i,lp)>=80.E+03.AND.plasma_grid_Z(i,lp)<=700.E+03) THEN
            ! rho[kg cm-3] = n[cm-3] * mass[kg]
            total_rho = ( &
              &     ON_m3(i,lp,mp) * mass_kg(1) &
              &  +  HN_m3(i,lp,mp) * mass_kg(2) &
              &  + N2N_m3(i,lp,mp) * mass_kg(3) &
              &  + O2N_m3(i,lp,mp) * mass_kg(4) &
              &  +  HE_m3(i,lp,mp) * mass_kg(5) &
              &  + N4S_m3(i,lp,mp) * mass_kg(6) &
              & ) * M3_to_CM3 * AMU

! convert unit: eV/cm3/s --> J/cm3/s --> J/kg/s
 !t     NHEAT_mks(1:NPTS2D,NMP) = NHEAT_cgs(1:NPTS2D,NMP) * eV2J / total_rho(1:NPTS2D,NMP) 
! each component of the heating rate
            IF ( total_rho>0.0 ) THEN
              hrate_mks3d(i,lp,mp,1:7) = hrate_cgs(1:7,i) * eV2J / total_rho
            ELSE

              print *,'!STOP! INVALID total_rho',total_rho,i,lp,mp
              STOP
!  hrate_mks(1:7,i,lp) = 0.0

            END IF ! ( total_rho>0.0 ) THEN

          ELSE    ! IF(plasma_grid_Z(i)<80.or.plasma_grid_Z(i)>700) THEN

            hrate_mks3d(i,lp,mp,1:7) = 0.0

          END IF  !(plasma_grid_Z(i)>=80.AND.plasma_grid_Z(i)<=700) THEN      

        END DO i_loop ! DO i=JMIN_IN(lp),JMAX_IS(lp)
!nm20121020      end do ! lp=1,NLP
!nm20121020 !SMS$PARALLEL END
    END SUBROUTINE get_neutral_heating_rate
  
END MODULE module_heating_rate

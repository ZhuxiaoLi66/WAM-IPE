!nm20170427: implement gradually shift from msis to wam
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

      MODULE module_NEUTRAL_MKS            

      USE module_precision
      USE module_IPE_dimension
      USE module_FIELD_LINE_GRID_MKS
      USE module_input_parameters
      USE module_unit_conversion
      USE module_physical_constants

      IMPLICIT NONE

      PRIVATE
      PUBLIC :: neutral


      CONTAINS
!---------------------------
      subroutine neutral (utime_local, start_time) 


      implicit none

      INTEGER(KIND=int_prec), INTENT(in) :: utime_local  !universal time[sec]
      INTEGER(KIND=int_prec), INTENT(in) :: start_time

      integer :: npts  !=IPDIM
      INTEGER(KIND=int_prec) :: IN, IS

      real(8) :: ut_hour 
      real(8) :: f107D_dum, f107A_dum 
      real(8),              dimension(IPDIM)   :: glon_deg, glat_deg, alt_km
      real(8),              dimension(7)       :: AP_dum
      REAL(KIND=real_prec), dimension(3,IPDIM) :: Vn_ms1

! Neutral Wind vectors:
! Vn_ms1 : Neutral Wind on geographic frame:  1:east; 2:north; 3:up
! Un_ms1 : Neutral Wind on geomagnetic frame: 1:east; 2:down/equatorward; 3:parallel (South to North)

      INTEGER(KIND=int_prec) :: i, lp, mp, midpoint, ipts, ipts1
      REAL (KIND=real_prec)  :: dotprod
      INTEGER(KIND=int_prec) :: ihTopN,ihTopS,ihTop,jth
      INTEGER(KIND=int_prec) :: ihem,iStep,midPoints
      REAL (KIND=real_prec)  :: dht, dist, dist1
      REAL (KIND=real_prec)  :: scale_height_O, scale_height_O2, scale_height_N2
      REAL (KIND=real_prec)  :: Hk_O, Hk_O2, Hk_N2, Hk1_O, Hk1_O2, Hk1_N2

! Dummy variables for when IPE is coupled to WAM.....

      REAL (KIND=real_prec),dimension(IPDIM) :: ON_m3_dummy, N2N_m3_dummy, O2N_m3_dummy, TN_k_dummy, Tinf_k_dummy
      REAL (KIND=real_prec),dimension(3,IPDIM) :: Vn_ms1_dummy


!------

      if ( sw_ctip_input ) then
        LPI = INT( ( utime_local - utime0LPI ) / real(input_params_interval) ) + 1 + input_params_begin
      else
        LPI=1
      end if
      f107D_dum  = F107_new(lpi)
      f107A_dum  = F107d_new(lpi)
      AP_dum(1) = apa_eld(lpi)
      AP_dum(2) =  ap_eld(lpi)
      AP_dum(3) =  ap_eld(lpi-INT( 3.*60*60/input_params_interval))
      AP_dum(4) =  ap_eld(lpi-INT( 6.*60*60/input_params_interval))
      AP_dum(5) =  ap_eld(lpi-INT( 9.*60*60/input_params_interval))
      AP_dum(6) = apa_eld(lpi-INT(12.*60*60/input_params_interval))
      AP_dum(7) = apa_eld(lpi-INT(36.*60*60/input_params_interval))
      ut_hour = REAL(utime_local)/3600. !convert from sec to hours
     


!sw_neutral
!1: IPE coupled to WAM
!3: IPE utilizing MSIS (ie, standalone mode).


if (sw_neutral.eq.1) then

!---------------------------------------------------------------------------
! This block of code if IPE is using WAM for neutral parameters.....
!---------------------------------------------------------------------------

!SMS$PARALLEL(dh, lp, mp) BEGIN
      DO mp = 1,mpstop
        DO lp = 1,NLP

          IN = JMIN_IN(lp)
          IS = JMAX_IS(lp)
          NPTS = IS - IN + 1

          glon_deg(1:NPTS) = plasma_grid_3d(IN:IS,lp,mp,IGLON)*180./pi
          glat_deg(1:NPTS) = 90. - plasma_grid_3d(IN:IS,lp,mp,IGCOLAT)*180./pi
          alt_km  (1:NPTS) = plasma_grid_Z(IN:IS,lp) * M_TO_KM  !/ 1000. 

!---------------------------------------------------------------------------
! When coupled to WAM, IPE still needs neutral hydrogen and helium parameters
! from MSIS. Other parameters are 'dummies' (not needed)
!---------------------------------------------------------------------------

call get_thermosphere (npts, nyear, nday, ut_hour, f107D_dum, f107A_dum, AP_dum &   
                 , glon_deg, glat_deg, alt_km &
                 , he_m3( IN:IS,lp,mp) &
                 , on_m3_dummy(IN:IS) &
                 , o2n_m3_dummy(IN:IS) &
                 , n2n_m3_dummy(IN:IS) &
                 , hn_m3(IN:IS,lp,mp) &
                 , n4s_m3(IN:IS,lp,mp) &
                 , tn_k_dummy(IN:IS) &
                 , tinf_k_dummy(IN:IS) &
                 , Vn_ms1(1:3,1:NPTS))

        if( utime_local == start_time )then
        
          ! At the first timestep - the geographic wind comes from the value read in from
          ! the plasma parameters file....

          do jth=1,3
            do i=IN,IS
              Vn_ms1(jth,i-IN+1) = Vn_ms1_4output(i-IN+1,lp,mp,jth)
            end do
          end do

        endif

        midpoint = IN + (IS-IN)/2

        do ipts=in,midpoint
           if ( plasma_grid_Z(ipts,lp)<=mesh_height_max .and. mesh_height_max<plasma_grid_Z(ipts+1,lp) ) then
              ihTopN=ipts
              exit
           endif
           ! if midpoint < mesh_height_max[km]
           if ( ipts==midpoint) ihTopN = midpoint
        end do 
        
        do ipts=is,midpoint, -1
           if ( plasma_grid_Z(ipts,lp)<=mesh_height_max .and. mesh_height_max<plasma_grid_Z(ipts-1,lp) ) then
              ihTopS=ipts
              exit
           endif
           ! if midpoint < mesh_height_max[km]
           if ( ipts==midpoint) ihTopS = midpoint
        end do       

        if( simulation_is_warm_start .or. utime_local > start_time )then

!---------------------------------------------------------------------------
! The WamField array has indices of:
!           1:Tn; 2:Vn_geographic_east; 3:Vn_geographic_north; 4:Vn_geographic_up; 
!           5:O_density; 6:O2_density; 7:N2_density
!
! NH and SH denote Northern and Southern Hemispheres respectively
!---------------------------------------------------------------------------

          tn_k(IN:ihTopN,lp,mp)           =  WamField(IN:ihTopN,lp,mp,1)          !Tn < 800km NH
          tn_k(ihTopS:IS,lp,mp)           =  WamField(ihTopS:IS,lp,mp,1)          !Tn < 800km SH
          tn_k(ihTopN+1:midpoint  ,lp,mp) =  WamField(ihTopN,lp,mp,1)             !Tn >800km NH
          tn_k(midpoint+1:ihTopS-1,lp,mp) =  WamField(ihTopS,lp,mp,1)             !Tn >800km SH
          tinf_k(IN:midpoint  ,lp,mp)     =  WamField(ihTopN,lp,mp,1)             !Tn Inf NH 
          tinf_k(midpoint+1:IS,lp,mp)     =  WamField(ihTopS,lp,mp,1)             !Tn Inf SH
          
          Vn_ms1(1,IN-IN+1:ihTopN-IN+1)   =  WamField(IN:ihTopN,lp,mp,2)      !Vn_geographic_east < 800km NH
          Vn_ms1(1,ihTopS-IN+1:IS-IN+1)   =  WamField(ihTopS:IS,lp,mp,2)      !Vn_geographic_east < 800km SH
          Vn_ms1(2,IN-IN+1:ihTopN-IN+1)   =  WamField(IN:ihTopN,lp,mp,3)      !Vn_geographic_north < 800km NH
          Vn_ms1(2,ihTopS-IN+1:IS-IN+1)   =  WamField(ihTopS:IS,lp,mp,3)      !Vn_geographic_north < 800km SH
          Vn_ms1(3,IN-IN+1:ihTopN-IN+1)   =  WamField(IN:ihTopN,lp,mp,4)      !Vn_geographic_up < 800km NH
          Vn_ms1(3,ihTopS-IN+1:IS-IN+1)   =  WamField(ihTopS:IS,lp,mp,4)      !Vn_geographic_up < 800km SH
          
          on_m3( IN:ihTopN,lp,mp)         =  WamField(IN:ihTopN,lp,mp,5)          !O < 800km NH           
          on_m3( ihTopS:IS,lp,mp)         =  WamField(ihTopS:IS,lp,mp,5)          !O < 800km SH           
          o2n_m3( IN:ihTopN,lp,mp)        =  WamField(IN:ihTopN,lp,mp,6)          !O2 < 800km NH        
          o2n_m3( ihTopS:IS,lp,mp)        =  WamField(ihTopS:IS,lp,mp,6)          !O2 < 800km SH       
          n2n_m3( IN:ihTopN,lp,mp)        =  WamField(IN:ihTopN,lp,mp,7)          !N2 < 800km NH       
          n2n_m3( ihTopS:IS,lp,mp)        =  WamField(ihTopS:IS,lp,mp,7)          !N2 < 800km SH       

         DO ihem=1,2

          if ( ihem==1 ) then 
             istep=+1
             ihTop=ihTopN
             midPoints=midPoint
          else if ( ihem==2 ) then 
             istep=-1
             ihTop=ihTopS
             midPoints=midPoint+1               
          end if
 
          DO ipts=ihTop+istep, midPoints, iStep 
                     
                    Vn_ms1(1,ipts)                  =  WamField(ihTop,lp,mp,2)          !Vn_geographic_east > 800km 
                    Vn_ms1(2,ipts)                  =  WamField(ihTop,lp,mp,3)          !Vn_geographic_north > 800km 
                    Vn_ms1(3,ipts)                  =  WamField(ihTop,lp,mp,4)          !Vn_geographic_up > 800km   
                        

!---------------------------------------------------------------------------
! O, O2, and N2 densities above 800km decrease exponentially with height.....
!---------------------------------------------------------------------------

                        ipts1=ipts-iStep 

                        dist = earth_radius/(earth_radius+plasma_grid_Z(ipts,lp))
                        dist1 = earth_radius/(earth_radius+plasma_grid_Z(ipts1,lp))

                        Hk_O  = GSCON * Tn_k(ipts ,lp,mp) / (massn_kg(1)*G0*dist*dist)
                        Hk1_O = GSCON * Tn_k(ipts1,lp,mp) / (massn_kg(1)*G0*dist1*dist1)
                        Hk_O2  = GSCON * Tn_k(ipts ,lp,mp) / (massn_kg(2)*G0*dist*dist)
                        Hk1_O2 = GSCON * Tn_k(ipts1,lp,mp) / (massn_kg(2)*G0*dist1*dist1)
                        Hk_N2  = GSCON * Tn_k(ipts ,lp,mp) / (massn_kg(3)*G0*dist*dist)
                        Hk1_N2 = GSCON * Tn_k(ipts1,lp,mp) / (massn_kg(3)*G0*dist1*dist1)

                        scale_height_O  = (Hk_O+Hk1_O) / 2.0
                        scale_height_O2 = (Hk_O2+Hk1_O2) / 2.0
                        scale_height_N2 = (Hk_N2+Hk1_N2) / 2.0

                        dht = plasma_grid_Z(ipts1,lp) - plasma_grid_Z(ipts,lp)


                        on_m3(ipts,lp,mp)               =  on_m3(ipts1,lp,mp) * exp(dht/scale_height_O)   !O  > 800km 
                        o2n_m3(ipts,lp,mp)              =  o2n_m3(ipts1,lp,mp) * exp(dht/scale_height_O2) !O2 > 800km        
                        n2n_m3(ipts,lp,mp)              =  n2n_m3(ipts1,lp,mp) * exp(dht/scale_height_N2) !N2 > 800km            
                        
                        
                        !---------------------------------------------------------------------------
                        ! Finally, Make sure these densities do not go to zero at high altitudes....
                        !---------------------------------------------------------------------------
                        
                        if (on_m3(ipts,lp,mp).le.1.0e-12) on_m3(ipts,lp,mp) = 1.0e-12
                        if (o2n_m3(ipts,lp,mp).le.1.0e-12) o2n_m3(ipts,lp,mp) = 1.0e-12
                        if (n2n_m3(ipts,lp,mp).le.1.0e-12) n2n_m3(ipts,lp,mp) = 1.0e-12
                        
                      end do  
                    end do        

                   ! Write the geographic wind across to a 3-D parameter for output....

                   do jth=1,3
                     do i=IN,IS
                       Vn_ms1_4output(i-IN+1,lp,mp,jth)=Vn_ms1(jth,i-IN+1)
                     end do
                   end do


          endif ! selection between cold and warm start


          DO i=IN,IS
            ipts = i-IN+1 !1:NPTS

! un(3)=Ue3=d3*U: positive parallel to a field line, Eq(5.6) 
               dotprod = apexD(i,lp,mp,east ,3)*apexD(i,lp,mp,east ,3)  &
                    &  + apexD(i,lp,mp,north,3)*apexD(i,lp,mp,north,3) &
                    &  + apexD(i,lp,mp,up   ,3)*apexD(i,lp,mp,up   ,3)

               IF ( dotprod > 0.0 ) THEN
                  Un_ms1(i,lp,mp,3) = & 
                       &     ( apexD(i,lp,mp,east ,3)*Vn_ms1(1,ipts)     &
                       &     + apexD(i,lp,mp,north,3)*Vn_ms1(2,ipts)     &
                       &     + apexD(i,lp,mp,up   ,3)*Vn_ms1(3,ipts) ) / &
                       &     SQRT(  dotprod   )
               ELSE
                  Un_ms1(i,lp,mp,3) = 0.0
               END IF

!dbg20110131: the midpoint values become NaN otherwise because of inappropriate D1/3 values...
               IF ( lp>=1 .AND. lp<=6 .AND. i==midpoint )   Un_ms1(i,lp,mp,:) = Un_ms1(i-1,lp,mp,:) 

          END DO  !: DO i=IN,IS
        END DO !: DO lp = 1,NLP
      END DO   !: DO mp = 1,NMP
!SMS$PARALLEL END

      else if (sw_neutral.eq.3) then 

!---------------------------------------------------------------------------
! This block of code if IPE is using MSIS/HWM for neutral parameters.....
!---------------------------------------------------------------------------

!SMS$PARALLEL(dh, lp, mp) BEGIN
      apex_longitude_loop2: DO mp = 1,mpstop
        apex_latitude_height_loop2: DO lp = 1,NLP

          IN = JMIN_IN(lp)
          IS = JMAX_IS(lp)
          NPTS = IS - IN + 1

          glon_deg(1:NPTS) = plasma_grid_3d(IN:IS,lp,mp,IGLON)*180./pi
          glat_deg(1:NPTS) = 90. - plasma_grid_3d(IN:IS,lp,mp,IGCOLAT)*180./pi
          alt_km  (1:NPTS) = plasma_grid_Z(IN:IS,lp) * M_TO_KM  !/ 1000.

call get_thermosphere (npts, nyear, nday, ut_hour, f107D_dum, f107A_dum, AP_dum &
                 , glon_deg, glat_deg, alt_km &
                 , he_m3(IN:IS,lp,mp) &
                 , on_m3(IN:IS,lp,mp) &
                 , o2n_m3(IN:IS,lp,mp) &
                 , n2n_m3(IN:IS,lp,mp) &
                 , hn_m3(IN:IS,lp,mp) &
                 , n4s_m3(IN:IS,lp,mp) &
                 , tn_k(IN:IS,lp,mp) &
                 , tinf_k(IN:IS,lp,mp) &
                 , Vn_ms1(1:3,1:NPTS))

          do jth=1,3
            do i=IN,IS
              Vn_ms1_4output(i-IN+1,lp,mp,jth)=Vn_ms1(jth,i-IN+1)
            end do
          end do

          DO i=IN,IS
            ipts = i-IN+1 !1:NPTS

               dotprod = apexD(i,lp,mp,east ,3)*apexD(i,lp,mp,east ,3)  &
                    &  + apexD(i,lp,mp,north,3)*apexD(i,lp,mp,north,3) &
                    &  + apexD(i,lp,mp,up   ,3)*apexD(i,lp,mp,up   ,3)

               IF ( dotprod > 0.0 ) THEN
                  Un_ms1(i,lp,mp,3) = & 
                       &     ( apexD(i,lp,mp,east ,3)*Vn_ms1(1,ipts)     &
                       &     + apexD(i,lp,mp,north,3)*Vn_ms1(2,ipts)     &
                       &     + apexD(i,lp,mp,up   ,3)*Vn_ms1(3,ipts) ) / &
                       &     SQRT(  dotprod   )
               ELSE
                  Un_ms1(i,lp,mp,3) = 0.0
               END IF

               IF ( lp>=1 .AND. lp<=6 .AND. i==midpoint )   Un_ms1(i,lp,mp,:) = Un_ms1(i-1,lp,mp,:) 

          END DO  !: DO i=IN,IS

        END DO  apex_latitude_height_loop2
      END DO  apex_longitude_loop2
!SMS$PARALLEL END

      endif  !sw_neutral switch

      end subroutine neutral

      END MODULE module_NEUTRAL_MKS            

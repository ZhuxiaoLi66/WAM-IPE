!nm20130201: separated from interpolate_flux_tube_v18.f90
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
      MODULE module_Qinterpolation
      IMPLICIT NONE
      PRIVATE
      PUBLIC :: Qinterpolation

      CONTAINS
      SUBROUTINE Qinterpolation (mp,lp &
&, lp0,mp0 &
&, iR, Qint_dum, TSP ) !nm20170328: 

      USE module_precision
!     plasma_grid_3d,plasma_grid_Z,plasma_grid_GL,plasma_3d_old are all IN arrays
      USE module_FIELD_LINE_GRID_MKS,ONLY:JMIN_IN,JMAX_IS,plasma_grid_3d,plasma_grid_Z,plasma_grid_GL,ht90,ISL,IBM,IGR,IQ,IGCOLAT,IGLON,plasma_3d_old
      USE module_input_parameters,ONLY:sw_debug,mype,lps,lpe,mps,mpe,sw_ihepls,sw_inpls
      USE module_IPE_dimension,ONLY: ISPEC,ISPET,IPDIM ,nlp
      USE module_physical_constants,ONLY: earth_radius,pi,zero
      USE module_PLASMA,ONLY:utime_save
      IMPLICIT NONE
!--- INPUT ---
      INTEGER (KIND=int_prec),INTENT(IN) :: mp
      INTEGER (KIND=int_prec),INTENT(IN) :: lp
      INTEGER (KIND=int_prec),INTENT(IN) :: lp0,mp0
      INTEGER (KIND=int_prec),INTENT(IN) :: iR !=iB+1 !add R
!nm20170328:
      INTEGER (KIND=int_prec),INTENT(IN) :: TSP     !N(1:4) perp.transport
!--- OUTPUT ---
      REAL(KIND=real_prec8),INTENT(OUT) :: Qint_dum(iR,IPDIM)  !1d:species


!---local
      INTEGER (KIND=int_prec) :: ip,i,jth
      INTEGER (KIND=int_prec) :: ispecial,isouth,inorth,i1d,ip1d,midpoint ,kk
      REAL(KIND=real_prec8) :: factor2

!nm20170328      INTEGER (KIND=int_prec),PARAMETER :: TSP=3      !N(1:3) perp.transport 
     INTEGER (KIND=int_prec),PARAMETER :: iT=ISPEC+3   !add T(1:3)


      INTEGER (KIND=int_prec),PARAMETER :: iB=iT+1 !add B

!---




!array initialization:solved lahey err: variable undefined value
Qint_dum(:,:)=zero

!NOTE: grid point distributions are identical in all mp: between mp_t0(ihem,1)v.s.mp_t0(ihem,2)
! however Q values are not identical in all mp!!! 




    flux_tube_loopT1_Q: DO ip=JMIN_IN(lp),JMAX_IS(lp)
      ip1d=ip-JMIN_IN(lp)+1

!check the foot point values
!ispecial=2: interpolation at/below IN

!JFM  Will this be in the halo? plasma_grid_3d has it's halo set in module_read_plasma_grid_global by the SERIAL directive.
!     And the distance from lp and mp is checked against the hal size in perpendicular_transport.
       IF ( plasma_grid_3d( JMIN_IN(lp0) , lp0,mp0,IQ) < plasma_grid_3d(ip ,lp,mp,IQ) ) THEN
          ispecial=2
          isouth=JMIN_IN(lp0)+1  !not used!!!
          inorth=JMIN_IN(lp0)
          i1d  =+1

!ispecial=1: interpolation at/below IS
       ELSE IF ( plasma_grid_3d( JMAX_IS(lp0) , lp0,mp0,IQ) >= plasma_grid_3d(ip,lp,mp,IQ) ) THEN
          ispecial=1
          isouth=JMAX_IS(lp0)   !not used!!!
          inorth=JMAX_IS(lp0)-1
          i1d  =JMAX_IS(lp0)-JMIN_IN(lp0)+1

       ELSE

!search for north & south grid point of i1
          flux_tube_loopT0: DO i=JMIN_IN(lp0),JMAX_IS(lp0)

             IF ( plasma_grid_3d(i,lp0,mp0,IQ) < plasma_grid_3d(ip,lp,mp,IQ) ) THEN
! ispecial=0: normal interpolation
                ispecial = 0
                inorth = i-1
                isouth = i
                i1d   = i-JMIN_IN(lp0)+1
         
! factor2 IS NOT equal for all mp
                factor2 = ( plasma_grid_3d(ip,lp,mp,IQ)       - plasma_grid_3d(isouth,lp0,mp0,IQ) ) /&
                        & ( plasma_grid_3d(inorth,lp0,mp0,IQ) - plasma_grid_3d(isouth,lp0,mp0,IQ) )


                if ( factor2<zero.or.factor2>1.0) then
                   WRITE(6,*) 'sub-Intrp:!STOP! INVALID factor2',factor2 ,ip,lp,mp,isouth,lp0,mp0
                   STOP
                endif

                EXIT flux_tube_loopT0
             END IF

          END DO flux_tube_loopT0 !: DO i=IN,IS

       END IF ! IF ( plasma_grid_3d( JMIN_IN(lp0) , lp0,mp0,IQ) < plasma_grid_3d(ip ,lp,mp,IQ) ) THEN


!ispecial=0: normal interpolation
       if(ispecial == 0) then
! calculate all the ionospheric parameters      
!        ni1_in(ip)=(factor2*(ni(inorth,mp0,1) - ni(isouth,mp0,1))) + ni(isouth,mp0,1)


!N 1:TSP: density
!not sure if LOG is necessary for densities???
          jth_loop0:         DO jth=1,iT !=TSP+3

             IF ( jth>TSP.AND.jth<=ISPEC )  CYCLE jth_loop0
             !nm20170328: he+ ihepls<=0
             if ( jth==3.and.sw_ihepls<=0 ) CYCLE jth_loop0 
             !nm20170328: n+ inpls<=0
             if ( jth==4.and.sw_inpls<=0 ) CYCLE jth_loop0

             !nm20160421 DOUBLE log interpolation
             if(jth<=TSP)then      !for densities
                Qint_dum(jth, ip1d) = 10**( factor2*(ALOG10(plasma_3d_old(iNorth,lp0,mp0,jth)) - ALOG10(plasma_3d_old(iSouth,lp0,mp0,jth))) + ALOG10(plasma_3d_old(iSouth,lp0,mp0,jth)) )
             else 
                Qint_dum(jth, ip1d) = (factor2*(plasma_3d_old(inorth,lp0,mp0,jth) - plasma_3d_old(isouth,lp0,mp0,jth))) + plasma_3d_old(isouth,lp0,mp0,jth)
             end if

!dbg20160417 error stop
             if ( Qint_dum(jth, ip1d)<=zero ) then 
!SMS$IGNORE begin
                print*,mype,utime_save,'subQint(0):INVALID Qint !STOP!',Qint_dum(jth, ip1d) &
                     &,plasma_3d_old(inorth,lp0,mp0,jth) &
                     &,plasma_3d_old(isouth,lp0,mp0,jth) &
                     &,lp,mp,lp0,mp0,ip1d,inorth,isouth,factor2,jth
                if(jth<=TSP) print*,jth,'after LOG!',( factor2*(ALOG10(plasma_3d_old(iNorth,lp0,mp0,jth)) - ALOG10(plasma_3d_old(iSouth,lp0,mp0,jth))) + ALOG10(plasma_3d_old(iSouth,lp0,mp0,jth)) )
!SMS$IGNORE end
                STOP
             end if !( jth==1.and.Qint_dum(jth, ip1d)<=0. ) then 

          END DO jth_loop0 !jth=1,iT !=TSP+3

!B iT+1=iB: B magnetic field intensity
          Qint_dum(iB, ip1d) = (factor2*(plasma_grid_3d(inorth,lp0,mp0,IBM) - plasma_grid_3d(isouth,lp0,mp0,IBM))) + plasma_grid_3d(isouth,lp0,mp0,IBM)

!NOTE: R should be the same for all mp!!!
!R iB+1=iR: R = RE + Z
          Qint_dum(iR, ip1d) =( (factor2*(plasma_grid_Z(inorth,lp0) - plasma_grid_Z(isouth,lp0))) + plasma_grid_Z(isouth,lp0) ) +earth_radius


!ispecial=1: interpolation at/below IS_t0
       ELSE if(ispecial == 1) then

          jth_loop1: DO jth=1,iT
             IF ( jth>TSP.AND.jth<=ISPEC )  CYCLE jth_loop1
             !nm20170328: he+ ihepls<=0
             if ( jth==3.and.sw_ihepls<=0 ) CYCLE jth_loop1 
             !nm20170328: n+ inpls<=0
             if ( jth==4.and.sw_inpls<=0 ) CYCLE jth_loop1
             Qint_dum(jth   ,ip1d) = plasma_3d_old(isouth,lp0,mp0,jth)


!dbg20160417
if ( jth==1.and.Qint_dum(jth, ip1d)<=zero ) then 
!SMS$IGNORE begin
print *,mype,'sub-Qint(1)',lp,mp,lp0,mp0,jth,isouth, plasma_3d_old(isouth,lp0,mp0,jth)
!SMS$IGNORE end
STOP
endif




          END DO jth_loop1 !jth
          !B:
          Qint_dum(iB      ,ip1d) = plasma_grid_3d(isouth,lp0,mp0,IBM)
          !R:
          Qint_dum(iR      ,ip1d) = plasma_grid_Z(isouth,lp0) +earth_radius

!ispecial=2: interpolation at/below IN_t0
       ELSE if(ispecial == 2) then
          jth_loop2: DO jth=1,iT
             IF ( jth>TSP.AND.jth<=ISPEC )  CYCLE jth_loop2
             !nm20170328: he+ ihepls<=0
             if ( jth==3.and.sw_ihepls<=0 ) CYCLE jth_loop2 
             !nm20170328: n+ inpls<=0
             if ( jth==4.and.sw_inpls<=0 ) CYCLE jth_loop2
             Qint_dum(jth   ,ip1d) = plasma_3d_old(inorth,lp0,mp0,jth)
             !dbg20160417
             if ( jth==1.and.Qint_dum(jth, ip1d)<=zero ) then 
!SMS$IGNORE begin
               print *,mype,'sub-Qint(2)',lp,mp,lp0,mp0,jth,inorth, plasma_3d_old(inorth,lp0,mp0,jth)
!SMS$IGNORE end
               STOP
             endif
          END DO jth_loop2!jth
          !B
          Qint_dum(iB      ,ip1d) = plasma_grid_3d(inorth,lp0,mp0,IBM)
          !R
          Qint_dum(iR      ,ip1d) = plasma_grid_Z(inorth,lp0) +earth_radius
          
       END IF !(ispecial == 1) then
    END DO flux_tube_loopT1_Q !: DO ip=IN,IS


  END SUBROUTINE Qinterpolation
END MODULE module_Qinterpolation


MODULE module_init_plasma_grid
IMPLICIT NONE

      CONTAINS
!---------------------------
! initialise plasma grids
SUBROUTINE init_plasma_grid ( )
USE module_read_plasma_grid_global,only: read_plasma_grid_global
USE module_precision
USE module_IPE_dimension      ,ONLY: NMP,NLP,ISTOT
USE module_physical_constants ,ONLY: earth_radius, pi, G0,zero
USE module_input_parameters   ,ONLY: sw_grid,parallelBuild,mpHaloSize,nprocs,mype
USE module_FIELD_LINE_GRID_MKS,ONLY: Pvalue,JMIN_IN,JMAX_IS, r_meter2D      &
&, plasma_grid_mag_colat,plasma_grid_3d,apexD,apexE,Be3,plasma_grid_Z,ISL,IBM,IGR  &
& ,IQ,IGCOLAT,IGLON,east,north,up,mlon_rad,dlonm90km,apexDscalar,l_mag
!sms$insert USE module_prepPoleVal,ONLY: prepPoleVal
!USE module_cal_apex_param,ONLY:cal_apex_param

INTEGER (KIND=int_prec)           :: i,mp,lp,in,is
REAL    (KIND=real_prec)          :: sinI
INTEGER (KIND=int_prec),parameter :: sw_sinI=0  !0:flip; 1:APEX
INTEGER (KIND=int_prec)           :: midpoint
!local
      REAL    (KIND=real_prec)              :: ufac
      REAL    (KIND=real_prec),DIMENSION(3) :: bhat !eq(3.14)
      REAL    (KIND=real_prec),DIMENSION(3) :: a,b,c      
!dbg20130814
      INTEGER (KIND=int_prec)               :: ii



!if the new GLOBAL 3D version
CALL read_plasma_grid_global

! make sure to use the MKS units.
print *,"Z_meter calculation completed"

Pvalue = zero
do lp=1,NLP
  IN = JMIN_IN(lp)
  CALL Get_Pvalue_Dipole ( r_meter2D(IN,lp), plasma_grid_mag_colat(IN,lp), Pvalue(lp) )
enddo

!SMS$PARALLEL(dh, lp, mp) BEGIN
apex_longitude_loop: DO mp = 1,NMP

!.. p coordinate (L-shell) is a single value along a flux tube 
!NOTE: in FLIP, PCO is only used in setting up the rough plasmasphere H+ initial profiles (See PROFIN). It does not have to be accurate.

!dbg20120112:      Pvalue(:) = zero
    apex_latitude_height_loop:   DO lp = 1,NLP

      IN = JMIN_IN(lp)
      IS = JMAX_IS(lp)

      midpoint = IN + ( IS - IN )/2

IF ( sw_grid==0 ) THEN  !APEX
! assuming Newtonian gravity: G0 is gravity at the sea level (z=0) 
!NOTE: positive in NORTHern hemisphere; negative in SOUTHern hemisphere
         flux_tube: DO i=IN,IS

!nm20130201: need to calculate some more apex parameters 
!           CALL cal_apex_param (i,lp,mp,sinI)
!---
!nm20130814: double check if apexD has the values...
            if ( apexD(i,lp,mp,east,1)==0.0.AND.apexD(i,lp,mp,east,2)==0.0 ) then

               if ( i==midpoint ) then
                  ii = i-1  !assign the northward neighboring value
               else
                  print *,'sub-init_plasma_grid: STOP! INVALID apexD!',i,lp,mp,apexD(i,lp,mp,:,1),apexD(i,lp,mp,:,2)
                  STOP 
               end if

            else
              ii = i
           end if

! calculate D from eq 3.15: | d1 X d2 |
           a(1) = apexD(ii,lp,mp,east,1)
           a(2) = apexD(ii,lp,mp,north,1)
           a(3) = apexD(ii,lp,mp,up,1)
           b(1) = apexD(ii,lp,mp,east,2)
           b(2) = apexD(ii,lp,mp,north,2)
           b(3) = apexD(ii,lp,mp,up,2)

!      call cross_product (a,b,c)
!---
! calculate cross product   a X b 


           c(1) = a(2)*b(3) - a(3)*b(2) 
           c(2) = a(3)*b(1) - a(1)*b(3) 
           c(3) = a(1)*b(2) - a(2)*b(1) 
!---

           apexDscalar(i,lp,mp) = &
                &     ABS ( &
                & c(1)*c(1) + c(2)*c(2) + c(3)*c(3) &
                & )


! calculate bhat
           bhat(east)  = apexD(ii,lp,mp,east, 3) * apexDscalar(i,lp,mp)
           bhat(north) = apexD(ii,lp,mp,north,3) * apexDscalar(i,lp,mp)
           bhat(up)    = apexD(ii,lp,mp,up,   3) * apexDscalar(i,lp,mp)

           sinI = -bhat(up) !output

!dbg           print *,i,lp,mp,'bhat',bhat
           ufac  = SQRT(bhat(north)**2 + bhat(east)**2)

!dbg           print *,'ufac',ufac

! l_mag: unit vector
!(1) magnetic eastward exactly horizontal
!    l_e = bhat x k /|bhat x k|
           if ( ufac > 0 ) then
              l_mag(i,lp,mp,east ,1) =  bhat(north)/ufac
              l_mag(i,lp,mp,north,1) = -bhat(east) /ufac
              l_mag(i,lp,mp,up   ,1) =  zero
           else 
              print *,'sub-init_plasma_grid: STOP! INVALID ufac!',ufac
              STOP
           endif
      
!i and j are longitude and latitude index, bhat is the unit vector
!in direction of the geomagnetic filed line from subroutine apxmall,
!k is unit vector in vertical

!unit vector in upward direction
!(2) magnetic upward exactly horizontal
! l_u = l_e x bhat
           l_mag(i,lp,mp,east, 2) = l_mag(i,lp,mp,north,1)*bhat(up)   - l_mag(i,lp,mp,up   ,1)*bhat(north)
           l_mag(i,lp,mp,north,2) = l_mag(i,lp,mp,up   ,1)*bhat(east) - l_mag(i,lp,mp,east ,1)*bhat(up)
           l_mag(i,lp,mp,up,   2) = l_mag(i,lp,mp,east ,1)*bhat(north)- l_mag(i,lp,mp,north,1)*bhat(east)
!with the vector l you can use the formula in Art's paper.

!inserted from
!      END SUBROUTINE  cal_apex_param
!      END MODULE module_cal_apex_param


!dbg20110831
!d print *,'calling sub-Get_sinI'&
!d &, i,lp,mp,sw_sinI, sinI&
!d &, plasma_grid_mag_colat(i), apexD(3,i,mp)%east, apexD(3,i,mp)%north, apexD(3,i,mp)%up

           CALL Get_sinI ( sw_sinI, sinI, plasma_grid_mag_colat(i,lp) &
     &, apexD(i,lp,mp,east,3), apexD(i,lp,mp,north,3), apexD(i,lp,mp,up,3) ) 
           plasma_grid_3d(i,lp,mp,IGR)  =  G0 * ( earth_radius * earth_radius ) / ( r_meter2D(i,lp) * r_meter2D(i,lp) ) * sinI * (-1.0)

         END DO flux_tube

!nm20120304: introducing the flip grid
!reminder:
!(1) neutral wind should be calculated using sinI from flip_grid: "SINDIP"
ELSE IF ( sw_grid==1 ) THEN  !FLIP
  if(parallelBuild) then
    print*,'sw_grid=1 does not work for a parallel run'
    print*,'Stopping in module_init_plasma_grid'
    STOP
  endif
  print *,'calling get_FLIP_grid',mp,lp
  CALL get_flip_grid (mp,lp)
END IF !( sw_grid==0 ) THEN  !APEX

       END DO apex_latitude_height_loop   !: DO lp = 1,NLP
     END DO apex_longitude_loop         !: DO mp = 1,NMP 

!SMS$PARALLEL END

!sms$insert call prepPoleVal


     mlon_rad(:) = zero
!nm20160419
!     DO mp = 1,NMP+1
     DO mp = 1-mpHaloSize,NMP+mpHaloSize
       mlon_rad(mp) = REAL( (mp-1),real_prec ) * dlonm90km *pi/180.00
 !print*,'mp=',mp,'mlon=',mlon_rad(mp)
     END DO

END SUBROUTINE init_plasma_grid

END MODULE module_init_plasma_grid

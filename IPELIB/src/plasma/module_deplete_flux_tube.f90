module module_deplete_flux_tube
PRIVATE
PUBLIC ::  deplete_flux_tube
CONTAINS
subroutine deplete_flux_tube ( utime, mp, lp, HPEQ )
      USE module_precision
      USE module_input_parameters,ONLY: HPEQ_flip,start_time, sw_depleted_flip, start_time_depleted
      USE module_FIELD_LINE_GRID_MKS,ONLY:JMIN_IN,plasma_grid_GL
      USE module_physical_constants,ONLY:pi,zero
IMPLICIT NONE
!
      INTEGER (KIND=int_prec), INTENT(IN) :: utime !universal time [sec]
      INTEGER (KIND=int_prec), INTENT(IN) :: mp    !longitude
      INTEGER (KIND=int_prec), INTENT(IN) :: lp    !latitude
      DOUBLE PRECISION HPEQ
!--local
      REAL (KIND=real_prec), PARAMETER :: mlat_boundary = 53.9 ![deg]
      REAL (KIND=real_prec) :: mlat !deg

!
!nm20140729 original
      IF ( sw_depleted_flip==0 ) THEN 

         IF ( utime == start_time ) THEN
            HPEQ      = HPEQ_flip
         ELSE
            HPEQ      = 0.0
         END IF
         RETURN

      ELSE IF ( sw_depleted_flip==1 ) THEN 

!nm20141027: deplete flux tube ONLY ONCE at start_time_depleted
!nm20141027         IF ( utime < start_time_depleted ) THEN
         IF ( utime /= start_time_depleted ) THEN
            HPEQ = zero

            IF ( utime == start_time ) THEN
               HPEQ      = HPEQ_flip
            END IF
           
            RETURN

!nm20141027         ELSE !       IF ( utime >= start_time_depleted ) THEN
         ELSE !       IF ( utime == start_time_depleted ) THEN

!(1) fixed boundary location=75deg (as an example)
! (2.1) HPEQ=0.8 !the rate of the depletion

            mlat = ( 90. - plasma_grid_GL(JMIN_IN(lp),lp) * 180./pi )
            IF ( mlat >= mlat_boundary ) THEN
               HPEQ      = -0.2
            ELSE
               HPEQ      = zero
            END IF
print *,'utime=',utime,' mp=', mp,' lp=',lp,' mlat=',mlat,' HPEQ=',HPEQ,sw_depleted_flip, start_time, start_time_depleted,mlat_boundary
            RETURN
         END IF !( utime >= start_time_depleted ) THEN
      END IF!( sw_depleted_flip==0 ) THEN 
!
end subroutine deplete_flux_tube
end module module_deplete_flux_tube

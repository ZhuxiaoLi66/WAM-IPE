MODULE module_calcPoleVal
USE module_precision
USE module_IPE_dimension,ONLY: NLP,NMP,ISTOT
IMPLICIT NONE
include "gptl.inc"

CONTAINS
!nm20160420: store special pole values in poleVal
      SUBROUTINE calcPoleVal
      USE module_input_parameters   ,ONLY: SMScomm,NumPolevalProcs,sendCount,barriersOn
      USE module_FIELD_LINE_GRID_MKS,ONLY: plasma_3d,plasma_mp,plasma_mpG,JMIN_IN,JMAX_IS,poleVal,MPends,DISPLS,recvCounts
      INCLUDE "mpif.h"
      REAL    (KIND=real_prec) :: sumpv
      INTEGER (KIND=int_prec)  :: jth,kount,n,status,ret,mp1,mp,lp,i

!SMS$PARALLEL(dh, lp, mp) BEGIN
      ret = gptlstart ('calcPoleVal_plasma_mp')
      do jth=1,ISTOT
        mp1 = 0
        do mp = 1,NMP
          mp1 = mp1 + 1
          do lp=1,1
            do i=JMIN_IN(1),JMAX_IS(1)
              plasma_mp(i,mp1,jth) = plasma_3d(i,lp,mp,jth)
            enddo
          enddo
        enddo
      enddo
      ret = gptlstop  ('calcPoleVal_plasma_mp')
!SMS$PARALLEL END

      if(barriersOn) then
        ret = gptlstart ('calcPoleVal_gatherv_barrier')
!sms$insert      call ppp_barrier(status)
        ret = gptlstop  ('calcPoleVal_gatherv_barrier')
      endif
      ret = gptlstart ('calcPoleVal_gatherv')
      call MPI_GATHERV(plasma_mp, sendCount, MPI_REAL, plasma_mpG, recvCounts, DISPLS, MPI_REAL, 0, SMScomm, status)
      if(status /= 0) then
!SMS$ignore begin
        print*,'MPI_GATHERV error in calcPoleVal',status
        stop
      endif
!SMS$ignore end
        
      ret = gptlstop ('calcPoleVal_gatherv')

      if(barriersOn) then
        ret = gptlstart ('calcPoleVal_serial_barrier')
!sms$insert      call ppp_barrier(status)
        ret = gptlstop  ('calcPoleVal_serial_barrier')
      endif
      ret = gptlstart ('calcPoleVal_serial')
!SMS$SERIAL(<poleVal,OUT>:default=ignore) BEGIN
      ret = gptlstart ('poleVal')
      do jth=1,ISTOT
        do i=JMIN_IN(1),JMAX_IS(1)
            kount = 0
            sumpv = 0.0
          do n=1,NumPolevalProcs
            do mp=1,MPends(n)
              sumpv = sumpv + plasma_mpG(i,mp,jth,n)
              kount = kount + 1
            enddo
          enddo
          if(kount /= NMP) then
            print*,'Error in module_sub_plasms kount /= NMP',kount,NMP
            stop
          endif
          poleVal(i,jth) = sumpv / REAL(NMP) !sum on the MP (2nd) dimension
        end do
      end do 
      ret = gptlstop ('poleVal')
!SMS$SERIAL END
      ret = gptlstop  ('calcPoleVal_serial')
      END SUBROUTINE calcPoleVal

END MODULE module_calcPoleVal

MODULE module_prepPoleVal
USE module_precision
USE module_IPE_dimension,ONLY: NLP,NMP,ISTOT
IMPLICIT NONE

CONTAINS
!Initialize for MPI_gatherV in module_sub_plasma.f90 to calculate poleVal
  SUBROUTINE prepPoleVal
  USE module_input_parameters   ,ONLY: nprocs,SMScomm,NumPolevalProcs,sendCount,NumPolevalProcs,mype
  USE module_FIELD_LINE_GRID_MKS,ONLY: JMIN_IN,JMAX_IS,DISPLS,MPends,recvCounts,plasma_mp,plasma_mpG,MaxFluxTube
!SMS$INSERT INCLUDE "mpif.h"
  INTEGER (KIND=int_prec) :: jth,NumMP,i,mp,lp,status,PE,NumMPs(nprocs),MPendMax

!SMS$INSERT call GET_SMS_MPI_COMMUNICATOR(SMScomm)
!SMS$PARALLEL(dh, lp, mp) BEGIN
      sendCount=0
      do jth=1,ISTOT  
        NumMP = 0
        do mp=1,NMP  
          do lp=1,1
            NumMP = NumMP + 1
            do i=JMIN_IN(1),JMAX_IS(1)
              sendCount = sendCount + 1
            enddo
          enddo
        enddo
      enddo
!SMS$INSERT call MPI_ALLREDUCE(NumMP,MPendMax,1,MPI_INTEGER,MPI_MAX,SMScomm,status)
      if(status /= 0) then
!SMS$ignore begin
        print*,'MPI_ALLREDUCE error in module_init_plasma_grid',mype,status
!SMS$ignore end
        stop
      endif
!SMS$INSERT call MPI_GATHER(sendCount,1,MPI_INTEGER,recvCounts,1,MPI_INTEGER,0,SMScomm,status)
      if(status /= 0) then
!SMS$ignore begin
        print*,'MPI_GATHER sendCount error in module_init_plasma_grid',mype,status
!SMS$ignore end
        stop
      endif
!SMS$INSERT call MPI_GATHER(NumMP,1,MPI_INTEGER,NumMPs,1,MPI_INTEGER,0,SMScomm,status)
      if(status /= 0) then
!SMS$ignore begin
        print*,'MPI_GATHER NumMP error in module_init_plasma_grid',mype,status
!SMS$ignore end
        stop
      endif

!SMS$SERIAL(default=ignore) BEGIN
      DISPLS(1) = 0
      do PE = 2,nprocs
        DISPLS  (PE) = DISPLS  (PE-1) + recvCounts(PE-1)
      enddo
      NumPolevalProcs = 0
      do PE = 1,nprocs
        if(recvCounts(PE) > 0) then
          NumPolevalProcs = NumPolevalProcs + 1
          MPends(NumPolevalProcs) = NumMPs(PE) 
        endif
      enddo
!SMS$SERIAL END
!SMS$PARALLEL END

!SMS$ignore begin
      allocate(plasma_mp   (MaxFluxTube,MPendMax,ISTOT                )) !Must be in ignore or else it has halos
      if(mype == 0) then
        allocate(plasma_mpG(MaxFluxTube,MPendMax,ISTOT,NumPolevalProcs)) !Must be in ignore or else it has halos
      else
        allocate(plasma_mpG(1,1,1,1))
      endif
!SMS$ignore end

  END SUBROUTINE prepPoleVal
END MODULE module_prepPoleVal

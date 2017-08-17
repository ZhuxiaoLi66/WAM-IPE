subroutine stop
USE module_input_parameters,ONLY: MaxLpHaloUsed,MaxMpHaloUsed,mype,lps,lpe,mps,mpe,parallelBuild
implicit none
include "gptl.inc"
integer       :: MAXlpHalo! Max (over all PEs) lp halo size used
integer       :: MAXmpHalo! Max (over all PEs) mp halo size used
integer       :: comm
character(80) :: string
integer       :: ret

if(parallelBuild) then
  print*
  MAXlpHalo = MaxLpHaloUsed
  MAXmpHalo = MaxMpHaloUsed
!SMS$REDUCE(MAXlpHalo,max)
!SMS$REDUCE(MAXmpHalo,max)
  print"(' Maximum of MAXlpHalo over all processors: ',i0)",MAXlpHalo
  print"(' Maximum of MAXmpHalo over all processors: ',i0)",MAXmpHalo
  print*
  print*,'   PE   lps   lpe   mps   mpe   MAXlpHalo   MAXmpHalo'
!SMS$IGNORE BEGIN
  write(string,100) mype,lps,lpe,mps,mpe,MaxLpHaloUsed,MaxMpHaloUsed
100 format(5i6,i8,i12)
!SMS$IGNORE END
!SMS$INSERT call SMSPrintModeOrdered(string)
endif

! Print timing results to file named timing.mype
!ret = gptlpr (mype)

if(parallelBuild) then
!SMS$INSERT call GET_SMS_MPI_COMMUNICATOR(COMM)
  ret = gptlpr_summary(COMM)
else
  ret = gptlpr_summary()
endif

!ret = gptlprint_memusage ('Memory usage:')

print*,'IPE completed successfully'

end subroutine stop

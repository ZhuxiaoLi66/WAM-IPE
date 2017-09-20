subroutine stop
USE module_input_parameters,ONLY: MaxLpHaloUsed,MaxMpHaloUsed,mype,nprocs,lps,lpe,mps,mpe,SMScomm
implicit none
include "gptl.inc"
integer       :: MAXlpHalo! Max (over all PEs) lp halo size used
integer       :: MAXmpHalo! Max (over all PEs) mp halo size used
integer       :: comm
character(80) :: string
integer       :: ret
real(8)       :: TOTALTIME ! returned from gptl

if(nprocs > 1) then
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
!if (mype == 0 .or. mype == nprocs-1) then
  ret = gptlpr (mype)
!endif

!SMS$INSERT ret = gptlpr_summary(SMScomm)

!ret = gptlprint_memusage ('Memory usage:')
ret = gptlget_wallclock ('Total', 0, TOTALTIME)  ! The "0" is thread number
print*,'Total time =' , TOTALTIME
print*,'IPE completed successfully'

end subroutine stop

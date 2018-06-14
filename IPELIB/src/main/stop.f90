subroutine stop
USE module_input_parameters,ONLY: MaxLpHaloUsed,MaxMpHaloUsed,nprocs
implicit none
integer       :: MAXlpHalo! Max (over all PEs) lp halo size used
integer       :: MAXmpHalo! Max (over all PEs) mp halo size used

if(nprocs > 1) then
  MAXlpHalo = MaxLpHaloUsed
  MAXmpHalo = MaxMpHaloUsed
!SMS$REDUCE(MAXlpHalo,max)
!SMS$REDUCE(MAXmpHalo,max)
endif

print*,'IPE completed successfully'

end subroutine stop

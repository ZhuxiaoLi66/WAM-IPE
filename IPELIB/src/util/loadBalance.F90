Program loadBalance
!This routine does load balancing for the IPE code.
!The IPE code has a flux tube for each lp point, but the tubes are of differing lengths causing load imbalance.
!This code combines some of the lp points to create groups that have similar numbers of flux tube points.
!For load balancing the groups, not the lp points, will be apportioned among the processors.

implicit none
integer,parameter   :: ngroups=40       !Number of groups NLP will be reduced to.
integer,parameter   :: Ntries=1000      !Number of tries before giving up.
integer,parameter   :: imax  =100       !Maximum number of lp points in a group.
real   ,parameter   :: targetRatio=1.05 !Target for maximux ratio: (target # group points)/MaxFLuxTube
integer,allocatable :: JMAX_IS(:)       !Flux tube length for each lp
integer             :: MaxFluxTube      !Maximum flux tube length
integer             :: NLP              !Number of lp points
integer             :: maxPts           !target_loop vairable
integer             :: maxPtsStart      !target loop start
integer             :: maxPtsEnd        !target loop end
integer             :: sumTubes (imax)  !Number of flux tube points in a group
integer             :: sumTubess(imax)  !Number of flux tube points in a group
integer             :: layout (0:imax)  !Starting lp in a group
integer             :: layouts(0:imax)  !Starting lp in a group
integer             :: i                !Index for the  groups
integer             :: d                !Index for the direction
integer             :: lp,lps,lpe,lpd   !Indexes for lp
real                :: finalRatio       !The final ratio: (final # group points)/MaxFLuxTube
character(13)       :: direction(2)     !The direction for the search: forward or backward
character(LEN=72)   :: fname
integer             :: istat

direction(1) = 'forward pass '
direction(2) = 'backward pass'
read(88,*) NLP
allocate(JMAX_IS(NLP))
read(88,*) (i,JMAX_IS(lp),lp=1,NLP)
MaxFluxTube = maxval(JMAX_IS)
print"(' NLP = ',i0)",NLP
print"(' Average flux tube length  ',i0)",sum(JMAX_IS)/NLP
print"(' Maximum flux tube length ' ,i0)",MaxFluxTube
print*

maxPtsStart = targetRatio*MaxFluxTube
maxPtsEnd   = maxPtsStart+Ntries
target_loop: do maxPts = maxPtsStart,maxPtsEnd
  lps=1
  lpe=NLP
  lpd=1
  do d=1,2 !First forward pass, then backward pass
    i=1
    sumTubes=0
    layout  =0
    do lp=lps,lpe,lpd
      if(sumTubes(i)+JMAX_IS(lp) > maxPts) then
        i=i+1
        if(i>imax) then
          print*,'i > imax',i,imax
          stop
        endif
        layout  (i) = lp
        sumTubes(i) = sumTubes(i)+JMAX_IS(lp)
        cycle
      endif
      layout  (i) = lp
      sumTubes(i) = sumTubes(i)+JMAX_IS(lp)
    enddo !lp
    if(i < ngroups) then
      print*,'i < ngroups',i,ngroups
      stop
    elseif(i  == ngroups) then
      finalRatio = float(maxPts)/float(MaxFluxTube)
      print"(' Success on ',a13,'after ',i0,' iterations')",direction(d),maxPts-maxPtsStart+1
      print"(' Target ratio     ',f6.2)",targetRatio
      print"(' Target MAX points',i6  )",maxPtsStart
      print"(' Final ratio      ',f6.2)",finalRatio
      print"(' MIN points       ',i6  )",minval(sumTubes(1:i))
      print"(' MAX points       ',i6  )",maxval(sumTubes(1:i))
      exit target_loop
    elseif(maxPts == maxPtsEnd .and. d==2) then
      print*,'End of target_loop',i,ngroups
      print*,'Could not find solution'
      stop
    endif
    lps=NLP
    lpe=1
    lpd=-1
  enddo !d
enddo target_loop
if(d==2) then
  layouts    = layout
  layouts(0) = NLP+1
  sumTubess  = sumTubes
  do i=1,ngroups
    layout  (i) = layouts  (ngroups-i)-1
    sumTubes(i) = sumTubess(ngroups-i+1)
  enddo
endif
    
print*
print*,'       bin        range     points'
layout(0) = 0
do i=1,ngroups
  print"(2i10,'-',i3,i10)",i,layout(i-1)+1,layout(i),sumTubes(i)
enddo

fname = './load_balance_groups1'
OPEN (100,FILE=fname,STATUS='NEW',IOSTAT=istat)
write(100,*) layout(1:ngroups)

end Program loadBalance



Program loadBalance
!This routine does LP load balancing for the IPE code.
!This code reads in timings (times) from a representative run of the IPE code parallelized on NLP points in the LP direction only.
!The representative run should br made with "barriersOn=f' in IPE.inp.
!Currently flux_tube_solver times are used which are grepped from timing.* by:
!grep '   flux_tube_solver   ' timing.* | cut -c 1-7 --complement | tr -d ':' | awk '{$2=$4=""}1' | sort -h > times
!This code combines some of the LP points to create groups (bins) that have similar flux_tube_solver times.
!To create the bins this code uses the the max time over all processors ( maxval(times) ) as the target time (targetTime).
!Ideally all the bins would have a total time of targetTime for a minimum number of bins and perfect load balancing.
!In practice the times don't add up that way and it's generally good to allow some bins to have a total time (maxTimeAllowed) of slightly more than targetTime.
!The code has one input called "tolerance" which controls the size of maxtimeAllowed by : maxTimeAllowed = tolerance*targetTime
!tolerance must be >= 1.0. A good value for tolerance is 1.02. 
!If you want to decrease the number of bins, say from 61 to 60 increase tolerance say from 1.02 to 1.03.
!Conversely, if you want to increase the number of bins, say from 59 to 60 decrease tolerance.
!For parallel distribution the groups, not the lp points, will be apportioned among the processors.

implicit none
integer,parameter :: NLP=170            !Number of lp points
integer,parameter :: maxBins=NLP        !Maximum number of lp points in a group.
integer           :: layout(0:maxBins,2)!Starting lp in a group
integer           :: layouts(0:maxBins) !Starting lp in a group
integer           :: i                  !Index for the  groups
integer           :: nt                 !Index for the time in a group
integer           :: d                  !Index for the direction
integer           :: lp,lps,lpe,lpd     !Indexes for lp
integer           :: p                  !Processor number starting from 0
integer           :: Ncalls             !Number of calls for a processor
integer           :: istat              !Status for open
integer           :: Ngroups(2),ng      !Number of groups
integer           :: status             !Status of the read statement.
real              :: minTime            !Minimum input time
real              :: meanTime           !Mean input time
real              :: maxTime            !Maximum input time
real              :: totalTime          !Sum of all the input times
real              :: targetTime         !The time we are shooting for
real              :: tolerance=1.02     !See explanation above.
real              :: maxdiff(2)         !Maximum difference between low and high for each direction
real              :: sumTime (maxBins,2)!Running sum of the times in a group
real              :: sumTimes(maxBins)  !Temporary storage for sumTimes
real              :: times(NLP)         !The run times which are read in from the 'ties' file
real              :: testTime           !The sum of times in a group for testing
real              :: maxTimeAllowed     !The maximum time allowed in a group.
real              :: dum1,dum2,dum3,dum4!Dummy values for reading file 'times'
character(LEN=13) :: direction(2)       !The direction for the load balancing: forward or backward
character(LEN=72) :: fname              !The name of the output file

direction(1) = 'forward pass '
direction(2) = 'backward pass'
open(88,file='times')
do lp=1,NLP
  read(88,*) p,Ncalls,times(lp),dum1,dum2,dum3,dum4
  if(p /= lp-1) then
    print*,'Error reading pe2s'
    print*,lp,p,Ncalls,times(lp),dum1,dum2,dum3,dum4
    stop
  endif
enddo
Print"('Enter tolerance >= 1.0 or press Ctrl-D for the default value 1.02: ',$)"
read(*,*,iostat=status) tolerance
if(tolerance < 1.0) then
  print*,'tolerance must be >= 1.0, you entered ',tolerance
  stop
endif
maxTime        = maxval(times)
minTime        = minval(times)
meanTime       = sum(times)/NLP
totalTime      = sum(times)
targetTime     = maxtime
maxTimeAllowed = tolerance*targetTime
maxdiff        = 0.0
lps            = 1
lpe            = NLP
lpd            = 1
do d=1,2 !First forward pass, then backward pass
  i=1
  sumTime(:,d) = 0.0
  layout (:,d) = 0
  LP_LOOP:  do lp=lps,lpe,lpd
    do nt=1,maxBins
      testTime = sumTime(i,d)+times(lp)
      if( testTime>= targetTime) then !Above or equal to the target
        if (testTime-targetTime <= targetTime-sumTime(i,d).and.testTime <= maxTimeAllowed) then !Less above target than below so use this one
          if(i>maxBins) then
            print*,'i > maxBins',i,maxBins
            stop
          endif
          sumTime(i,d) = testTime
          maxDiff(d) = max(maxDiff(d),abs(sumTime(i,d)-targetTime))
          layout (i,d) = lp
          i=i+1
          cycle LP_LOOP
        else !more above target than below so use the last one
          layout (i,d) = lp -1 + 2*(d/2)
          maxDiff(  d) = max(maxDiff(d),abs(sumTime(i,d)-targetTime))
          i=i+1
          cycle
        endif
      else ! Below target so continue
        if(lp == lpe) then !At the end of the times
          if(sumTime(i-1,d)+testTime-targetTime <= maxDiff(d)) then !Add testTime to the previous bin
            i = i-1
          endif
          sumTime(i,d) = sumTime(i,d)+testTime
          maxDiff(  d) = max(maxDiff(d),abs(sumTime(i,d)-targetTime))
          layout (i,d) = lp
          exit LP_LOOP
        endif
      endif
      sumTime(i,d) = sumTime(i,d)+times(lp)
      cycle LP_LOOP
    enddo ! nt loop
    print*,'Exited nt loop without filling all the bins'
    print*,d,lp,i
    stop
  enddo LP_LOOP
  lps=NLP
  lpe=1
  lpd=-1
  Ngroups(d) = i
  !print*,direction(d)
  !print*,layout(1:Ngroups(d),d)
  !print*,'Points =',Ngroups(d)
  !print*,'maxDiff=', maxDiff(d)
  !print*,maxval(abs(sumTime(1:Ngroups(d),d)-targetTime))
  !print*,'Min    =', minval(sumTime(1:Ngroups(d),d))
  !print*,'Max    =', maxval(sumTime(1:Ngroups(d),d))
  !print*,'Spread =', maxval(sumTime(1:Ngroups(d),d)) - minval(sumTime(1:Ngroups(d),d))
  !print*,sumTime(1:Ngroups(d),d)
enddo !d

if(maxDiff(1) < maxDIff(2)) then
  if(maxval(sumTime(1:Ngroups(1),1)) <= maxTimeAllowed) then
    d = 1
  else
    d = 2
  endif
else
  d = 2
endif

ng=Ngroups(d)
if(d==2) then
  layouts    = layout(:,2)
  layouts(0) = NLP+1
  sumTimes  = sumTime(:,d)
  do i=1,ng
    layout (i,d) = layouts (ng-i)-1
    sumTime(i,d) = sumTimes(ng-i+1)
  enddo
else
  layout(0,d) = 0
endif

print*
print"(' NLP             ',i10  )",NLP
print"(' Min    time     ',f10.2)",minTime
print"(' Mean   time     ',f10.2)",meanTime
print"(' Max    time     ',f10.2)",maxTime
print"(' Total  time     ',f10.2)",totalTime
print"(' Target time     ',f10.2)",targetTime
print"(' tolerance       ',f10.5)",tolerance
print"(' Max time allowed',f10.2)",maxTimeAllowed
print*
print*,direction(d)
print"(' Number of Groups',i10  )", ng
print"(' Max Difference  ',f10.2)", maxDiff(d)
print"(' Min             ',f10.2)", minval(sumTime(1:ng,d))
print"(' Max             ',f10.2)", maxval(sumTime(1:ng,d))
print"(' Spread          ',f10.2)", maxval(sumTime(1:ng,d)) - minval(sumTime(1:ng,d))
print*
print*,'       bin     LP range     time       spread'
do i=1,ng
  print"(i10,i9'-',i3,2f12.2)",i-1,layout(i-1,d)+1,layout(i,d),sumTime(i,d),sumTime(i,d)-targetTime
enddo

fname = './load_balance_groups1'
OPEN (100,FILE=fname,STATUS='REPLACE',IOSTAT=istat)
write(100,*) layout(1:ng,d)

end Program loadBalance


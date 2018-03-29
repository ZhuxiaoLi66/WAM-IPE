program loadBalance
implicit none
integer :: NLP,n,i,k,j,goal,sum,nlb,MAXabove,MINbelow,bias
integer,allocatable :: Ntubes(:),balance(:),Lbalance(:)

read(88,*) NLP
print*,'NLP      =',NLP
allocate(Ntubes(NLP),balance(NLP),Lbalance(NLP))
read(88,*) (n,Ntubes(i),i=1,NLP)
goal     = Ntubes(1)
bias     = .12*goal
i        = NLP
k        = 0
sum      = 0
MINbelow = goal
MAXabove = goal
do while (i>1)
  k=k+1
  balance(k) = i
  do j=i,1,-1
    sum = sum + Ntubes(j)
    if(sum+Ntubes(j-1)> goal) then
      if(goal-sum < sum+Ntubes(j-1)-goal .or. sum+Ntubes(j-1)-goal > bias) then
        MINbelow = min(MINbelow,sum)
        i = j-1
        sum = 0
        exit
      else
        MAXabove = max(MAXabove,sum+Ntubes(j-1))
        i = j-2
        sum = 0
        exit
      endif
    endif
  enddo
enddo
nlb = k+1
print*,'nlb      =',nlb
print*,'Target   =',goal
print*,'bias     =',bias
print*,'MINbelow =',goal-MINbelow
print*,'MAXabove =',MAXabove-goal
balance(nlb) = 1
do i=1,nlb
  Lbalance(i) = balance(nlb-i+1)
enddo
print*,(Lbalance(i),i=1,nlb)
end program loadBalance

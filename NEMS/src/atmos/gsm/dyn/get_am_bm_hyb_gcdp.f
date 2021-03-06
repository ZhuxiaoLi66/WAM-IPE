      subroutine am_bm_hyb_gcdp
!
! hmhj : this is modified hybrid by finite difference from henry juang
!        work for temperature and enthalpy with thickness in pressure
!        as the prognostic variables
! program log
! 2011 02 20 : henry juang, coded to have mass_dp options
! 2013 09 30 : henry juang, remove linear divergence for vertical motion
!                           in thermodynamic linear term for ndslfv
!
      USE gfs_dyn_MACHINE , ONLY : kind_grid
      use gfs_dyn_resol_def
      use namelist_dynamics_def , only : ref_temp,ndslfv, 
     &                                   semi_implicit_temp_profile
      use gfs_dyn_coordinate_def
      use gfs_dyn_physcons, rd => con_rd
     &                    , cp => con_cp, rearth => con_rerth

      IMPLICIT NONE 

      REAL(KIND=KIND_EVOD)                                              
     &     PK5REF(LEVP1),DPKREF(LEVS),RPKREF(LEVS),rp2ref(levs),
     &     DBKREF(LEVS), BBKREF(LEVS),rdlref(levs),
     &     HV0(LEVS),HV(LEVS),HVK(LEVS),
     &     C0KREF(LEVP1),dprp,hrp,tcprp2,tcpcprp2,
     &     TREF,PSREF,RKAA,KAPPA,RKAPPA
      real alpha(levs),betta(levs),gamma(levs),delta(levs)
      real zm(levs-1,levs-1),tm(levs-1,levs)
      real pm(levs,levs-1),wtm(levs-1,levs)
      real dup(levs),dum(levs),hmm(levs,levs)
      real det,wmkm1,wmkp1
      integer lll(levs-1),mmm(levs-1)
      INTEGER K,I,J,N
!                                                                       
! idea change
      real test1(150)
      data test1/300.,300.,300., 300.,300.,300.,
     &300.,  300.,  300.,  300.,  300.,  300.,
     &300.,  300.,  300.,  300.,  300.,  300.,
     &300.,  300.,  300.,  300.,  300.,  300.,
     &300.,  300.,  300.,  300.,  300.,  300.,
     &300.,  300.,  300.,  300.,  300.,  300.,
     &300.,  300.,  300.,  300.,  300.,  300.,
     &300.,  300.,  300.,  300.,  300.,  300.,
     &300.,  300.,  300.,  300.,  300.,  300.,
     &300.,  300.,  300.,  300.,  300.,  300.,
     &300.,  300.,  300.,  300.,  300.,  300.,
     &300.,  300.,  300.,  300.,  300.,  300.,
     &300.,  300.,  300.,  300.,  300.,  300.,
     &300.,  300.,  300.,  300.,  300.,  300.,
     &300.,  300.,  300.,  300.,  300.,  300.,
     &300.,  300.,  300.,  300.,  300.,  300.,
     &300.,  300.,  300.,  300.,  300.,  300.,
     &300.,  300.,  300.,  300.,  303.,  314.,
     &337.,  376.,  434.,  495.,  590.,  698.,
     &815.,  938., 1066., 1196., 1326., 1455.,
     &1579., 1699., 1812., 1919., 2017., 2106.,
     &2185., 2256., 2316., 2367., 2408., 2440.,
     &2457., 2475., 2487., 2494., 2498., 2499.,
     &2500., 2500., 2500., 2500., 2500., 2500.,
     &2500., 2500., 2500., 2500., 2500., 2500./

!     print *,' enter  get_am_bm_hyb_gc_fd'

!     PSREF=80.
      PSREF=101.316                                                         
      KAPPA=RD/CP                                                       
      RKAPPA=1./KAPPA

      if( thermodyn_id.eq.3 ) then
        TREF=ref_temp*cp
        RKAA=KAPPA/(REARTH*REARTH)
      else
        TREF=ref_temp                                                         
        RKAA=RD/(REARTH*REARTH)
      endif
!
! idea change
      if( semi_implicit_temp_profile ) then
        print *,' use layer mean temperature for semi-implicit '
      else
        print *,' use constant temperature for semi-implicit '
        if( thermodyn_id.eq.3 ) then
          do k=1,levs
            thref (K)=test1(K)*cp
          enddo
        else
          do k=1,levs
            thref (K)=test1(K)
          enddo
        endif
      endif
                                                                        
      DO K=1,LEVP1                                                      
       PK5REF(K)=AK5(K)+BK5(K)*PSREF+CK5(K)                                  
      ENDDO                                                             
                                                                        
      DO K=1,LEVS                                                       
!      THREF (K)=TREF
       DPKREF(K)=PK5REF(K)-PK5REF(K+1)                                  
       rdlref(k)=0.5/dpkref(k)
       RPKREF(K)=1.0/(PK5REF(K)+PK5REF(K+1)) 
       RP2REF(K)=rpkref(k)*rpkref(k)
       DBKREF(K)=BK5(K)-BK5(K+1)                                  
       BBKREF(K)=BK5(K)+BK5(K+1)                                  
      ENDDO                                                             
 
      c0kref=0
      do k=2,levs
       c0kref(k)=ck5(k)*rkappa/(thref(k-1)+thref(k))
      enddo
                                                                        
! -----------------------------------------------------
c hm
      HMHYB = 0.0
c hm1
      do i=1,levs
        hrp=thref(i)*rpkref(i)
        do j=i,levs
          hmhyb(i,j)=2.*hrp
        enddo
      enddo
      do j=1,levs-1
        hrp=thref(j)*rpkref(j)
        do i=j+1,levs
          hmhyb(i,j)=hmhyb(i,j)+2.*hrp
        enddo
      enddo
c hm2
      hmm = 0.0
      do i=1,levs
        hrp = thref(i)*dpkref(i)*rp2ref(i)
        hmm(i,i) = hrp
        do j=i+1,levs
          hmm(i,j) = 2*hrp
        enddo
      enddo
      do j=1,levs
        do i=1,levs
          hmhyb(i,j)=hmhyb(i,j)-hmm(i,j)
        enddo
      enddo
c hm3
      hmm = 0.0
      dum = 0.0
      do i=2,levs
        hrp = 2.*thref(i-1)*dpkref(i-1)*rp2ref(i-1)
        dum(i-1) = dum(i-1) + hrp
        do j=i,levs
          dum(j) = dum(j) + 2.*hrp
        enddo
        do j=1,levs
          hmm(i,j) = dum(j)
        enddo
      enddo
      do j=1,levs
        do i=1,levs
          hmhyb(i,j)=hmhyb(i,j)-hmm(i,j)
        enddo
      enddo
! apply RKAA to the sum
      do j=1,levs
        do i=1,levs
          hmhyb(i,j)=RKAA*hmhyb(i,j)
        enddo
      enddo
! -----------------------------------------------------
c am
      AMHYB = 0.0
      do j=1,levs
        dprp=dpkref(j)*rpkref(j)
        amhyb(j,j)=amhyb(j,j)+dprp
        do i=j+1,levs
          amhyb(i,j)=amhyb(i,j)+2.*dprp
        enddo
      enddo
! apply RKAA to the sum
      do j=1,levs
        do i=1,levs
          AMHYB(i,j)=RKAA*amhyb(i,j)
        enddo
      enddo
! -----------------------------------------------------
c bm
      BMHYB = 0.0
      do i=1,levs
        BMHYB(i,i)=kappa*thref(i)*rpkref(i)*dpkref(i)
      enddo
      do j=2,levs
        do i=1,j-1
          BMHYB(i,j)=2.*kappa*thref(i)*rpkref(i)*dpkref(j)
        enddo
      enddo

c need zm, tm and pm for bm+
! alpha, betta, gamma
      alpha(levs)=0.0
      betta(   1)=0.0
      do k=2,levs
        alpha(k-1)=(pk5ref(k)+pk5ref(k+1))/(pk5ref(k-1)+pk5ref(k))
        alpha(k-1)=alpha(k-1)**kappa
      enddo
      do k=1,levs
        gamma(k)=1.0 - kappa*DPKREF(k)*RPKREF(k)*2.0
        delta(k)=1.0 + kappa*DPKREF(k)*RPKREF(k)*2.0
      enddo
      do k=1,levs-1
        betta(k+1)=(pk5ref(k)+pk5ref(k+1))/(pk5ref(k+1)+pk5ref(k+2))
        betta(k+1)=betta(k+1)**kappa
      enddo
! zm
! zm [levs-1,levs-1]
      dup(levs)=0.0
      dum(1 )=0.0
      do k=1,levs-1
        dup(k  )=delta(k)*thref(k)-betta(k+1)*thref(k+1)
        dum(k+1)=alpha(k)*thref(k)-gamma(k+1)*thref(k+1)
      enddo
 
      zm=0.0		! (levs-1,levs-1)
      k=2
        wmkm1=c0kref(k)*rdlref(k-1)
        wmkp1=c0kref(k)*rdlref(  k)
        zm(k-1,k-1)=wmkm1*dup(k-1)+wmkp1*dum(k)-1.0
        zm(k-1,k  )=wmkp1*dup(k)
      do k=3,levs-1
        wmkm1=c0kref(k)*rdlref(k-1)
        wmkp1=c0kref(k)*rdlref(  k)
        zm(k-1,k-2)=wmkm1*dum(k-1)
        zm(k-1,k-1)=wmkm1*dup(k-1)+wmkp1*dum(k)-1.0
        zm(k-1,k  )=wmkp1*dup(k)
      enddo
      k=levs
        wmkm1=c0kref(k)*rdlref(k-1)
        wmkp1=c0kref(k)*rdlref(  k)
        zm(k-1,k-2)=wmkm1*dum(k-1)
        zm(k-1,k-1)=wmkm1*dup(k-1)+wmkp1*dum(k)-1.0
      call iminv(zm,levs-1,det,lll,mmm)
 
! tm
! [levs-1,levs]
      tm=0.0
      do k=2,levs
        tm(k-1,k-1)=-c0kref(k)*kappa*thref(k-1)*dpkref(k-1)*rpkref(k-1)
        tm(k-1,k  )= c0kref(k)*kappa*thref(k  )*dpkref(k  )*rpkref(k  )
      enddo
      do k=2,levs
        do n=1,levs
          tm(k-1,n)=tm(k-1,n)-bk5(k)*dpkref(n)
        enddo
      enddo
      do k=2,levs
        do n=k,levs
          tm(k-1,n)=tm(k-1,n)+(1.-2.*c0kref(k)*kappa*
     &          (thref(k-1)*rpkref(k-1)+thref(k)*rpkref(k)))*dpkref(n)
        enddo
      enddo
! wtm
! = zm * tm
! [levs-1,levs]
      wtm=0.0
      do i=1,levs
        do n=1,levs-1
          do j=1,levs-1
            wtm(j,i)=wtm(j,i)+zm(j,n)*tm(n,i)
          enddo
        enddo
      enddo
! pm
! [levs, levs-1]
      pm=0.0
      k=1
        pm(k,k  )=(delta(k)*thref(k)-betta(k+1)*thref(k+1))*
     &            rdlref(k)
      do k=2,levs-1
        pm(k,k-1)=(alpha(k-1)*thref(k-1)-gamma(k)*thref(k))*
     &            rdlref(k)
        pm(k,k  )=(delta(k)*thref(k)-betta(k+1)*thref(k+1))*
     &            rdlref(k)
      enddo
      k=levs
        pm(k,k-1)=(alpha(k-1)*thref(k-1)-gamma(k)*thref(k))*
     &            rdlref(k)
!
! bm+ = pm * wtm
!
      do i=1,levs
        do k=1,levs
          do j=1,levs-1
            bmhyb(k,i)=bmhyb(k,i)+pm(k,j)*wtm(j,i)
          enddo
        enddo
      enddo
 
!   
! sm = sv + [0,wtm] - [wtm 0]
!
      smhyb = 0.0
      do k=1,levs                                                       
       smhyb(k,k)=dpkref(k)
      enddo                                                             
      do j=1,levs
          do i=2,levs
            smhyb(i,j) = smhyb(i,j) + wtm(i-1,j)
          enddo
          do i=1,levs-1
            smhyb(i,j) = smhyb(i,j) - wtm(i,j)
          enddo
      enddo

!
!     print *,' end of get_am_bm_hyb_gc_fd. '
!!
      RETURN                                                            
      END                                                               

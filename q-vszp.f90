
!! ------------------------------------------------------------------------
!  setup AO basis with q/CN  AO dep. i.e. q-vSZP for individual atoms
!  part of GP3 pilot code
!! ------------------------------------------------------------------------

subroutine setupbas(prt)
      use bascom
      implicit none
      logical, intent(in)  :: prt

! local variables
      integer i,j,k,l,iat,ish,iish,iao,nao,lshell,pr,lladr(0:3),lladr2(0:3)
      data lladr  /1,3,6,10/
      data lladr2 /1,3,5, 7/
      real*8 DEX(-1:96), XNORM,DEX2,PI,EXPO,IAM,alp,tmp,x1,f2
      real*8,allocatable :: cn(:)

      if (prt) then
      write(*,*) ' AO modification setup (EEQ, (q+aq^2+bCN^0.5) :'
      write(*,*) ' at        q          CN       q change    CN change     total '
      endif

      allocate(cn(n))

      call ncoord_erf(n,at,rab,-3.75d0,cn) ! slightly better than original kn=-7.5 is optimum for EEQ->basis

      call eeq(n,at,rab,chrg,cn,q_eeq)     ! version with some EEQ parameters opt. for H-Kr total energies

      call ncoord_basq(n,at,rab,cn)        ! mod CN for AO compression/field term

      PI=4.D0*ATAN(1.D0)
      DO I=-1,10
         DEX(I)=DEX2(I)
      ENDDO   

      k  = 0
      j  = 0
      iish=0
      npr= 0
      ncao=0
      nao =0
      do i=1,n   
         iat = at(i)
         if (prt)  &
   &     write(*,'(i4,7F12.5)') iat,q_eeq(i),cn(i),& 
   &                            q_eeq(i)-bas_par(2,iat)*q_eeq(i)**2, bas_par(1,iat) * sqrt(cn(i)), &
   &                            q_eeq(i)-bas_par(2,iat)*q_eeq(i)**2+ bas_par(1,iat) * sqrt(cn(i))
         do ish=1,bas_nsh(iat)
            iish = iish + 1
            shmap(ish,i)   =iish  ! needed for shell ES
            caoshell(ish,i)=k
             aoshell(ish,i)=j
            lshell=bas_lsh(ish,iat)
            k=k+lladr (lshell)          
            j=j+lladr2(lshell)               
            do iao=1,lladr2(lshell)              
               nao = nao + 1
               aoat(nao) = i
            enddo
            x1 = q_eeq(i)-bas_par(2,iat)*q_eeq(i)**2 &
   &                    + bas_par(1,iat) * sqrt(cn(i)) ! effective charge = charge + CN dep. (element parameters)
            do iao=1,lladr(lshell)              
               ncao = ncao + 1
               prim_npr(ncao)=bas_npr(ish,iat)
               alp = 0
               tmp = 0
               do pr=1,bas_npr(ish,iat)
                  npr= npr+ 1
                  expo = bas_ec(1,pr,ish,iat)
                  iam = lshell
                  XNORM=(2.D0*EXPO/PI)**0.75D0*(4.D0*EXPO)**(IAM/2.D0)/SQRT(DEX(2*IAM-1))
                  prim_exp(npr) =  bas_ec(1,pr,ish,iat)                                       ! primitive data
!                 prim_cnt(npr) =  bas_ec(2,pr,ish,iat) * xnorm
                  prim_cnt(npr) = (bas_ec(2,pr,ish,iat) + bas_ec(3,pr,ish,iat) * x1) * xnorm  ! cont. dep.
                  if(lshell.eq.2.and.iao.gt.3) prim_cnt(npr)=prim_cnt(npr)*sqrt(3.0d0)
               enddo
            enddo
         enddo
      enddo

      if(iish.ne.nsh) stop 'weird iish'

      nsao=0
      do i=1,n   
         iat = at(i)
         do ish=1,bas_nsh(iat)
            lshell=bas_lsh(ish,iat)
            do iao=1,lladr2(lshell)              
               nsao = nsao + 1
               shell2ao(nsao) = ish
            enddo
         enddo
      enddo

      k=0  
      do i=1,ncao
         prim_count(i)=k
         k=k+prim_npr(i)
      enddo   

      if (prt) then
      write(*,*) 'ncao ',ncao
      write(*,*) 'nsao ',nsao
      write(*,*) 'npr  ',npr 
      endif

end

! ------------------------------------------------------------------------
! slightly modified EEQ, SG 2/23
! ------------------------------------------------------------------------

subroutine eeq(n,at,rab,chrg,cn,q)
      use iso_fortran_env, id => output_unit, wp => real64
      implicit none
      integer, intent(in)  :: n                ! number of atoms     
      integer, intent(in)  :: at(n)            ! ordinal numbers            
      real(wp),intent(in)  :: rab(n*(n+1)/2)   ! dist
      real(wp),intent(in)  :: chrg             ! total charge on system
      real(wp),intent(in)  :: cn(n)            ! CN                       
      real(wp),intent(out) :: q(n)             ! output charges

!  local variables
      integer  :: m,i,j,k,ii,ij,nfrag
      integer  :: info,lwork
      integer,allocatable :: ipiv(:)
      real(wp) :: gammij,tsqrt2pi,r2,tmp,rij
      real(wp) :: chieeq(86), gameeq(86), cnfeeq(86), alpeeq(86)
      real(wp),allocatable :: A (:,:),x (:),work (:), alp(:), gam(:)
!  parameter
      parameter (tsqrt2pi = 0.797884560802866_wp)
! ------------------------------------------------------------------------
!  PARAMETRISATION BY S. SPICHER (NEW)
! ------------------------------------------------------------------------
! ------------------------------------------------------------------------
   parameter( chieeq =(/& ! mod. by SG for wB97X-D3/q-vSZP energy see below
    1.27000000_wp, 1.25000000_wp, 0.44000000_wp, 0.99000000_wp, 1.30000000_wp, &
    1.48000000_wp, 1.63000000_wp, 1.66000000_wp, 1.66000000_wp, 1.05000000_wp, &
    0.44000000_wp, 0.68000000_wp, 1.16000000_wp, 1.15000000_wp, 1.67000000_wp, &
    1.52000000_wp, 1.60000000_wp, 1.20000000_wp, 0.35000000_wp, 0.70000000_wp, & ! Ca
    0.76482096_wp, 0.98457281_wp, 0.96702598_wp, 1.05266584_wp, 0.93274875_wp, &
    1.04025281_wp, 0.92738624_wp, 1.07419210_wp, 1.07900668_wp, 1.04712861_wp, & ! Zn
    1.11000000_wp, 1.14000000_wp, 1.45000000_wp, 1.49000000_wp, 1.48000000_wp, &
    1.11000000_wp, 0.36273870_wp, 0.58797255_wp, 0.71961946_wp, 0.96158233_wp, &
    0.89585296_wp, 0.81360499_wp, 1.00794665_wp, 0.92613682_wp, 1.09152285_wp, &
    1.14907070_wp, 1.13508911_wp, 1.08853785_wp, 1.11005982_wp, 1.12452195_wp, &
    1.21642129_wp, 1.36000000_wp, 1.42000000_wp, 1.33000000_wp, 0.34125098_wp, & ! Cs
    0.58884173_wp, 0.68441115_wp, 0.56999999_wp, 0.56999999_wp, 0.56999999_wp, &
    0.56999999_wp, 0.56999999_wp, 0.56999999_wp, 0.56999999_wp, 0.56999999_wp, &
    0.56999999_wp, 0.56999999_wp, 0.56999999_wp, 0.56999999_wp, 0.56999999_wp, &
    0.56999999_wp, 0.87936784_wp, 1.02761808_wp, 0.93297476_wp, 1.10172128_wp, &
    0.97350071_wp, 1.16695666_wp, 1.23997927_wp, 1.18464453_wp, 1.14191734_wp, &
    1.12334192_wp, 1.01485321_wp, 1.12950808_wp, 1.30804834_wp, 1.33689961_wp, &
    1.37000000_wp /))
   parameter( gameeq =(/&
   -0.35015861_wp, 1.04121227_wp, 0.09281243_wp, 0.09412380_wp, 0.26629137_wp, &  
    0.19408787_wp, 0.05317918_wp, 0.03151644_wp, 0.32275132_wp, 1.30996037_wp, &
    0.24206510_wp, 0.04147733_wp, 0.11634126_wp, 0.13155266_wp, 0.15350650_wp, &
    0.15250997_wp, 0.17523529_wp, 0.28774450_wp, 0.42937314_wp, 0.01896455_wp, &
    0.07179178_wp,-0.01121381_wp,-0.03093370_wp, 0.02716319_wp,-0.01843812_wp, &
   -0.15270393_wp,-0.09192645_wp,-0.13418723_wp,-0.09861139_wp, 0.18338109_wp, & ! 30
    0.08299615_wp, 0.11370033_wp, 0.19005278_wp, 0.10980677_wp, 0.12327841_wp, &
    0.25345554_wp, 0.58615231_wp, 0.16093861_wp, 0.04548530_wp,-0.02478645_wp, &
    0.01909943_wp, 0.01402541_wp,-0.03595279_wp, 0.01137752_wp,-0.03697213_wp, & ! 45
    0.08009416_wp, 0.02274892_wp, 0.12801822_wp,-0.02078702_wp, 0.05284319_wp, &
    0.07581190_wp, 0.09663758_wp, 0.09547417_wp, 0.07803344_wp, 0.64913257_wp, & ! 55
    0.15348654_wp, 0.05054344_wp, 0.11000000_wp, 0.11000000_wp, 0.11000000_wp, &
    0.11000000_wp, 0.11000000_wp, 0.11000000_wp, 0.11000000_wp, 0.11000000_wp, &
    0.11000000_wp, 0.11000000_wp, 0.11000000_wp, 0.11000000_wp, 0.11000000_wp, &
    0.11000000_wp,-0.02786741_wp, 0.01057858_wp,-0.03892226_wp,-0.04574364_wp, &
   -0.03874080_wp,-0.03782372_wp,-0.07046855_wp, 0.09546597_wp, 0.21953269_wp, &
    0.02522348_wp, 0.15263050_wp, 0.08042611_wp, 0.01878626_wp, 0.08715453_wp, &
    0.10500484_wp /))
   parameter( cnfeeq =(/& ! changed handish for RGs 
    0.04916110_wp, 0.25000000_wp,-0.12349591_wp,-0.02665108_wp,-0.02631658_wp, &
    0.06005196_wp, 0.09279548_wp, 0.11689703_wp, 0.15704746_wp, 0.21000000_wp, &
   -0.10002962_wp,-0.07712863_wp,-0.02170561_wp,-0.04964052_wp, 0.14250599_wp, &
    0.07126660_wp, 0.13682750_wp, 0.19000000_wp,-0.10219289_wp,-0.08979338_wp, &
   -0.08273597_wp,-0.01754829_wp,-0.02765460_wp,-0.02558926_wp,-0.08010286_wp, &
   -0.04163215_wp,-0.09369631_wp,-0.03774117_wp,-0.05759708_wp, 0.02431998_wp, &
   -0.01056270_wp,-0.02692862_wp, 0.07657769_wp, 0.06561608_wp, 0.08006749_wp, &
    0.14139200_wp,-0.05351029_wp,-0.06701705_wp,-0.07377246_wp,-0.02927768_wp, &
   -0.03867291_wp,-0.06929825_wp,-0.04485293_wp,-0.04800824_wp,-0.01484022_wp, &
    0.07917502_wp, 0.06619243_wp, 0.02434095_wp,-0.01505548_wp,-0.03030768_wp, &
    0.01418235_wp, 0.08953411_wp, 0.08967527_wp, 0.07277771_wp,-0.02129476_wp, &
   -0.06188828_wp,-0.06568203_wp,-0.11000000_wp,-0.11000000_wp,-0.11000000_wp, &
   -0.11000000_wp,-0.11000000_wp,-0.11000000_wp,-0.11000000_wp,-0.11000000_wp, &
   -0.11000000_wp,-0.11000000_wp,-0.11000000_wp,-0.11000000_wp,-0.11000000_wp, &
   -0.11000000_wp,-0.03585873_wp,-0.03132400_wp,-0.05902379_wp,-0.02827592_wp, &
   -0.07606260_wp,-0.02123839_wp, 0.03814822_wp, 0.02146834_wp, 0.01580538_wp, &
   -0.00894298_wp,-0.05864876_wp,-0.01817842_wp, 0.07721851_wp, 0.07936083_wp, &
    0.05849285_wp /))
   parameter( alpeeq =(/&
    0.55159092_wp, 0.66205886_wp, 0.90529132_wp, 1.51710827_wp, 2.86070364_wp, &
    1.88862966_wp, 1.32250290_wp, 1.23166285_wp, 1.77503721_wp, 1.11955204_wp, &
    1.28263182_wp, 1.22344336_wp, 1.70936266_wp, 1.54075036_wp, 1.38200579_wp, &
    2.18849322_wp, 1.36779065_wp, 1.27039703_wp, 1.64466502_wp, 1.58859404_wp, &
    1.65357953_wp, 1.50021521_wp, 1.30104175_wp, 1.46301827_wp, 1.32928147_wp, &
    1.02766713_wp, 1.02291377_wp, 0.94343886_wp, 1.14881311_wp, 1.47080755_wp, &
    1.76901636_wp, 1.98724061_wp, 2.41244711_wp, 2.26739524_wp, 2.95378999_wp, &
    1.20807752_wp, 1.65941046_wp, 1.62733880_wp, 1.61344972_wp, 1.63220728_wp, &
    1.60899928_wp, 1.43501286_wp, 1.54559205_wp, 1.32663678_wp, 1.37644152_wp, &
    1.36051851_wp, 1.23395526_wp, 1.65734544_wp, 1.53895240_wp, 1.97542736_wp, &
    1.97636542_wp, 2.05432381_wp, 3.80138135_wp, 1.43893803_wp, 1.75505957_wp, &
    1.59815118_wp, 1.76401732_wp, 1.63999999_wp, 1.63999999_wp, 1.63999999_wp, &
    1.63999999_wp, 1.63999999_wp, 1.63999999_wp, 1.63999999_wp, 1.63999999_wp, &
    1.63999999_wp, 1.63999999_wp, 1.63999999_wp, 1.63999999_wp, 1.63999999_wp, &
    1.63999999_wp, 1.47055223_wp, 1.81127084_wp, 1.40189963_wp, 1.54015481_wp, &
    1.33721475_wp, 1.57165422_wp, 1.04815857_wp, 1.78342098_wp, 2.79106396_wp, &
    1.78160840_wp, 2.47588882_wp, 2.37670734_wp, 1.76613217_wp, 2.66172302_wp, &
    2.82773085_wp /))

      nfrag = 1
      m=n+nfrag ! # atoms + chrg constrain + frag constrain

      allocate(A(m,m),x(m),work(m*m),ipiv(m),alp(n),gam(n))
!  setup RHS
      do i=1,n
         x(i) =-chieeq(at(i)) + cnfeeq(at(i))*sqrt(cn(i)) 
       alp(i) = alpeeq(at(i)) *0.85d0 ! changed for q-vSZP (S.G.), opt. for sum of Etot of H-Kr (MG) sets
       gam(i) = gameeq(at(i)) 
      enddo

      A = 0
!  setup A matrix  
      do i=1,n
      A(i,i)=tsqrt2pi/alp(i)+gam(i)
      k = i*(i-1)/2
      do j=1,i-1
         ij = k+j
         rij=rab(ij)
         gammij=1./sqrt(alp(i)**2+alp(j)**2) 
         tmp = erf(gammij*rij)
         A(j,i) = tmp/rij
         A(i,j) = A(j,i)
      enddo
      enddo

!  fragment charge constrain
      do i=1,nfrag 
        x(n+i)=chrg    
        do j=1,n
!        if(fraglist(j).eq.i) then
            A(n+i,j)=1       
            A(j,n+i)=1      
!        endif
        enddo
      enddo

!     call prmat(6,A,m,m,'A eg')
      lwork=m*m
      call DSYSV('U', m, 1, A, m, IPIV, x, m, WORK, LWORK, INFO)
      q(1:n)=x(1:n)

      if(info.ne.0) stop 'EEQ *SYSV failed'

      if(n.eq.1) q(1)=chrg

end subroutine eeq          

! ------------------------------------------------------------------------
! standard erf CN
! ------------------------------------------------------------------------

subroutine ncoord_erf(nat,at,rab,kn,cn)
   use iso_fortran_env, only : wp => real64
   use com, only: rcov
   implicit none

   integer, intent(in)  :: nat
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: rab(nat*(nat+1)/2)
   real(wp),intent(in)  :: kn         
   real(wp),intent(out) :: cn(nat)

   integer  :: i, j, k
   real(wp) :: r, rcovij, tmp, arg

   cn  = 0._wp

   k = 0
   do i = 1, nat
      do j = 1, i-1
         k = k + 1
         rcovij=(4./3.)*(rcov(at(i))+rcov(at(j)))
         arg = (rab(k)-rcovij)/rcovij
         tmp = 0.5d0 * (1d0 + erf(kn*arg)) 
         cn(i) = cn(i) + tmp      
         cn(j) = cn(j) + tmp 
      enddo
      k = k + 1
   enddo

end subroutine ncoord_erf

! ------------------------------------------------------------------------
! standard erf CN but with changed radius for H and fixed erf exponent
! ------------------------------------------------------------------------

subroutine ncoord_basq(nat,at,rab,cn)        
   use iso_fortran_env, only : wp => real64

   implicit none

   integer, intent(in)  :: nat
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: rab(nat*(nat+1)/2)
   real(wp),intent(out) :: cn(nat)

   integer  :: i, j, k
   real(wp) :: r, rcovij, tmp, arg, kn
   real(wp) :: rcov(118) = 1.889725949_wp * [ &
!  & 0.32_wp,0.46_wp, & ! H,He
   & 0.47_wp,0.46_wp, & ! H,He
   & 1.20_wp,0.94_wp,0.77_wp,0.75_wp,0.71_wp,0.63_wp,0.64_wp,0.67_wp, & ! Li-Ne
   & 1.40_wp,1.25_wp,1.13_wp,1.04_wp,1.10_wp,1.02_wp,0.99_wp,0.96_wp, & ! Na-Ar
   & 1.76_wp,1.54_wp, & ! K,Ca
   &                 1.33_wp,1.22_wp,1.21_wp,1.10_wp,1.07_wp, & ! Sc-
   &                 1.04_wp,1.00_wp,0.99_wp,1.01_wp,1.09_wp, & ! -Zn
   &                 1.12_wp,1.09_wp,1.15_wp,1.10_wp,1.14_wp,1.17_wp, & ! Ga-Kr
   & 1.89_wp,1.67_wp, & ! Rb,Sr
   &                 1.47_wp,1.39_wp,1.32_wp,1.24_wp,1.15_wp, & ! Y-
   &                 1.13_wp,1.13_wp,1.08_wp,1.15_wp,1.23_wp, & ! -Cd
   &                 1.28_wp,1.26_wp,1.26_wp,1.23_wp,1.32_wp,1.31_wp, & ! In-Xe
   & 2.09_wp,1.76_wp, & ! Cs,Ba
   &         1.62_wp,1.47_wp,1.58_wp,1.57_wp,1.56_wp,1.55_wp,1.51_wp, & ! La-Eu
   &         1.52_wp,1.51_wp,1.50_wp,1.49_wp,1.49_wp,1.48_wp,1.53_wp, & ! Gd-Yb
   &                 1.46_wp,1.37_wp,1.31_wp,1.23_wp,1.18_wp, & ! Lu-
   &                 1.16_wp,1.11_wp,1.12_wp,1.13_wp,1.32_wp, & ! -Hg
   &                 1.30_wp,1.30_wp,1.36_wp,1.31_wp,1.38_wp,1.42_wp, & ! Tl-Rn
   & 2.01_wp,1.81_wp, & ! Fr,Ra
   &         1.67_wp,1.58_wp,1.52_wp,1.53_wp,1.54_wp,1.55_wp,1.49_wp, & ! Ac-Am
   &         1.49_wp,1.51_wp,1.51_wp,1.48_wp,1.50_wp,1.56_wp,1.58_wp, & ! Cm-No
   &                 1.45_wp,1.41_wp,1.34_wp,1.29_wp,1.27_wp, & ! Lr-
   &                 1.21_wp,1.16_wp,1.15_wp,1.09_wp,1.22_wp, & ! -Cn
   &                 1.36_wp,1.43_wp,1.46_wp,1.58_wp,1.48_wp,1.57_wp ] ! Nh-Og

   cn  = 0._wp

   kn = -3.75d0 ! opt. 

   k = 0
   do i = 1, nat
      do j = 1, i-1
         k = k + 1
         rcovij=1.00*(rcov(at(i))+rcov(at(j))) ! opt. (instead for 4/3 for standard CN)
         arg = (rab(k)-rcovij)/rcovij
         tmp = 0.5d0 * (1d0 + erf(kn*arg))         
         cn(i) = cn(i) + tmp      
         cn(j) = cn(j) + tmp 
      enddo
      k = k + 1
   enddo

end subroutine ncoord_basq

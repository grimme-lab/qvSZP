module chargscfcts
   use mctc_io
   use mctc_env, only: i4, wp
   use tblite_ceh_ceh, only: ceh_guess
   use multicharge_lapack, only : sytrf, sytrs
   implicit none
   private

   public :: eeq,calcrab,ncoord_basq,extcharges
contains
   subroutine eeq(tmpmol,rab,chrg,cn,orig,scal,scal2,scal3,scal4,q,efield)
      implicit none
      type(structure_type),intent(in)  :: tmpmol       ! Molecular structure
      real(wp),intent(in)              :: rab(tmpmol%nat*(tmpmol%nat+1)/2)   ! dist
      real(wp),intent(in)              :: chrg         ! total charge on system
      real(wp),intent(in)              :: cn(tmpmol%nat)        ! CN
      logical ,intent(in)              :: orig         ! logical. if true use original eeq model
      real(wp),intent(in)              :: scal(86)     ! scale orig xi parameters
      real(wp),intent(in)              :: scal2(86)    ! scale gamma paramters
      real(wp),intent(in)              :: scal3(86)    ! scale CN fac
      real(wp),intent(in)              :: scal4(86)    ! scale alpha
      real(wp),intent(in)              :: efield(3)    ! electric field
      real(wp),intent(out)             :: q(tmpmol%nat)         ! output charges

      integer(i4) :: info
!  local variables
      integer  :: m,i,j,k,ij,nfrag,msq
      integer,allocatable :: ipiv(:)
      real(wp) :: gammij,tmp,rij
      real(wp), allocatable :: A(:,:), x(:), work(:), alp(:), gam(:)
!  parameter
      real(wp), parameter :: tsqrt2pi = 0.797884560802866_wp
! ------------------------------------------------------------------------
!  PARAMETRISATION BY S. SPICHER (NEW)
! ------------------------------------------------------------------------
      real(wp), parameter :: chieeq(86) =(/& ! mod. by SG for wB97X-D3/q-vSZP energy see below
         1.27000000_wp, 1.25000000_wp, 0.44000000_wp, 0.99000000_wp, 1.30000000_wp, &
         1.50000000_wp, 1.63000000_wp, 1.66000000_wp, 1.66000000_wp, 1.05000000_wp, &
         0.44000000_wp, 0.68000000_wp, 1.16000000_wp, 1.16000000_wp, 1.67000000_wp, &
         1.53000000_wp, 1.60000000_wp, 1.20000000_wp, 0.34000000_wp, 0.67000000_wp, & ! Ca
         0.72000000_wp, 0.95000000_wp, 0.96000000_wp, 1.00000000_wp, 0.92000000_wp, &
         1.05000000_wp, 0.93000000_wp, 1.08000000_wp, 0.96000000_wp, 0.96000000_wp, & ! Zn
         1.11000000_wp, 1.14000000_wp, 1.45000000_wp, 1.49000000_wp, 1.48000000_wp, &
         1.11000000_wp, 0.30000000_wp, 0.54000000_wp, 0.73000000_wp, 0.91000000_wp, & ! Zr
         0.88000000_wp, 0.81000000_wp, 0.98000000_wp, 0.94000000_wp, 1.12000000_wp, &
         1.23000000_wp, 1.08000000_wp, 0.90000000_wp, 1.13000000_wp, 1.12000000_wp, &
         1.35000000_wp, 1.46000000_wp, 1.47000000_wp, 1.37000000_wp, 0.40000000_wp, & ! Cs
         0.50000000_wp, 0.68441115_wp, 0.56999999_wp, 0.56999999_wp, 0.56999999_wp, &
         0.56999999_wp, 0.56999999_wp, 0.56999999_wp, 0.56999999_wp, 0.56999999_wp, &
         0.56999999_wp, 0.56999999_wp, 0.56999999_wp, 0.56999999_wp, 0.56999999_wp, &
         0.56999999_wp, 0.87936784_wp, 1.02761808_wp, 0.93297476_wp, 1.10172128_wp, &
         0.97350071_wp, 1.16695666_wp, 1.23997927_wp, 1.18464453_wp, 1.14191734_wp, & ! Hg
         1.17000000_wp, 1.02000000_wp, 1.21000000_wp, 1.33000000_wp, 1.43000000_wp, &
         1.43000000_wp /)
      real(wp), parameter :: gameeq(86) =(/&
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
         0.10500484_wp /)
      real(wp), parameter :: cnfeeq(86) =(/& ! changed handish for RGs
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
         0.05849285_wp /)
      real(wp), parameter :: alpeeq(86) =(/&
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
         2.82773085_wp /)

      nfrag = 1
      m = tmpmol%nat + nfrag ! # atoms + chrg constrain + frag constrain
      msq = m**2
      allocate(alp(tmpmol%nat), x(m), ipiv(m), gam(tmpmol%nat), work(msq))
      allocate(A(m,m), source=0.0_wp)

!  setup RHS
      if(orig)then
         do i=1,tmpmol%nat
            x(i) =-chieeq(tmpmol%num(tmpmol%id(i))) + cnfeeq(tmpmol%num(tmpmol%id(i)))*sqrt(cn(i)) &
            & + sum( efield(1:3)*tmpmol%xyz(1:3,i) )
            alp(i) = alpeeq(tmpmol%num(tmpmol%id(i)))
            gam(i) = gameeq(tmpmol%num(tmpmol%id(i)))
         enddo
      else
         do i=1,tmpmol%nat
            x(i) =-chieeq(tmpmol%num(tmpmol%id(i))) * scal(tmpmol%num(tmpmol%id(i))) + &
            & scal3(tmpmol%num(tmpmol%id(i)))*cnfeeq(tmpmol%num(tmpmol%id(i)))*sqrt(cn(i)) &
            & + sum( efield(1:3)*tmpmol%xyz(1:3,i) )
            alp(i) = alpeeq(tmpmol%num(tmpmol%id(i))) * scal4(tmpmol%num(tmpmol%id(i)))
            gam(i) = gameeq(tmpmol%num(tmpmol%id(i))) * scal2(tmpmol%num(tmpmol%id(i)))
         enddo
      endif

!  setup A matrix
      do i=1,tmpmol%nat
         A(i,i)=tsqrt2pi/alp(i)+gam(i)
         k = i*(i-1)/2
         do j=1,i-1
            ij = k+j
            rij=rab(ij)
            gammij= 1.0_wp / sqrt( alp(i)**2 + alp(j)**2 )
            tmp = erf( gammij*rij )
            A(j,i) = tmp/rij
            A(i,j) = A(j,i)
         enddo
      enddo

!  fragment charge constrain
      do i=1,nfrag
         x(tmpmol%nat + i)=chrg
         do j=1,tmpmol%nat
!        if(fraglist(j).eq.i) then
            A(tmpmol%nat+i,j)=1
            A(j,tmpmol%nat+i)=1
!        endif
         enddo
      enddo

      call sytrf(A, ipiv, info=info, uplo='l')
      if (info == 0) then
         call sytrs(A, x, ipiv, info=info, uplo='l')
      else
         error stop 'EEQ *SYSV failed'
      endif
      q(1:tmpmol%nat)=x(1:tmpmol%nat)

      if(tmpmol%nat .eq. 1) q(1)=chrg

   end subroutine eeq

   subroutine calcrab(tmpmol,rab)
      implicit none
      type(structure_type),intent(in)  :: tmpmol
      real(wp), intent(inout) :: rab(tmpmol%nat*(tmpmol%nat+1)/2)

      integer                 ::    i,j,k
      real(wp)                ::    rab2, r

      k=0
      do i=1,tmpmol%nat
         do j=1,i-1
            k = k + 1
            rab2= &
            &              (tmpmol%xyz(1,i)-tmpmol%xyz(1,j))**2+&
            &              (tmpmol%xyz(2,i)-tmpmol%xyz(2,j))**2+&
            &              (tmpmol%xyz(3,i)-tmpmol%xyz(3,j))**2
            r = sqrt(rab2)
            rab(k) = r
            if(r.lt.1d-6) stop 'cold fusion!'
         enddo
         k = k +1
      enddo

   end subroutine calcrab

   subroutine ncoord_basq(tmpmol,rab,kn,cn)

      implicit none

      type(structure_type),intent(in)  :: tmpmol
      real(wp),intent(in)  :: rab(tmpmol%nat*(tmpmol%nat+1)/2)
      real(wp),intent(in)  :: kn
      real(wp),intent(inout) :: cn(tmpmol%nat)

      integer  :: i, j, k
      real(wp) :: rcovij, tmp, arg

      !  covalent radii (taken from Pyykko and Atsumi, Chem. Eur. J. 15, 2009,
!  188-197), values for metals decreased by 10 %
      real(wp),parameter :: rcov(118) = 1.889725949_wp * [ &
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

      cn = 0.0_wp
      k = 0
      do i = 1, tmpmol%nat
         do j = 1, i-1
            k = k + 1
            rcovij= 1.0_wp * ( rcov(tmpmol%num(tmpmol%id(i))) + rcov(tmpmol%num(tmpmol%id(j))) )
            arg = (rab(k)-rcovij)/rcovij
            tmp = 0.5_wp * (1.0_wp + erf(kn*arg))
            cn(i) = cn(i) + tmp
            cn(j) = cn(j) + tmp
         enddo
         k = k + 1
      enddo

   end subroutine ncoord_basq

   subroutine extcharges(tmpmol,fname,q,cn)
      implicit none
      type(structure_type),intent(in)  :: tmpmol         ! Molecular structure
      character(len=*), intent(in)     :: fname          ! file name of external data
      real(wp),intent(out)             :: q(tmpmol%nat)  ! output charges
      real(wp),intent(out), optional   :: cn(tmpmol%nat) ! CN
      character(len=256)               :: tmpstr
      integer                          :: i,myunit,ierr
      logical                          :: ex

      q = 0.0_wp
      if (present(cn)) cn = 0.0_wp
      inquire(file=fname,exist=ex)
      if (.not. ex) then
         write(*,*) 'WARNING: file ',trim(fname),' does not exist'
         error stop 'ERROR: file with external charges does not exist'
      else
         open(newunit=myunit,file=trim(adjustl(fname)),status='old',action='read')
         do i=1,size(q)
            read(myunit,'(a)') tmpstr
            if (present(cn)) then
               read(tmpstr,*,iostat=ierr) q(i), cn(i)
            else
               read(tmpstr,*,iostat=ierr) q(i)
            endif
         enddo
         close(myunit)
      endif

      if (present(cn)) then
         write(*,'(a,a)') 'External charges and CN read from file: ',trim(adjustl(fname))
         do i = 1,size(q)
            write(*,'(a,a2,i0,a,f12.8,a,a,i0,a,f12.8)') 'q( ',tmpmol%sym(tmpmol%id(i)),i,' ) = ',q(i), &
            & ' CN( ',tmpmol%sym(tmpmol%id(i)),i,' ) = ',cn(i)
         enddo
      else
         write(*,'(a,a)') 'External charges read from file: ',trim(adjustl(fname))
         do i = 1,size(q)
            write(*,'(a,a2,i0,a,f12.8)') 'q( ',tmpmol%sym(tmpmol%id(i)),i,' ) = ',q(i)
         enddo
      endif
      write(*,*)  ' '
   end subroutine extcharges

end module chargscfcts

module write_output
   use mctc_env, only: wp, error_type, fatal_error
   use mctc_io, only: structure_type
   use mctc_io_convert, only : autoaa
   use ioroutines, only: basis_type,ecp_type
   implicit none
   private

   public :: wrorca, orcaconfig

   type :: orcaconfig
      character(:), allocatable :: outn
      character(:), allocatable :: guess
      character(:), allocatable :: scfconv
      character(:), allocatable :: dfa

      logical          :: geoopt        = .false. ! turn on !OPT keyword
      logical          :: nocosx        = .false. ! turns off RIJCOSX, seminumerical exchange
      logical          :: notrah        = .false. ! turns off TRAH
      logical          :: nososcf       = .false. ! turns off SOSCF
      logical          :: dipgrad       = .false. ! dipole moment gradients
      logical          :: polgrad       = .false. ! polarizability derivatives
      logical          :: polar         = .false. ! polarizability calc
      logical          :: noelprop      = .false.
      logical          :: nouseshark    = .false. ! use different integral library
      logical          :: verbose       = .false.
      logical          :: indefile      = .false.
      logical          :: hfref         = .false.

      logical, allocatable :: ghostatoms(:)

      integer          :: defgrid = 2
      integer          :: coremem = 5000
      integer          :: mpi = 4

      real(wp)         :: d4_s6 = 1.00_wp
      real(wp)         :: d4_s8 = 1.00_wp
      real(wp)         :: d4_a1 = 0.35_wp
      real(wp)         :: d4_a2 = 5.60_wp
      real(wp)         :: d4_s9 = 0.00_wp
      real(wp)         :: efield(3) = 0.0_wp

   end type orcaconfig

contains

   subroutine wrorca(orcainp, mol, bas, ecp, q, cn, verbose, &
   & error, q_short, cn_short)
      type(orcaconfig), intent(in)      :: orcainp
      type(structure_type), intent(in)   :: mol
      type(basis_type), intent(in)      :: bas
      type(ecp_type), intent(in)        :: ecp
      real(wp), intent(in)              :: q(mol%nat)
      real(wp), intent(in)              :: cn(mol%nat)
      logical, intent(in)               :: verbose
      type(error_type), allocatable, intent(out)     :: error
      real(wp), intent(in), optional    :: q_short(:)
      real(wp), intent(in), optional    :: cn_short(:)

      integer :: myunit,i,n,j,k

      character(len=1)     :: ltmp

      logical :: ecpex

      real(wp) :: scalfac
      real(wp),allocatable             :: sccoeff(:,:,:)

      allocate(sccoeff(mol%nat,20,20))

      ! Start writing the ORCA input file
      open(newunit=myunit,file=orcainp%outn)

      if (orcainp%hfref) then
         write(myunit,'(a)')     "! HF NORI ExtremeSCF"
         write(myunit,'(a)')     "! NOCOSX"
         write(myunit,'(a,/)')   "! LargePrint"
         write(myunit,'(a)')     "%scf"
         write(myunit,'(a,a)')   "  guess ", trim(adjustl(orcainp%guess))
         write(myunit,'(a)')     "  SCFMode Conventional"
         write(myunit,'(a,/)')   "end"
      else
         write(myunit,'(a,a,a)') "! ",orcainp%dfa, " def2/J PrintBasis"
         write(myunit,'(a,a)',advance='NO') "! ",orcainp%scfconv
         write(myunit,'(1x,a,i1)') "DEFGRID", orcainp%defgrid
         if(orcainp%notrah) write(myunit,'(''! NoTRAH'')')
         if(orcainp%nososcf) write(myunit,'(''! NoSOSCF'')')
         if(orcainp%geoopt) write(myunit,'(''! Opt'')')
         if(orcainp%nocosx) write(myunit,'(''! NOCOSX'')')
         if(orcainp%dipgrad) write(myunit,'(''! Freq'')')
         if(orcainp%polgrad) write(myunit,'(''! NumFreq'')')
         if(orcainp%nouseshark) write(myunit,'(''! NoUseShark'',/)')

         if(orcainp%mpi.gt.0) write(myunit,'(''%pal'',/,''  nprocs'',i4,/,''end'',/)') orcainp%mpi
         write(myunit,'(''%MaxCore '',i6,/)') orcainp%coremem

         write(myunit,'(a)')    "%method"
         write(myunit,'(a,f9.4)')    "  D4S6  ", orcainp%d4_s6
         write(myunit,'(a,f9.4)')    "  D4S8  ", orcainp%d4_s8
         write(myunit,'(a,f9.4)')    "  D4A1  ", orcainp%d4_a1
         write(myunit,'(a,f9.4)')    "  D4A2  ", orcainp%d4_a2
         write(myunit,'(a,f9.4)')    "  D4S9  ", orcainp%d4_s9
         write(myunit,'(a,/)')  "end"

         write(myunit, '(a)') "%scf"
         write(myunit,'(''  guess '',a)') trim(adjustl(orcainp%guess))
         if (sum(abs(orcainp%efield)) > 0.0_wp) then
            write(myunit,'(a,f12.8,a,f12.8,a,f12.8)') "  efield", &
            & orcainp%efield(1),", ", orcainp%efield(2),", ", orcainp%efield(3)
         endif
         write(myunit, '(a,/)') "end"
      endif

      if (orcainp%noelprop) then
         write(myunit,'(a)') "%elprop"
         write(myunit,'(a)') "  Dipole false"
         write(myunit,'(a)') "  Quadrupole false"
         write(myunit,'(a)') "end"
      endif

      if(orcainp%polar.or.orcainp%polgrad) write(myunit,'(''%elprop polar 1'',/,''end'',/)')

      n = 0
      sccoeff = 0.0_wp
      do i = 1, mol%nat
         if (.not. bas%sccoeff(mol%num(mol%id(i)))) cycle
         if (present(q_short)) then
            if (.not.orcainp%ghostatoms(i)) then
               n = n + 1
               scalfac = q_short(n) - ( bas%scalparam(mol%num(mol%id(i)),2) * q_short(n)**2 ) & ! Dependency on q and q^2
                  + ( bas%scalparam(mol%num(mol%id(i)),1) * sqrt(cn_short(n)) ) &                  ! Dependency on square root of CN
                  + ( bas%scalparam(mol%num(mol%id(i)),3) * cn_short(n) * q_short(n) )          ! Cross term dependency of CN and q
            else
               scalfac = q(i) - ( bas%scalparam(mol%num(mol%id(i)),2) * q(i)**2 ) & ! Dependency on q and q^2
                  + ( bas%scalparam(mol%num(mol%id(i)),1) * sqrt(cn(i)) ) &            ! Dependency on square root of CN
                  + ( bas%scalparam(mol%num(mol%id(i)),3) * cn(i) * q(i) )          ! Cross term dependency of CN and q
            endif
         else
            scalfac = q(i) - ( bas%scalparam(mol%num(mol%id(i)),2) * q(i)**2 ) & ! Dependency on q and q^2
               + ( bas%scalparam(mol%num(mol%id(i)),1) * sqrt(cn(i)) ) &            ! Dependency on square root of CN
               + ( bas%scalparam(mol%num(mol%id(i)),3) * cn(i) * q(i) )          ! Cross term dependency of CN and q
         endif
         do j = 1, bas%nbf(mol%num(mol%id(i)))
            do k = 1, bas%npr(mol%num(mol%id(i)),j)
               sccoeff(i,j,k) = bas%coeff(mol%num(mol%id(i)),j,k) + &
                  scalfac * bas%qcoeff(mol%num(mol%id(i)),j,k)
            enddo
         enddo
      enddo
      write(*,'(a)') ""

      ecpex = .false.
      do i=1,maxval(mol%id,1)
         if ( sum(abs(ecp%exp(mol%num(i),:,:))) > 0.0_wp .and. mol%num(i) >= ecp%atmin ) then
            if ( .not. ecpex ) then
               write(myunit,'(a)') "%basis"
               ! Check if an atom with ordinal number > 86 is in the molecule
               if ( any(mol%num > 86) ) then
                  ! Print the following line to the output file
                  !   AuxJ  "AutoAux"      # Use AutoAux to generate the AuxJ fitting basis set
                  write(myunit,'(a)') '  AuxJ  "AutoAux" # Use AutoAux to generate the AuxJ fitting basis set'
               end if
               ecpex = .true.
            endif
            write(myunit,'(a,a2)') "  NewECP ", mol%sym(i)
            write(myunit,'(a,i2)') "  N_core ", ecp%ncore(mol%num(i))
            if (ecp%lmax(mol%num(i)) == 0) ltmp = 's'
            if (ecp%lmax(mol%num(i)) == 1) ltmp = 'p'
            if (ecp%lmax(mol%num(i)) == 2) ltmp = 'd'
            if (ecp%lmax(mol%num(i)) == 3) ltmp = 'f'
            if (ecp%lmax(mol%num(i)) == 4) ltmp = 'g'
            if (ecp%lmax(mol%num(i)) == 5) ltmp = 'h'
            if (ecp%lmax(mol%num(i)) == 6) ltmp = 'i'
            write(myunit,'(a,a2)') "  lmax ", ltmp
            do j=1,ecp%nbf(mol%num(i))
               if (ecp%angmom(mol%num(i),j) == 0) ltmp = 's'
               if (ecp%angmom(mol%num(i),j) == 1) ltmp = 'p'
               if (ecp%angmom(mol%num(i),j) == 2) ltmp = 'd'
               if (ecp%angmom(mol%num(i),j) == 3) ltmp = 'f'
               if (ecp%angmom(mol%num(i),j) == 4) ltmp = 'g'
               if (ecp%angmom(mol%num(i),j) == 5) ltmp = 'h'
               if (ecp%angmom(mol%num(i),j) == 6) ltmp = 'i'
               write(myunit,'(3x,a1,2x,i3)') ltmp, ecp%npr(mol%num(i),ecp%sindex(mol%num(i),j))
               do k=1,ecp%npr(mol%num(i),ecp%sindex(mol%num(i),j))
                  write(myunit,'(i5,2x,2f14.8,2x,i3)') k, ecp%exp(mol%num(i),ecp%sindex(mol%num(i),j),k), &
                     ecp%coeff(mol%num(i),ecp%sindex(mol%num(i),j),k), &
                     ecp%nfactor(mol%num(i),ecp%sindex(mol%num(i),j),k)
               enddo
            enddo
            write(myunit,'(2x,a)') "end"
         else
            if (mol%num(i) <= 2) then
               if (verbose) write(*,'(a,a,a,/)') "No ECP assigned for element ",trim(mol%sym(i)),"."
            else
               write(*,'(a,a,a)') "ERROR: No ECP assigned for element ",trim(mol%sym(i)),"."
               error stop
            endif
         endif
      enddo
      if (ecpex) write(myunit,'(a,/)') "end"

      if (allocated(error)) then
         print '(a)', error%message
         error stop
      end if

      write(myunit,'(a,2x,2i3)') "* xyz", nint(mol%charge), mol%uhf+1
      do i=1,mol%nat
         if (present(q_short)) then
            if (orcainp%ghostatoms(i)) then
               write(myunit,'(a2,2x,a,2x,3F20.14)') mol%sym(mol%id(i)),":",mol%xyz(:,i)*autoaa
            else
               write(myunit,'(a2,2x,3F20.14)') mol%sym(mol%id(i)),mol%xyz(:,i)*autoaa
            endif
         else
            write(myunit,'(a2,2x,3F20.14)') mol%sym(mol%id(i)),mol%xyz(:,i)*autoaa
         endif
         if ( bas%exp(mol%num(mol%id(i)),1,1) > 0.0_wp ) then
            write(myunit,'(2x,a,1x,a2)') "NewGTO"
            do j=1,bas%nbf(mol%num(mol%id(i)))
               if (bas%angmom(mol%num(mol%id(i)),j) == 0) ltmp = 'S'
               if (bas%angmom(mol%num(mol%id(i)),j) == 1) ltmp = 'P'
               if (bas%angmom(mol%num(mol%id(i)),j) == 2) ltmp = 'D'
               if (bas%angmom(mol%num(mol%id(i)),j) == 3) ltmp = 'F'
               if (bas%angmom(mol%num(mol%id(i)),j) == 4) ltmp = 'G'
               if (bas%angmom(mol%num(mol%id(i)),j) == 5) ltmp = 'H'
               if (bas%angmom(mol%num(mol%id(i)),j) == 6) ltmp = 'I'
               write(myunit,'(3x,a1,2x,i3)') ltmp,bas%npr(mol%num(mol%id(i)),j)
               do k=1,bas%npr(mol%num(mol%id(i)),j)
                  write(myunit,'(i5,2x,2f20.14)') k, &
                     bas%exp(mol%num(mol%id(i)),j,k),sccoeff(i,j,k)
               enddo
            enddo
            write(myunit,'(2x,a)') "end"
            if (mol%num(i) > 86) then
               write(myunit,'(2x,a,/,5x,a,/,2x,a)') 'NewAuxJGTO', '"AutoAux"', 'end'
            end if
         else
            call fatal_error(error,"No basis set for desired element available.")
         endif
      enddo
      write(myunit,'(a)') "*"

      close(myunit)

      write(*,'(a,a,a)') "Successfully wrote input file: '",orcainp%outn,"'"

   end subroutine wrorca

end module write_output

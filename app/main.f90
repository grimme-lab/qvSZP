program main
   use mctc_io
   use mctc_io_convert, only : autoaa
   use mctc_env
   use ioroutines, only: rdfile,rdbas,rdecp
   use basistype, only: basis_type,ecp_type
   use miscellaneous, only: helpf
   implicit none
   integer              :: narg,i,myunit,j,k
   integer              :: mpi      = 4
   integer              :: defgrid  = 2
   integer              :: charge   = 0
   integer              :: nopen    = 0
   integer              :: chrg     = 0
   integer              :: coremem  = 5000

   character(len=120)   :: atmp,guess,filen,outn,bfilen,efilen
   character(len=:),allocatable :: scfconv
   character(len=1)     :: ltmp

   logical              :: indguess, polar, beta, polgrad, dipgrad, geoopt, nocosx
   logical              :: tightscf, strongscf, verbose, suborca, nouseshark,sugg_guess
   logical              :: ploteldens, help, uhfgiven, da,indbfile,indefile,indcharge

   type(structure_type)             :: mol
   type(error_type), allocatable    :: error
   type(basis_type)                 :: bas
   type(ecp_type)                   :: ecp

   filen       = 'coord' ! input  filename
   outn        = 'wb97x3c.inp'   ! output filename
   scfconv     = 'NormalSCF'
   polar       = .false. ! polarizability calc
   beta        = .false. ! hyperpolarizabilities
   polgrad     = .false. ! polarizability derivatives
   dipgrad     = .false. ! dipole moment gradients
   geoopt      = .false. ! turn on !OPT keyword
   nocosx      = .false. ! turns on RIJCOSX, seminumerical exchange
   verbose     = .false. ! verbose output
   nouseshark  = .false. ! use different integral library
   tightscf    = .false. ! SCF conv criterium
   strongscf   = .false. ! SCF conv criterium
   indguess    = .false.
   uhfgiven    = .false.
   help        = .false.
   ploteldens  = .false.
   indbfile    = .false.
   indefile    = .false.
   indcharge   = .false.
   sugg_guess  = .false.

   ! get number of arguments
   narg = command_argument_count()

   do i=1,narg
      call get_command_argument(i,atmp)
      if(index(atmp,'--mpi').ne.0) then
         call get_command_argument(i+1,atmp)
         read(atmp,*) mpi
      endif
      call get_command_argument(i,atmp)
      if(index(atmp,'--struc').ne.0) then
         call get_command_argument(i+1,atmp)
         filen = trim(atmp)
      endif
      if(index(atmp,'--basisfile').ne.0) then
         call get_command_argument(i+1,atmp)
         indbfile=.true.
         bfilen = trim(atmp)
      endif
      if(index(atmp,'--ecpfile').ne.0) then
         call get_command_argument(i+1,atmp)
         indefile=.true.
         efilen = trim(atmp)
      endif
      if(index(atmp,'--chrg').ne.0) then
         call get_command_argument(i+1,atmp)
         indcharge=.true.
         read(atmp,*) charge
      endif
      if(index(atmp,'--uhf').ne.0) then
         call get_command_argument(i+1,atmp)
         uhfgiven=.true.
         read(atmp,*) nopen
      endif
      if(index(atmp,'--memory').ne.0) then
         call get_command_argument(i+1,atmp)
         read(atmp,*) coremem
      endif
      if(index(atmp,'--defgrid').ne.0) then
         call get_command_argument(i+1,atmp)
         read(atmp,*) defgrid
      endif
      if(index(atmp,'--guess').ne.0) then
         indguess=.true.
         call get_command_argument(i+1,atmp)
         guess=trim(atmp)
      endif
      if(index(atmp,'--polar').ne.0) polar=.true.
      ! if(index(atmp,'-hyppol').ne.0) beta=.true.
      if(index(atmp,'--polgrad').ne.0) polgrad=.true.
      if(index(atmp,'--dipgrad').ne.0) dipgrad=.true.
      if(index(atmp,'--geoopt').ne.0) geoopt=.true.
      if(index(atmp,'--nocosx').ne.0) nocosx=.true.
      if(index(atmp,'--tightscf').ne.0)   scfconv='TightSCF'
      if(index(atmp,'--strongscf').ne.0)  scfconv='StrongSCF'
      if(index(atmp,'--v').ne.0) verbose=.true.
      if(index(atmp,'--suborca').ne.0) suborca=.true.
      if(index(atmp,'--nouseshark').ne.0) nouseshark=.true.
      if(index(atmp,'--help').ne.0) help=.true.
      if(index(atmp,'--plot').ne.0) ploteldens=.true.
      if(index(atmp,'--suggestedguess').ne.0) sugg_guess=.true.
   enddo
   if (.not. uhfgiven) then
      inquire(file='.UHF',exist=da)
      if(da)then
         open(newunit=myunit,file='.UHF')
         read(myunit,*) nopen
         uhfgiven=.true.
         close(myunit)
      endif
   endif

   if (.not. indcharge) then
      inquire(file='.CHRG',exist=da)
      if(da)then
         open(newunit=myunit,file='.CHRG')
         read(myunit,*) charge
         close(myunit)
      endif
   endif

   if (help) then
      call helpf()
      stop
   endif

   call rdfile(trim(filen),mol,chrg)
   chrg = chrg - charge
   if (.not. uhfgiven .and. (chrg .eq. 1 .or. (mod(chrg,2) .ne. 0))) then
      write(*,'(a)') "Use a .UHF file or '--uhf <int>' to indicate the number of unpaired electrons."
      call fatal_error(error, "Odd number of electrons for closed-shell calculation. ")
      if (allocated(error)) then
         print '(a)', error%message
         error stop
      end if
   endif

! start writing
   open(newunit=myunit,file=outn)
   write(myunit,'(''! RKS WB97X-D4 def2/J'')')
   if (verbose) then
      write(myunit,'(''! PRINTBASIS LARGEPRINT'')')
   endif
   write(myunit,'(a,a)',advance='NO') "! ",scfconv
   write(myunit,'(1x,a,i1,/)') "DEFGRID", defgrid
   if(geoopt) write(myunit,'(''! Opt'')')
   if(nocosx) write(myunit,'(''! NOCOSX'')')
   if(dipgrad) write(myunit,'(''! Freq'')')
   if(polgrad) write(myunit,'(''! NumFreq'')')
   if(nouseshark) write(myunit,'(''! NoUseShark'',/)')

   if(mpi.gt.0) write(myunit,'(''%pal'',/,''  nprocs'',i4,/,''end'',/)') mpi
   write(myunit,'(''%MaxCore '',i6,/)') coremem

   write(myunit,'(a)')    "%method"
   write(myunit,'(a)')    "  D4A1    0.2464"
   write(myunit,'(a)')    "  D4A2    4.737"
   write(myunit,'(a)')    "  D4S6    1.00"
   write(myunit,'(a)')    "  D4S8    0.00"
   write(myunit,'(a)')    "  D4S9    1.00"
   write(myunit,'(a,/)')  "end"


   if(.not.indguess) then
      if (sugg_guess) then
         if (any((mol%num >= 37 .and. mol%num <= 45) .or. (mol%num >= 48 .and. mol%num <= 54))) then
            indguess=.true.
            guess='hueckel'
         endif
         if (any(mol%num == 21 .or. mol%num == 47 .or. &
            mol%num == 74 .or. mol%num == 82 .or. mol%num == 83)) then
            indguess=.true.
            guess='hcore'
         endif
      endif
   endif

   if (indguess .or. beta) then
      write(myunit, '(a)') "%scf"
      if (indguess) write(myunit,'(''  guess '',a20)') guess
      ! if (beta) write(myunit, '(2x,a)') "efield XXX,YYY,ZZZ"
      write(myunit, '(a,/)') "end"
   endif

   if(polar.or.polgrad.or.beta) write(myunit,'(''%elprop polar 1'',/,''end'',/)')
   if (verbose) then
      write(myunit,'(''%output'')')
      write(myunit,'(''       print[P_Hirshfeld] 1'')')
      write(myunit,'(''       print[P_BondOrder_M] 1'')')
      write(myunit,'(''       print[P_basis] 2'')')
      write(myunit,'(''end'')')
   endif

   if (indbfile) then
      call rdbas(bas,verbose,trim(bfilen))
   else
      call rdbas(bas,verbose)
   endif

   write(myunit,'(a)') "%basis"
   do i=1,maxval(mol%id,1)
      if ( bas%exp(mol%num(i),1,1) > 0.0_wp ) then
         write(myunit,'(2x,a,1x,a2)') "NewGTO",mol%sym(i)
         do j=1,bas%nbf(mol%num(i))
            if (bas%angmom(mol%num(i),j) == 0) ltmp = 'S'
            if (bas%angmom(mol%num(i),j) == 1) ltmp = 'P'
            if (bas%angmom(mol%num(i),j) == 2) ltmp = 'D'
            if (bas%angmom(mol%num(i),j) == 3) ltmp = 'F'
            if (bas%angmom(mol%num(i),j) == 4) ltmp = 'G'
            if (bas%angmom(mol%num(i),j) == 5) ltmp = 'H'
            if (bas%angmom(mol%num(i),j) == 6) ltmp = 'I'
            write(myunit,'(3x,a1,2x,i3)') ltmp,bas%npr(mol%num(i),j)
            do k=1,bas%npr(mol%num(i),j)
               write(myunit,'(i5,2x,2f14.8)') k,bas%exp(mol%num(i),j,k),bas%coeff(mol%num(i),j,k)
            enddo
         enddo
         write(myunit,'(2x,a)') "end"
      else
         call fatal_error(error,"No basis set for atoms with Z > 86")
      endif
   enddo

   if (indefile) then
      call rdecp(ecp,verbose,efilen)
   else
      call rdecp(ecp,verbose)
   endif

   do i=1,maxval(mol%id,1)
      if ( sum(abs(ecp%exp(mol%num(i),:,:))) > 0.0_wp .and. mol%num(i) >= ecp%atmin ) then
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
            if (verbose) then
               write(*,*) "sortindex of element i and bf j: ",mol%num(i), j, ecp%sindex(mol%num(i),j)
            endif
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
         write(*,'(a,a,a,/)') "No ECP assigned for element ",trim(mol%sym(i)),"."
      endif
   enddo
   write(myunit,'(a,/)') "end"

   if (allocated(error)) then
      print '(a)', error%message
      error stop
   end if

   if (ploteldens) then
      write(myunit,'(a)') "%plots"
      write(myunit,'(a)') "  dim1 60"
      write(myunit,'(a)') "  dim2 60"
      write(myunit,'(a)') "  dim3 60"
      write(myunit,'(a)') "  min1 -15"
      write(myunit,'(a)') "  max1  15"
      write(myunit,'(a)') "  min2 -15"
      write(myunit,'(a)') "  max2  15"
      write(myunit,'(a)') "  min3 -15"
      write(myunit,'(a)') "  max3  15"
      write(myunit,'(a)') "  Format Gaussian_Cube"
      write(myunit,'(a)') '  Eldens("eldens.cube");'
      write(myunit,'(a,/)') "end"
   endif

   write(myunit,'(a,2x,2i3)') "* xyz", charge, nopen+1
   do i=1,mol%nat
      write(myunit,'(a2,2x,3F22.14)') mol%sym(mol%id(i)),mol%xyz(:,i)*autoaa
   enddo
   write(myunit,'(a)') "*"

   close(myunit)

   write(*,'(a,a)') "Successfully wrote input file: ",outn

end program main

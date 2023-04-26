program main
   use mctc_io
   use mctc_io_convert, only : autoaa
   use mctc_env
   use ioroutines, only: rdfile,rdbas,rdecp,check_ghost_atoms,search_ghost_atoms
   use basistype, only: basis_type,ecp_type
   use chargscfcts, only: eeq,calcrab,ncoord_basq
   use miscellaneous, only: helpf
   implicit none
   integer              :: narg,i,myunit,j,k,n
   integer              :: mpi      = 4
   integer              :: defgrid  = 2
   integer              :: charge   = 0
   integer              :: nopen    = 0
   integer              :: chrg     = 0
   integer              :: coremem  = 5000
   integer,allocatable  :: tmpids(:)
   integer,allocatable  :: tmpids_short(:)

   character(len=120)   :: atmp,guess,filen,outn,bfilen,efilen
   character(len=:),allocatable :: scfconv
   character(len=1)     :: ltmp

   logical, allocatable :: ghostatoms(:)
   logical              :: dummy

   logical              :: indguess, polar, polgrad, dipgrad, geoopt, nocosx
   logical              :: tightscf, strongscf, verbose, suborca, nouseshark,sugg_guess
   logical              :: help, uhfgiven, da,indbfile,indefile,indcharge,indd4param
   logical              :: ecpex

   type(structure_type)             :: mol,molshort
   type(error_type), allocatable    :: error
   type(basis_type)                 :: bas
   type(ecp_type)                   :: ecp

   real(wp),allocatable             :: distvec(:),cnvec(:),qeeq(:), sccoeff(:,:,:)
   real(wp),allocatable             :: distvec_short(:),cnvec_short(:),qeeq_short(:)

   real(wp)                         :: scalfac
   real(wp)                         :: d4_s6 = 1.00_wp
   real(wp)                         :: d4_s8 = 1.00_wp
   real(wp)                         :: d4_a1 = 0.35_wp
   real(wp)                         :: d4_a2 = 5.60_wp
   real(wp)                         :: d4_s9 = 0.00_wp

   ! ####### Set up EEQ parameters for coefficient scaling #######
   real(wp),parameter               :: unity(86)      = 1.0_wp
   real(wp),parameter               :: alphascal(86)  = 0.75_wp
   real(wp),parameter               :: chiscal(86)    = 1.20_wp
   real(wp),parameter               :: gamscal(86)    = 0.90_wp
   ! ###########################################################

   filen       = 'coord' ! input  filename
   outn        = 'wb97xd4-qvszp.inp'   ! output filename
   scfconv     = 'NormalSCF'
   polar       = .false. ! polarizability calc
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
         filen=trim(adjustl(atmp))
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
      if(index(atmp,'--d4par').ne.0) then
         indd4param=.true.
         call get_command_argument(i+1,atmp)
         read(atmp,*) d4_s6
         call get_command_argument(i+2,atmp)
         read(atmp,*) d4_s8
         call get_command_argument(i+3,atmp)
         read(atmp,*) d4_a1
         call get_command_argument(i+4,atmp)
         read(atmp,*) d4_a2
         call get_command_argument(i+5,atmp)
         read(atmp,*) d4_s9
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

   !####### END OF INITIALIZATION PART ######

   if (index((filen),'.xyz').ne.0) then
      call check_ghost_atoms(filen,dummy,error)
      if (allocated(error)) then
         print '(a)', error%message
         error stop
      end if
   end if

   if (dummy) then
      call search_ghost_atoms(filen,ghostatoms,error)
      if (allocated(error)) then
         print '(a)', error%message
         error stop
      end if
   end if

   if (dummy) then
      call rdfile("struc_raw.xyz",mol,chrg)
      call rdfile("struc_short.xyz",molshort)
   else
      call rdfile(filen,mol,chrg)
   end if

   allocate(tmpids(mol%nat))
   allocate(distvec(mol%nat*(mol%nat+1)/2))
   allocate(cnvec(mol%nat),qeeq(mol%nat))
   allocate(sccoeff(mol%nat,20,20))

   if (dummy) then
      allocate(cnvec_short(molshort%nat),qeeq_short(molshort%nat) &
         ,distvec_short(molshort%nat*(molshort%nat+1)/2))
      allocate(tmpids_short(molshort%nat))
   endif

   do i = 1, mol%nat
      tmpids(i) = mol%num(mol%id(i))
   enddo
   if (dummy) then
      do i = 1, molshort%nat
         tmpids_short(i) = molshort%num(molshort%id(i))
      enddo
   endif

   distvec  =  0.0_wp
   cnvec    =  0.0_wp
   qeeq     =  0.0_wp
   if (dummy) then
      distvec_short  =  0.0_wp
      cnvec_short    =  0.0_wp
      qeeq_short     =  0.0_wp
   endif

   call calcrab(mol,distvec)
   call ncoord_basq(mol,distvec,-3.75_wp,cnvec)
   if (dummy) then
      call calcrab(molshort,distvec_short)
      call ncoord_basq(molshort,distvec_short,-3.75_wp,cnvec_short)
   endif

   chrg = chrg - charge
   if (.not. uhfgiven .and. (chrg .eq. 1 .or. (mod(chrg,2) .ne. 0))) then
      write(*,'(a)') "Use a .UHF file or '--uhf <int>' to indicate the number of unpaired electrons."
      call fatal_error(error, "Odd number of electrons for closed-shell calculation. ")
      if (allocated(error)) then
         print '(a)', error%message
         error stop
      end if
   endif

   call eeq(mol%nat,tmpids,distvec,real(charge,wp),cnvec,.False.,unity,gamscal,chiscal,alphascal,qeeq)
   if (dummy) then
      if (charge /= 0) then
         error stop "Non-zero charge not supported for dummy atoms."
      endif
      call eeq(molshort%nat,tmpids_short,distvec_short, &
         real(charge,wp),cnvec_short,.False.,unity,gamscal,chiscal,alphascal,qeeq_short)
   endif

   if (indbfile) then
      call rdbas(bas,verbose,trim(bfilen))
   else
      call rdbas(bas,verbose)
   endif

   ! Start writing the ORCA input file
   open(newunit=myunit,file=outn)
   write(myunit,'(a)') "! WB97X-D4 def2/J PrintBasis"
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
   write(myunit,'(a,f5.2)')    "  D4S6  ", d4_s6
   write(myunit,'(a,f5.2)')    "  D4S8  ", d4_s8
   write(myunit,'(a,f5.2)')    "  D4A1  ", d4_a1
   write(myunit,'(a,f5.2)')    "  D4A2  ", d4_a2
   write(myunit,'(a,f5.2)')    "  D4S9  ", d4_s9
   write(myunit,'(a,/)')  "end"

   if (sugg_guess) then
      guess='hueckel'
   endif

   if (sugg_guess .or. indguess) then
      write(myunit, '(a)') "%scf"
      write(myunit,'(''  guess '',a20)') guess
      write(myunit, '(a,/)') "end"
   endif

   if(polar.or.polgrad) write(myunit,'(''%elprop polar 1'',/,''end'',/)')

   if (verbose) then
      if (dummy) then
         write(*,'(a)') "WARNING!! Verbose output not optimized for presence of dummy atoms."
      endif
      do i = 1, mol%nat
         write(*,'(1x,a,i3,a,i3)')               "Z for atom               ",i,":",mol%num(mol%id(i))
         write(*,'(1x,a,i3,a,f9.5,f9.5)')               "scaling factors for atom ",i,":", &
            bas%scalparam(mol%num(mol%id(i)),1),bas%scalparam(mol%num(mol%id(i)),2)
         write(*,'(1x,a,i3,a,f9.5)')               "qEEQ of atom             ",i,":",qeeq(i)
         write(*,'(1x,a,i3,a,f9.5,/)') "CN of atom               ",i,":",cnvec(i)
      enddo
   endif

   n = 0
   sccoeff = 0.0_wp
   if (verbose) then
      write(*,'(a)') "Scaling prefactor for: #atom, #bf, #primitive, final coeff., scaling factor"
   endif
   do i = 1, mol%nat
      if (.not. bas%sccoeff(mol%num(mol%id(i)))) cycle
      if (dummy) then
         if (.not.ghostatoms(i)) then
            n = n + 1
            scalfac = qeeq_short(n) - ( bas%scalparam(mol%num(mol%id(i)),2) * qeeq_short(n)**2 ) &
               + ( bas%scalparam(mol%num(mol%id(i)),1) * sqrt(cnvec_short(n)) )
         else
            scalfac = qeeq(i) - ( bas%scalparam(mol%num(mol%id(i)),2) * qeeq(i)**2 ) &
               + ( bas%scalparam(mol%num(mol%id(i)),1) * sqrt(cnvec(i)) )
         endif
      else
         scalfac = qeeq(i) - ( bas%scalparam(mol%num(mol%id(i)),2) * qeeq(i)**2 ) &
            + ( bas%scalparam(mol%num(mol%id(i)),1) * sqrt(cnvec(i)) )
      endif
      do j = 1, bas%nbf(mol%num(mol%id(i)))
         do k = 1, bas%npr(mol%num(mol%id(i)),j)
            sccoeff(i,j,k) = bas%coeff(mol%num(mol%id(i)),j,k) + &
               scalfac * bas%qcoeff(mol%num(mol%id(i)),j,k)
            if (verbose) then
               write(*,'(a,i3,i3,i3,f14.8,f14.8)') "Scaling prefactor for: ",i,j,k,sccoeff(i,j,k), scalfac
            endif
         enddo
      enddo
   enddo
   write(*,'(a)') ""

   if (indefile) then
      call rdecp(ecp,verbose,efilen)
   else
      call rdecp(ecp,verbose)
   endif

   ecpex = .false.
   do i=1,maxval(mol%id,1)
      if ( sum(abs(ecp%exp(mol%num(i),:,:))) > 0.0_wp .and. mol%num(i) >= ecp%atmin ) then
         if ( .not. ecpex ) then
            write(myunit,'(a)') "%basis"
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
   if (ecpex) write(myunit,'(a,/)') "end"

   if (allocated(error)) then
      print '(a)', error%message
      error stop
   end if

   write(myunit,'(a,2x,2i3)') "* xyz", charge, nopen+1
   do i=1,mol%nat
      if (dummy) then
         if (ghostatoms(i)) then
            write(myunit,'(a2,2x,a,2x,3F22.14)') mol%sym(mol%id(i)),":",mol%xyz(:,i)*autoaa
         else
            write(myunit,'(a2,2x,3F22.14)') mol%sym(mol%id(i)),mol%xyz(:,i)*autoaa
         endif
      else
         write(myunit,'(a2,2x,3F22.14)') mol%sym(mol%id(i)),mol%xyz(:,i)*autoaa
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
               write(myunit,'(i5,2x,2f14.8)') k, &
                  bas%exp(mol%num(mol%id(i)),j,k),sccoeff(i,j,k)
            enddo
         enddo
         write(myunit,'(2x,a)') "end"
      else
         call fatal_error(error,"No basis set for atoms with Z > 86")
      endif
   enddo
   write(myunit,'(a)') "*"

   close(myunit)

   write(*,'(a,a)') "Successfully wrote input file: ",outn

end program main
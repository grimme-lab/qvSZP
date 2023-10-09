program main
   use mctc_io, only: structure_type, write_structure
   use mctc_env, only: error_type, wp, fatal_error
   use ioroutines, only: rdfile, rdbas, rdecp_default, rdecp_qvSZPs, check_ghost_atoms, &
   & search_ghost_atoms, basis_type, ecp_type
   use chargscfcts, only: eeq,calcrab,ncoord_basq,extcharges,ceh
   use miscellaneous, only: helpf
   use write_output, only: orcaconfig, wrorca
   implicit none

   integer              :: narg,i,myunit,errint
   integer              :: charge   = 0
   integer              :: nopen    = 0
   integer              :: chrg     = 0
   integer              :: verbosity = 0
   integer,allocatable  :: tmpids(:)
   integer,allocatable  :: tmpids_short(:)

   character(len=180)   :: atmp
   character(len=:), allocatable :: cm,version,filen,bfilen,efilen,extcall

   logical              :: dummy = .false.
   logical              :: qvSZPs = .false.
   logical              :: tightscf, strongscf, verbose
   logical              :: help, uhfgiven, da,indbfile,indefile,indcharge,indd4param
   logical              :: prversion

   type(structure_type)             :: mol,molshort
   type(error_type), allocatable    :: error
   type(basis_type)                 :: bas
   type(ecp_type)                   :: ecp
   type(orcaconfig)                 :: orcainp

   real(wp),allocatable             :: distvec(:),cn(:),q(:), sccoeff(:,:,:), increment(:)
   real(wp),allocatable             :: distvec_short(:),cn_short(:),q_short(:)

   real(wp)                         :: energy_shift = 0.0_wp

   ! ####### Set up EEQ parameters for coefficient scaling #######
   real(wp),parameter               :: unity(86)      = 1.0_wp
   real(wp),parameter               :: alphascal(86)  = 0.75_wp
   real(wp),parameter               :: chiscal(86)    = 1.20_wp
   real(wp),parameter               :: gamscal(86)    = 0.90_wp
   ! ###########################################################

   verbose     = .false. ! verbose output
   tightscf    = .false. ! SCF conv criterium
   strongscf   = .false. ! SCF conv criterium
   uhfgiven    = .false.
   help        = .false.
   indbfile    = .false.
   indefile    = .false.
   indcharge   = .false.
   prversion   = .false.

   filen             = 'coord' ! input  filename
   orcainp%outn      = 'wb97xd4-qvszp.inp'   ! output filename
   orcainp%scfconv   = 'NormalSCF'
   orcainp%guess     = "PModel"
   orcainp%dfa       = "wB97X-D4"
   cm                = "ceh"

   !##### VERSION STRING #####
   version = "2.1"
   !##########################

   ! get number of arguments
   narg = command_argument_count()

   do i=1,narg
      call get_command_argument(i,atmp)
      if(index(atmp,'--mpi').ne.0) then
         call get_command_argument(i+1,atmp)
         read(atmp,*) orcainp%mpi
      endif
      call get_command_argument(i,atmp)
      if(index(atmp,'--struc').ne.0) then
         call get_command_argument(i+1,atmp)
         filen=trim(adjustl(atmp))
      endif
      if(index(atmp,'--bfile').ne.0) then
         call get_command_argument(i+1,atmp)
         indbfile=.true.
         bfilen = trim(adjustl(atmp))
      endif
      if(index(atmp,'--efile').ne.0) then
         call get_command_argument(i+1,atmp)
         indefile=.true.
         efilen = trim(adjustl(atmp))
      endif
      !####### SET ORCA q-vSZPs DEFAULT PARAMETERS ######
      if(index(atmp,'--qvSZPs').ne.0) then 
         qvSZPs = .true.
         orcainp%dfa = "wB97M-D4"

         orcainp%d4_s6 = 1.0000_wp
         orcainp%d4_s8 = 0.4769_wp
         orcainp%d4_a1 = 0.2074_wp
         orcainp%d4_a2 = 5.2000_wp
         orcainp%d4_s9 = 1.0000_wp

         orcainp%outn = "wb97m-3c.inp"
      endif
      !##################################################
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
         read(atmp,*) orcainp%coremem
      endif
      if(index(atmp,'--defgrid').ne.0) then
         call get_command_argument(i+1,atmp)
         read(atmp,*) orcainp%defgrid
      endif
      if(index(atmp,'--guess').ne.0) then
         call get_command_argument(i+1,atmp)
         orcainp%guess=trim(adjustl(atmp))
      endif
      if(index(atmp,'--dfa').ne.0) then
         call get_command_argument(i+1,atmp)
         orcainp%dfa=trim(adjustl(atmp))
      endif
      if(index(atmp,'--cm').ne.0) then
         call get_command_argument(i+1,atmp)
         cm=trim(adjustl(atmp))
      endif
      if(index(atmp,'--conv').ne.0) then
         call get_command_argument(i+1,atmp)
         orcainp%scfconv=trim(adjustl(atmp))
      endif
      if(index(atmp,'--d4par').ne.0) then
         indd4param=.true.
         call get_command_argument(i+1,atmp)
         read(atmp,*) orcainp%d4_s6
         call get_command_argument(i+2,atmp)
         read(atmp,*) orcainp%d4_s8
         call get_command_argument(i+3,atmp)
         read(atmp,*) orcainp%d4_a1
         call get_command_argument(i+4,atmp)
         read(atmp,*) orcainp%d4_a2
         call get_command_argument(i+5,atmp)
         read(atmp,*) orcainp%d4_s9
      endif
      if(index(atmp,'--efield').ne.0) then
         call get_command_argument(i+1,atmp)
         read(atmp,*) orcainp%efield(1)
         call get_command_argument(i+2,atmp)
         read(atmp,*) orcainp%efield(2)
         call get_command_argument(i+3,atmp)
         read(atmp,*) orcainp%efield(3)
      endif
      if(index(atmp,'--polar').ne.0) orcainp%polar=.true.
      if(index(atmp,'--polgrad').ne.0) orcainp%polgrad=.true.
      if(index(atmp,'--dipgrad').ne.0) orcainp%dipgrad=.true.
      if(index(atmp,'--geoopt').ne.0) orcainp%geoopt=.true.
      if(index(atmp,'--nocosx').ne.0) orcainp%nocosx=.true.
      if(index(atmp,'--tightscf').ne.0)   orcainp%scfconv='TightSCF'
      if(index(atmp,'--strongscf').ne.0)  orcainp%scfconv='StrongSCF'
      if(index(atmp,'--nouseshark').ne.0) orcainp%nouseshark=.true.
      if(index(atmp,'--noelprop').ne.0) orcainp%noelprop=.true.
      if(index(atmp,'--suggestedguess').ne.0) orcainp%guess='hueckel'
      !> HF reference calculation
      if(index(atmp,'--hfref').ne.0) then
         orcainp%hfref=.true.
         orcainp%outn = "hf_q-vSZP.inp"
      endif
      !> ORCA input file name
      if(index(atmp,'--outname').ne.0) then
         call get_command_argument(i+1,atmp)
         orcainp%outn= trim(adjustl(atmp)) // ".inp"
      endif
      !> Verbosity
      if(index(atmp,'--ov').ne.0 .or. &
      & index(atmp,'--orcaverbose').ne.0) orcainp%verbose=.true.
      if(index(atmp,'--v').ne.0 .or. &
      & index(atmp,'--verbose').ne.0) verbose=.true.
      !> Help
      if(index(atmp,'--help').ne.0) help=.true.
      !> Version
      if(index(atmp,'--version').ne.0) prversion=.true.
   enddo

   !> Print version
   if (prversion) then
      write(*,'(a,a)') "qvSZP version ", version
      stop
   endif

   !> Check input consistency
   if (orcainp%noelprop .and. orcainp%polar) then
      print '(a)', "Error: --noelprop and --polar cannot be used together"
      error stop
   endif

   !> Print header - actual program starts here
   write(*,'(a)')       "     -------------------------------------"
   write(*,'(a)')       "     |    q-vSZP ORCA INPUT GENERATOR    |"
   write(*,'(a,a,a)')   "     |               v",version,"                |"
   write(*,'(a)')       "     |        M. Müller, S. Grimme       |"
   write(*,'(a,/)')     "     -------------------------------------"
   if (help) then
      call helpf()
      stop
   endif

   if (.not. uhfgiven) then
      inquire(file='.UHF',exist=da)
      if(da)then
         open(newunit=myunit,file='.UHF')
         read(myunit,*) nopen
         uhfgiven=.true.
         write(*,'(a,i0,/)') "Number of unpaired electrons read from .UHF file: ", nopen
         close(myunit)
      endif
   endif

   if (.not. indcharge) then
      inquire(file='.CHRG',exist=da)
      if(da)then
         open(newunit=myunit,file='.CHRG')
         read(myunit,*) charge
         close(myunit)
         write(*,'(a,i0,/)') "Charge read from .CHRG file: ", charge
      endif
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
      call search_ghost_atoms(filen,orcainp%ghostatoms,error)
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

   mol%uhf   = nopen
   allocate(tmpids(mol%nat))
   allocate(distvec(mol%nat*(mol%nat+1)/2))
   allocate(cn(mol%nat),q(mol%nat))
   allocate(sccoeff(mol%nat,20,20))

   if (dummy) then
      allocate(cn_short(molshort%nat),q_short(molshort%nat) &
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
   cn    =  0.0_wp
   q     =  0.0_wp
   if (dummy) then
      distvec_short  =  0.0_wp
      cn_short    =  0.0_wp
      q_short     =  0.0_wp
   endif

   call calcrab(mol,distvec)
   call ncoord_basq(mol,distvec,-4.00_wp,cn)
   if (dummy) then
      call calcrab(molshort,distvec_short)
      call ncoord_basq(molshort,distvec_short,-4.00_wp,cn_short)
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
   mol%charge = real(charge,wp)

   write(*,'(a,a,a)') "---- CHARGE MODEL: ", cm," ----"
   select case (cm)
    case default
      call fatal_error(error, "Unknown charge model.")
    case('eeq')
      call eeq(mol,distvec,real(charge,wp),cn,.False., &
      & unity,gamscal,chiscal,alphascal,q,orcainp%efield)
      if (dummy) then
         if (charge /= 0) then
            error stop "Non-zero charge not supported for dummy atoms."
         endif
         call eeq(molshort,distvec_short,real(charge,wp),cn_short,.False., &
         & unity,gamscal,chiscal,alphascal,q_short,orcainp%efield)
      endif
    case('extq')
      if (index((filen),'coord').eq.0) then
         call write_structure(mol, 'coord', error)
         if (allocated(error)) then
            print '(a)', error%message
            error stop
         end if
      endif
      ! Call gp3 per system call to calculate the external charges and store them in the file
      ! 'ceh.charges'. The file is then read in the next step.
      if (sum(abs(orcainp%efield)) > 0.0_wp) then
         write(extcall,'(f12.8,1x,f12.8,1x,f12.8)') orcainp%efield(1), orcainp%efield(2), orcainp%efield(3)
         extcall = "gp3 -ceh -efield " // trim(adjustl(extcall)) // " > gp3.out 2> gp3.err"
         write(*,'(a,f12.8,f12.8,f12.8)') "External electric field is initialized: ", &
         & orcainp%efield(1), orcainp%efield(2), orcainp%efield(3)
         call execute_command_line(extcall, exitstat=errint)
      else
         call execute_command_line('gp3 -ceh > gp3.out 2> gp3.err', exitstat=errint)
      endif
      if (errint /= 0) then
         call fatal_error(error, "Error in GP3 call! ")
      endif
      call extcharges(mol,'ceh.charges',q,cn)
      if (dummy) then
         ! Copy coord file to coord.bak file
         call execute_command_line('rm ceh.charges')
         call execute_command_line('cp coord coord.bak', exitstat=errint)
         call write_structure(molshort, 'coord', error)
         if (allocated(error)) then
            print '(a)', error%message
            error stop
         end if
         call execute_command_line('gp3 -ceh > gp3.out 2> gp3.err', exitstat=errint)
         if (errint /= 0) then
            call fatal_error(error, "Error in GP3 call! ")
         endif
         call extcharges(molshort,'ceh.charges',q_short,cn_short)
      endif
    case('ceh')
      if (verbose) verbosity = 2
      call ceh(mol,orcainp%efield,q,error,verbosity)
      if (dummy) then
         call ceh(molshort,orcainp%efield,q_short,error,verbosity)
      endif
   end select
   if(verbose) write(*,'(a)') "---------------------------"
   if (allocated(error)) then
      print '(a)', error%message
      error stop
   end if

   if (indbfile) then
      call rdbas(bas,verbose,trim(bfilen))
   else
      call rdbas(bas,verbose)
   endif

   if (.not. qvSZPs) then
      if (indefile) then
         call rdecp_default(ecp,verbose,efilen)
      else
         call rdecp_default(ecp,verbose)
      endif
   else
      if (.not. indefile) error stop "ECP file location must be given with --efile <filename>."
      call rdecp_qvSZPs(ecp,verbose,efilen,increment)
      do i = 1, mol%nat
         if ( .not. sum(abs(ecp%exp(mol%num(mol%id(i)),:,:))) > 0.0_wp ) then
            error stop "ERROR. No ECP/ACP found for atom " // trim(adjustl(mol%sym(mol%id(i)))) // "."
         endif
      enddo
   endif

   if (verbose) then
      if (dummy) write(*,'(a)') "MOLECULE WITH GHOST ATOMS"
      write(*,'(a)') "----  ------  -------"
      write(*,'(a)') "Atom  Symbol  Charge"
      write(*,'(a)') "----  ------  -------"
      do i = 1, mol%nat
         write(*,'(i4,2x,a,2x,f12.8)') i, mol%sym(mol%id(i)), q(i)
      enddo
      write(*,'(a)') "----  ------ --------------------"
      write(*,'(a)') "Atom  Symbol  Coordination number"
      write(*,'(a)') "----  ------  -------------------"
      do i = 1, mol%nat
         write(*,'(i4,2x,a,2x,f12.8)') i, mol%sym(mol%id(i)), cn(i)
      enddo
      if (dummy) then
         write(*,'(a)') "----------------------------"
         write(*,'(a)') "MOLECULE WITHOUT GHOST ATOMS"
         write(*,'(a)') "----  ------  -------"
         write(*,'(a)') "Atom  Symbol  Charge"
         write(*,'(a)') "----  ------  -------"
         do i = 1, molshort%nat
            write(*,'(i4,2x,a,2x,f12.8)') i, molshort%sym(molshort%id(i)), q_short(i)
         enddo
         write(*,'(a)') "----  ------ --------------------"
         write(*,'(a)') "Atom  Symbol  Coordination number"
         write(*,'(a)') "----  ------  -------------------"
         do i = 1, molshort%nat
            write(*,'(i4,2x,a,2x,f12.8)') i, molshort%sym(molshort%id(i)), cn_short(i)
         enddo
      endif
   endif

   if (dummy) then
      call wrorca(orcainp, mol, bas, ecp, q, cn, verbose, &
      & error, q_short, cn_short)
   else
      call wrorca(orcainp, mol, bas, ecp, q, cn, verbose, &
      & error)
   endif

   if (qvSZPs) then
      do i = 1, mol%nat
         energy_shift = energy_shift + increment(mol%num(mol%id(i)))
      enddo
      write(*,'(/,a)') "CAUTION: Sum of q-vSZPs atomic increments has to be added!"
      write(*,'(a)') " -----  THIS IS NOT DONE AUTOMATICALLY!  -----"
      open(newunit=myunit,file='eshift.qvSZPs', status='replace', action='write')
      write(myunit,'(f16.10)') energy_shift
      close(myunit)
      write(*,'(a)') " -----  written to file 'eshift.qvSZPs'  -----"
      write(*,'(a,f14.9)') "TOTAL ENERGY SHIFT: ", energy_shift
   endif


end program main

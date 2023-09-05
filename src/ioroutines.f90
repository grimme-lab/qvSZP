module ioroutines
   use mctc_io
   use mctc_env
   use stdlib_sorting
   implicit none
   private

   public :: rdfile,rdbas,rdecp,check_ghost_atoms,search_ghost_atoms,basis_type,ecp_type

   type :: basis_type
      integer               :: atmax ! atom index
      integer, allocatable  :: npr(:,:) ! number of primitive functions per atom
      integer, allocatable  :: nbf(:) ! total number of primitive functions
      integer, allocatable  :: angmom(:,:) ! angular momentum of each primitive function
      real(wp), allocatable :: exp(:,:,:) ! exponent of each primitive function
      real(wp), allocatable :: coeff(:,:,:) ! contraction coefficient of each primitive function
      logical, allocatable  :: sccoeff(:)
      real(wp), allocatable :: qcoeff(:,:,:)
      real(wp), allocatable :: scalparam(:,:)   ! scaling parameters for the basis functions

   end type basis_type

   type :: ecp_type
      integer               :: atmax ! atom index
      integer               :: atmin ! atom index
      integer, allocatable  :: npr(:,:) ! number of primitive functions per atom
      integer, allocatable  :: nbf(:) ! total number of primitive functions
      integer, allocatable  :: ncore(:) ! total number of primitive functions
      integer, allocatable  :: lmax(:) ! total number of primitive functions
      integer, allocatable  :: nfactor(:,:,:) ! total number of primitive functions
      integer, allocatable  :: angmom(:,:) ! angular momentum of each primitive function
      integer(int_size), allocatable :: sindex(:,:)
      real(wp), allocatable :: exp(:,:,:) ! exponent of each primitive function
      real(wp), allocatable :: coeff(:,:,:) ! contraction coefficient of each primitive function
   end type ecp_type

contains
   subroutine rdfile(filename,imol,nel)

      character(len=*), intent(in)      :: filename
      type(structure_type),intent(out)  :: imol
      integer,optional,intent(out)      :: nel
      type(error_type), allocatable     :: error
      integer                           :: i

      call read_structure(imol, filename, error)
      if (allocated(error)) then
         print '(a)', error%message
         error stop
      end if

      if (present(nel)) then
         nel = 0
         do i=1,imol%nat
            nel = nel + imol%num(imol%id(i))
         end do
      endif

      imol%sym = to_symbol(imol%num)

   end subroutine rdfile

   subroutine rdbas(ibasis,verb,basisfilename)

      type(basis_type),intent(out)           :: ibasis
      logical, intent(in), optional          :: verb
      character(len=*), intent(in),optional  :: basisfilename
      type(error_type), allocatable          :: error

      character(len=:), allocatable          :: fname,homedir
      character(len=120)                     :: atmp
      character(len=1)                       :: ltmp

      logical                                :: da,checked

      integer, allocatable                   :: nbf(:),npr(:,:)
      character(len=1), allocatable          :: angmom(:,:)

      integer                                :: myunit,char_length,ierr,iread
      integer                                :: iat,tmpnpr,imax
      integer                                :: i,j,l,k

      real(wp)                               :: tmpvec(4) = 1000.00_wp
      real(wp)                               :: iattmp(3)

      allocate(nbf(100),npr(100,20),angmom(100,20),ibasis%sccoeff(118),ibasis%scalparam(118,3))
      ibasis%sccoeff = .false.
      ibasis%scalparam = 0.0_wp

      if (.not. present(basisfilename)) then
         call get_environment_variable("HOME", length=char_length)
         allocate(character(len=char_length) :: homedir)
         CALL get_environment_variable("HOME", value=homedir, status=ierr)
         select case (ierr)
          case (0)
            ! do nothing
          case (1)
            call fatal_error(error, "HOME environment variable is not set!")
            if (allocated(error)) then
               print '(a)', error%message
               error stop "I/O error stop."
            end if
          case (2)
            call fatal_error(error, "this compiler does not support environment variables")
            if (allocated(error)) then
               print '(a)', error%message
               error stop "I/O error stop."
            end if
         end select
         fname=trim(homedir)//'/.basis_vDZP'
      else
         fname=basisfilename
      end if

      write(*,'(a,a)') "Used basis set file: ",fname

      inquire(file=fname,exist=da)
      if(da)then
         open(newunit=myunit,file=fname,action='read',status='old',iostat=ierr)
      else
         call fatal_error(error, "basis set file cannot be opened!")
         if (allocated(error)) then
            print '(a)', error%message
            error stop "I/O error stop."
         end if
      endif
      iread =  0
      l     =  0
      npr   =  0
      nbf   =  0
      imax  =  0
      do while (iread >= 0)
         read(myunit,'(a)',iostat=iread) atmp
         l = l + 1
         if (iread < 0) then
            exit
         end if
         if (iread > 0) then
            write(*,*) "Current line number: ",l
            error stop "I/O error in basis set read in."
         end if
         if(index(atmp,'*').ne.0) then
            cycle
         else
            iattmp=1000.00_wp
            atmp=trim(atmp)
            read(atmp,*,iostat=iread) iat, iattmp(:)
            if (iread /= 0) then
               write(*,*) "Current line number: ",l
               error stop "I/O error in reading element-wise coefficients from line."
            end if
            if (maxval(iattmp) <= 999.0_wp) then
               ibasis%scalparam(iat,1) = iattmp(1)
               ibasis%scalparam(iat,2) = iattmp(2)
               ibasis%scalparam(iat,3) = iattmp(3)
            end if
            checked = .false.
            tmpvec = 1000.00_wp
            ! ATTENTION: IOSTAT error only detects read-in errors, not end of file as this is is
            ! what we want to test for here.
            if (iread > 0) then
               write(*,*) "Current line number: ",l
               error stop "I/O error in basis set read in."
            end if
            nbf(iat) = 0
            npr(iat,:) = 0


            if (iat > imax) then
               imax = iat
            end if

            read(myunit,'(a)',iostat=iread) atmp
            l = l + 1
            if (iread /= 0) then
               write(*,*) "Current line number: ",l
               error stop "I/O error in basis set read in."
            end if
            k = 0
            do while (index(atmp,'*').eq.0)
               read(atmp,*) tmpnpr,ltmp
               k = k + 1
               nbf(iat) = k
               npr(iat,nbf(iat)) = tmpnpr
               angmom(iat,nbf(iat)) = ltmp
               do j=1,tmpnpr
                  read(myunit,'(a)',iostat=iread) atmp
                  l = l + 1
                  if (.not. checked) then
                     read(atmp,*,iostat=iread) tmpvec(:)
                     if (tmpvec(3) < 999.0_wp) then
                        ibasis%sccoeff(iat) = .true.
                     end if
                     checked = .true.
                  endif
               end do
               read(myunit,'(a)',iostat=iread) atmp
               l = l + 1
               if (iread /= 0) then
                  write(*,*) "Current line number: ",l
                  error stop "I/O error in basis set read in."
               end if
            enddo
         endif
      enddo

      ibasis%atmax = imax
      allocate(ibasis%npr(imax,20),ibasis%angmom(imax,20),ibasis%nbf(imax))
      allocate(ibasis%exp(imax,20,20),ibasis%coeff(imax,20,20),ibasis%qcoeff(imax,20,20))
      close(myunit)

      ibasis%npr = 0
      ibasis%nbf = 0
      ibasis%angmom = 0
      ibasis%exp = 0.0_wp
      ibasis%coeff = 0.0_wp
      ibasis%qcoeff = 0.0_wp

      do i=1,ibasis%atmax
         ibasis%nbf(i) = nbf(i)
         do j=1,ibasis%nbf(i)
            ibasis%npr(i,j) = npr(i,j)
            if (angmom(i,j) == "s") ibasis%angmom(i,j) = 0
            if (angmom(i,j) == "p") ibasis%angmom(i,j) = 1
            if (angmom(i,j) == "d") ibasis%angmom(i,j) = 2
            if (angmom(i,j) == "f") ibasis%angmom(i,j) = 3
            if (angmom(i,j) == "g") ibasis%angmom(i,j) = 4
            if (angmom(i,j) == "h") ibasis%angmom(i,j) = 5
            if (angmom(i,j) == "i") ibasis%angmom(i,j) = 6
            if (angmom(i,j) == "j") error stop "too high primitives"
            ibasis%npr(iat,:) = npr(iat,:)
         enddo
      enddo


      open(newunit=myunit,file=fname,action='read',status='old',iostat=ierr)
      iread =  0
      l     =  0
      do while (iread >= 0)
         read(myunit,'(a)',iostat=iread) atmp
         l = l + 1
         if (iread < 0) then
            write(*,'(a,1x,i3,/)') "End of file reached after reading basis set for element:", iat
            exit
         end if
         if (iread > 0) then
            write(*,*) "Current line number: ",l
            error stop "I/O error in basis set read in."
         end if
         if(index(atmp,'*').ne.0) then
            cycle
         else
            read(atmp,*,iostat=iread) iat
            if (iread /= 0) then
               write(*,*) "Current line number: ",l
               error stop "I/O error in basis set read in."
            end if
            do i=1,ibasis%nbf(iat)
               read(myunit,*,iostat=iread) atmp
               l = l + 1
               do j=1,ibasis%npr(iat,i)
                  if (ibasis%sccoeff(iat)) then
                     read(myunit,*,iostat=iread) ibasis%exp(iat,i,j),ibasis%coeff(iat,i,j), &
                        ibasis%qcoeff(iat,i,j)
                  else
                     read(myunit,*,iostat=iread) ibasis%exp(iat,i,j),ibasis%coeff(iat,i,j)
                  end if
                  l = l + 1
               end do
            enddo
         endif
      enddo

   end subroutine rdbas

   subroutine rdecp(iecp,verb,ecpfilename)

      type(ecp_type),intent(out)           :: iecp
      logical, intent(in), optional          :: verb
      character(len=*), intent(in),optional  :: ecpfilename
      type(error_type), allocatable          :: error

      character(len=:), allocatable          :: fname,homedir
      character(len=120)                     :: atmp,btmp,ctmp,dtmp,etmp

      logical                                :: da,angfound

      integer, allocatable                   :: nbf(:),npr(:,:),ncore(:),lmax(:)
      integer, allocatable                   :: angmom(:,:)

      integer                                :: myunit,char_length,ierr,iread
      integer                                :: iat,imax,ltmp
      integer                                :: i,j,l,imin

      integer(int_size), allocatable        :: sortindex(:)

      allocate(nbf(118),npr(118,20),angmom(118,20),ncore(118),lmax(118))

      if (.not. present(ecpfilename)) then
         call get_environment_variable("HOME", length=char_length)
         allocate(character(len=char_length) :: homedir)
         CALL get_environment_variable("HOME", value=homedir, status=ierr)
         select case (ierr)
          case (0)
            ! do nothing
          case (1)
            call fatal_error(error, "HOME environment variable is not set!")
            if (allocated(error)) then
               print '(a)', error%message
               error stop "I/O error stop."
            end if
          case (2)
            call fatal_error(error, "this compiler does not support environment variables")
            if (allocated(error)) then
               print '(a)', error%message
               error stop "I/O error stop."
            end if
         end select
         fname=trim(homedir)//'/.ecp_vDZP'
      else
         fname=ecpfilename
      end if

      write(*,'(a,a)') "Used ECP set file: ",fname

      inquire(file=fname,exist=da)
      if(da)then
         open(newunit=myunit,file=fname,action='read',status='old',iostat=ierr)
      else
         call fatal_error(error, "ECP set file cannot be opened!")
         if (allocated(error)) then
            print '(a)', error%message
            error stop "I/O error stop."
         end if
      endif
      iread =  0
      l     =  0
      npr   =  0
      nbf   =  0
      imax  =  0
      angmom=  0
      imin  = 118
      do while (iread >= 0)
         read(myunit,'(a)',iostat=iread) atmp
         l = l + 1
         if (iread < 0) then
            exit
         end if
         if (iread > 0) then
            write(*,*) "Current line number: ",l
            error stop "I/O error in basis set read in. Pt. 1."
         end if
         if( (index(atmp,'*').ne.0) .or. (index(atmp,'#').ne.0) ) then
            cycle
         else
            read(atmp,*,iostat=iread) iat
            if (iread /= 0) then
               write(*,*) "Current line number: ",l
               error stop "I/O error in basis set read in. Pt. 2."
            end if
            nbf(iat) = 0
            npr(iat,:) = 0

            if (iat > imax) then
               imax = iat
            end if

            if (iat < imin) then
               imin = iat
            end if

            read(myunit,*,iostat=iread) btmp, ctmp, ncore(iat), dtmp, etmp, lmax(iat)
            l = l + 1
            if (iread /= 0) then
               write(*,*) "Current line number: ",l
               error stop "I/O error in basis set read in. Pt. 3."
            end if

            angfound = .false.
            do while ( (index(atmp,'*').eq.0) .and. (index(atmp,'#').eq.0) )
               read(myunit,'(a)',iostat=iread) atmp
               l = l + 1
               if (iread /= 0) then
                  write(*,*) "Current line number: ",l
                  error stop "I/O error in basis set read in. Pt. 4."
               end if

               if(index(atmp(1:1),'s').eq.1) then
                  ltmp=0
                  angfound=.true.
               else if(index(atmp(1:1),'p').eq.1) then
                  ltmp=1
                  angfound=.true.
               else if(index(atmp(1:1),'d').eq.1) then
                  ltmp=2
                  angfound=.true.
               else if(index(atmp(1:1),'f').eq.1) then
                  ltmp=3
                  angfound=.true.
               else if(index(atmp(1:1),'g').eq.1) then
                  ltmp=4
                  angfound=.true.
               else if(index(atmp(1:1),'h').eq.1) then
                  ltmp=5
                  angfound=.true.
               else if(index(atmp,'*').eq.1 .or. index(atmp,'#').eq.1) then
                  exit
               endif

               if (angfound) then
                  nbf(iat) = nbf(iat) + 1
                  angmom(iat,nbf(iat)) = ltmp
                  npr(iat,nbf(iat)) = 0
                  angfound = .false.
                  cycle
               end if

               npr(iat,nbf(iat)) = npr(iat,nbf(iat)) + 1

            enddo
         endif
      enddo

      iecp%atmax = imax
      allocate(iecp%npr(imax,20),iecp%angmom(imax,20),iecp%nbf(imax),iecp%ncore(imax),iecp%lmax(imax))
      allocate(iecp%exp(imax,20,20),iecp%coeff(imax,20,20),iecp%nfactor(imax,20,20))
      allocate(iecp%sindex(imax,20))
      close(myunit)

      iecp%npr = 0
      iecp%nbf = 0
      iecp%angmom = 0
      iecp%exp = 0.0_wp
      iecp%coeff = 0.0_wp
      iecp%ncore = 0
      iecp%lmax = 0
      iecp%sindex = 0

      do i=1,iecp%atmax
         iecp%nbf(i) = nbf(i)
         do j=1,iecp%nbf(i)
            iecp%npr(i,j) = npr(i,j)
            iecp%angmom(i,j) = angmom(i,j)
            iecp%npr(iat,:) = npr(iat,:)
            iecp%ncore(i) = ncore(i)
            iecp%lmax(i) = lmax(i)
         enddo
      enddo
      iecp%atmin = imin

      open(newunit=myunit,file=fname,action='read',status='old',iostat=ierr)
      iread =  0
      l     =  0
      do while (iread >= 0)
         read(myunit,'(a)',iostat=iread) atmp
         l = l + 1
         if (iread < 0) then
            write(*,'(a,1x,i3,/)') "End of file reached after reading ECP set for element:", iat
            exit
         end if
         if (iread > 0) then
            write(*,*) "Current line number: ",l
            error stop "I/O error in basis set read in."
         end if
         if( (index(atmp,'*').ne.0) .or. (index(atmp,'#').ne.0) ) then
            cycle
         else
            read(atmp,*,iostat=iread) iat
            read(myunit,*,iostat=iread) btmp
            l = l + 1
            if (iread /= 0) then
               write(*,*) "Current line number: ",l
               error stop "I/O error in basis set read in."
            end if
            do i=1,iecp%nbf(iat)
               read(myunit,*,iostat=iread) btmp
               l = l + 1
               if (iread /= 0) then
                  write(*,*) "Current line number: ",l
                  error stop "I/O error in basis set read in."
               end if
               do j=1,iecp%npr(iat,i)
                  read(myunit,*,iostat=iread) iecp%coeff(iat,i,j), iecp%nfactor(iat,i,j), iecp%exp(iat,i,j)
                  l = l + 1
                  if (iread /= 0) then
                     write(*,*) "Current line number: ",l
                     error stop "I/O error in basis set read in."
                  end if
               end do
            enddo
         endif
      enddo

      do i=iecp%atmin,iecp%atmax
         allocate(sortindex(iecp%nbf(i)))
         call sort_index(iecp%angmom(i,1:iecp%nbf(i)),sortindex(1:iecp%nbf(i)))
         iecp%sindex(i,1:iecp%nbf(i)) = sortindex(1:iecp%nbf(i))
         deallocate(sortindex)
      enddo

      close(myunit)

   end subroutine rdecp

   subroutine check_ghost_atoms(filename,ghost,error)
      character(len=*), intent(in)      :: filename
      logical, intent(out)              :: ghost
      type(error_type), allocatable     :: error
      character(len=120)                :: atmp
      character(len=3)                  :: delete_char = " : "
      integer                           :: ierr
      integer                           :: myunit

      ghost = .false.
      open(newunit=myunit,file=filename,action='read',status='old',iostat=ierr)
      ! check for the existence of the character delete_char in the file "filename"
      if (ierr /= 0) then
         call fatal_error(error, "coordinate file cannot be opened!")
         if (allocated(error)) then
            print '(a)', error%message
            error stop "I/O error stop."
         end if
      end if
      do while (ierr == 0)
         read(myunit,'(a)',iostat=ierr) atmp
         if (ierr /= 0) then
            exit
         end if
         if (index(atmp,delete_char).ne.0) then
            write(*,'(a,/,a)') "Ghost atom(s) detected in coordinate file.", &
               "Writing dummy coordinate files 'struc_raw.xyz' and 'struc_short.xyz'."
            ghost = .true.
            exit
         end if
      enddo
      close(myunit)

   end subroutine check_ghost_atoms

   subroutine search_ghost_atoms(filename,ghostat,error)
      character(len=*), intent(in)      :: filename
      logical, allocatable, intent(out) :: ghostat(:)
      type(error_type), allocatable     :: error
      character(len=120)                :: atmp
      character(len=1)                  :: delete_char = ':'
      integer                           :: i,j,ierr
      integer                           :: myunit,myunit2,myunit3,nat,natshort

      natshort = 0
      open(newunit=myunit,file=filename,action='read',status='old',iostat=ierr)
      open(newunit=myunit2,file="struc_raw.xyz",action='write',status='replace',iostat=ierr)

      i=0
      if (ierr /= 0) then
         call fatal_error(error, "coordinate file cannot be opened!")
         if (allocated(error)) then
            print '(a)', error%message
            error stop "I/O error stop."
         end if
      end if
      read(myunit,'(a)',iostat=ierr) atmp
      if (ierr /= 0) then
         call fatal_error(error, "coordinate file cannot be read!")
         if (allocated(error)) then
            print '(a)', error%message
            error stop "I/O error stop."
         end if
      end if
      read(atmp,*,iostat=ierr) nat
      if (ierr /= 0) then
         call fatal_error(error, "coordinate file cannot be read!")
         if (allocated(error)) then
            print '(a)', error%message
            error stop "I/O error stop."
         end if
      end if
      write(myunit2,'(i0)',iostat=ierr) nat
      read(myunit,'(a)',iostat=ierr) atmp
      write(myunit2,'(a)',iostat=ierr) trim(adjustl(atmp))
      allocate(ghostat(nat))
      ghostat = .false.
      do while (ierr == 0)
         read(myunit,'(a)',iostat=ierr) atmp
         if (ierr /= 0) then
            exit
         end if
         i = i + 1
         if (index(atmp,delete_char).ne.0) then
            ghostat(i) = .true.
            natshort = natshort + 1
            atmp = trim(adjustl(atmp))
            ! remove "delete_char" from string
            do j = 1, len(atmp)
               if (atmp(j:j) == delete_char) then
                  atmp = atmp(1:j-1) // " " // atmp(j+1:len(atmp))
               end if
            end do
            write(myunit2,'(a)',iostat=ierr) trim(atmp)
            cycle
         else
            write(myunit2,'(a)',iostat=ierr) trim(adjustl(atmp))
         end if
      end do

      close(myunit2)

      natshort = nat - natshort
      rewind(myunit)
      open(newunit=myunit3,file="struc_short.xyz",action='write',status='replace',iostat=ierr)

      i=0
      if (ierr /= 0) then
         call fatal_error(error, "coordinate file cannot be opened!")
         if (allocated(error)) then
            print '(a)', error%message
            error stop "I/O error stop."
         end if
      end if
      read(myunit,'(a)',iostat=ierr) atmp
      if (ierr /= 0) then
         call fatal_error(error, "coordinate file cannot be read!")
         if (allocated(error)) then
            print '(a)', error%message
            error stop "I/O error stop."
         end if
      end if
      read(atmp,*,iostat=ierr) nat
      if (ierr /= 0) then
         call fatal_error(error, "coordinate file cannot be read!")
         if (allocated(error)) then
            print '(a)', error%message
            error stop "I/O error stop."
         end if
      end if
      write(myunit3,'(i0)',iostat=ierr) natshort
      read(myunit,'(a)',iostat=ierr) atmp
      write(myunit3,'(a)',iostat=ierr) trim(adjustl(atmp))
      do i = 1, nat
         read(myunit,'(a)',iostat=ierr) atmp
         if (ierr /= 0) then
            exit
         end if
         if (.not. index(atmp,delete_char).ne.0) then
            write(myunit3,'(a)',iostat=ierr) trim(adjustl(atmp))
         end if
      end do

      close(myunit)
      close(myunit3)

   end subroutine search_ghost_atoms

end module ioroutines

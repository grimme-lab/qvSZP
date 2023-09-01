module test_qvSZP
   use chargscfcts, only : eeq
   use stdlib_io, only : getline
   use testdrive, only : error_type, unittest_type, new_unittest, check
   use mstore, only : get_structure
   use mctc_io, only : structure_type, new
   use mctc_env, only : wp
   implicit none
   private

   public :: collect_charges

contains

   !> Collect all exported unit tests
   subroutine collect_charges(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [new_unittest("eeq-charges", test_eeq_neutral)]
   end subroutine collect_charges

   !> Check substitution of a single line
   subroutine test_eeq_neutral(error)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      character(len=:), allocatable :: line
      real(wp),allocatable             :: distvec(:),cnvec(:),qeeq(:)
      real(wp) :: efield(3) = 0.0_wp
      integer :: i
      ! ####### Set up EEQ parameters for coefficient scaling #######
      real(wp),parameter               :: unity(86)      = 1.0_wp
      real(wp),parameter               :: alphascal(86)  = 0.75_wp
      real(wp),parameter               :: chiscal(86)    = 1.20_wp
      real(wp),parameter               :: gamscal(86)    = 0.90_wp
      ! ###########################################################
      real(wp), parameter :: ref(16) = reshape([ &
      0.30785024_wp, &
     -0.02792822_wp, &
     -0.11358521_wp, &
      0.03271103_wp, &
      0.10217515_wp, &
      0.25825112_wp, &
     -0.20769464_wp, &
     -0.03447769_wp, &
     -0.25817610_wp, &
      0.24460001_wp, &
      0.26974150_wp, &
     -0.01050015_wp, &
     -0.28007123_wp, &
     -0.42862573_wp, &
      0.03898681_wp, &
      0.10674311_wp], shape(ref))

      allocate(distvec(mol%nat*(mol%nat+1)/2))
      allocate(cnvec(mol%nat),qeeq(mol%nat))
      distvec  =  0.0_wp
      cnvec    =  0.0_wp
      qeeq     =  0.0_wp
      call get_structure(mol, "MB16-43", "12")
      call calcrab(mol,distvec)
      call ncoord_basq(mol,distvec,-3.75_wp,cnvec)

      call eeq(mol,distvec,mol%charge,cnvec,.False., &
      & unity,gamscal,chiscal,alphascal,qeeq,efield)
      ! line = "This is a valid example"
      ! call check(error, line, "This is a valid example")
      do i = 1, mol%nat
         call check(error, qeeq(i), ref(i), thr=1.0e-6_wp)
      end do

   end subroutine test_eeq_neutral
end module test_qvSZP

program tester
   use, intrinsic :: iso_fortran_env, only : error_unit
   use testdrive, only : run_testsuite
   use test_qvSZP, only : collect_charges
   implicit none
   integer :: stat

   stat = 0
   call run_testsuite(collect_charges, error_unit, stat)

   if (stat > 0) then
      write(error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
      error stop
   end if

end program tester

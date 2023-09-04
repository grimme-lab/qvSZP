module test_qvSZP
   use chargscfcts, only : eeq
   use stdlib_io, only : getline
   use testdrive, only : error_type, unittest_type, new_unittest, check
   use mstore, only : get_structure
   use mctc_io, only : structure_type, new
   use mctc_env, only : wp
   use chargscfcts, only : calcrab, ncoord_basq, eeq
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
       0.14519974_wp, &
      -0.13570696_wp, &
      -0.13513992_wp, &
      -0.21127141_wp, &
       0.16304431_wp, &
       0.22455622_wp, &
      -0.25454477_wp, &
       0.05184459_wp, &
      -0.29019112_wp, &
       0.38616228_wp, &
       0.34946124_wp, &
       0.09644326_wp, &
      -0.26588399_wp, &
      -0.40907177_wp, &
       0.12349849_wp, &
       0.16159982_wp], shape(ref))

      call get_structure(mol, "MB16-43", "12")
      allocate(distvec(mol%nat*(mol%nat+1)/2))
      allocate(cnvec(mol%nat),qeeq(mol%nat))
      distvec  =  0.0_wp
      cnvec    =  0.0_wp
      qeeq     =  0.0_wp
      call calcrab(mol,distvec)
      call ncoord_basq(mol,distvec,-3.75_wp,cnvec)

      call eeq(mol,distvec,mol%charge,cnvec,.False., &
      & unity,gamscal,chiscal,alphascal,qeeq,efield)
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

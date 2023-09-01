module test_qvSZP
   use chargscfcts, only : eeq
   use stdlib_io, only : getline
   use testdrive, only : error_type, unittest_type, new_unittest, check
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
     character(len=:), allocatable :: line
     line = "This is a valid example"
     call check(error, line, "This is a valid example")
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
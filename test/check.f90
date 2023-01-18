program check
   use mctc_io
   use mctc_env
   use ioroutines
   implicit none
   type(structure_type)             :: mol
   type(error_type), allocatable    :: error
   integer                          :: nel


   call rdfile("test/aucl4.xyz",mol,nel)
   if (allocated(error)) then
      print '(a)', error%message
      error stop
   end if

   if (mol%nat .ne. 5) then
    call fatal_error(error,"wrong number of atoms")
   endif
   if (allocated(error)) then
      print '(a)', error%message
      error stop
   end if

end program check

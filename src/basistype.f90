module basistype
   use mctc_io
   use mctc_env
   use stdlib_sorting, only: int_size
   implicit none
   private

   public :: basis_type, ecp_type

   type :: basis_type
      integer               :: atmax ! atom index
      integer, allocatable  :: npr(:,:) ! number of primitive functions per atom
      integer, allocatable  :: nbf(:) ! total number of primitive functions
      integer, allocatable  :: angmom(:,:) ! angular momentum of each primitive function
      real(wp), allocatable :: exp(:,:,:) ! exponent of each primitive function
      real(wp), allocatable :: coeff(:,:,:) ! contraction coefficient of each primitive function
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

end module basistype

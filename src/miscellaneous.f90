module miscellaneous
  implicit none
  private

  public :: helpf
contains
  subroutine helpf()
   write(*,*) "Possible commands:"
   write(*,*) "--mpi <int>   # set number of MPI processes"
   write(*,*) "--defgrid <int> # set grid size"
   write(*,*) "--guess <guess options> # SCF guess options, see ORCA manual"
   write(*,*) "--polar"
   write(*,*) "--hyppol"
   write(*,*) "--polgrad"
   write(*,*) "--dipgrad"
   write(*,*) "--geoopt"
   write(*,*) "--nocosx"
   write(*,*) "--tightscf    # (tight convergence)"
   write(*,*) "--strongscf   # (strong convergence)"
   write(*,*) "--v           # verbose mode with extended printout of O4wB97X3c and for ORCA itself"
   write(*,*) "--nouseshark  # (use different integral library)"
   write(*,*) "--plot        # (Plot the electron density with the following settings)"
  end subroutine helpf
end module miscellaneous
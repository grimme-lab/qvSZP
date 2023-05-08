module miscellaneous
  implicit none
  private

  public :: helpf
contains
  subroutine helpf()
   write(*,*) "Possible commands:"
   write(*,*) "--struc <filename> # set structure file. DEFAULT 'coord'"
   write(*,*) "--basisfile <filename> # set basis set file"
   write(*,*) "--ecpfile <filename> # set ECP file"
   write(*,*) "--mpi <int>   # set number of MPI processes"
   write(*,*) "--memory <int> # set memory in MB"
   write(*,*) "--chrg <int>  # set charge"
   write(*,*) "--uhf <int>  # set number of unpaired electrons"
   write(*,*) "--defgrid <int> # set grid size"
   write(*,*) "--guess <guess options> # SCF guess options, see ORCA manual"
   write(*,*) "--efield <x> <y> <z> # set electric field vector"
   write(*,*) "--d4par <s6> <s8> <a1> <a2> <s9> # set D4 dispersion parameters"
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
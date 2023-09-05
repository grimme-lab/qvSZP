module miscellaneous
   implicit none
   private

   public :: helpf
contains
   subroutine helpf()
      write(*,'(a)') "Input generator-related commands:"
      write(*,'(a25,10x,a)') "--version",         "# print version"
      write(*,'(a25,10x,a)') "--help",            "# print this help"
      write(*,'(a25,10x,a)') "--v",               "# verbose mode with extended printout of input generator"
      write(*,'(a25,10x,a)') "--struc <filename>","# set structure file. DEFAULT 'coord'"
      write(*,'(a25,10x,a)') "--bfile <filename>","# set basis set file location"
      write(*,'(a25,10x,a)') "--efile <filename>","# set ECP file location"
      write(*,'(a25,10x,a)') "--chrg <int>",      "# set charge"
      write(*,'(a25,10x,a)') "--uhf <int>",       "# set number of unpaired electrons"
      write(*,'(a25,10x,a)') "--efield <x> <y> <z>","# set electric field vector"
      write(*,'(a25,10x,a)') "--outname <filename>","# set output file name. NOTE: An'.inp' suffix will be added."
      write(*,'(a25,10x,a,/)') "--cm <name of desired charge model>","# set charge model. DEFAULT: 'ceh'"
      write(*,'(a)') "ORCA settings-related commands:"
      write(*,'(a25,10x,a)') "--mpi <int>",       "# set number of MPI processes"
      write(*,'(a25,10x,a)') "--memory <int>",    "# set memory in MB"
      write(*,'(a25,10x,a)') "--defgrid <int>",   "# set grid size"
      write(*,'(a25,10x,a)') "--guess <guess options>","# SCF guess options, see ORCA manual"
      write(*,'(a25,10x,a)') "--d4par <s6> <s8> <a1> <a2> <s9>","# set D4 dispersion parameters"
      write(*,'(a25,10x,a)') "--polar",           "# calculate polarizability"
      write(*,'(a25,10x,a)') "--polgrad",         "# calculate polarizability gradient"
      write(*,'(a25,10x,a)') "--dipgrad",         "# calculate dipole moment gradient"
      write(*,'(a25,10x,a)') "--geoopt",          "# geometry optimization"
      write(*,'(a25,10x,a)') "--nocosx",          "# do not use COSX approximation"
      write(*,'(a25,10x,a)') "--tightscf",        "# (tight convergence)"
      write(*,'(a25,10x,a)') "--strongscf",       "# (strong convergence)"
      write(*,'(a25,10x,a)') "--nouseshark",      "# (use different integral library)"
   end subroutine helpf
end module miscellaneous
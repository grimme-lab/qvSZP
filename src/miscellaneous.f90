module miscellaneous
   implicit none
   private

   public :: helpf
contains
   subroutine helpf()
      write(*,'(a)') "Input generator-related commands:"
      write(*,'(a35,10x,a)') "--version",         "# print version"
      write(*,'(a35,10x,a)') "--help",            "# print this help"
      write(*,'(a35,10x,a)') "--v",               "# verbose mode with extended printout of input generator"
      write(*,'(a35,10x,a)') "--struc <filename>","# set structure file. DEFAULT 'coord'"
      write(*,'(a35,10x,a)') "--bfile <filename>","# set basis set file path. DEFAULT '~/.basisq'"
      write(*,'(a35,10x,a)') "--efile <filename>","# set ECP file path. DEFAULT '~/.ecpq'"
      write(*,'(a35,10x,a)') "--chrg <int>",      "# set charge"
      write(*,'(a35,10x,a)') "--uhf <int>",       "# set number of unpaired electrons"
      write(*,'(a35,10x,a)') "--efield <x> <y> <z>","# set electric field vector"
      write(*,'(a35,10x,a)') "--outname <filename>","# set output file name. NOTE: An'.inp' suffix will be added."
      write(*,'(a35,10x,a)') "--cm <charge model>","# set charge model. DEFAULT: 'ceh'"
      write(*,'(a35,10x,a)') "--qvSZPs", &
      & "# q-vSZPs basis set mode for different ECP layout and atomic increments"
      write(*,'(/,a)') "ORCA settings-related commands:"
      write(*,'(a35,10x,a)') "--dfa <density functional>", "# set density functional"
      write(*,'(a35,10x,a)') "--mpi <int>",       "# set number of parallel MPI processes"
      write(*,'(a35,10x,a)') "--memory <int>",    "# set memory per core in MB"
      write(*,'(a35,10x,a)') "--defgrid <int>",   "# set grid size <int, [1,3]>"
      write(*,'(a35,10x,a)') "--guess <guess options>","# SCF guess options, see ORCA manual"
      write(*,'(a35,10x,a)') "--d4par <s6> <s8> <a1> <a2> <s9>","# set D4 dispersion parameters"
      write(*,'(a35,10x,a)') "--polar",           "# calculate polarizability"
      write(*,'(a35,10x,a)') "--polgrad",         "# calculate polarizability gradient"
      write(*,'(a35,10x,a)') "--dipgrad",         "# calculate dipole moment gradient"
      write(*,'(a35,10x,a)') "--geoopt",          "# geometry optimization"
      write(*,'(a35,10x,a)') "--nocosx",          "# do not use COSX approximation"
      write(*,'(a35,10x,a)') "--conv <convergence setting>","# set convergence setting, see ORCA manual"
      write(*,'(a35,10x,a)') "--tightscf",        "# (tight convergence)"
      write(*,'(a35,10x,a)') "--strongscf",       "# (strong convergence)"
      write(*,'(a35,10x,a)') "--nouseshark",      "# (use different integral library)"
      write(*,'(a35,10x,a)') "--hfref",           "# prepare HF reference calculation with super strong settings"
   end subroutine helpf
end module miscellaneous
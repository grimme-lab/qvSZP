# qvSZP
This tool sets up an ORCA calculation using the q-vSZP basis set (submitted to _The Journal of Chemical Physics_). It is intended to work with ORCA version 5.0.4 and higher.

## Installation
### Building with Fortran package Manager
To build the project, you can use the Fortran Package Manager (https://github.com/fortran-lang/fpm) in version 0.9 or newer.
To install the project in your preferred path, just use 
```
fpm install -profile release -prefix [path]
```
More information about FPM can be found in the respective documentation.
The project depends on other subprojects, the most important being `tblite` (https://github.com/tblite/tblite) and `stdlib` (https://github.com/fortran-lang/stdlib).

To build this project from the source code in this repository, you need to have
- a Fortran compiler supporting Fortran 2008
- a LAPACK / BLAS provider, like MKL or OpenBLAS

## Usage
To set up an input file using the q-vSZP basis set for a given geometry, you have to execute `qvSZP` in a directory with a molecular structure file (can be either `.xyz`, `coord`, or other common formats (see `mctc-lib` (https://github.com/grimme-lab/mctc-lib) for possible formats).

You need the files:
- `.basisq` and `.ecpq` in a known location (default: `$HOME/<file>`. The individual location of the files can be provided with the sample input below (or press `--help`).

Example program command-line calls:

```
qvSZP --struc coord.benzene
qvSZP --struc ch3.xyz --bfile /home/$USER/basissets/basisq --efile /home/$USER/basissets/ecpq --chrg -1
```
If no `--struc` file is explicitly given, `qvSZP` assumes a `coord` file.

See the `-help` flag for further input possibilities.
If you should observe instabilities with the `PModel` guess in ORCA, try to use `qvSZP` together with the `--suggestedguess` flag or provide an individual guess option with `--guess`.

# qvSZP
This tool sets up an ORCA calculation using the q-vSZP basis set ([_J. Chem. Phys. 159_, 164108 (**2023**)](https://doi.org/10.1063/5.0172373)). It is intended to work with ORCA version 5.0.4 and higher. The project depends on other subprojects, the most important being `tblite` (https://github.com/tblite/tblite) and `stdlib` (https://github.com/fortran-lang/stdlib).

The basis sets itself is located in `q-vSZP_basis/`. Besides the full q-vSZP basis set, also versions without polarization functions and without core electrons for *f* elements are available.

## Installation
### Building with Fortran package Manager
To build the project, you can use the Fortran Package Manager (https://github.com/fortran-lang/fpm) in version `0.9.0` or newer.
To install the project in your preferred path, just use 
```
fpm install -profile release -prefix [path]
```
Depending on your system, you might need to adapt the BLAS/LAPACK library in the `fpm.toml` file:
```
link=["lapack"] # or "openblas"
```
On `macOS-arm64`, the respective library is `openblas`. The build command is then:
```
fpm build --profile release --compiler gfortran-14 --c-compiler gcc-14 --link-flag "-L/opt/homebrew/opt/openblas/lib" --c-flag "-I/opt/homebrew/opt/openblas/include"
```
More information about FPM can be found in the respective documentation.


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

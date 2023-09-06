# qvSZP
This tool sets up an ORCA calculation using the q-vSZP basis set (submitted to _The Journal of Chemical Physics_). It is intended to work with ORCA version 5.0.4 and higher. The project depends on other subprojects, the most important being `tblite` (https://github.com/tblite/tblite) and `stdlib` (https://github.com/fortran-lang/stdlib).

## Installation
### Release version (recommended)
The use of the statically-linked release binary [`qvSZP`] is recommended. The binary has to be added to a location belonging to your `$PATH` variable.

### Requirements for building from source
To build this project from the source code in this repository, you need to have one of the following two build systems:
- [meson](https://mesonbuild.com) version 0.57.2 or newer, with a build-system backend, *i.e.* [ninja](https://ninja-build.org) version 1.10 or newer _OR_
- [fpm](https://github.com/fortran-lang/fpm) version 0.9.0 or newer

Moreover, you need to have:
- a Fortran compiler supporting Fortran 2008
- a LAPACK / BLAS provider, like MKL or OpenBLAS
- for `meson` _required_, for `fpm` _not required_: Installed Fortran standard library (https://github.com/fortran-lang/stdlib):
  - `git clone git@github.com:fortran-lang/stdlib.git` and move into the new directory.
  - `FC=ifort CC=icc CXX=icpc cmake -B build -DBUILD_SHARED_LIBS=on -DCMAKE_INSTALL_PREFIX=$HOME/.local` (if you have problems, try to omit `-DBUILD_SHARED_LIBS=on`).
  - `cmake --build build`
  - `cmake --install build`

### Building with Fortran package Manager
To build the project, you can use the Fortran Package Manager (https://github.com/fortran-lang/fpm) in version 0.9 or newer.
To install the project in your preferred path, just use 
```
fpm install -profile release -prefix [path]
```
More information about FPM can be found in the respective documentation.

### Building with Meson
You can use meson to build the project from source. Caution: For building with meson, the Fortran Standard Library (`stdlib`) has to be installed beforehand (see below).
To install the `qvSZP` project in your preferred path (assuming the use of Intel compilers), setup a build in the (new) directory `_build` with: 
```
FC=ifort CC=icc CXX=icpc meson setup _build --buildtype release --prefix=[path] -Dfortran_link_args=-qopenmp
```
You can compile (and install) the build with the following commands:
```
meson compile -C _build
meson install -C _build
```

If you want to setup a static compilation, replace the above meson setup with the following:
```
FC=ifort CC=icc CXX=icpc meson setup _build --buildtype release -Dfortran_link_args="-static -qopenmp -lifcoremt" --default-library=static --prefix=[path]
```

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

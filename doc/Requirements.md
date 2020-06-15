## OpenFOAM&reg; System Requirements

OpenFOAM requires a functioning C++11 compiler and GNU `make` build toolchain.

### Minimum recommended versions

- gcc : 4.8.5
- cmake: 3.8 (required for ParaView and CGAL build)
- boost: 1.48 (required for CGAL build and some functionality)
- fftw: 3.3.7 (recommended - required for FFT-related functionality)
- paraview: 5.5.2 (for visualization)

If using the Intel&reg; compiler, `17.0.1 20161005` is the minimum
usable version.

To check the installed versions

| Program       | To check the version  |
|---------------|-----------------------|
| gcc           | `gcc --version`       |
| icc           | `icc --version`       |
| cmake         | `cmake --version`     |
| openmpi       | `orterun --version`   |


#### Cautionary note for system openmpi

When using system openmpi, some caution is required due to the
[openmpi issue 5375](https://github.com/open-mpi/ompi/issues/5375) that
affected several versions of the openmpi2 or openmpi3 series.

| series    | major/minor | Minimum version
|-----------|-------------|-------------------|
| openmpi1  | all         | not affected      |
| openmpi2  | 2.0         | avoid this series |
| openmpi2  | 2.1         | 2.1.6 |
| openmpi3  | 3.0         | 3.0.4 |
| openmpi3  | 3.1         | 3.1.4 |
| openmpi4  | all         | not affected |



### Additional utilities

- flex : ***not 2.6.2, 2.6.3*** (fails for building scotch)
- m4 : no known minimum level
- QT : 5.9 (optional - for ParaView build)


### Ubuntu (eg, 20.04)

Install dependencies by executing the following lines on the command line:
```
sudo apt-get update
sudo apt-get install build-essential autoconf autotools-dev cmake gawk gnuplot
sudo apt-get install flex libfl-dev libreadline-dev zlib1g-dev openmpi-bin libopenmpi-dev mpi-default-bin mpi-default-dev
sudo apt-get install libgmp-dev libmpfr-dev libmpc-dev
```
If you intend to use system components, you can also install the following:
```
apt-get install libscotch-dev libptscotch-dev libfftw3-dev libboost-system-dev libboost-thread-dev libcgal-dev
```

Additional libraries will be required if compiling ParaView from
source, however it is recommended to skipped this initially since
it generally represent the main compilation difficulty.
In many cases, a system installation of paraview or a
[precompiled binary][download ParaView]
will be much easier.

Some libraries can be installed from the operating system, or
compiled from the ThirdParty directory.
The default configuration for OpenFOAM assumes OpenMPI from the system
and ThirdParty installations for most others (boost, CGAL, FFTW,
Kahip, Scotch). This is generally the most portable configuration
across various Linux distributions, but it may be desirable to use
more system libraries on Ubuntu.

To inspect the available system versions, use the `apt-cache show`
command. For example,
```
sudo apt-cache show libboost-dev
sudo apt-cache show libfftw3-dev
...
```

| Program   | apt-cache show  | Ubuntu  | Version |
|-----------|-----------------|---------|---------|
| boost     | libboost-dev    | 20.04   | 1.71.0  |
| CGAL      | libcgal-dev     | 20.04   | 5.0.2   |
| FFTW      | libfftw3-dev    | 20.04   | 3.3.8   |
| scotch    | libscotch-dev   | 20.04   | 6.0.9   |


| Program   | Ubuntu    | Program version |
|-----------|-----------|-----------------|
| gcc       | 20.04     | 9.3.0           |
| openmpi   | 20.04     | 4.0.3           |
| cmake     | 20.04     | 3.16.3          |
| flex      | 20.04     | 2.6.4           |
| m4        | 20.04     | 1.4.18          |


### openSUSE (eg, Leap-15.1)

Install the dependencies by copying and pasting the following lines to
the command line:

```
sudo zypper install -t pattern devel_C_C++
sudo zypper install cmake gnuplot flex libfl-devel readline-devel zlib-devel openmpi-devel
sudo zypper install libgmp-devel libmpfr-devel libmpc-devel
```
If you intend to use system components, you can also install the following:
```
sudo zypper install fftw3-devel libboost_system-devel libboost_thread-devel
```
but note that scotch and cgal are only available via the science repository.

This installs

| Program   | openSUSE  | Program version |
|-----------|-----------|-----------------|
| gcc       | 15.1      | 7.4.3           |
| openmpi   | 15.1      | 1.10.7          |
| cmake     | 15.1      | 3.10.2          |
| flex      | 15.1      | 2.6.4           |
| m4        | 15.1      | 1.4.18          |


#### OpenMPI

Check that the openmpi installation can be found:
```
orterun --version
```
And the command `mpicc --show` should display a complete compilation
line. This information is used in OpenFOAM to obtain the
compilation and link options.

For openSUSE it is common that this command cannot be found.
The reason being that the operating system can have several different
MPI vendors and versions installed and the user or sysadmin needs to
defined the preferred MPI. For this task, the `mpi-selector` and
`mpi-selector-menu` programs are used. These are automatically installed
by the system when openmpi has been installed.

If openmpi has been selected, the output will resemble the following:
```
$ mpi-selector --query

default:openmpi
level:user
```
If this is not the case, the `mpi-selector-menu` can be used to define
the preferred default.

***Note that changes in the preferred MPI do not take effect until the
next login.***

If all else fails, a brute force solution can also be used by simply
adding the corresponding lines to the `$HOME/.bashrc` file:
```
export PATH=$PATH:/usr/lib64/mpi/gcc/openmpi/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib64/mpi/gcc/openmpi/lib
```
This solution is not particularly nice, but may also be necessary if
some other part of the operating system installation is incomplete.


### Using IntelMPI

To run with IntelMPI, the LD_PRELOAD environment variable can be used
to pre-load the appropriate Intel&reg; MPI binding library. For more
details, see INTELMPI release note page 13.
The following line can be the `$HOME/.bashrc` file:
```
export LD_PRELOAD="libmpi.so"
```


### Additional libraries

When compiled ParaView from source additional dependencies will be
required.
A partial list is given in the [ThirdParty requirements][link third-require].
**Please help us with keeping that information up-to-date and accurate.**

<!-- Links -->

[page ParaView]:  http://www.paraview.org/
[download ParaView]: https://www.paraview.org/download/


<!-- OpenFOAM -->

[link openfoam-readme]: https://develop.openfoam.com/Development/openfoam/blob/develop/README.md
[link openfoam-config]: https://develop.openfoam.com/Development/openfoam/blob/develop/doc/Config.md
[link openfoam-build]: https://develop.openfoam.com/Development/openfoam/blob/develop/doc/Build.md
[link openfoam-require]: https://develop.openfoam.com/Development/openfoam/blob/develop/doc/Requirements.md
[link third-readme]: https://develop.openfoam.com/Development/ThirdParty-common/blob/develop/README.md
[link third-build]: https://develop.openfoam.com/Development/ThirdParty-common/blob/develop/BUILD.md
[link third-require]: https://develop.openfoam.com/Development/ThirdParty-common/blob/develop/Requirements.md

---
Copyright 2019-2020 OpenCFD Ltd

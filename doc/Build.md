## OpenFOAM&reg; Quick Build Guide

Prior to building, ensure that the [system requirements][link openfoam-require]
are satisfied (including any special [cross-compiling][link openfoam-cross]
considerations), and source the correct OpenFOAM environment.
For example, for the OpenFOAM-v1906 version:
```
source /installation/path/OpenFOAM-v1906/etc/bashrc
```

## Preliminaries

The [third-party][repo third] directory includes a
[readme][link third-readme] and additional information about
[requirements][link third-require] and a more detailed
[build guide][link third-build].

Some known build issues related to specific compiler and VTK library versions
can be found in the [$WM_PROJECT_DIR/doc/BuildIssues.md][link openfoam-issues]
file.

If you need to change the default versions for third-party libraries,
or use system libraries for some components, please some additional
information about the [config structure][link openfoam-config].

## Compile OpenFOAM

The compilation process is self-contained and will compile and install
all OpenFOAM code and dependencies.

- Test the system readiness (optional, not supported for cross-compilation)
```
foamSystemCheck
```
- Change to the main OpenFOAM directory ($WM_PROJECT_DIR).
  If this fails, the environment is not correctly configured.
```
foam
```
- Compile OpenFOAM
```
./Allwmake -s -l
```
- In case you need to stop the compilation, continue later by running
`./Allwmake` again.

## Compile OpenFOAM faster

For faster compilation, users should take advantage of multi-processor
machines when building the code. This is supported directly by `wmake`
and the `Allwmake` scripts. For example,
```
wmake -j               # Use all cores
wmake -j 8             # Use specified number of cores
```
It can also be helpful to use the builtin queuing (the `-queue`
option), which collects subdirectories and dispatches to make in
larger chunks.

The following compilation sequence can be useful:
```
./Allwmake -j -s -q -l
```
This compiles with all cores (-j), reduced output (-s, -silent), with
queuing (-q, -queue) and logs (-l, -log) the output to a file such as
`log.linux64GccDPInt32Opt` for later inspection.

If you use the `-k` option (`-keep-going` = ignore errors) to compile
as much as possible on the first pass, be certain to follow that with
second pass (without the `-k` option) at a later stage to ensure that
you haven't missed any error messages.

## Post-compilation steps

- Open a new shell and source the OpenFOAM environment to see all
  changes (refer to top of page).
- Validate the build (not supported for cross-compilation) by running
```
foamInstallationTest
```
- Create the user `run` directory:
```
mkdir -p $FOAM_RUN
```
- Test the installation with a simple tutorial:
```
run
cp -r $FOAM_TUTORIALS/incompressible/simpleFoam/pitzDaily ./
cd pitzDaily
blockMesh
simpleFoam
```

### ParaView

OpenFOAM ships with ParaView sources for post-processing OpenFOAM
field results. However, the paraview version distributed with
the operating system or a [binary package][download ParaView]
will often be sufficient, and avoids additional compilation complexity.

If you do wish to compile ParaView from source, it is recommended
that you do so ***after*** completing an initial compilation of OpenFOAM.
This gets the process started much more quickly. At a later stage,
OpenFOAM can be updated to compile with paraview. Only the affected
applications will be compiled (eg, the blockMesh reader module) and the
balance of the OpenFOAM installation will not affected.

If you decide to compiling in two passes, you only need to execute the
top-level `Allwmake` a second time. Do **not** use `wclean` to force a
complete rebuild! This is unnecessary.

More details in the [ThirdParty build guide][link third-build].


------------

<!-- Links -->

[page ParaView]:  http://www.paraview.org/
[download ParaView]: https://www.paraview.org/download/


<!-- OpenFOAM -->

[repo openfoam]: https://develop.openfoam.com/Development/OpenFOAM-plus/
[repo third]: https://develop.openfoam.com/Development/ThirdParty-plus/

[link openfoam-readme]: https://develop.openfoam.com/Development/OpenFOAM-plus/blob/develop/README.md
[link openfoam-issues]: https://develop.openfoam.com/Development/OpenFOAM-plus/blob/develop/doc/BuildIssues.md
[link openfoam-config]: https://develop.openfoam.com/Development/OpenFOAM-plus/blob/develop/doc/Config.md
[link openfoam-build]: https://develop.openfoam.com/Development/OpenFOAM-plus/blob/develop/doc/Build.md
[link openfoam-cross]: https://develop.openfoam.com/Development/OpenFOAM-plus/blob/develop/doc/Cross-Compile-mingw.md
[link openfoam-require]: https://develop.openfoam.com/Development/OpenFOAM-plus/blob/develop/doc/Requirements.md
[link third-readme]: https://develop.openfoam.com/Development/ThirdParty-plus/blob/develop/README.md
[link third-build]: https://develop.openfoam.com/Development/ThirdParty-plus/blob/develop/BUILD.md
[link third-require]: https://develop.openfoam.com/Development/ThirdParty-plus/blob/develop/Requirements.md

---
Copyright 2019 OpenCFD Ltd

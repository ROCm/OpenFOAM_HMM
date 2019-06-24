# Notes for cross-compiling with mingw

## Minimum version

The mingw cross-compiler should be at least version 8.2.0 (tested) or
slightly older. Versions that are much older may have faultly regex
implementations.

## Host setup

On openSUSE use the packages for compilation:
```
mingw64-cross-binutils
mingw64-cross-cpp
mingw64-cross-gcc
mingw64-cross-gcc-c++
mingw64-filesystem
mingw64-headers
mingw64-runtime

mingw64-libwinpthread1
mingw64-winpthreads-devel

mingw64-libfftw3
mingw64-fftw3-devel
```
This setup is missing `zlib`, so download that manually and compile as a
*static* library.
```
CC="$(wmake -show-c)" CFLAGS="$(wmake -show-cflags)" ./configure --static
make
```

The resulting output files (zconf.h, zlib.h) and (libz.a) either need
to be installed in system locations where OpenFOAM can find them, or if
they are to be shipped directly with OpenFOAM, they can also be placed
in the `src/OpenFOAM/include` and `platforms/XXX/lib` paths.

If the header files are only needed during compilation, it can be a
fairly convenient hack to simply place copies of them in the
`src/OSspecific/MSwindows` directory.

Flex is used in a few locations within OpenFOAM for generating code.
The generated C++ code requires the `FlexLexer.h` header file, but
its `/usr/include` location will be ignored by the cross-compiler.

As another ugly hack, a copy of this file can be made in a standard
project include location. For example,
```
ln -s /usr/include/FlexLexer.h src/OSspecific/MSwindows
```

The last point to consider when cross-compiling is the behaviour of
the OpenFOAM wmake toolchain used during compilation. These are found
under `wmake/src`. If the `Makefile` is used directly, executables
will be created that work on the target platform (Windows), but *not*
on the host platform (which is what is required). This is addressed
directly by the `wmake/src/Allmake` script, which will use the system
gcc to create host binaries for the wmake toolchain. If, for some
reason, you also require target wmake toolchain binaries, you will
need invoke make manually within the `wmake/src` directory.


The settings for cross-compilation are normally defined in the
`etc/pref.sh` file with contents like this:
```
# For mingw cross-compile

export WM_COMPILER=Mingw
export WM_MPLIB=MSMPI

export WM_LABEL_SIZE=32
# other settings...
```

Additional adjustments may be required in some other places. For example
in `etc/config.sh/FFTW`
```
fftw_version=fftw-system
export FFTW_ARCH_PATH=/usr/x86_64-w64-mingw32/sys-root/mingw
```


## Run-time setup

When using the cross-compiled executables and libraries, the
corresponding runtime libraries will be required.
These will need to be copied across from the Linux host system to the
target machine.
On openSUSE these runtime libraries are provided by the packages:
```
mingw64-libgcc_s_seh1
mingw64-libstdc++6
```

When running, the `WM_PROJECT_DIR` environment must be set.
OpenFOAM will otherwise not be able to locate its files.


## Known limitations (2019-06-24)

- kahip does not build
- boost should build ok, but no CGAL support (ie, no foamyHexMesh)
- no ParaView plugin, runTimePostProcessing
- reacting EulerFoam solvers have too many interdependencies and do
  not yet compile cleanly.
  It is advisable to compile with the wmake `-k` option to keep going
  even when the EulerFoam solvers fail to compile.

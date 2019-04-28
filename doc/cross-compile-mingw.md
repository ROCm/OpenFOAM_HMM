Notes for cross-compiling with mingw.

On openSUSE use the packages:
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

This setup is missing zlib, so download that manually and compile as a
*static* library.

```
CC="$(wmake -show-c) CFLAGS="$(wmake -show-cflags) ./configure --static
```

The resulting output files (zconf.h zlib.h) and (libz.a) either need
to be installed in system locations where OpenFOAM can find it, or if
they are to be shipped directly with OpenFOAM, they can also be placed
in the `src/OpenFOAM/include` and `platforms/XXX/lib` paths.

When using the cross-compiled executables and libraries, additional
runtime libraries are required. On openSUSE these are provided by the
packages:
```
mingw64-libgcc_s_seh1
mingw64-libstdc++6
```

In a few locations within OpenFOAM, flex is used to generate code. The
generated C++ code requires the `FlexLexer.h` header file.
This is normally located under `/usr/include/FlexLexer.h`, which will
be ignored by the cross-compiler.

As a fairly ugly hack, a copy of this file can be made in a standard
project include location. For example,

```
mkdir src/OpenFOAM/lnInclude
cp /usr/include/FlexLexer.h src/OpenFOAM/lnInclude
```


The last point to consider when cross-compiling is the behaviour of
the OpenFOAM tools used during compilation. These are found under
`wmake/src`. When these are compiled they will create executables that
work on the target platform (Windows), but *not* on the host platform.

The workaround:

1. Activate the native OpenFOAM environment (Linux with system gcc)
and use that to compile the build tools
```
wmake/src/Allmake
```
This can be skipped if you already have an existing OpenFOAM build.

2. Activate the OpenFOAM for cross-compiling (Linux with mingw)
and use that to compile the build tools
```
wmake/src/Allmake
```

3. Copy the contents of the native platform build tools into the
cross-compilation platform directory. For example,
```
cp wmake/platforms/linux64Gcc/* wmake/platforms/linux64Mingw/
```

The `wmake/platforms/linux64Mingw/` directory should now contain
native and cross-compiled tools:
```
dirToString
wmkdep
wmkdepend

dirToString.exe
wmkdep.exe
wmkdepend.exe
```

The native tools are the one that will (automatically) be used
throughout the balance of the cross-compilation process.


Adjust paths for third-party sources. Eg, in `etc/config.sh/FFTW` you
may wish to have the following:
```
fftw_version=fftw-system
export FFTW_ARCH_PATH=/usr/x86_64-w64-mingw32/sys-root/mingw
```


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


Known limitations (2019-05-01)

- No CGAL support (ie, no foamyHexMesh)
- No ParaView plugin, runTimepostProcessing
- reacting EulerFoam solvers have too many interdependencies and do
  not yet compile.

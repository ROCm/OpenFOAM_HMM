## Known Build Issues (OpenFOAM-v1906)

### Thermo problems with Clang

Clang builds required updates to the thermophysical libraries to prevent
optimised builds from generating sigFpe's.  The changes are wrapped in `#ifdef`
`__clang__` statements to not affect other compilers.

The following tutorials experience known failures:

- combustion/XiFoam/RAS/moriyoshiHomogeneous
- multiphase/reactingTwoPhaseEulerFoam/laminar/bubbleColumnEvaporatingDissolving


This will be further investigated to identify the root cause.


### Intel MPI with Gcc/Clang

Either `I_MPI_ROOT` (preferred) or `MPI_ROOT` can be used to specify
the Intel-MPI installation directory path.
The ThirdParty build of ptscotch uses `mpiicc` for Intel-MPI instead
of the usual `mpicc`. When gcc or clang are used, it is quite likely
that the `I_MPI_CC` environment variable also needs to be set
accordingly.
See `mpiicc -help` for more information about environment variables.


### VTK

If using the runTimePostProcessing to create on-the-fly images, you
can simply just compile ParaView and these libraries will be used.
If you elect to use a separate VTK compilation (for example for
off-screen rendering), it is advisable to reuse the VTK libraries that
are provided with ParaView by making an appropriate symlink
prior to using makeVTK. This doesn't just reduce disk-space, but works
much better than using the VTK tar file.

Using runTimePostProcessing with the *plain* VTK libraries does
generally work, but may not exit cleanly:
```
symbol lookup error: .../linux64Gcc/VTK-7.1.0/lib/libvtkCommonExecutionModel-7.1.so.1:
undefined symbol: _ZN33vtkFilteringInformationKeyManager13ClassFinalizeEv

symbol lookup error: .../linux64Gcc/VTK-7.1.0/lib/libvtkCommonDataModel-7.1.so.1:
undefined symbol: _ZN49vtkInformationQuadratureSchemeDefinitionVectorKeyD1Ev
```

This error appears to be suppressed if VTK is compiled with a `Debug` build-type.


### Building on older systems

If the system gcc is too old for building OpenFOAM, a third-party gcc or
clang/llvm installation can be used. If building clang/llvm, note that
there are also minimum gcc/g++ requirements as listed in the
detailed [build guide][link third-build].

If your system compiler is too old to build the minimum required gcc or
clang/llvm, it is just simply too old.


### ThirdParty clang without gmp/mpfr

If using ThirdParty clang without gmp/mpfr, the ThirdParty makeCGAL
script will need to be run manually and specify that there is no
gmp/mpfr. Eg,
```
cd $WM_THIRD_PARTY_DIR
./makeCGAL gmp-none mpfr-none
```

Subequent compilation with Allwmake will now run largely without any
problems, except that the components linking against CGAL
(foamyMesh and surfaceBooleanFeatures) will also try to link against
a nonexistent mpfr library. As a workaround, the link-dependency can
be removed in wmake/rules/General/CGAL :
```
CGAL_LIBS = \
    -L$(BOOST_ARCH_PATH)/lib \
    -L$(BOOST_ARCH_PATH)/lib$(WM_COMPILER_LIB_ARCH) \
    -L$(CGAL_ARCH_PATH)/lib \
    -L$(CGAL_ARCH_PATH)/lib$(WM_COMPILER_LIB_ARCH) \
    -lCGAL
```

A robuster solution is still being sought.


### Building with spack

If you are building with spack, note that the `depends_on` for paraview
resolves poorly. The `+qt` dependency (for building the reader module)
may need to be specified as a preference by including the following in
your `~/.spack/packages.yaml` file:
```
packages:
    paraview:
        variants: +qt
```
It appears that spack will otherwise ignore any `paraview+qt` version
and attempt to install a `paraview~qt` version instead.


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
[link openfoam-require]: https://develop.openfoam.com/Development/OpenFOAM-plus/blob/develop/doc/Requirements.md
[link third-readme]: https://develop.openfoam.com/Development/ThirdParty-plus/blob/develop/README.md
[link third-build]: https://develop.openfoam.com/Development/ThirdParty-plus/blob/develop/BUILD.md
[link third-require]: https://develop.openfoam.com/Development/ThirdParty-plus/blob/develop/Requirements.md

---
Copyright 2019 OpenCFD Ltd

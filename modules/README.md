[[_TOC_]]

## OpenFOAM Modules

This directory is a location for additional OpenFOAM components or
tools to placed and have them built as part of the normal OpenFOAM
build process. It is assumed that each subdirectory contain an
appropriate `Allwmake` (or `Allwmake.override`) file.


### Build locations

Any individual _module_ will normally also be able to exist outside of
the module directory structure and will typically build into user
locations (`$FOAM_USER_APPBIN`, `$FOAM_USER_LIBBIN`).

When compiled from the top-level OpenFOAM `Allwmake` or the
`modules/Allwmake`, they should build into OpenFOAM project locations
(`$FOAM_APPBIN`, `$FOAM_LIBBIN`). This can be adjusted by
supplying an alternative `-prefix=` to the corresponding Allwmake
command.

| Command    | Install location |
|------------|------------------|
| ./Allwmake -prefix=user | `$FOAM_USER_APPBIN`, `$FOAM_USER_LIBBIN` |
| ./Allwmake -prefix=group | `$FOAM_SITE_APPBIN`, `$FOAM_SITE_LIBBIN` |
| ./Allwmake -prefix=openfoam | `$FOAM_APPBIN`, `$FOAM_LIBBIN` |
| ./Allwmake -prefix=/some/pathname | `/some/pathname/bin`, `/some/pathname/lib` |


### Adding additional components

These additional components may be added as [git submodules][man git-submodule],
by script or by hand.


#### git

On the first use, it will be necessary to register the submodules:
```
git submodule init
```

This will clone the relevant submodules from their respective
repositories.

The following will indicate the current state:
```
git submodule status
```

On the first use, or after merging upstream changes in the OpenFOAM
repository, it will be necessary to update the submodules:
```
git submodule update
```

A quick overview of `git submodule` can be in this
[*blog*][blog git-submodule] with full details in the
[*manpage*][man git-submodule].


An easy way to see which submodules are actually in use:
```
cat .gitmodules
```

Which will reveal content resembling the following:
```
[submodule "avalanche"]
    path = modules/avalanche
    url = https://develop.openfoam.com/Community/avalanche.git
[submodule "cfmesh"]
    path = modules/cfmesh
    url = https://develop.openfoam.com/Community/integration-cfmesh.git
...
```

### Documentation (doxygen)

To build the doxygen information for the components, it is also
necessary to link the directories to the doc/ subdirectory.
This is a purely manual operation.


### Developer Information

#### Build locations

To accomodate building into various locations, the module code should
be adapted with the following changes:

- ***Make/files***
   ```
   ...
   EXE = $(FOAM_MODULE_APPBIN)/someExecutable

   LIB = $(FOAM_MODULE_LIBBIN)/libSomeLibrary
   ```

- `Make/options` should include this
  ```
  include $(GENERAL_RULES)/module-path-user
  ...
  ```

The following changes to `Make/options` are universally applicable
(ie, work with older or other versions of OpenFOAM), but more verbose.

- `Make/options` with the following
  ```
  sinclude $(GENERAL_RULES)/module-path-user

  /* Failsafe - user locations */
  ifeq (,$(FOAM_MODULE_APPBIN))
  FOAM_MODULE_APPBIN = $(FOAM_USER_APPBIN)
  endif
  ifeq (,$(FOAM_MODULE_LIBBIN))
  FOAM_MODULE_LIBBIN = $(FOAM_USER_LIBBIN)
  endif
  ...
  ```

<!-- General Information -->

[man git-submodule]:  https://git-scm.com/docs/git-submodule
[blog git-submodule]: http://blog.joncairns.com/2011/10/how-to-use-git-submodules/

---

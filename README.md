# About OpenFOAM
OpenFOAM is a free, open source CFD software [released and developed primarily by OpenCFD Ltd](http://www.openfoam.com) since 2004. It has a large user base across most areas of engineering and science, from both commercial and academic organisations. OpenFOAM has an extensive range of features to solve anything from complex fluid flows involving chemical reactions, turbulence and heat transfer, to acoustics, solid mechanics and electromagnetics.  [More...](http://www.openfoam.com/documentation)


OpenFOAM is professionally released every six months to include
customer sponsored developments and contributions from the community -
individual and group contributors, fork re-integrations (including from
FOAM-extend and OpenFOAM Foundation Ltd) - in this Official Release
sanctioned by the OpenFOAM Worldwide Trademark Owner aiming towards
one OpenFOAM.


# Copyright
OpenFOAM is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.  See the file COPYING in this directory or [http://www.gnu.org/licenses/](http://www.gnu.org/licenses), for a description of the GNU General Public License terms under which you can copy the files.


# OpenFOAM Trademark
OpenCFD Ltd grants use of its OpenFOAM trademark by Third Parties on a licence basis. ESI Group and OpenFOAM Foundation Ltd are currently permitted to use the Name and agreed Domain Name. For information on trademark use, please refer to the [trademark policy guidelines](http://www.openfoam.com/legal/trademark-policy.php).

Please [contact OpenCFD](http://www.openfoam.com/contact) if you have any questions on the use of the OpenFOAM trademark.

Violations of the Trademark are continuously monitored, and will be duly prosecuted.


# Using OpenFOAM

If OpenFOAM has already been compiled on your system, simply source
the appropriate `etc/bashrc` or `etc/cshrc` file and get started.
For example, for the OpenFOAM-v1906 version:
```
source /installation/path/OpenFOAM-v1906/etc/bashrc
```

# Compiling OpenFOAM

If you are compiling OpenFOAM from source, please see the relevant
guides:

| Location    | Readme    | Requirements | Build |
|-------------|-----------|--------------|-------|
| [OpenFOAM][repo openfoam] | [readme][link openfoam-readme] | [system requirements][link openfoam-require] | [build][link openfoam-build] |
| [ThirdParty][repo third] | [readme][link third-readme] | [system requirements][link third-require] | [build][link third-build] |


# How do I know which version I am currently using?

The value of the `$WM_PROJECT_DIR` or even `$WM_PROJECT_VERSION` are
not guaranteed to have any correspondence to the OpenFOAM release
(API) value. If OpenFOAM has already been compiled, the build-time
information is embedded into each application. For example, as
displayed from `blockMesh -help`:
```
Using: OpenFOAM-v1812.local (1812) (see www.OpenFOAM.com)
Build: 65d6551ff7-20190530 (patch=190531)
Arch:  LSB;label=32;scalar=64
```
This output contains all of the more interesting information that we need:

| item          | value         |
|---------------|---------------|
| version       | v1812.local   |
| api           | 1812          |
| commit        | 65d6551ff7    |
| author date   | 20190530      |
| patch-level   | (20)190531    |

As can be seen in this example, the git build information is
supplemented by the date when the last change was authored, which can
be helpful when the repository contains local changes. If you simply
wish to know the current API and patch levels directly, the
`wmakeBuildInfo` script provides the relevant information even
when OpenFOAM has not yet been compiled:
```
$ wmakeBuildInfo
make
    api = 1812
    patch = 190531
    branch = master
    build = 65d6551ff7-20190530
```
Similar information is available with `foamEtcFile`, using the
`-show-api` or `-show-patch` options. For example,
```
$ foamEtcFile -show-api
1812

$ foamEtcFile -show-patch
190531
```
This output will generally be the easiest to parse for scripts.
The `$FOAM_API` convenience environment variable may not reflect the
patching changes made within the currently active environment and
should be used with caution.


# ThirdParty directory

OpenFOAM normally ships with a directory of 3rd-party software and
build scripts for some 3rd-party software that is either necessary or
at least highly useful for OpenFOAM, but which are not necessarily
readily available on every operating system or cluster installation.

These 3rd-party sources are normally located in a directory parallel
to the OpenFOAM directory. For example,
```
/path/parent
|-- OpenFOAM-v1906
\-- ThirdParty-v1906
```
There are, however, many cases where this simple convention is inadequate:

* When no additional 3rd party software is actually required (ie, the
  operating system or cluster installation provides it)

* When we have changed the OpenFOAM directory name to some arbitrary
  directory name, e.g. openfoam-sandbox1906, etc..

* When we would like any additional 3rd party software to be located
  inside of the OpenFOAM directory to ensure that the installation is
  encapsulated within a single directory structure. This can be
  necessary for cluster installations, or may simply be a convenient
  means of performing a software rollout for individual workstations.

* When we have many different OpenFOAM directories for testing or
  developing various different features but wish to use or reuse the
  same 3rd party software for them all.

The solution for these problems is a newer, more intelligent discovery when locating the ThirdParty directory with the following precedence:

1. PROJECT/ThirdParty
   * for single-directory installations
2. PREFIX/ThirdParty-VERSION
   * this corresponds to the traditional approach
3. PREFIX/ThirdParty-vAPI
   * allows for an updated value of VERSION, *eg*, `v1906-myCustom`,
     without requiring a renamed ThirdParty. The API value would still
     be `1906` and the original `ThirdParty-v1906/` would be found.
4. PREFIX/ThirdParty-API
   * this is the same as the previous example, but using an unadorned
     API value. This also makes sense if the chosen version name also
     uses the unadorned API value in its naming, *eg*,
     `1906-patch190131`, `1906.19W03`
5. PREFIX/ThirdParty-common
   * permits maximum reuse for various versions, but only for
     experienced user who are aware of potential version
     incompatibilities

If none of these directories are found to be suitable, it reverts to using PROJECT/ThirdParty as a dummy location (even if the directory does not exist). This is a safe fallback value since it is within the OpenFOAM directory structure and can be trusted to have no negative side-effects.
In the above, the following notation has been used:

| name          | value         | meaning       |
|---------------|---------------|---------------|
| PROJECT       | `$WM_PROJECT_DIR`     | The OpenFOAM directory |
| PREFIX        | `dirname $WM_PROJECT_DIR` | The OpenFOAM parent directory |
| API           | `foamEtcFiles -show-api` |  The api or release version |
| VERSION       | `$WM_PROJECT_VERSION` | The version we've chosen |

To reduce the potential of false positive matches (perhaps some other
software also uses ThirdParty-xxx for its naming), the directory test
is accompanied by a OpenFOAM-specific sanity test. The OpenFOAM
ThirdParty directory will contain either an `Allwmake` file or a
`platforms/` directory.


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


# Useful Links
- [Download and installation instructions](http://www.openfoam.com/download/)
- [Documentation](http://www.openfoam.com/documentation)
- [Reporting bugs/issues/feature requests](http://www.openfoam.com/code/bug-reporting.php)
- [OpenFOAM Community](http://www.openfoam.com/community/)
- [Contacting OpenCFD](http://www.openfoam.com/contact/)

Copyright 2016-2019 OpenCFD Ltd

OpenFOAM Modules
================

This directory is a location for additional OpenFOAM components or tools
to placed and have them built as part of the normal OpenFOAM build
process. It is assumed that each subdirectory contain an appropriate
Allwmake file, and that they in all likelihood also build into
`$FOAM_APPBIN` and `$FOAM_LIBBIN` instead of
`$FOAM_USER_APPBIN` and `$FOAM_USER_LIBBIN`.

These additional components may be added as [git submodules][man git-submodule],
by script or by hand.


### git

On the first use, it will be necessary to register the submodules:

    git submodule init


This will clone the relevant submodules from their respective
repositories.


The following will indicate the current state:

    git submodule status


On the first use, or after merging upstream changes in the OpenFOAM
repository, it will be necessary to update the submodules:

    git submodule update


A quick overview of `git submodule` can be in this
[*blog*][blog git-submodule] with full details in the
[*manpage*][man git-submodule].


An easy way to see which submodules are actually in use:

    `cat .gitmodules`

Which will reveal content resembling the following:

    [submodule "cfmesh"]
        path = modules/cfmesh
        url = https://develop.openfoam.com/Community/integration-cfmesh.git


### doxygen

To build the doxygen information for the components, it is also
necessary to link the directories to the doc/ subdirectory.
This is a purely manual operation.

<!-- General Information -->

[man git-submodule]:  https://git-scm.com/docs/git-submodule
[blog git-submodule]: http://blog.joncairns.com/2011/10/how-to-use-git-submodules/

---

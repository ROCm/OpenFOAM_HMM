OpenFOAM Modules
================

This directory is a location to place additional OpenFOAM components
or tools and have them built as part of the normal OpenFOAM build
process. It is assumed that each subdirectory contain an appropriate
Allwmake file, and that they should in all likelihood also build into
`$FOAM_APPBIN` and `$FOAM_LIBBIN` instead of
`$FOAM_USER_APPBIN` and `$FOAM_USER_LIBBIN`.

These additional components may be added as git submodules, by script
or by hand.

To build the doxygen information for the components, it is also
necessary to link the directories to the doc/ subdirectory.

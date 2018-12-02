OpenFOAM Configuration
----------------------

The main OpenFOAM settings are located in the parent `etc/` directory.
Both POSIX (bash, dash,...) and csh shells are supported.
To configure OpenFOAM, source either the `etc/bashrc` or the
`etc/cshrc` file, as appropriate for your shell.

These source the following files in the `config.sh/` or
`config.csh/` directories:

* `setup` : finalize setup of OpenFOAM environment (called by bashrc,cshrc)
* `settings` : core settings
* `aliases` : aliases for interactive shells
* `unset` : sourced to clear as many OpenFOAM environment settings as possible
* `mpi` : MPI communications library settings
* `ensight` : application settings for EnSight
* `paraview` : application settings for ParaView
* `scotch` : application settings for compiling against scotch
* `metis` : application settings for compiling against metis

The `config.*/example` directories contain additional example configuration
files for the corresponding shell:

* `compiler` : an example of fine tuning ThirdParty compiler settings
* `openmpi` : an example of fine tuning openmpi settings for OpenFOAM
* `paraview` : an example of chaining to the standard config/paraview
   with a different ParaView_VERSION
* `prefs`: an example of supplying alternative site-defined settings

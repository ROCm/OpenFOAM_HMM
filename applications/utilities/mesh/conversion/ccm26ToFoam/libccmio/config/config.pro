#  Settings in this file will override those in the default configuration.
#

######################### Compilation Options #####################
# Default release/debug options.  
isEmpty(RELEASEMODE) {
   RELEASEMODE = debug
}

#########################  MeshKernel Options #####################
# The base directory to install includes, libs, and objects 
# Defaults to PATH_TO_SRC/..
KERNELBASE=

# Directory to install include (header) files.
# Defaults to KERNELBASE/include
KERNELINC=

# Directory to install library files.
# Defaults to KERNELBASE/lib/BUILDNAME/RELEASEMODE
KERNELLIB=

# Directory to store object files.
# Defaults to KERNELBASE/obj/BUILDNAME/RELEASEMODE
KERNELOBJ=

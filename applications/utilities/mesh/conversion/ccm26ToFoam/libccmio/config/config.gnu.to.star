#!/bin/sh

# $Id: config.gnu.to.star,v 1.4 2006/06/05 21:12:16 geoffp Exp $

# Convert the modified gnu machine type used in ammbatch from the config.system
# script to the star machine type.  Return unknown if there is no known match.

# Also process MACHMOD ($2) to see if we want to compile with different than
# the default options.

# echo Input is ${1:-null}-${2:-null}

case ${1:-null}-${2:-null} in
    rs6000-ibm-aix4.3.0-32)            echo aix_4.3-com                     ;;
    rs6000-ibm-aix4.3.0-null)          echo aix64_4.3-pwr3                  ;;
    rs6000-ibm-aix5.1.0-32)            echo aix_5.1-pwr4                    ;;
    rs6000-ibm-aix5.1.0-null)          echo aix64_5.1-pwr4                  ;;
    hppa1.1-hp-hpux-10.20-null)        echo hpux_10.20-com                  ;;
    hppa2.0-hp-hpux-11.00-32)          echo hpux_11.00-pa8000               ;;
    hppa2.0-hp-hpux-11.00-null)        echo hpux64_11.00-pa8000             ;;
    ia64-hp-hpux11.22-null)            echo hpux64_11.22-itanium2           ;;
    mips-sgi-irix6-32)                 echo irix_6.5-mips3                  ;;
    mips-sgi-irix64_6-32)              echo irix_6.5-mips3                  ;;
    mips-sgi-irix64_6-null)            echo irix64_6.5-mips4                ;;
    i586-unknown-linux_2.2-glibc-null) echo linux_2.2-x86-glibc_2.2.0       ;;
    i586-unknown-linux_2.4-glibc-null) echo linux_2.4-x86-glibc_2.3.2       ;;
    ia64-unknown-linux-null)           echo linux64_2.4-itanium-glibc_2.2.4 ;;
    alpha-dec-osf3-null)               echo osf1_4.0-com                    ;;
    alpha-dec-osf5-null)               echo osf1_5.1-com                    ;;
    sparc-sun-sunos-5.8-32)            echo sunos_5.8-ultra                 ;;
    sparc-sun-sunos-5.8-null)          echo sunos64_5.8-ultra               ;;
    i686-pc-mingw32-null)              echo windows-x86-gcc                 ;;
    i686-pc-mingw32-vc)                echo windows-x86                     ;;
    x86_64-unknown-linux-gnu-null)     echo linux64_2.4-x86-glibc_2.2.5     ;;
    ppc64-unknown-linux-gnu-null)      echo linux64_2.6-pwr4-glibc_2.3.3    ;;
    *)                                 echo unknown                         ;;
esac

# Automatic setting of emacs local variables.
# Local Variables:
# tab-width: 8
# End:

#!/bin/tcsh -x

find . -type f \( -name '*.h' -o -name '*.c' -o -name '*.cpp' -o -name '*.f' -o -name '*.F' \) -print | grep -iv dynamic | grep -iv pro_tcl | sort | xargs -t etags

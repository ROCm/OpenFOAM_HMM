#!/bin/sh
grep "CGNS_VERSION" | cut -f3 -d" " | awk '{ major = substr($1,1,1); minor = substr($1,2,2); revision = substr($1,3,4); printf "VERSION=%d.%02d.%03d\n", major, minor, revision; }'

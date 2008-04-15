#!/bin/sh
grep "@(#)ADF Library Version" | cut -f2 -d"=" | cut -f1 -d">" | awk '{ major = index("ABCDEFGHIJKLMNOPQRSTUVWXYZ",substr($4,1,1)); if (major == 0) index("abcdefghijklmnopqrstuvwxyz",substr($4,1,1)); minor = substr($4,2,2); revision = 0; printf "VERSION=%d.%02d.%03d\n", major, minor, revision; }'

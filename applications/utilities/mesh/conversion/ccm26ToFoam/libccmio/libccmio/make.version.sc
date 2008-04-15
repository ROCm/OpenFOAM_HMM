#$!/bin/sh

awk '/^#define kCCMIOMajorVersion/ { major = $3; } /^#define kCCMIOMinorVersion/ { minor = $3; } /^#define kCCMIORevision/ { revision = $3 } END { printf "VERSION=%d.%02d.%03d\n", major, minor, revision; }'

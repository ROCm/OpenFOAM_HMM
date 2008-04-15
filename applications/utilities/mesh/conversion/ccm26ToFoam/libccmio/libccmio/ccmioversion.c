#include "ccmioversion.h"

/* ## does not evaluate macros, so we need to have several layers of indirection
   so that the kCCMIO* macros get evaluated before VersionConcatenate() is
   evaluated.*/
#define kVersionIdentifier VersionIdentifierImp(kCCMIOMajorVersion, kCCMIOMinorVersion, kCCMIORevision)
#define VersionIdentifierImp(a, b, c) VersionConcatenate(a,b,c)

#define VersionConcatenate(major, minor, revision) \
	version_##major##_##minor##_##revision

char kVersionIdentifier[] = kCCMIOVersionStr;

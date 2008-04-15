#define kCCMIOMajorVersion	2	/* Requires large-scale re-write */
#define kCCMIOMinorVersion	6	/* Changes to existing prototypes */
#define kCCMIORevision		1	/* Bug fixes, new functions */

#define kCCMIOVersion 20601	/* mmnnrr (major, minor, rev) */

#define Stringify(x)	StringifyImp(x)
#define StringifyImp(x)	#x

#define kCCMIOVersionStr	Stringify(kCCMIOMajorVersion)"."Stringify(kCCMIOMinorVersion)"."Stringify(kCCMIORevision)

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include "sff.h"
#include "sffbuffer.h"

const char kFilename[] = "/tmp/sff.adf";

const int kVersion = 1234;
const char kTitle[] = "Test title";
const char kI32ScalarName[] = "scalar_i32";
const int kI32ScalarValue = 1234567890;
const char kF64ScalarName[] = "scalar_f64";
const double kF64ScalarValue = 1234567890.0987654321;
const char kArrayNodeName[] = "arrays";
const char kBufferedArrayNodeName[] = "buffered arrays";
const char k1DArrayNodeName[] = "1D array";
const char k1DArrayPath[] = "/arrays/1D array";
const char k1DBufferedArrayPath[] = "/buffered arrays/1D array";
const char k2DArrayNodeName[] = "2D array";
#define k1DArraySize	10
const int k1DArrayValue[k1DArraySize] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
const char k2DArrayPath[] = "/arrays/2D array";
const char k2DBufferedArrayPath[] = "/buffered arrays/2D array";
const int k2DArraySize[2] = { 4, 5 };
#define k2DArraySize_Dim0	4
#define k2DArraySize_Dim1	5
const int k2DArrayValue[k2DArraySize_Dim0][k2DArraySize_Dim1] =
    { { 1, 2, 3, 4, 5 },
      { 6, 7, 8, 9, 10},
      { 11, 12, 13, 14, 15 },
      { 16, 17, 18, 19, 20 } };
const char kHierarchyNodeName[] = "hierarchy";
#define kHierarchyWidth		3
#define kHierarchyDepth		4
const char kHierarchyValue[kHierarchyDepth][kHierarchyWidth][16] = 
      { { "",		"grandfather",		"" },
	{ "father1",	"father2",		"father3" },
	{ "dead son",	"son2.1",		"son3.1" },
	{ "grandson",		"",		"" } };
struct {
    char name[16];
    int nChildren;
} kHierarchyTestTable[] =
{ { "grandfather", 3 }, { "father1", 0 }, { "father2", 1 }, { "father3", 1 },
  { "", 0 } };
const char kDeletedNodePath[] = "/hierarchy/grandfather/father1/dead son";
const char kDeletedChildPath[] = "/hierarchy/grandfather/father1/dead son/grandson";
const char kMovedNodeOriginalPath[] = "/hierarchy/grandfather/father2/father3";
const char kMovedNodePath[] = "/hierarchy/grandfather/father3";
const char kMovedChildPath[] = "/hierarchy/grandfather/father3/son3.1";
const int kBufferSize = 3;

int gErr = 0;

void PrintError( SFFError err, char *fmt, ... )
{
    va_list args;

    va_start(args, fmt);
    vprintf(fmt, args);
    va_end(args);
    printf("  (error %d)\n", (int)err);
}

void Test( SFFError err, int result, char *title )
{
    printf("Test '%s' ", title);
    if (err == kSFFNoErr && result)
	printf ("passed");
    else if (err == kSFFNoErr)
    {
	printf("FAILED (bad value)");
	gErr++;
    }
    else
    {
	printf("FAILED (error %d)", (int)err);
	gErr++;
    }
    printf(".\n");
}

int InHierarchy( char *name, int level )
{
    int i;

    for (i = 0;  i < kHierarchyWidth;  ++i)
    {
	if (strcmp(kHierarchyValue[level][i], name) == 0)
	    return(TRUE);
    }
    return(FALSE);
}

int HasRightChildren( char *name, SFFNode node )
{
    int i, n;

    for (i = 0;  kHierarchyTestTable[i].name[0] != '\0';  ++i)
    {
	if (strcmp(kHierarchyTestTable[i].name, name) == 0)
	{
	    n = -1;	/* In case an error occurs and n isn't set, we'll know*/
	    SFFGetNumberOfChildren(NULL, node, &n);
	    if (n != kHierarchyTestTable[i].nChildren)
		return(FALSE);
	    return(TRUE);
	}
	    
    }
    return(FALSE);
}

void ReadTest( const char *filename )
{
    char *title, *name;
    int x, y, version, int32, val, badVals, nth;
    double float64;
    SFFNode root, node, grandfather;
    SFFBuffer buff;
    SFFError err = kSFFNoErr;

    if ((err = SFFOpen(filename, &root)))
	Test(err, TRUE, "SFFOpenFile:  existing file");
    if (err)
	return;

    /* Test title and version.  (Title also tests 1D arrays) */
    SFFGetVersion(&err, root, &version);
    Test(err, version == kVersion, "SFFGetVersion");
    err=0; SFFGetTitle(&err, root, &title);
    Test(err, strcmp(kTitle, title) == 0, "SFFGetTitle");
    free(title);

    /* Test scalar data reads */
    err=0; SFFGetNode(&err, root, kI32ScalarName, &node);
    Test(err, TRUE, "GetNode kI32ScalarName");
    err=0; SFFReadData(&err, node, (void *)&int32, kSFFInt32, 1);
    Test(err, int32 == kI32ScalarValue, "Reading Int32 value");
    
    err=0; SFFGetNode(&err, root, kF64ScalarName, &node);
    Test(err, TRUE, "GetNode kF64ScalarName");
    err=0; SFFReadData(&err, node, (void *)&float64, kSFFFloat64, 1);
    Test(err, float64 == kF64ScalarValue, "Reading Float64 value");

    /* Test 1D array read */
    err=0; SFFGetNode(&err, root, k1DArrayPath, &node);
    printf("Reading 1D array:  ");
    badVals = 0;
    for (x = 0;  x < k1DArraySize;  ++x)
    {
	SFFReadDataPoint(&err, node, (void *)&val, x, kSFFEnd);
	if (val == k1DArrayValue[x])
	    printf(". ");
	else
	    printf("X ");
	badVals += val != k1DArrayValue[x];
    }
    printf("\n");
    Test(err, badVals == 0, "Reading 1D array");

    /* Test 2D array read */
    err=0; SFFGetNode(&err, root, k2DArrayPath, &node);
    printf("Reading 2D array:\n\t");
    badVals = 0;
    for (y = 0;  y < k2DArraySize_Dim1;  ++y)
    {
	for (x = 0;  x < k2DArraySize_Dim0;  ++x)
	{
	    SFFReadDataPoint(&err, node, (void *)&val, x, y, kSFFEnd);
	    if (val == k2DArrayValue[x][y])
		printf(". ");
	    else
		printf("X ");
	    badVals += val != k2DArrayValue[x][y];
	}
	printf("\n\t");
    }
    printf("\n");
    Test(err, badVals == 0, "Reading 2D array");

    /* Test 1D buffered array read */
    SFFSetBufferSize(3);
    err=0; SFFGetNode(&err, root, k1DBufferedArrayPath, &node);
    SFFCreateBuffer(&err, node, kSFFRead, &buff);
    printf("Reading buffered 1D array:  ");
    badVals = 0;
    for (x = 0;  x < k1DArraySize;  ++x)
    {
	SFFBufferReadDataPoint(&err, buff, (void *)&val, x, kSFFEnd);
	if (val == k1DArrayValue[x])
	    printf(". ");
	else
	    printf("X ");
	badVals += val != k1DArrayValue[x];
    }
    printf("\n");
    SFFDestroyBuffer(&err, buff);
    Test(err, badVals == 0, "Reading 1D array");

    /* Test 2D buffered array read */
    err=0; SFFGetNode(&err, root, k2DBufferedArrayPath, &node);
    SFFCreateBuffer(&err, node, kSFFRead, &buff);
    printf("Reading 2D array:\n\t");
    badVals = 0;
    for (y = 0;  y < k2DArraySize_Dim1;  ++y)
    {
	for (x = 0;  x < k2DArraySize_Dim0;  ++x)
	{
	    SFFBufferReadDataPoint(&err, buff, (void *)&val, x, y, kSFFEnd);
	    if (val == k2DArrayValue[x][y])
		printf(". ");
	    else
		printf("X ");
	    badVals += val != k2DArrayValue[x][y];
	}
	printf("\n\t");
    }
    printf("\n");
    SFFDestroyBuffer(&err, buff);
    Test(err, badVals == 0, "Reading 2D array");

    /* Test deletion */
    err=0;  SFFGetNode(&err, root, kDeletedNodePath, &node);
    Test(kSFFNoErr, err == kSFFNoNodeErr, "SFFDeleteNode verification");
    err=0;  SFFGetNode(&err, root, kDeletedChildPath, &node);
    Test(kSFFNoErr, err == kSFFNoNodeErr, "SFFDeleteNode hierarchy deletion verification");

    /* Test moved node */
    err=0;  SFFGetNode(&err, root, kMovedNodePath, &node);
    Test(err, TRUE, "SFFMoveNode verification");
    err=0;  SFFGetNode(&err, root, kMovedChildPath, &node);
    Test(err, TRUE, "SFFMoveNode hierarchy move verification");

    /* Test hierarchy structure */
    err = 0;  badVals = 0;
    SFFGetNode(&err, root, kHierarchyNodeName, &node);
    SFFGetNumberOfChildren(&err, node, &val);
    badVals += !(val == 1);	/* Top-level hierarchy is only grandfather */
    nth = 0;
    SFFGetNextChild(&err, node, &nth, &grandfather);
    SFFGetName(&err, grandfather, &name);
    badVals += !InHierarchy(name, 0);
    SFFGetNumberOfChildren(&err, grandfather, &val);
    badVals += !HasRightChildren(name, grandfather);
    free(name);
    for (x = 0; SFFGetNextChild(NULL, grandfather, &x, &node) == kSFFNoErr; ++x)
    {
	SFFGetName(&err, node, &name);
	badVals += !InHierarchy(name, 1);
	badVals += !HasRightChildren(name, node);
	free(name);
    }
    Test(err, badVals == 0, "Hierarchy verification");
    
    err = SFFClose(root);
    Test(err, TRUE, "SFFCloseFile:  existing file");
}

void WriteTest( const char *filename )
{
    int x, y;
    SFFNode root, node, aryParent, parent, son, grandson;
    SFFBuffer buff;
    SFFError err = kSFFNoErr;

    err = SFFOpen(filename, &root);
    Test(err, TRUE, "SFFOpenFile");

    err=0; SFFSetVersion(&err, root, kVersion);
    Test(err, TRUE, "SFFSetVersion");
    err=0; SFFSetTitle(&err, root, kTitle);
    Test(err, TRUE, "SFFSetTitle");

    /* Write scalar values */
    err=0; SFFCreateNode(&err, root, FALSE, kI32ScalarName, kI32ScalarName,
			 &node);
    Test(err, TRUE, "SFFCreateNode(kI32ScalarName)");
    SFFSetDataType(&err, node, kSFFInt32, 1, kSFFEnd);
    Test(err, TRUE, "SFFSetDataType(kSFFInt32, 1)");
    SFFWriteData(&err, node, (void *)&kI32ScalarValue);
    Test(err, TRUE, "SFFWriteDataPoint(kI32ScalarValue, 1)");

    err=0; SFFCreateNode(&err, root, FALSE, kF64ScalarName, kF64ScalarName,
			 &node);
    Test(err, TRUE, "SFFCreateNode(kF64ScalarName)");
    SFFSetDataType(&err, node, kSFFFloat64, 1, kSFFEnd);
    Test(err, TRUE, "SFFSetDataType(kSFFFloat64, 1)");
    SFFWriteData(&err, node, (void *)&kF64ScalarValue);
    Test(err, TRUE, "SFFWriteDataPoint(kF64ScalarValue, 1)");

    /* Write array nodes */
    err=0; SFFCreateNode(&err, root, FALSE, kArrayNodeName, kArrayNodeName,
			 &aryParent);
    Test(err, TRUE, "SFFCreateNode(kArrayNodeName)");

    /* 1D array */
    err=0; SFFCreateNode(&err, aryParent, FALSE, k1DArrayNodeName,
			 k1DArrayNodeName, &node);
    SFFSetDataType(&err, node, kSFFInt32, k1DArraySize, kSFFEnd);
    for (x = 0;  x < k1DArraySize;  ++x)
	SFFWriteDataPoint(&err, node, (void *)&k1DArrayValue[x], x, kSFFEnd);
    Test(err, TRUE, "Create 1D array");

    /* 2D array */
    err=0; SFFCreateNode(&err, aryParent, FALSE, k2DArrayNodeName,
			 k2DArrayNodeName, &node);
    SFFSetDataType(&err, node, kSFFInt32, k2DArraySize_Dim0, k2DArraySize_Dim1,
		   kSFFEnd);
    for (x = 0;  x < k2DArraySize_Dim0;  ++x)
    {
	for (y = 0;  y < k2DArraySize_Dim1;  ++y)
	    SFFWriteDataPoint(&err, node, (void *)&k2DArrayValue[x][y],
			      x, y, kSFFEnd);
    }
    Test(err, TRUE, "Create 2D array");

    /* 1D buffered array */
    err=0; SFFCreateNode(&err, root, FALSE, kBufferedArrayNodeName,
			 kBufferedArrayNodeName, &aryParent);
    SFFCreateNode(&err, aryParent, FALSE, k1DArrayNodeName, k1DArrayNodeName,
		  &node);
    SFFSetDataType(&err, node, kSFFInt32, k1DArraySize, kSFFEnd);
    SFFCreateBuffer(&err, node, kSFFWrite, &buff);
    for (x = 0;  x < k1DArraySize;  ++x)
	SFFBufferWriteDataPoint(&err, buff, (void *)&k1DArrayValue[x],
				x, kSFFEnd);
    SFFDestroyBuffer(&err, buff);
    Test(err, TRUE, "Create 1D buffered array");

    /* 2D buffered array */
    err=0; SFFCreateNode(&err, aryParent, FALSE, k2DArrayNodeName,
			 k2DArrayNodeName, &node);
    SFFSetDataType(&err, node, kSFFInt32, k2DArraySize_Dim0, k2DArraySize_Dim1,
		   kSFFEnd);
    SFFCreateBuffer(&err, node, kSFFWrite, &buff);
    for (x = 0;  x < k2DArraySize_Dim0;  ++x)
    {
	for (y = 0;  y < k2DArraySize_Dim1;  ++y)
	    SFFBufferWriteDataPoint(&err, buff, (void *)&k2DArrayValue[x][y],
			      x, y, kSFFEnd);
    }
    SFFDestroyBuffer(&err, buff);
    Test(err, TRUE, "Create 2D buffered array");

    /* Create hierarchy test */
    err=0;
    SFFCreateNode(&err, root, FALSE, kHierarchyNodeName, NULL, &node);
    SFFCreateNode(&err, node, FALSE, kHierarchyValue[0][1], NULL, &parent);
    SFFCreateNode(&err, parent, FALSE, kHierarchyValue[1][0], NULL, &son);
    SFFCreateNode(&err, son, FALSE, kHierarchyValue[2][0], NULL, &grandson);
    SFFCreateNode(&err, grandson, FALSE, kHierarchyValue[3][0], NULL, NULL);
    SFFCreateNode(&err, parent, FALSE, kHierarchyValue[1][1], NULL, &son);
    SFFCreateNode(&err, son, FALSE, kHierarchyValue[2][1], NULL, NULL);
    SFFCreateNode(&err, son, FALSE, kHierarchyValue[1][2], NULL, &grandson);
    SFFCreateNode(&err, grandson, FALSE, kHierarchyValue[2][2], NULL, NULL);
    Test(err, TRUE, "Hierarchy creation (before move and delete)");
    SFFGetNode(&err, root, kDeletedNodePath, &node);
    err=0;  SFFDeleteNode(&err, node);
    Test(err, TRUE, "SFFDeleteNode");
    err=0;  SFFGetNode(&err, root, kMovedNodeOriginalPath, &node);
    SFFMoveNode(&err, node, parent);
    Test(err, TRUE, "SFFMoveNode");

    err = SFFClose(root);
    Test(err, TRUE, "SFFCloseFile");
}

int main(int argc, char *argv[])
{
    if (remove(kFilename) != 0)
    	{ printf("Could not remove '%s'!\n", kFilename);  return(1); }

    WriteTest(kFilename);
    ReadTest(kFilename);

    if (gErr != 1)
	printf("There were %d errors.", gErr);
    else
	printf("There was 1 error.");
    printf("  Tests ");
    if (!gErr)
	printf("passed");
    else
	printf("FAILED");
    printf(".\n");

    return(gErr);
}


/* Automatic setting of emacs local variables. */
/* Local Variables: */
/* mode: C++ */
/* tab-width: 8 */
/* End: */

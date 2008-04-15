/*@@
 *  Program: Star File Format Library  - $RCSfile: ccmioread.c,v $
 *  Author:  Geoff Prewett
 *  Date:    January 5, 2004
 *
 *
 *  Star File Format Library - Copyright (C) 2004 by adapco, Ltd.
 *
 *  This program is the property of adapco, Ltd. and contains
 *  confidential and proprietary information.  The unauthorized use,
 *  distribution, or duplication of this program is prohibited.
 *  All rights reserved.
 *
 *  $Id: ccmioread.c,v 1.18 2006/06/30 16:58:44 geoffp Exp $
 */

/* This file is not intended to be compiled directly.  Instead it
   implements a "template" function. The macros TYPE, TYPE_ABBREV, and CCMIOTYPE
   must be defined before inclusion and this file must be included for
   each type of "template" that will be instantiated.  A poor man's template
   function. */

/* ------------ Global section (not templated) ------------------------------ */
#ifndef CCMIO_READ_GLOBAL_DEF
#define CCMIO_READ_GLOBAL_DEF
typedef void (*ADF_RW_All_Data_Func)(const double, char*, int *);
typedef void (*ADF_RW_Block_Data_Func)(const double, const long, const long,
				       char*, int *);

CCMIOError CCMIOExtendedADFIO( CCMIOError *err, CCMIONode node, 
			       CCMIOIOType ioType, CCMIODataType dataType,
			       int nDims, CCMIOIndex const *dims,
			       char *data, CCMIOIndex start, CCMIOIndex end )
{
    int i, thisNodeNDims, thisNodeDims[ADF_MAX_DIMENSIONS];
    CCMIOIndex nItems, itemSize = 1, idx2bytes, itemsPerNode;
    CCMIOSize dataTypeSize;
    ADFError adfErr = NO_ERROR;
    ADF_RW_All_Data_Func ADF_RW_All_Data = &ADF_Read_All_Data;
    ADF_RW_Block_Data_Func ADF_RW_Block_Data = &ADF_Read_Block_Data;

    CHECK_ERROR(err);

    if (start == kCCMIOStart && ioType == kCCMIOWrite)
	if (CCMIOSetDataTypev(err, node, dataType, nDims, dims) != kCCMIONoErr)
	    return(*err);

    /* start and end are in conceptual units.  For example, for a 2D array
       verts[100][3], there are 100 units (each of size three), but map[100]
       also has 100 units (size one).  This is complicated a bit because
       ADF is in FORTRAN order, so the order of the dimensions is reversed
       from C. */
    nItems = dims[nDims - 1];
    for (i = 0;  i < nDims - 1;  ++i)
	itemSize *= dims[i];

    /* start and end are 0-based indices.  ADF uses 1-based indices
       (i.e. FORTRAN indexing), so they will need to be converted before
       we actually do I/O. */
    if (end == 0 || end > nItems)
	end = nItems;
    if (start > end)
	return(*err = kCCMIOBadParameterErr);

    if (ioType == kCCMIOWrite)
    {
	ADF_RW_All_Data = (ADF_RW_All_Data_Func)&ADF_Write_All_Data;
	ADF_RW_Block_Data = &ADF_Write_Block_Data;
    }
    
    dataTypeSize = CCMIOGetDataTypeSize(dataType);
    idx2bytes = dataTypeSize;
    if (nDims > 1)
	idx2bytes *= dims[0];
    if (ioType == kCCMIOWrite)
	itemsPerNode = kCCMIOMaxADFNodeSize / idx2bytes;
    else
    {
	/* We need to read the dimensions of the main node to find
	   the items per node, in case a different version of the library
	   used a different value for kCCMIOMaxADFNodeSize */
	if (GetADFNodeDimensions(NULL, node, &thisNodeNDims, thisNodeDims)
	    							 != kCCMIONoErr)
	    return(*err = kCCMIOCorruptFileErr);
	itemsPerNode = thisNodeDims[nDims - 1];
    }

    if (end <= itemsPerNode)
    {
	if (start == 0 && (end == nItems || end == itemsPerNode))
	    (*ADF_RW_All_Data)(node.node, data, &adfErr);
	else
	    (*ADF_RW_Block_Data)(node.node, (long)start * itemSize + 1,
				 (long)end * itemSize,
				 data, &adfErr);
    }
    else
    {
	int i;
	char nodeExName[kCCMIOMaxStringLength + 1];
	CCMIOIndex dimsCopy[ADF_MAX_DIMENSIONS];
	CCMIONode nodeEx;

	/* Note:  this check must go here instead of after the assignment to
	   itemsPerNode because an empty string needs to have nDims = 1 and
	   dims[0] = 0, which would return an incorrect error.  If we put
	   it here, an empty string cannot reach the check. */
	if (itemsPerNode == 0)
	    return (*err = kCCMIOArrayDimensionToLargeErr);

	for (i = 0;  i < nDims;  ++i)
	    dimsCopy[i] = dims[i];

	while (start < end)
	{
	    CCMIOIndex startNodeNum = start / itemsPerNode;
	    CCMIOIndex endNodeNum = (end - 1) / itemsPerNode;
	    CCMIOIndex whereDoesThisNodeStart =
		    			  (start / itemsPerNode) * itemsPerNode;

	    /* Figure out which node we need to be reading from */
	    if (start < itemsPerNode)
		nodeEx = node;
	    else
	    {
		snprintf(nodeExName, kCCMIOMaxStringLength, "%s-%u",
			 kCCMIOExtendedDataName, startNodeNum);
		if (ioType == kCCMIOWrite)
		{
		    CCMIOCreateNode(err, node, TRUE, nodeExName,
				    kCCMIOExtendedDataLabel, &nodeEx);

		    if (dims[nDims-1] - whereDoesThisNodeStart < itemsPerNode)
			dimsCopy[nDims - 1] = dims[nDims-1] - whereDoesThisNodeStart;
		    else
			dimsCopy[nDims - 1] = itemsPerNode;
		    CCMIOSetDataTypev(err, nodeEx, dataType, nDims, dimsCopy);
		    if (*err !=kCCMIONoErr)
			return(*err);
		}
		else
		{
		    if (CCMIOGetNode(NULL, node, nodeExName, &nodeEx)
								  !=kCCMIONoErr)
			return(*err = kCCMIOCorruptFileErr);
		}
	    }

	    /* Read/write the next section of data */
	    if (startNodeNum == endNodeNum && end - start >= itemsPerNode)
	    {
		(*ADF_RW_All_Data)(nodeEx.node, data, &adfErr);
		data += (end - start + 1) * idx2bytes;  /*data += size of node*/
		start += (end - start + 1) * dims[nDims - 1];
	    }
	    else
	    {
		CCMIOIndex thisNodeOffsetIdx = start - whereDoesThisNodeStart+1;
		CCMIOIndex thisNodeEndIdx;

		if (end - start >= itemsPerNode || startNodeNum != endNodeNum)
		    /* Note that the starting index is 1 (ADF uses Fortran idx)
		       so we don't want to subtract 1 here... */
		    thisNodeEndIdx = (startNodeNum + 1) * itemsPerNode;
		else
		    thisNodeEndIdx = end;
		thisNodeEndIdx = thisNodeEndIdx - whereDoesThisNodeStart;

		if (thisNodeOffsetIdx <= 1 && thisNodeEndIdx == itemsPerNode)
		    (*ADF_RW_All_Data)(nodeEx.node, data, &adfErr);
		else
		    (*ADF_RW_Block_Data)(nodeEx.node,
				((long)(thisNodeOffsetIdx-1)*itemSize)+1,
				(long)thisNodeEndIdx * itemSize,
					 data, &adfErr);
		if (IsADFError(adfErr))
		    return(*err = ADFToCCMIOError(adfErr));

		start += thisNodeEndIdx - thisNodeOffsetIdx + 1;
		data += (thisNodeEndIdx - thisNodeOffsetIdx + 1) * idx2bytes;
	    }
	}
    }
    return(*err = ADFToCCMIOError(adfErr));
}
#endif /* CCMIO_READ_GLOBAL_DEF */
/* -------------------- end global section ---------------------------------- */

/* Need to use nested macros because concatenation (and stringification)
   do not expand their arguments.  But STRCAT does, so then it can pass it
   onto MAKE_NAME to get concatenated. */
#define STRCAT(X, Y)		MAKE_NAME(X, Y)
#define MAKE_NAME(X, Y)		X##Y

#define CONVERT		STRCAT(ConvertTo,TYPE)
TYPE CONVERT( const void *data, CCMIODataType type )
{
    if (type == kCCMIOInt32)
	return((TYPE)(*(int *)data));
    if (type == kCCMIOFloat32)
	return((TYPE)(*(float *)data));
    if (type == kCCMIOFloat64)
	return((TYPE)(*(double *)data));
    return((TYPE)0xbad0bad0);
}

#define OLDREAD STRCAT(CCMIOOldRead,TYPE_ABBREV)
CCMIOError OLDREAD( CCMIOError *err, CCMIONode node, int dimension,
		    int swapDims, TYPE *data, CCMIOIndex start,
		    CCMIOIndex end )
{
    CCMIOIndex size;
    int nDims, adfErr;
    CCMIOIndex i, j, k, *dims;
    CCMIODataType nodeType;
    CHECK_ERROR(err);

    CCMIOGetDimensions(err, node, &nDims, &dims);
    if (nDims != dimension)
    {
	free(dims);
	return(*err = kCCMIOWrongDataTypeErr);
    }
    CCMIOGetDataType(err, node, &nodeType);
    CCMIOGetDataSize(err, node, &size);

    if (swapDims)
    {
	if (end == 0 || end > (unsigned int)dims[nDims - 1])
	    end = dims[0];
    }
    else
    {
	if (end == 0 || end > (unsigned int)dims[nDims - 1])
	    end = dims[nDims - 1];
    }

    /* This needs to be fixed if old files are going to be supported. */
    {
        char *buffer = (char*) malloc(size);
        unsigned int typeSize = CCMIOGetDataTypeSize(nodeType);
        TYPE *cBuffer;
        if (!buffer)
        {
	    free(dims);
	    return(*err = kCCMIONoMemoryErr);
        }
        /*You'd think it wouldn't be necessary to allocate an entire new buffer.
          Unfortunately, the original algorithm is completely messed up and
          requires it. */
        cBuffer = (TYPE *)malloc(size / typeSize * sizeof(TYPE));
        if (!cBuffer)
        {
	    free(dims);  free(buffer);
	    return(*err = kCCMIONoMemoryErr);
        }

        ADF_Read_All_Data(node.node, (char *) buffer, &adfErr);
        if (!IsADFError(adfErr))
        {
	    int l;  /* This needs to be signed so that the loop terminates */
	    unsigned int delta, itemSize = 1;

	    if (swapDims)
	    {
	        for (i = 1;  i < (unsigned int)nDims;  ++i)
		    itemSize *= dims[i];
	    }
	    else
	    {
	        for (l = (unsigned int)nDims - 2;  l >= 0;  --l)
		    itemSize *= dims[l];
	    }
            delta = itemSize * (end - start);
            
	    if (nDims == 2)
	    {
	        for (i = 0;  i < dims[0];  ++i)
		    for (j = 0;  j < dims[1]; ++j)
		        cBuffer[i*dims[1] + j] = CONVERT(buffer + typeSize * (dims[0] * j + i), nodeType);
	    }
	    else if (nDims == 3)
	    {
	        for (i = 0;  i < dims[0];  ++i)
		    for (j = 0;  j < dims[1]; ++j)
		        for (k = 0;  k < dims[2]; ++k)
			    cBuffer[i + j * dims[0] + k * dims[1] * dims[2]] = CONVERT(buffer + typeSize*(i*dims[1]*dims[2] + j*dims[0] + k), nodeType);
	    }

	    for (i = 0, j = itemSize*start;  i < delta;  ++i, ++j)
	        data[i] = cBuffer[j];

        }
        free(buffer);
        free(cBuffer);
    }
    free(dims);
    return(*err = ADFToCCMIOError(adfErr));
}

#define READ STRCAT(CCMIORead,TYPE_ABBREV)
CCMIOError READ( CCMIOError *err, CCMIONode node, int dimension, TYPE *data,
		 CCMIOIndex start, CCMIOIndex end )
{
    unsigned int nItems, itemSize = 1;
    int nDims;
    CCMIOIndex *dims;
    CCMIODataType nodeType;

    CHECK_ERROR(err);

    CCMIOGetDimensions(err, node, &nDims, &dims);
    if (nDims != dimension)
    {
	free(dims);
	return(*err = kCCMIOWrongDataTypeErr);
    }

    CCMIOGetDataType(err, node, &nodeType);
    if (nodeType == CCMIOTYPE)
	CCMIOExtendedADFIO(err, node, kCCMIORead, CCMIOTYPE, nDims, dims,
			   (char *)data, start, end);
    else
    {
	void *buffer;
	unsigned int i, n, typeSize = CCMIOGetDataTypeSize(nodeType);

	/* start and end are in conceptual units.  For example, for a 2D array
	   verts[100][3], there are 100 units (each of size three), but map[100]
	   also has 100 units (size one).  This is complicated a bit because
	   ADF is in FORTRAN order, so the order of the dimensions is reversed
	   from C. */
	nItems = dims[nDims - 1];
	for (i = 0; (int)i < nDims - 1;  ++i)
	    itemSize *= dims[i];

	if (end == 0 || end > nItems)
	    end = nItems;

	if (start > end)
	{
	    free(dims);
	    return(*err = kCCMIOBadParameterErr);
	}

	buffer = malloc((end - start + 2) * itemSize * typeSize);
	if (!buffer)
	{
	    free(dims);
	    return(*err = kCCMIONoMemoryErr);
	}
	if (CCMIOExtendedADFIO(err, node, kCCMIORead, nodeType, nDims, dims,
			       (char *)buffer, start, end) == kCCMIONoErr)
	{
	    n = (end - start) * itemSize;
	    for (i = 0;  i < n;  ++i)
		data[i] = CONVERT((char *)buffer + i * typeSize, nodeType);
	}
	free(buffer);
    }

    free(dims);
    return(*err);
}

#define WRITE STRCAT(CCMIOWrite,TYPE_ABBREV)
CCMIOError WRITE( CCMIOError *err, CCMIONode node, int nDims,
		  const CCMIOIndex *dims, const TYPE *data,
		  CCMIOIndex start, CCMIOIndex end )
{
    CHECK_ERROR(err);

    CCMIOExtendedADFIO(err, node, kCCMIOWrite, CCMIOTYPE, nDims, dims,
		       (char *)data, start, end);

    return(*err);
}

#define READ_1D STRCAT(CCMIORead1,TYPE_ABBREV)
CCMIOError READ_1D( CCMIOError *err, CCMIONode node, TYPE *data,
		    CCMIOIndex start, CCMIOIndex end )
{
    return(READ(err, node, 1, data, start, end));
}

#define WRITE_1D STRCAT(CCMIOWrite1,TYPE_ABBREV)
CCMIOError WRITE_1D( CCMIOError *err, CCMIONode node, CCMIOIndex n,
		     const TYPE *data, CCMIOIndex start, CCMIOIndex end )
{
/*  long dims[] = { n };     Sun compilers don't like this code */
    CCMIOIndex dims[1];
    dims[0] = n;

    return(WRITE(err, node, 1, dims, data, start, end));
}

#define READ_2D STRCAT(CCMIORead2,TYPE_ABBREV)
CCMIOError READ_2D( CCMIOError *err, CCMIONode node, TYPE *data,
		    CCMIOIndex start, CCMIOIndex end )
{
    return(READ(err, node, 2, data, start, end));
}

#define WRITE_2D STRCAT(CCMIOWrite2,TYPE_ABBREV)
CCMIOError WRITE_2D( CCMIOError *err, CCMIONode node,
		     CCMIOIndex x, CCMIOIndex y,
		     const TYPE *data, CCMIOIndex start, CCMIOIndex end )
{
/*  long dims[] = { x, y };	Sun compilers don't like this */
    CCMIOIndex dims[2];
    dims[0] = x;
    dims[1] = y;

    return(WRITE(err, node, 2, dims, data, start, end));
}

#define READ_3D STRCAT(CCMIORead3,TYPE_ABBREV)
CCMIOError READ_3D( CCMIOError *err, CCMIONode node, TYPE *data,
		    CCMIOIndex start, CCMIOIndex end )
{
    return(READ(err, node, 3, data, start, end));
}

#define WRITE_3D STRCAT(CCMIOWrite3,TYPE_ABBREV)
CCMIOError WRITE_3D( CCMIOError *err, CCMIONode node,
		     CCMIOIndex x, CCMIOIndex y, CCMIOIndex z,
		     const TYPE *data, CCMIOIndex start, CCMIOIndex end )
{
/*  long dims[] = { x, y, z };		Sun compilers do not like this */
    CCMIOIndex dims[3];
    dims[0] = x;
    dims[1] = y;
    dims[2] = z;

    return(WRITE(err, node, 3, dims, data, start, end));
}

#define READOPT STRCAT(CCMIOReadOpt,TYPE_ABBREV)
CCMIOError READOPT( CCMIOError *err, CCMIOID parent, const char *name,
		    TYPE *data )
{
    char buffer[64];	/* This ought to be big enough for any type */
    int adfErr;
    CCMIONode node;
    CCMIODataType type;
    CHECK_ERROR_AND_CLEAR_PTR(err, data, 0);

    CCMIOGetNode(err, parent.node, name, &node);
    CCMIOGetDataType(err, node, &type);
    ADF_Read_Block_Data(node.node, 1, 1, buffer, &adfErr);
    *err = ADFToCCMIOError(adfErr);
    if (*err == kCCMIONoErr)
	*data = CONVERT(buffer, type);

    return(*err);
}

#define READOPT_1D STRCAT(CCMIOReadOpt1,TYPE_ABBREV)
CCMIOError READOPT_1D( CCMIOError *err, CCMIOID parent, const char *name,
		       TYPE *data, CCMIOIndex start, CCMIOIndex end )
{
    CCMIONode node;
    CHECK_ERROR(err);
    if (!data)  return(*err = kCCMIOBadParameterErr);

    CCMIOGetNode(err, parent.node, name, &node);
    return(READ_1D(err, node, data, start, end));
}

#define READOPT_2D STRCAT(CCMIOReadOpt2,TYPE_ABBREV)
CCMIOError READOPT_2D( CCMIOError *err, CCMIOID parent, const char *name,
		       TYPE *data, CCMIOIndex start, CCMIOIndex end )
{
    CCMIONode node;
    CHECK_ERROR(err);
    if (!data)  return(*err = kCCMIOBadParameterErr);

    CCMIOGetNode(err, parent.node, name, &node);
    if (parent.version < 20300)
    {
	int nDims, swap = FALSE;
	CCMIOIndex *dims;
	CCMIOGetDimensions(err, node, &nDims, &dims);
	if (nDims == 2 && dims[0] > dims[1])
	    swap = TRUE;
	free(dims);
	return(OLDREAD(err, node, 2, swap, data, start, end));
    }
    return(READ_2D(err, node, data, start, end));
}

#define READOPT_3D STRCAT(CCMIOReadOpt3,TYPE_ABBREV)
CCMIOError READOPT_3D( CCMIOError *err, CCMIOID parent, const char *name,
		       TYPE *data, CCMIOIndex start, CCMIOIndex end )
{
    CCMIONode node;
    CHECK_ERROR(err);
    if (!data)  return(*err = kCCMIOBadParameterErr);

    CCMIOGetNode(err, parent.node, name, &node);
    if (parent.version < 20300)
    {
	int nDims, swap = FALSE;
	CCMIOIndex *dims;
	CCMIOGetDimensions(err, node, &nDims, &dims);
	if (nDims == 2 && (dims[0] > dims[1] || dims[0] > dims[2]))
	    swap = TRUE;
	free(dims);
	return(OLDREAD(err, node, 2, swap, data, start, end));
    }
    return(READ_3D(err, node, data, start, end));
}

#define WRITEOPT_1D STRCAT(CCMIOWriteOpt1,TYPE_ABBREV)
CCMIOError WRITEOPT_1D( CCMIOError *err, const CCMIOID parent,
                        const char *name, const CCMIOIndex n, const TYPE *data,
			const CCMIOIndex start, const CCMIOIndex end )
{
    static CCMIONode node;
    CHECK_ERROR(err);
    if (!data)  return(*err = kCCMIOBadParameterErr);

    if (start == 0)
        CCMIOCreateNode(err, parent.node, TRUE, name, name, &node);

    return(WRITE_1D(err, node, n, data, start, end));
}
#define WRITEOPT_2D STRCAT(CCMIOWriteOpt2,TYPE_ABBREV)
CCMIOError WRITEOPT_2D( CCMIOError *err, const CCMIOID parent,
                        const char *name, 
			const CCMIOIndex x, const CCMIOIndex y,
			const TYPE *data,
			const CCMIOIndex start, const CCMIOIndex end )
{
    static CCMIONode node;
    CHECK_ERROR(err);
    if (!data)  return(*err = kCCMIOBadParameterErr);

    if (start == 0)
        CCMIOCreateNode(err, parent.node, TRUE, name, name, &node);

    return(WRITE_2D(err, node, y, x, data, start, end));
}

#define WRITEOPT_3D STRCAT(CCMIOWriteOpt3,TYPE_ABBREV)
CCMIOError WRITEOPT_3D( CCMIOError *err, CCMIOID parent, const char *name,
			CCMIOIndex x, CCMIOIndex y, CCMIOIndex z,
			const TYPE *data,
			CCMIOIndex start, CCMIOIndex end )
{
    CCMIONode node;
    CHECK_ERROR(err);
    if (!data)  return(*err = kCCMIOBadParameterErr);

    CCMIOCreateNode(err, parent.node, TRUE, name, name, &node);
    return(WRITE_3D(err, node, z, y, x, data, start, end));
}


/* Automatic setting of emacs local variables. */
/* Local Variables: */
/* mode: C++ */
/* tab-width: 8 */
/* End: */

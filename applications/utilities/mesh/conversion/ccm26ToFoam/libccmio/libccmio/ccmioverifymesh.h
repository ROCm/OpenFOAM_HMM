#ifndef CCMIOMESH_VERIFY_H
#define CCMIOMESH_VERIFY_H

/*@@
 *  Program: Star File Format Library  - $RCSfile: ccmioverifymesh.h,v $
 *  Author:  Geoff Prewett
 *  Date:    Sept 18, 2003
 *
 *
 *  Star File Format Library - Copyright (C) 2003 by adapco, Ltd.
 *
 *  This program is the property of adapco, Ltd. and contains
 *  confidential and proprietary information.  The unauthorized use,
 *  distribution, or duplication of this program is prohibited.
 *  All rights reserved.
 *
 *  $Id: ccmioverifymesh.h,v 1.3 2004/10/25 19:11:55 wayne Exp $
 */

#ifdef __cplusplus
extern "C" {
#endif

/** \ingroup mesh
    Should be called only after all data has been written to the processor.
    Sanity checks the mesh for consistency.  Displays all errors in the mesh.*/
extern CCMIOError CCMIOCheckProcessorMesh( CCMIOError *err, CCMIOProcessor proc );

#ifdef __cplusplus
}
#endif
#endif /* CCMIOMESH_VERIFY_H */


/* Automatic setting of emacs local variables. */
/* Local Variables: */
/* mode: C++ */
/* tab-width: 8 */
/* End: */

#ifndef KAHIP_H
#define KAHIP_H

/* *** DUMMY VERSION of kaHIP_interface.h - this file should not be included if you have KaHIP
 *     installed in the correct position in $WM_THIRD_PARTY_DIR - see
 *     decompositionMethods/kahipDecomp/Make/options
 */

#warning "Dummy kahip.h - included since it cannot find KaHIP installation."

#ifdef __cplusplus
extern "C" {
#endif

void kaffpa
(
    // [inputs]
    int* n,
    int* vwgt,
    int* xadj,
    int* adjcwgt,
    int* adjncy,
    int* nparts,
    double* imbalance,
    bool suppress_output,
    int seed,
    int mode,
    // [outputs]
    int* edgecut,
    int* part
);

#ifdef __cplusplus
}
#endif

#endif

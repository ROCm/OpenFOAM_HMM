#include "readTimeControls.H"

correctPhi = pimple.dict().lookupOrDefault("correctPhi", false);

checkMeshCourantNo = pimple.dict().lookupOrDefault("checkMeshCourantNo", false);

massFluxInterpolation =
    pimple.dict().lookupOrDefault("massFluxInterpolation", false);

ddtCorr = pimple.dict().lookupOrDefault("ddtCorr", true);

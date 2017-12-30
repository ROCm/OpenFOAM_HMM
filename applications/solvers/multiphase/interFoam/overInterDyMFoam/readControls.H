#include "readTimeControls.H"

correctPhi = pimple.dict().lookupOrDefault<Switch>("correctPhi", false);

checkMeshCourantNo =
    pimple.dict().lookupOrDefault<Switch>("checkMeshCourantNo", false);

moveMeshOuterCorrectors =
    pimple.dict().lookupOrDefault<Switch>("moveMeshOuterCorrectors", false);

massFluxInterpolation =
    pimple.dict().lookupOrDefault("massFluxInterpolation", false);

ddtCorr = pimple.dict().lookupOrDefault("ddtCorr", true);

adjustFringe =  pimple.dict().lookupOrDefault("oversetAdjustPhi", false);

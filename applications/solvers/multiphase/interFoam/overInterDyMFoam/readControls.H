#include "readTimeControls.H"

correctPhi = pimple.dict().getOrDefault("correctPhi", false);

checkMeshCourantNo =
    pimple.dict().getOrDefault("checkMeshCourantNo", false);

moveMeshOuterCorrectors =
    pimple.dict().getOrDefault("moveMeshOuterCorrectors", false);

massFluxInterpolation =
    pimple.dict().getOrDefault("massFluxInterpolation", false);

ddtCorr = pimple.dict().getOrDefault("ddtCorr", true);

adjustFringe =  pimple.dict().getOrDefault("oversetAdjustPhi", false);

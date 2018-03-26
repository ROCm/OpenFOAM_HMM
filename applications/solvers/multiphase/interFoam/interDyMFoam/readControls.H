#include "readTimeControls.H"

correctPhi = pimple.dict().lookupOrDefault("correctPhi", true);

checkMeshCourantNo =
    pimple.dict().lookupOrDefault("checkMeshCourantNo", false);

moveMeshOuterCorrectors =
    pimple.dict().lookupOrDefault("moveMeshOuterCorrectors", false);

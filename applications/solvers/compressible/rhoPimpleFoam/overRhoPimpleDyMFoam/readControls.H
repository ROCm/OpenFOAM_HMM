#include "readTimeControls.H"

correctPhi = pimple.dict().getOrDefault("correctPhi", false);

checkMeshCourantNo =
    pimple.dict().getOrDefault("checkMeshCourantNo", false);


ddtCorr = pimple.dict().getOrDefault("ddtCorr", true);

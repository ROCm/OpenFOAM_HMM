#include "createTimeControls.H"

bool correctPhi
(
    pimple.dict().lookupOrDefault("correctPhi", false)
);

bool checkMeshCourantNo
(
    pimple.dict().lookupOrDefault("checkMeshCourantNo", false)
);

bool massFluxInterpolation
(
    pimple.dict().lookupOrDefault("massFluxInterpolation", false)
);

bool adjustFringe
(
    pimple.dict().lookupOrDefault("oversetAdjustPhi", false)
);

bool ddtCorr
(
    pimple.dict().lookupOrDefault("ddtCorr", true)
);

#include "createTimeControls.H"

bool correctPhi
(
    pimple.dict().getOrDefault("correctPhi", false)
);

bool checkMeshCourantNo
(
    pimple.dict().getOrDefault("checkMeshCourantNo", false)
);

bool massFluxInterpolation
(
    pimple.dict().getOrDefault("massFluxInterpolation", false)
);

bool adjustFringe
(
    pimple.dict().getOrDefault("oversetAdjustPhi", false)
);

bool ddtCorr
(
    pimple.dict().getOrDefault("ddtCorr", true)
);

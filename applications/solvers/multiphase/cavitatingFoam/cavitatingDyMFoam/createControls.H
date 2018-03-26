#include "createTimeControls.H"

bool correctPhi
(
    pimple.dict().lookupOrDefault("correctPhi", true)
);

scalar maxAcousticCo
(
    readScalar(runTime.controlDict().lookup("maxAcousticCo"))
);

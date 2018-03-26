#include "readTimeControls.H"

correctPhi = pimple.dict().lookupOrDefault("correctPhi", true);
maxAcousticCo = readScalar(runTime.controlDict().lookup("maxAcousticCo"));

#include "readDyMControls.H"

correctPhi = pimple.dict().lookupOrDefault("correctPhi", true);
maxAcousticCo = readScalar(runTime.controlDict().lookup("maxAcousticCo"));

#include "readDyMControls.H"

correctPhi = pimple.dict().lookupOrDefault("correctPhi", true);
maxAcousticCo = runTime.controlDict().get<scalar>("maxAcousticCo");

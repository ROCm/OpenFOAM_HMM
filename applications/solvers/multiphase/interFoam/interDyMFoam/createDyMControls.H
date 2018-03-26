bool correctPhi
(
    pimple.dict().lookupOrDefault("correctPhi", true)
);

bool checkMeshCourantNo
(
    pimple.dict().lookupOrDefault("checkMeshCourantNo", false)
);

bool moveMeshOuterCorrectors
(
    pimple.dict().lookupOrDefault("moveMeshOuterCorrectors", false)
);

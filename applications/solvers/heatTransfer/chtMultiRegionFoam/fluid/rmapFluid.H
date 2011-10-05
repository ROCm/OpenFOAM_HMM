// Note:Map rho and rho.oldTime() since fluid rho assigned to at
// end of iteration.
rmap
(
    allRho(),
    rho,
    faceProcAddressing[procI],
    cellProcAddressing[procI],
    boundaryProcAddressing[procI]
);
rmap
(
    allRho().oldTime(),
    rho.oldTime(),
    faceProcAddressing[procI],
    cellProcAddressing[procI],
    boundaryProcAddressing[procI]
);

// Necessary? Probably only for boundary values since bcs on
// h are not the same as those on allh

rmap
(
    allh(),
    h,
    faceProcAddressing[procI],
    cellProcAddressing[procI],
    boundaryProcAddressing[procI]
);

procAlpha.set(procI, fvc::interpolate(turb.alphaEff()));

rmap
(
    allVolAlpha(),
    turb.alphaEff()(),
    faceProcAddressing[procI],
    cellProcAddressing[procI],
    boundaryProcAddressing[procI]
);

rmap
(
    allSource(),
    (DpDt + rad.Sh(thermo))(),
    faceProcAddressing[procI],
    cellProcAddressing[procI],
    boundaryProcAddressing[procI]
);

procPhi.set
(
    procI,
    new surfaceScalarField(phi)
);

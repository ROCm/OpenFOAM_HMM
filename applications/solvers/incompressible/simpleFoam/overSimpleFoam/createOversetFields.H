//- Overset specific

// Add solver-specific interpolations
{
    dictionary oversetDict;
    oversetDict.add("U", true);
    oversetDict.add("p", true);
    oversetDict.add("HbyA", true);
    oversetDict.add("grad(p)", true);

    const_cast<dictionary&>
    (
        mesh.schemesDict()
    ).add
    (
        "oversetInterpolationRequired",
        oversetDict,
        true
    );
}

// Mask field for zeroing out contributions on hole cells
#include "createCellMask.H"

#include "createInterpolatedCells.H"

bool adjustFringe
(
    simple.dict().lookupOrDefault("oversetAdjustPhi", false)
);
bool massFluxInterpolation
(
    simple.dict().lookupOrDefault("massFluxInterpolation", false)
);

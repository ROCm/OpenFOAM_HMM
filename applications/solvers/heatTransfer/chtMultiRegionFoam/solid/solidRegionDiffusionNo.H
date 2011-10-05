    scalar DiNum = -GREAT;

    forAll(solidToProc, i)
    {
#       include "setRegionSolidFields.H"

        DiNum = max
        (
            solidRegionDiffNo
            (
                procMeshes[solidToProc[i]],
                runTime,
                rho*cp,
                K
            ),
            DiNum
        );
    }

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    Test-deltaCoeffs

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "Distribution.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void assessDeltaCoeffs
(
    const fvMesh& mesh,
    const surfaceScalarField& d,
    scalar dCRes,
    fileName prefix
)
{
    Distribution<scalar> dCDist(dCRes);

    label minInternalFaceI = findMin(d.internalField());
    label maxInternalFaceI = findMax(d.internalField());

    if (minInternalFaceI >= 0 && maxInternalFaceI >= 0)
    {
        Pout<< "Internal field" << endl;

        Pout<< "    min deltaCoeff "
            << d.internalField()[minInternalFaceI]
            << " on face " << minInternalFaceI  << nl
            << "    max deltaCoeff "
            << d.internalField()[maxInternalFaceI]
            << " on face " << maxInternalFaceI << nl
            << endl;

        forAll(d.internalField(), fI)
        {
            dCDist.add(d.internalField()[fI]);
        }
    }

    forAll(mesh.boundaryMesh(), patchI)
    {
        label minPatchFaceI = findMin(d.boundaryField()[patchI]);
        label maxPatchFaceI = findMax(d.boundaryField()[patchI]);

        Pout<< "Patch " << mesh.boundaryMesh()[patchI].name() << nl
            << "    type " << mesh.boundaryMesh()[patchI].type() << endl;

        if (minPatchFaceI >= 0 && maxPatchFaceI >= 0)
        {
            Pout<< "    min deltaCoeff "
                << d.boundaryField()[patchI][minPatchFaceI]
                << " on face "
                << minPatchFaceI + mesh.boundaryMesh()[patchI].start() << nl
                << "    max deltaCoeff "
                << d.boundaryField()[patchI][maxPatchFaceI]
                << " on face "
                << maxPatchFaceI + mesh.boundaryMesh()[patchI].start() << nl
                << endl;

            forAll(d.boundaryField()[patchI], fI)
            {
                dCDist.add(d.boundaryField()[patchI][fI]);
            }
        }
        else
        {
            Pout<< endl;
        }
    }

    reduce(dCDist, sumOp< Distribution<scalar> >());

    dCDist.write
    (
        prefix + "_time_" + mesh.time().timeName()
    );

}

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    argList::addOption
    (
        "deltaCoeffResolution",
        "scalar",
        "deltaCoeff distribution resolution - default is 1"
    );
    argList::addBoolOption
    (
        "readCellCentres",
        "read and override the cell centres"
    );
    argList::addOption
    (
        "deltaCoeffRatio",
        "scalar",
        "max deltaCoeff ratio, will write new cellCentres if violated."
    );
    argList::addOption
    (
        "ccMotionCoeff",
        "scalar",
        "fraction of distance to move read cell centre towards centroid"
        "is deltaCoeffRatio is violated, default to 0.1."
    );

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    const scalar dCRes = args.optionLookupOrDefault
    (
        "deltaCoeffResolution",
        1.0
    );

    const bool readCellCentres = args.optionFound("readCellCentres");

    instantList timeDirs = timeSelector::select0(runTime, args);

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        Info<< "Time = " << runTime.timeName() << nl << endl;

        mesh.readUpdate();
        mesh.clearOut();
        mesh.globalData();

        surfaceScalarField dC_centroids = mesh.deltaCoeffs();

        Info<< "Assessing centroid delta coeffs" << endl;

        assessDeltaCoeffs
        (
            mesh,
            dC_centroids,
            dCRes,
            "centroid_deltaCoeffDistribuition"
        );

        if (readCellCentres)
        {
            pointIOField overrideCCs
            (
                IOobject
                (
                    "cellCentres",
                    mesh.pointsInstance(),
                    polyMesh::meshSubDir,
                    runTime,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            );

            Pout<< "Read " << overrideCCs.size() << " cell centres" << endl;

            pointField centroids = mesh.cellCentres();

            mesh.clearOut();

            mesh.overrideCellCentres(overrideCCs);

            surfaceScalarField dC_readCentre = mesh.deltaCoeffs();

            Info<< "Assessing read cell centre delta coeffs" << endl;

            assessDeltaCoeffs
            (
                mesh,
                dC_readCentre,
                dCRes,
                "readCentre_deltaCoeffDistribuition"
            );

            scalarField deltaCoeffRatio =
                dC_readCentre.internalField()
               /dC_centroids.internalField();

            Pout<< "deltaCoeffRatio min " << min(deltaCoeffRatio) << nl
                << "deltaCoeffRatio max " << max(deltaCoeffRatio)
                << endl;

            scalar dCRatio = -GREAT;

            if (args.optionReadIfPresent("deltaCoeffRatio", dCRatio))
            {
                scalar ccMotionCoeff = args.optionLookupOrDefault
                (
                    "ccMotionCoeff",
                    0.1
                );

                if (dCRatio <= 1.0)
                {
                    FatalErrorIn(args.executable())
                        << "deltaCoeffRatio must be > 1.0"
                        << exit(FatalError);
                }

                scalar rdCRatio = 1/dCRatio;

                pointIOField newCCs
                (
                    IOobject
                    (
                        "cellCentres_deltaCoeffMoved",
                        mesh.pointsInstance(),
                        polyMesh::meshSubDir,
                        runTime,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    overrideCCs
                );

                forAll(mesh.cells(), cI)
                {
                    cell cFaces = mesh.cells()[cI];

                    scalar dCRatioMin = 1;
                    scalar dCRatioMax = 1;

                    forAll(cFaces, cFI)
                    {
                        label fI = cFaces[cFI];

                        if (mesh.isInternalFace(fI))
                        {
                            dCRatioMin = min
                            (
                                dCRatioMin,
                                deltaCoeffRatio[fI]
                            );

                            dCRatioMax = max
                            (
                                dCRatioMax,
                                deltaCoeffRatio[fI]
                            );
                        }
                    }

                    if (dCRatioMax > dCRatio || dCRatioMin < rdCRatio)
                    {
                        Pout<< "cell " << cI << " has min/max dCRatios: "
                            << dCRatioMin << " " << dCRatioMax << endl;

                        newCCs[cI] +=
                            ccMotionCoeff*(centroids[cI] - overrideCCs[cI]);
                    }
                }

                Pout<< "Writing " << newCCs.objectPath() << endl;

                newCCs.write();
            }
        }
    }

    Info<< nl << "End" << nl << endl;

    return 0;
}


// ************************************************************************* //

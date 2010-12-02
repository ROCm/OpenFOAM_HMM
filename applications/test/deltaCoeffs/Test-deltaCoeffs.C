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

int main(int argc, char *argv[])
{
    argList::addOption
    (
        "deltaCoeffResolution",
        "scalar",
        "deltaCoeff distribution resolution - default is 1"
    );

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    const scalar dCRes = args.optionLookupOrDefault
    (
        "deltaCoeffResolution",
        1
    );

    Distribution<scalar> dCDist(dCRes);

    const surfaceScalarField& d = mesh.deltaCoeffs();

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
        "deltaCoeffDistribuition_time_" + mesh.time().timeName()
    );

    Info<< nl << "End" << nl << endl;

    return 0;
}


// ************************************************************************* //

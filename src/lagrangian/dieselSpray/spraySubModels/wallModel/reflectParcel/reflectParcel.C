/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "reflectParcel.H"
#include "addToRunTimeSelectionTable.H"
#include "wallPolyPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(reflectParcel, 0);

    addToRunTimeSelectionTable
    (
        wallModel,
        reflectParcel,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::reflectParcel::reflectParcel
(
    const dictionary& dict,
    const volVectorField& U,
    spray& sm
)
:
    wallModel(dict, U, sm),
    U_(U),
    coeffsDict_(dict.subDict(typeName + "Coeffs")),
    elasticity_(readScalar(coeffsDict_.lookup("elasticity")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::reflectParcel::~reflectParcel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Return 'keepParcel'
bool Foam::reflectParcel::wallTreatment
(
    parcel& p,
    const label globalFaceI
) const
{
    label patchI = p.patch(globalFaceI);
    label faceI = p.patchFace(patchI, globalFaceI);

    const polyMesh& mesh = spray_.mesh();

    if (isA<wallPolyPatch>(mesh_.boundaryMesh()[patchI]))
    {
        // wallNormal defined to point outwards of domain
        vector Sf = mesh_.Sf().boundaryField()[patchI][faceI];
        Sf /= mag(Sf);

        if (!mesh.moving())
        {
            // static mesh
            scalar Un = p.U() & Sf;

            if (Un > 0)
            {
                p.U() -= (1.0 + elasticity_)*Un*Sf;
            }
        }
        else
        {
            // moving mesh
            vector Ub1 = U_.boundaryField()[patchI][faceI];
            vector Ub0 = U_.oldTime().boundaryField()[patchI][faceI];

            scalar dt = spray_.runTime().deltaTValue();
            const vectorField& oldPoints = mesh.oldPoints();

            const vector& Cf1 = mesh.faceCentres()[globalFaceI];

            vector Cf0 = mesh.faces()[globalFaceI].centre(oldPoints);
            vector Cf = Cf0 + p.stepFraction()*(Cf1 - Cf0);
            vector Sf0 = mesh.faces()[globalFaceI].normal(oldPoints);

            // for layer addition Sf0 = vector::zero and we use Sf
            if (mag(Sf0) > SMALL)
            {
                Sf0 /= mag(Sf0);
            }
            else
            {
                Sf0 = Sf;
            }

            scalar magSfDiff = mag(Sf - Sf0);

            vector Ub = Ub0 + p.stepFraction()*(Ub1 - Ub0);

            if (magSfDiff > SMALL)
            {
                // rotation + translation
                vector Sfp = Sf0 + p.stepFraction()*(Sf - Sf0);

                vector omega = Sf0 ^ Sf;
                scalar magOmega = mag(omega);
                omega /= magOmega+SMALL;

                scalar phiVel = ::asin(magOmega)/dt;

                scalar dist = (p.position() - Cf) & Sfp;
                vector pos = p.position() - dist*Sfp;
                vector vrot = phiVel*(omega ^ (pos - Cf));

                vector v = Ub + vrot;

                scalar Un = ((p.U() - v) & Sfp);

                if (Un > 0.0)
                {
                    p.U() -= (1.0 + elasticity_)*Un*Sfp;
                }
            }
            else
            {
                // translation
                vector Ur = p.U() - Ub;
                scalar Urn = Ur & Sf;
            /*
                if (mag(Ub-v) > SMALL)
                {
                    Info<< "reflectParcel:: v = " << v
                        << ", Ub = " << Ub
                        << ", faceI = " << faceI
                        << ", patchI = " << patchI
                        << ", globalFaceI = " << globalFaceI
                        << ", name = " << mesh_.boundaryMesh()[patchI].name()
                        << endl;
                }
            */
                if (Urn > 0.0)
                {
                    p.U() -= (1.0 + elasticity_)*Urn*Sf;
                }
            }
        }

    }
    else
    {
        FatalErrorIn("bool reflectParcel::wallTreatment(parcel& parcel) const")
            << " parcel has hit a boundary " << mesh_.boundary()[patchI].type()
            << " which not yet has been implemented." << nl
            << abort(FatalError);
    }
    return true;
}


// ************************************************************************* //

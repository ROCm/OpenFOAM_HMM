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

#include "rotorDiskSource.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"
#include "unitConversion.H"

using namespace Foam::constant;

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::rotorDiskSource::writeField
(
    const word& name,
    const List<Type>& values,
    const bool writeNow
) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    if (mesh_.time().outputTime() || writeNow)
    {
        tmp<fieldType> tfld
        (
            new fieldType
            (
                IOobject
                (
                    name,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensioned<Type>("zero", dimless, pTraits<Type>::zero)
            )
        );

        Field<Type>& fld = tfld().internalField();

        if (cells_.size() != values.size())
        {
            FatalErrorIn("") << "cells_.size() != values_.size()"
                << abort(FatalError);
        }

        forAll(cells_, i)
        {
            const label cellI = cells_[i];
            fld[cellI] = values[i];
        }

        tfld().write();
    }
}


template<class RhoType>
Foam::tmp<Foam::volVectorField> Foam::rotorDiskSource::calculateForces
(
    const RhoType& rho,
    const vectorField& U,
    const dimensionSet& dims
)
{
    tmp<volVectorField> tForce
    (
        new volVectorField
        (
            IOobject
            (
                "rotorForce",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector("zero", dims, vector::zero)
        )
    );

    vectorField& force = tForce().internalField();
    const scalarField& V = mesh_.V();


    // logging info
    scalar dragEff = 0.0;
    scalar liftEff = 0.0;
    scalar AOAmin = GREAT;
    scalar AOAmax = -GREAT;

    forAll(cells_, i)
    {
        if (area_[i] > ROOTVSMALL)
        {
            const label cellI = cells_[i];

            const scalar radius = x_[i].x();

            // apply correction due to flap in cartesian frame
            vector Uc = R_[i] & U[cellI];

            // velocity in local reference frame
            Uc = coordSys_.localVector(Uc);

            // set radial component of velocity to zero
            Uc.x() = 0.0;

            // remove blade linear velocity from blade normal component
            Uc.y() -= radius*omega_;

            // velocity magnitude
            scalar magUc = mag(Uc);

            // determine blade data for this radius
            // i1 = index of upper bound data point in blade list
            scalar twist = 0.0;
            scalar chord = 0.0;
            label i1 = -1;
            label i2 = -1;
            scalar invDr = 0.0;
            blade_.interpolate(radius, twist, chord, i1, i2, invDr);

            // effective angle of attack
            scalar alphaEff = alphag_[i] + twist - atan(Uc.z()/Uc.y());
            AOAmin = min(AOAmin, alphaEff);
            AOAmax = max(AOAmax, alphaEff);

            // determine profile data for this radius and angle of attack
            const label profile1 = blade_.profileID()[i1];
            const label profile2 = blade_.profileID()[i2];

            scalar Cd1 = 0.0;
            scalar Cl1 = 0.0;
            profiles_[profile1].Cdl(alphaEff, Cd1, Cl1);

            scalar Cd2 = 0.0;
            scalar Cl2 = 0.0;
            profiles_[profile2].Cdl(alphaEff, Cd2, Cl2);

            scalar Cd = invDr*(Cd2 - Cd1) + Cd1;
            scalar Cl = invDr*(Cl2 - Cl1) + Cl1;

            // apply tip effect for blade lift
            scalar tipFactor = 1.0;
            if (radius/rMax_ > tipEffect_)
            {
                tipFactor = 0.0;
            }

            // calculate forces
            const scalar pDyn = 0.5*rho[cellI]*sqr(magUc);
            const scalar f = pDyn*chord*nBlades_*area_[i]/(mathematical::twoPi);
            const vector localForce(0.0, f*Cd, tipFactor*f*Cl);

            // accumulate forces
            dragEff += localForce.y();
            liftEff += localForce.z();

            // convert force to global cartesian co-ordinate system
            force[cellI] = coordSys_.globalVector(localForce);

            force[cellI] /= V[cellI];
        }
    }


    Info<< type() << " output:" << nl
        << "    min/max(AOA)   = " << radToDeg(AOAmin) << ", "
        << radToDeg(AOAmax) << nl
        << "    Effective drag = " << dragEff << nl
        << "    Effective lift = " << liftEff << endl;


    return tForce;
}


// ************************************************************************* //

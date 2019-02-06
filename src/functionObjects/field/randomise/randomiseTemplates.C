/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2016 OpenFOAM Foundation
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

#include "volFields.H"
#include "Random.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::functionObjects::randomise::calcRandomised()
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;

    if (foundObject<VolFieldType>(fieldName_, false))
    {
        const VolFieldType& field = lookupObject<VolFieldType>(fieldName_);

        resultName_ = fieldName_ & "Random";

        auto trfield = tmp<VolFieldType>::New(field);
        auto& rfield = trfield.ref();

        Random rand(1234567);

        for (Type& cellval : rfield)
        {
            Type rndPert;
            rand.randomise01(rndPert);
            rndPert = 2.0*rndPert - pTraits<Type>::one;
            rndPert /= mag(rndPert);

            cellval += magPerturbation_*rndPert;
        }

        return store(resultName_, trfield);
    }

    return false;
}


// ************************************************************************* //

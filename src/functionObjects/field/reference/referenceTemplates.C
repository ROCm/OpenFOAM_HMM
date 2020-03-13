/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2020 OpenCFD Ltd.
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

#include "interpolation.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::functionObjects::reference::calcType()
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;

    const VolFieldType* vfPtr = findObject<VolFieldType>(fieldName_);

    if (vfPtr)
    {
        const VolFieldType& vf = *vfPtr;

        dimensioned<Type> offset("offset", vf.dimensions(), Zero, localDict_);

        dimensioned<Type> cellValue("value", vf.dimensions(), Zero);

        if (positionIsSet_)
        {
            cellValue.value() = -pTraits<Type>::one*GREAT;

            // Might trigger parallel comms (e.g. volPointInterpolation, if
            // result is not yet cached) so have all processors do it
            autoPtr<interpolation<Type>> interpolator
            (
                interpolation<Type>::New(interpolationScheme_, vf)
            );

            if (celli_ != -1)
            {
                cellValue.value() =
                    interpolator().interpolate(position_, celli_, -1);
            }

            reduce(cellValue.value(), maxOp<Type>());

            Log << "    sampled value: " << cellValue.value() << endl;
        }

        return store
        (
            resultName_,
            scale_*(vf - cellValue + offset)
        );
    }

    return false;
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "minMaxFields.H"
#include "volFields.H"
#include "dictionary.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::minMaxFields::calcMinMaxFields(const word& fieldName)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    if (obr_.foundObject<fieldType>(fieldName))
    {
        const fieldType& field = obr_.lookupObject<fieldType>(fieldName);
        scalar minValue = min(mag(field)).value();
        scalar maxValue = max(mag(field)).value();

        reduce(minValue, minOp<scalar>());
        reduce(maxValue, maxOp<scalar>());

        if (Pstream::master())
        {
            minMaxFieldsFilePtr_() << obr_.time().value() << tab
                << fieldName << tab << minValue << tab << maxValue << endl;

            if (log_)
            {
                Info<< "minMaxFields output:" << nl
                    << "    min(mag(" << fieldName << ")) = " << minValue << nl
                    << "    max(mag(" << fieldName << ")) = " << maxValue << nl
                    << endl;
            }
        }
    }
}


// ************************************************************************* //

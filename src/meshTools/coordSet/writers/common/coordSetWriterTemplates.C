/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::UPtrList<const Foam::Field<Type>>
Foam::coordSetWriter::repackageFields(const Field<Type>& field)
{
    UPtrList<const Field<Type>> fieldPtrs(1);
    fieldPtrs.set(0, &field);

    return fieldPtrs;
}


template<class Type>
Foam::UPtrList<const Foam::Field<Type>>
Foam::coordSetWriter::repackageFields(const UList<Field<Type>>& fieldValues)
{
    UPtrList<const Field<Type>> fieldPtrs(fieldValues.size());
    forAll(fieldValues, i)
    {
        fieldPtrs.set(i, &(fieldValues[i]));
    }

    return fieldPtrs;
}


template<class Type>
void Foam::coordSetWriter::writeTable
(
    Ostream& os,
    const coordSet& coords,
    const UList<Type>& values,
    const char* sep
)
{
    forAll(coords, pointi)
    {
        // Output coordinate (point or scalar) with separator
        if (coords.hasVectorAxis())
        {
            const vector& p = coords.vectorCoord(pointi);
            os << p.x() << sep << p.y() << sep << p.z();
        }
        else
        {
            os << coords.scalarCoord(pointi);
        }

        // Output component values with separator
        const auto& val = values[pointi];
        for (direction d=0; d < pTraits<Type>::nComponents; ++d)
        {
            os << sep << component(val, d);
        }
        os << nl;
    }
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenCFD Ltd.
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

#include "IOField.H"
#include "OFstream.H"
#include "pointField.H"
#include "vectorField.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::functionObjects::dataCloud::writeField
(
    Ostream& os,
    const vectorField& points,
    const Field<Type>& field
)
{
    const label len = field.size();

    for (label pointi=0; pointi<len; ++pointi)
    {
        const point& pt = points[pointi];
        const Type& val = field[pointi];

        os << pt.x() << ' ' << pt.y() << ' ' << pt.z();

        for (direction cmpt=0; cmpt < pTraits<Type>::nComponents; ++cmpt)
        {
            os << ' ' << component(val, cmpt);
        }
        os << nl;
    }
}


template<class Type>
bool Foam::functionObjects::dataCloud::writeField
(
    const fileName& outputName,
    const objectRegistry& obrTmp
) const
{
    // Fields are not always on all processors (eg, multi-component parcels).
    // Thus need to resolve between all processors.

    const auto* fldPtr = obrTmp.findObject<IOField<Type>>(fieldName_);

    if (!returnReduce((fldPtr != nullptr), orOp<bool>()))
    {
        return false;
    }

    const auto* pointsPtr = obrTmp.findObject<vectorField>("position");

    if (!pointsPtr)
    {
        // This should be impossible
        return false;
    }

    if (Pstream::master())
    {
        OFstream os(outputName);

        os << "# x y z " << fieldName_ << nl;

        // Master
        if (fldPtr)
        {
            writeField(os, *pointsPtr, *fldPtr);
        }

        // Slaves - recv
        for (int slave=1; slave<Pstream::nProcs(); ++slave)
        {
            IPstream fromSlave(Pstream::commsTypes::blocking, slave);
            vectorField points(fromSlave);
            Field<Type> fld(fromSlave);

            writeField(os, points, fld);
        }
    }
    else
    {
        // Slaves

        OPstream toMaster(Pstream::commsTypes::blocking, Pstream::masterNo());

        if (fldPtr)
        {
            toMaster
                << *pointsPtr
                << *fldPtr;
        }
        else
        {
            toMaster
                << vectorField()
                << Field<Type>();
        }
    }

    return true;
}


// ************************************************************************* //

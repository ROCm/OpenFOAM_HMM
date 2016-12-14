/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenCFD Ltd.
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

template<class Type>
bool Foam::functionObjects::particleDistribution::processField
(
    const objectRegistry& obr,
    const word& fieldName,
    const List<DynamicList<label>>& addr
)
{
    if (obr.foundObject<IOField<Type>>(fieldName))
    {
        const IOField<Type>& field =
            obr.lookupObject<IOField<Type>>(fieldName);

        if (addr.size())
        {
            forAll(addr, i)
            {
                const Field<Type> subField(field, addr[i]);
                for (direction d = 0; d < pTraits<Type>::nComponents; ++d)
                {
                    generateDistribution(fieldName, subField.component(d), i);
                }
            }
        }
        else
        {
            for (direction d = 0; d < pTraits<Type>::nComponents; ++d)
            {
                generateDistribution(fieldName, field.component(d));
            }
        }

        return true;
    }

    return false;
}


// ************************************************************************* //

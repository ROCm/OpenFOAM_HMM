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

#include "fvMesh.H"
#include "Time.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::basicSourceList::apply(fvMatrix<Type>& eqn)
{
    if (mesh_.time().timeIndex() == mesh_.time().startTimeIndex() + 1)
    {
        checkApplied();
    }

    const word& fieldName = eqn.psi().name();

    forAll(*this, i)
    {
        basicSource& source = this->operator[](i);

        label fieldI = source.applyToField(fieldName);

        if (source.isActive() && (fieldI != -1))
        {
            if (debug)
            {
                Info<< "Applying source " << source.name() << " to field "
                    << fieldName << endl;
            }

            source.addSup(eqn, fieldI);
            source.setValue(eqn, fieldI);
        }
    }
}


// ************************************************************************* //

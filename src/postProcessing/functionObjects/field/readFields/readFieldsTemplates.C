/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2010-2010 OpenCFD Ltd.
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

#include "readFields.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::readFields::setField(const word& fieldName)
{
    typedef GeometricField<Type, fvPatchField, volMesh> vfType;
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> sfType;

    if (obr_.foundObject<vfType>(fieldName))
    {
        if (debug)
        {
            Info<< "Field " << fieldName << " already in database" << endl;
        }

        vfType& vf =  const_cast<vfType&>(obr_.lookupObject<vfType>(fieldName));
        vf.checkOut();
    }
    if (obr_.foundObject<sfType>(fieldName))
    {
        if (debug)
        {
            Info<< "Field " << fieldName << " already in database" << endl;
        }

        sfType& sf =  const_cast<sfType&>(obr_.lookupObject<sfType>(fieldName));
        sf.checkOut();
    }

    const fvMesh& mesh = refCast<const fvMesh>(obr_);

    IOobject fieldHeader
    (
        fieldName,
        mesh.time().timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if
    (
        fieldHeader.headerOk()
     && fieldHeader.headerClassName() == vfType::typeName
    )
    {
        // store field on the mesh database
        Info<< "    Reading " << fieldName << endl;
        obr_.store(new vfType(fieldHeader, mesh));
    }
    else if
    (
        fieldHeader.headerOk()
     && fieldHeader.headerClassName() == sfType::typeName
    )
    {
        // store field on the mesh database
        Info<< "    Reading " << fieldName << endl;
        obr_.store(new sfType(fieldHeader, mesh));
    }
}


// ************************************************************************* //

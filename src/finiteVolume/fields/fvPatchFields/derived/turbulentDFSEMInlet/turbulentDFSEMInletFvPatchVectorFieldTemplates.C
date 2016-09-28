/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016 OpenCFD Ltd
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

#include "AverageIOField.H"
#include "pointToPointPlanarInterpolation.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::turbulentDFSEMInletFvPatchVectorField::interpolateOrRead
(
    const word& fieldName,
    const dictionary& dict,
    bool& interpolateField
) const
{
    if (dict.found(fieldName))
    {
        tmp<Field<Type> > tFld
        (
            new Field<Type>
            (
                fieldName,
                dict,
                this->patch().size()
            )
        );

        interpolateField = false;
        return tFld;
    }
    else
    {
        interpolateField = true;
        return interpolateBoundaryData<Type>(fieldName);
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::turbulentDFSEMInletFvPatchVectorField::interpolateBoundaryData
(
    const word& fieldName
) const
{
    const word& patchName = this->patch().name();

    // Note: reading from the '0' directory only
    IOobject io
    (
        fieldName,
        this->db().time().caseConstant(),
        "boundaryData"/patchName/"0",
        this->db(),
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE,
        false
    );

    Info<< "Turbulent DFSEM patch " << this->patch().name()
        << ": interpolating field " << fieldName
        << " from " << io.path() << endl;

    AverageIOField<Type> aFld(io);

    return patchMapper().interpolate(aFld);
}


// ************************************************************************* //

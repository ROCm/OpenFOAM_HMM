/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

#include "pointToPointPlanarInterpolation.H"
#include "Time.H"
#include "rawIOField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::turbulentDigitalFilterInletFvPatchVectorField::interpolateOrRead
(
    const word& fieldName,
    const dictionary& dict,
    bool& interpolateField
) const
{
    if (dict.found(fieldName))
    {
        tmp<Field<Type>> tFld
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
Foam::tmp<Foam::Field<Type>>
Foam::turbulentDigitalFilterInletFvPatchVectorField::interpolateBoundaryData
(
    const word& fieldName
) const
{
    const word& patchName = this->patch().name();

    // Reread values and interpolate
    const fileName valsFile
    (
        this->db().time().globalPath()
       /this->db().time().constant()
       /"boundaryData"
       /patchName
       /"0"
       /fieldName
    );

    IOobject io
    (
        valsFile,           // absolute path
        this->db().time(),
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false,              // no need to register
        true                // is global object (currently not used)
    );

    const rawIOField<Type> vals(io, false);


    Info<< "turbulentDigitalFilterInlet patch " << patchName
        << ": Interpolating field " << fieldName
        << " from " << valsFile << endl;

    return patchMapper().interpolate(vals);
}


template<class Form, class Type>
Form Foam::turbulentDigitalFilterInletFvPatchVectorField::randomSet
(
    const label len
)
{
    Form randomSet(len);

    std::generate
    (
        randomSet.begin(),
        randomSet.end(),
        [&]{return rndGen_.GaussNormal<Type>();}
    );

    return randomSet;
}


// ************************************************************************* //

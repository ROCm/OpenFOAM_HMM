/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd.
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

#include "fieldCoordinateSystemTransform.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(fieldCoordinateSystemTransform, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        fieldCoordinateSystemTransform,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fieldCoordinateSystemTransform::
fieldCoordinateSystemTransform
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    fieldSet_(mesh_),
    coordSys_(mesh_, dict.subDict("coordinateSystem"))
{
    read(dict);

    Info<< type() << " " << name << ":" << nl
        << "   Applying transformation from global Cartesian to local "
        << coordSys_ << nl << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::fieldCoordinateSystemTransform::
~fieldCoordinateSystemTransform()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::word
Foam::functionObjects::fieldCoordinateSystemTransform::transformFieldName
(
    const word& fieldName
) const
{
    return fieldName + ":Transformed";
}


bool Foam::functionObjects::fieldCoordinateSystemTransform::read
(
    const dictionary& dict
)
{
    if (fvMeshFunctionObject::read(dict))
    {
        fieldSet_.read(dict);
        return true;
    }

    return false;
}


bool Foam::functionObjects::fieldCoordinateSystemTransform::execute()
{
    fieldSet_.updateSelection();

    for (const word& fieldName : fieldSet_.selection())
    {
        transform<scalar>(fieldName);
        transform<vector>(fieldName);
        transform<sphericalTensor>(fieldName);
        transform<symmTensor>(fieldName);
        transform<tensor>(fieldName);
    }

    return true;
}


bool Foam::functionObjects::fieldCoordinateSystemTransform::write()
{
    forAllConstIters(fieldSet_, iter)
    {
        writeObject(transformFieldName(iter()));
    }

    return true;
}


// ************************************************************************* //

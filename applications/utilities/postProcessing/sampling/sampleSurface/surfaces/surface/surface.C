/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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

Description

\*---------------------------------------------------------------------------*/

#include "surface.H"
#include "polyMesh.H"
#include "primitiveMesh.H"
#include "meshSearch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(surface, 0);
defineRunTimeSelectionTable(surface, word);

autoPtr<surface> surface::New
(
    const word& sampleType,
    const polyMesh& mesh,
    meshSearch& searchEngine,
    const dictionary& dict
)
{
    wordConstructorTable::iterator cstrIter =
        wordConstructorTablePtr_
            ->find(sampleType);

    if (cstrIter == wordConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "surface::New(const word&, "
            "const polyMesh&, meshSearch&, const dictionary&)"
        )   << "Unknown sample type " << sampleType
            << endl << endl
            << "Valid sample types : " << endl
            << wordConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<surface>
    (
        cstrIter()
        (
            mesh,
            searchEngine,
            dict
        )
    );
}

} // End namespace Foam


bool Foam::surface::getBool
(
    const dictionary& dict,
    const word& key,
    const bool defaultVal
)
{
    if (dict.found(key))
    {
        return readBool(dict.lookup(key));
    }
    else
    {
        return defaultVal;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from mesh, name
Foam::surface::surface
(
    const polyMesh& mesh,
    meshSearch& searchEngine,
    const word& name
)
:
    mesh_(mesh),
    searchEngine_(searchEngine),
    name_(name)
{}


// Construct from dictionary
Foam::surface::surface
(
    const polyMesh& mesh,
    meshSearch& searchEngine,
    const dictionary& dict          
)
:
    mesh_(mesh),
    searchEngine_(searchEngine),
    name_(dict.lookup("name"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::surface::~surface()
{}


// ************************************************************************* //

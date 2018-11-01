/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2018 OpenCFD Ltd.
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

#include "surfaceWriter.H"

#include "MeshedSurfaceProxy.H"
#include "proxySurfaceWriter.H"

#include "HashTable.H"
#include "word.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(surfaceWriter, 0);
    defineRunTimeSelectionTable(surfaceWriter, word);
    defineRunTimeSelectionTable(surfaceWriter, wordDict);
    addNamedToRunTimeSelectionTable
    (
        surfaceWriter,
        surfaceWriter,
        word,
        null
    );
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::autoPtr<Foam::surfaceWriter>
Foam::surfaceWriter::New(const word& writeType)
{
    const auto cstrIter = wordConstructorTablePtr_->cfind(writeType);

    if (cstrIter.found())
    {
        return autoPtr<surfaceWriter>(cstrIter()());
    }
    else if (MeshedSurfaceProxy<face>::canWriteType(writeType))
    {
        // Unknown, but proxy handler can manage it
        return autoPtr<surfaceWriter>(new proxySurfaceWriter(writeType));
    }

    FatalErrorInFunction
        << "Unknown write type \"" << writeType << "\"\n\n"
        << "Valid types : "
        << wordConstructorTablePtr_->sortedToc() << nl
        << "Valid proxy types : "
        << MeshedSurfaceProxy<face>::writeTypes() << endl
        << exit(FatalError);

    return nullptr;
}


Foam::autoPtr<Foam::surfaceWriter>
Foam::surfaceWriter::New(const word& writeType, const dictionary& options)
{
    {
        // Constructor with options (dictionary)
        const auto cstrIter = wordDictConstructorTablePtr_->cfind(writeType);

        if (cstrIter.found())
        {
            return autoPtr<surfaceWriter>(cstrIter()(options));
        }
    }


    // Drop through to version without options

    const auto cstrIter = wordConstructorTablePtr_->cfind(writeType);

    if (cstrIter.found())
    {
        return autoPtr<surfaceWriter>(cstrIter()());
    }
    else if (MeshedSurfaceProxy<face>::canWriteType(writeType))
    {
        // Unknown, but proxy handler can manage it
        return autoPtr<surfaceWriter>
        (
            new proxySurfaceWriter(writeType, options)
        );
    }

    FatalErrorInFunction
        << "Unknown write type \"" << writeType << "\"\n\n"
        << "Valid types : "
        << wordConstructorTablePtr_->sortedToc() << nl
        << "Valid proxy types : "
        << MeshedSurfaceProxy<face>::writeTypes() << endl
        << exit(FatalError);

    return nullptr;
}


// ************************************************************************* //

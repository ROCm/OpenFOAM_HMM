/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
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

#include "kahipDecomp.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"

static const char* notImplementedMessage =
"You are trying to use kahip but do not have the kahipDecomp library loaded."
"\nThis message is from the dummy kahipDecomp stub library instead.\n"
"\n"
"Please install kahip and make sure that libkahip.so is in your "
"LD_LIBRARY_PATH.\n"
"The kahipDecomp library can then be built from "
"src/parallel/decompose/kahipDecomp and dynamically loading or linking"
" this library will add kahip as a decomposition method.\n";

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(kahipDecomp, 0);

    addToRunTimeSelectionTable
    (
        decompositionMethod,
        kahipDecomp,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::kahipDecomp::decomposeSerial
(
    const labelUList& adjncy,
    const labelUList& xadj,
    const UList<scalar>& cellWeights,
    List<label>& decomp
)
{
    FatalErrorInFunction
        << notImplementedMessage << exit(FatalError);

    return -1;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kahipDecomp::kahipDecomp
(
    const dictionary& decompDict
)
:
    metisLikeDecomp("kahip", decompDict)
{}


Foam::kahipDecomp::kahipDecomp
(
    const dictionary& decompDict,
    const word& regionName
)
:
    metisLikeDecomp("kahip", decompDict, regionName)
{}


// ************************************************************************* //

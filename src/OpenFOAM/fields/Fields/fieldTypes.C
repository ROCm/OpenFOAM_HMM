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

#include "fieldTypes.H"
#include "labelIOField.H"
#include "scalarIOField.H"
#include "vectorIOField.H"
#include "sphericalTensorIOField.H"
#include "symmTensorIOField.H"
#include "tensorIOField.H"

// * * * * * * * * * * * * * * * * Global Data * * * * * * * * * * * * * * * //

const Foam::wordList Foam::fieldTypes::basic
({
    Foam::labelIOField::typeName,
    Foam::scalarIOField::typeName,
    Foam::vectorIOField::typeName,
    Foam::sphericalTensorIOField::typeName,
    Foam::symmTensorIOField::typeName,
    Foam::tensorIOField::typeName
});


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

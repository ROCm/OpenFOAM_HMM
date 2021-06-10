/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "volumeExprDriver.H"
#include "volumeExprScanner.H"
#include "error.H"
#include "fvPatch.H"
#include "fvMesh.H"
#include "className.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace expressions
{
namespace volumeExpr
{

defineTypeNameAndDebug(parseDriver, 0);

addNamedToRunTimeSelectionTable
(
    fvExprDriver,
    parseDriver,
    dictionary,
    volume
);

addNamedToRunTimeSelectionTable
(
    fvExprDriver,
    parseDriver,
    idName,
    volume
);

addNamedToRunTimeSelectionTable
(
    fvExprDriver,
    parseDriver,
    dictionary,
    internalField
);

addNamedToRunTimeSelectionTable
(
    fvExprDriver,
    parseDriver,
    idName,
    internalField
);

} // End namespace volumeExpr
} // End namespace expressions
} // End namespace Foam


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

static inline const TimeState* lookupTimeState
(
    const polyMesh& m
)
{
    return &static_cast<const TimeState&>(m.time());
}

} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::expressions::volumeExpr::parseDriver::parseDriver
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    parsing::genericRagelLemonDriver(),
    expressions::fvExprDriver(dict, lookupTimeState(mesh)),
    mesh_(mesh),
    resultType_(),
    isLogical_(false),
    fieldGeoType_(NO_DATA),
    resultDimension_()
{}


Foam::expressions::volumeExpr::parseDriver::parseDriver
(
    const fvMesh& mesh,
    const parseDriver& driver
)
:
    parsing::genericRagelLemonDriver(),
    expressions::fvExprDriver(driver),
    mesh_(mesh),
    resultType_(),
    isLogical_(false),
    fieldGeoType_(NO_DATA),
    resultDimension_()
{
    resetTimeReference(mesh_.time());  // Extra safety
}


Foam::expressions::volumeExpr::parseDriver::parseDriver
(
    const word& meshName,
    const fvMesh& mesh
)
:
    parseDriver(mesh)
{
    //?? Info<< "Warn that meshName is ignored?" << nl;
}


Foam::expressions::volumeExpr::parseDriver::parseDriver
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    parsing::genericRagelLemonDriver(),
    expressions::fvExprDriver(dict, lookupTimeState(mesh)),
    mesh_(mesh),
    resultType_(),
    isLogical_(false),
    fieldGeoType_(NO_DATA),
    resultDimension_()
{}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::expressions::volumeExpr::parseDriver::readDict
(
    const dictionary& dict
)
{
    expressions::fvExprDriver::readDict(dict);
    dict.readIfPresent("dimensions", resultDimension_);

    return true;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

unsigned Foam::expressions::volumeExpr::parseDriver::parse
(
    const std::string& expr,
    size_t pos,
    size_t len
)
{
    scanner scan(this->debugScanner());

    scan.process(expr, pos, len, *this);

    return 0;
}


// ************************************************************************* //

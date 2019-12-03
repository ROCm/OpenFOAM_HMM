/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
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

#include "patchExprDriver.H"
#include "patchExprScanner.H"
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
namespace patchExpr
{

defineTypeNameAndDebug(parseDriver, 0);

addNamedToRunTimeSelectionTable
(
    fvExprDriver,
    parseDriver,
    dictionary,
    patch
);

addNamedToRunTimeSelectionTable
(
    fvExprDriver,
    parseDriver,
    idName,
    patch
);

} // End namespace patchExpr
} // End namespace expressions
} // End namespace Foam


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{
    static label getPatchID(const fvMesh& mesh, const word& patchName)
    {
        const auto& bMesh = mesh.boundaryMesh();

        const label patchId = bMesh.findPatchID(patchName);

        if (patchId < 0)
        {
            FatalErrorInFunction
                << "No patch " << patchName << " found in "
                << flatOutput(bMesh.names()) << nl
                << exit(FatalError);
        }
        return patchId;
    }


    static inline const fvPatch& findFvPatch
    (
        const fvMesh& mesh,
        const word& patchName
    )
    {
        return mesh.boundary()[getPatchID(mesh, patchName)];
    }

} // End namespace Foam


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

const Foam::fvPatch& Foam::expressions::patchExpr::parseDriver::getFvPatch
(
    const fvMesh& fvm,
    const dictionary& dict
)
{
    return findFvPatch
    (
        regionMesh(dict, fvm, true),
        dict.get<word>("patch")
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::expressions::patchExpr::parseDriver::parseDriver(const fvPatch& p)
:
    parsing::genericRagelLemonDriver(),
    expressions::fvExprDriver(),
    patch_(p)
{}


Foam::expressions::patchExpr::parseDriver::parseDriver
(
    const fvPatch& p,
    const dictionary& dict
)
:
    parsing::genericRagelLemonDriver(),
    expressions::fvExprDriver(dict),
    patch_(p)
{}


Foam::expressions::patchExpr::parseDriver::parseDriver
(
    const fvPatch& p,
    const parseDriver& driver_
)
:
    parsing::genericRagelLemonDriver(),
    expressions::fvExprDriver(driver_),
    patch_(p)
{}


Foam::expressions::patchExpr::parseDriver::parseDriver
(
    const word& patchName,
    const fvMesh& mesh
)
:
    parseDriver(findFvPatch(mesh, patchName))
{}


Foam::expressions::patchExpr::parseDriver::parseDriver
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    parseDriver(getFvPatch(mesh, dict), dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

unsigned Foam::expressions::patchExpr::parseDriver::parse
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

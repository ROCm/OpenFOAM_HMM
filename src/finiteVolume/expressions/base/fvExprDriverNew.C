/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2010-2018 Bernhard Gschaider
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

#include "fvExprDriver.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::expressions::fvExprDriver>
Foam::expressions::fvExprDriver::New
(
    const dictionary& dict
)
{
    return fvExprDriver::New
    (
        dict,
        defaultMesh()
    );
}


Foam::autoPtr<Foam::expressions::fvExprDriver>
Foam::expressions::fvExprDriver::New
(
    const dictionary& dict,
    const fvMesh& mesh
)
{
    const word driverType(dict.get<word>("valueType"));

    auto cstrIter = dictionaryConstructorTablePtr_->cfind(driverType);

    if (!cstrIter.found())
    {
        FatalIOErrorInLookup
        (
            dict,
            "valueType",
            driverType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    DebugInFunction << "Creating driver of type " << driverType << endl;

    resetDefaultMesh(mesh, false); // lazy

    return autoPtr<fvExprDriver>(cstrIter()(dict, mesh));
}


Foam::autoPtr<Foam::expressions::fvExprDriver>
Foam::expressions::fvExprDriver::New
(
    const word& driverType,
    const word& id,
    const fvMesh& mesh
)
{
    auto cstrIter = idNameConstructorTablePtr_->cfind(driverType);

    if (!cstrIter.found())
    {
        FatalErrorInLookup
        (
            "valueType",
            driverType,
            *idNameConstructorTablePtr_
        ) << exit(FatalError);
    }

    DebugInFunction << "Creating driver of type " << driverType << endl;

    resetDefaultMesh(mesh, false); // lazy

    return autoPtr<fvExprDriver>(cstrIter()(id, mesh));
}


// ************************************************************************* //

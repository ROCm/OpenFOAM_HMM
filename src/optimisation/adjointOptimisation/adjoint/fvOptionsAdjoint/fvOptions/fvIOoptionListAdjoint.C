/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2007-2019 PCOpt/NTUA
                            | Copyright (C) 2013-2019 FOSS GP
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

#include "fvIOoptionListAdjoint.H"
#include "fvMesh.H"
#include "Time.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::IOobject Foam::fv::IOoptionListAdjoint::createIOobject
(
    const fvMesh& mesh
) const
{
    IOobject io
    (
        "fvOptionsAdjoint",
        mesh.time().system(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (io.typeHeaderOk<IOdictionary>(true))
    {
        Info<< "Creating fintite volume adjoint options from " << io.name()
            << nl << endl;

        io.readOpt() = IOobject::MUST_READ_IF_MODIFIED;
        return io;
    }
    else
    {
        Info<< "No finite volume adjoint options present" << nl << endl;

        io.readOpt() = IOobject::NO_READ;
        return io;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::IOoptionListAdjoint::IOoptionListAdjoint
(
    const fvMesh& mesh
)
:
    IOdictionary(createIOobject(mesh)),
    optionList(mesh, *this)
{}


bool Foam::fv::IOoptionListAdjoint::read()
{
    if (regIOobject::read())
    {
        optionList::read(*this);
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //

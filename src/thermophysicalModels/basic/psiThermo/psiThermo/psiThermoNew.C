/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "psiThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::psiThermo> Foam::psiThermo::New
(
    const fvMesh& mesh
)
{
    IOdictionary thermoDict
    (
        IOobject
        (
            "thermophysicalProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            false
        )
    );

    word thermoTypeName;

    if (thermoDict.isDict("thermoType"))
    {
        const dictionary& thermoTypeDict(thermoDict.subDict("thermoType"));

        word type(thermoTypeDict.lookup("type"));
        word mixture(thermoTypeDict.lookup("mixture"));
        word transport(thermoTypeDict.lookup("transport"));
        word thermo(thermoTypeDict.lookup("thermo"));
        word equationOfState(thermoTypeDict.lookup("equationOfState"));
        word energy(thermoTypeDict.lookup("energy"));

        thermoTypeName =
            type + '<'
          + mixture + '<'
          + transport + '<'
          + thermo + '<'
          + equationOfState + ">,"
          + energy + ">>>";
    }
    else
    {
        thermoTypeName = word(thermoDict.lookup("thermoType"));
    }

    Info<< "Selecting thermodynamics package " << thermoTypeName << endl;

    fvMeshConstructorTable::iterator cstrIter =
        fvMeshConstructorTablePtr_->find(thermoTypeName);

    if (cstrIter == fvMeshConstructorTablePtr_->end())
    {
        FatalErrorIn("psiThermo::New(const fvMesh&)")
            << "Unknown psiThermo type " << thermoTypeName << nl << nl
            << "Valid psiThermo types are:" << nl
            << fvMeshConstructorTablePtr_->sortedToc() << nl
            << exit(FatalError);
    }

    return autoPtr<psiThermo>(cstrIter()(mesh));
}


// ************************************************************************* //

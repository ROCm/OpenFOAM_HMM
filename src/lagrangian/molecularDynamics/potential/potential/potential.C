/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "potential.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::potential::potential::readPotentialDict()
{
    Info<< nl <<  "Reading potential dictionary:" << endl;

    IOdictionary idListDict
    (
        IOobject
        (
            "idList",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    idList_ = List<word>(idListDict.lookup("idList"));

    IOdictionary moleculePropertiesDict
    (
        IOobject
        (
            "moleculeProperties",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    DynamicList<word> allSiteIdNames;

    forAll(idList_, i)
    {
        const word& id(idList_[i]);

        if (!moleculePropertiesDict.found(id))
        {
            FatalErrorIn("potential.C") << nl
                << id << " molecule subDict not found"
                << nl << abort(FatalError);
        }

        const dictionary& molDict(moleculePropertiesDict.subDict(id));

        List<word> siteIdNames = molDict.lookup("siteIds");

        forAll(siteIdNames, sI)
        {
            const word& siteId = siteIdNames[sI];

            if(findIndex(allSiteIdNames, siteId) == -1)
            {
                allSiteIdNames.append(siteId);
            }
        }
    }

    allSiteIdNames_.transfer(allSiteIdNames.shrink());

    Info<< nl << "Unique site ids found: " << allSiteIdNames_ << endl;

    List<word> tetherIdList(0);

    if (idListDict.found("tetherIdList"))
    {
        tetherIdList = List<word>(idListDict.lookup("tetherIdList"));
    }

    IOdictionary potentialDict
    (
        IOobject
        (
            "potentialDict",
            mesh_.time().system(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    potentialEnergyLimit_ = readScalar
    (
        potentialDict.lookup("potentialEnergyLimit")
    );

    if (potentialDict.found("removalOrder"))
    {
        List<word> remOrd = potentialDict.lookup("removalOrder");

        removalOrder_.setSize(remOrd.size());

        forAll(removalOrder_, rO)
        {
            removalOrder_[rO] = findIndex(idList_, remOrd[rO]);

            if (removalOrder_[rO] == -1)
            {
                FatalErrorIn("potentials.C") << nl
                    << "removalOrder entry: " << remOrd[rO]
                    << " not found in idList."
                    << nl << abort(FatalError);
            }
        }
    }

    // *************************************************************************
    // Pair potentials

    if (!potentialDict.found("pair"))
    {
        FatalErrorIn("potentials.C") << nl
            << "pair potential specification subDict not found"
            << abort(FatalError);
    }

    const dictionary& pairDict = potentialDict.subDict("pair");

    pairPotentials_.buildPotentials
    (
        allSiteIdNames_,
        pairDict,
        mesh_
    );

    // *************************************************************************
    // Tether potentials

    if (tetherIdList.size())
    {
        if (!potentialDict.found("tether"))
        {
            FatalErrorIn("potential.C") << nl
                << "tether potential specification subDict not found"
                << abort(FatalError);
        }

        const dictionary& tetherDict = potentialDict.subDict("tether");

        tetherPotentials_.buildPotentials
        (
            allSiteIdNames_,
            tetherDict,
            tetherIdList
        );
    }

    // *************************************************************************
    // External Forces

    gravity_ = vector::zero;

    if (potentialDict.found("external"))
    {

        Info << nl << "Reading external forces:" << endl;

        const dictionary& externalDict = potentialDict.subDict("external");

        // *********************************************************************
        // gravity

        if (externalDict.found("gravity"))
        {
            gravity_ = externalDict.lookup("gravity");
        }
    }

    Info << nl << tab << "gravity = " << gravity_ << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::potential::potential(const polyMesh& mesh)
:
    mesh_(mesh)
{
    readPotentialDict();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::potential::~potential()
{}


// ************************************************************************* //

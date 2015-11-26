/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015 OpenCFD Ltd.
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

#include "streamLine.H"
#include "fvMesh.H"
#include "streamLineParticleCloud.H"
#include "sampledSet.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(streamLine, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::streamLine::track()
{
    const fvMesh& mesh = dynamic_cast<const fvMesh&>(obr_);

    IDLList<streamLineParticle> initialParticles;
    streamLineParticleCloud particles
    (
        mesh,
        cloudName_,
        initialParticles
    );

    const sampledSet& seedPoints = sampledSetPtr_();

    forAll(seedPoints, i)
    {
        particles.addParticle
        (
            new streamLineParticle
            (
                mesh,
                seedPoints[i],
                seedPoints.cells()[i],
                lifeTime_               // lifetime
            )
        );
    }

    label nSeeds = returnReduce(particles.size(), sumOp<label>());

    if (log_) Info<< "    seeded " << nSeeds << " particles" << endl;

    // Read or lookup fields
    PtrList<volScalarField> vsFlds;
    PtrList<interpolation<scalar> > vsInterp;
    PtrList<volVectorField> vvFlds;
    PtrList<interpolation<vector> > vvInterp;
    label UIndex = -1;

    initInterpolations
    (
        nSeeds,
        UIndex,
        vsFlds,
        vsInterp,
        vvFlds,
        vvInterp
    );

    // Additional particle info
    streamLineParticle::trackingData td
    (
        particles,
        vsInterp,
        vvInterp,
        UIndex,         // index of U in vvInterp
        trackForward_,  // track in +u direction?
        nSubCycle_,     // automatic track control:step through cells in steps?
        trackLength_,   // fixed track length

        allTracks_,
        allScalars_,
        allVectors_
    );


    // Set very large dt. Note: cannot use GREAT since 1/GREAT is SMALL
    // which is a trigger value for the tracking...
    const scalar trackTime = Foam::sqrt(GREAT);

    // Track
    particles.move(td, trackTime);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::streamLine::streamLine
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    streamLineBase(name, obr, dict, loadFromFiles)
{
    // Check if the available mesh is an fvMesh otherise deactivate
    if (setActive<fvMesh>())
    {
        read(dict_);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::streamLine::~streamLine()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::streamLine::read(const dictionary& dict)
{
    if (active_)
    {
        streamLineBase::read(dict);

        bool subCycling = dict.found("nSubCycle");
        bool fixedLength = dict.found("trackLength");

        if (subCycling && fixedLength)
        {
            FatalIOErrorIn("streamLine::read(const dictionary&)", dict)
                << "Cannot both specify automatic time stepping (through '"
                << "nSubCycle' specification) and fixed track length (through '"
                << "trackLength')"
                << exit(FatalIOError);
        }

        nSubCycle_ = 1;
        if (dict.readIfPresent("nSubCycle", nSubCycle_))
        {
            trackLength_ = VGREAT;
            if (nSubCycle_ < 1)
            {
                nSubCycle_ = 1;
            }

            if (log_) Info
                << "    automatic track length specified through"
                << " number of sub cycles : " << nSubCycle_ << nl
                << endl;
        }
    }
}


// ************************************************************************* //

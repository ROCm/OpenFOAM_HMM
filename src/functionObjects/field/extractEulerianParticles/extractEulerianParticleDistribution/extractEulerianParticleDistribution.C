/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2016 OpenCFD Ltd.
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

#include "extractEulerianParticleDistribution.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(extractEulerianParticleDistribution, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        extractEulerianParticleDistribution,
        dictionary
    );
}
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void
Foam::functionObjects::extractEulerianParticleDistribution::initialiseBins()
{
    DebugInFunction << endl;

    const faceZone& fz = mesh_.faceZones()[zoneID_];

    // Agglomerate faceZone faces into nInjectorLocations_ global locations
    const indirectPrimitivePatch patch
    (
        IndirectList<face>(mesh_.faces(), fz),
        mesh_.points()
    );

    const label nFaces = fz.size();
    label nLocations = nInjectorLocations_;

    if (Pstream::parRun())
    {
        label nGlobalFaces = returnReduce(nFaces, sumOp<label>());
        scalar fraction = scalar(nFaces)/scalar(nGlobalFaces);
        nLocations = ceil(fraction*nInjectorLocations_);
        if (debug)
        {
            Pout<< "nFaces:" << nFaces
                << ", nGlobalFaces:" << nGlobalFaces
                << ", fraction:" << fraction
                << ", nLocations:" << nLocations
                << endl;
        }
    }

    pairPatchAgglomeration ppa(patch, 10, 50, nLocations, labelMax, 180);
    ppa.agglomerate();

    label nCoarseFaces = 0;
    if (nFaces != 0)
    {
        fineToCoarseAddr_ = ppa.restrictTopBottomAddressing();
        nCoarseFaces = max(fineToCoarseAddr_) + 1;
        coarseToFineAddr_ = invertOneToMany(nCoarseFaces, fineToCoarseAddr_);

        // Set coarse face centres as area average of fine face centres
        const vectorField& faceCentres = mesh_.faceCentres();
        const vectorField& faceAreas = mesh_.faceAreas();
        coarsePosition_.setSize(coarseToFineAddr_.size());
        forAll(coarseToFineAddr_, coarsei)
        {
            const labelList& fineFaces = coarseToFineAddr_[coarsei];

            scalar sumArea = 0;
            vector averagePosition(vector::zero);
            forAll(fineFaces, i)
            {
                label facei = fz[fineFaces[i]];
                scalar magSf = mag(faceAreas[facei]);
                sumArea += magSf;
                averagePosition += magSf*faceCentres[facei];
            }
            coarsePosition_[coarsei] = averagePosition/sumArea;
        }
    }

    // Create global addressing for coarse face addressing
    globalCoarseFaces_ = globalIndex(nCoarseFaces);

    Info<< "Created " << returnReduce(nCoarseFaces, sumOp<label>())
        << " coarse faces" << endl;
}


void
Foam::functionObjects::extractEulerianParticleDistribution::
writeBinnedParticleData()
{
    DebugInFunction << endl;

    // Gather particles ready for collection from all procs
    List<List<eulerianParticle> > allProcParticles(Pstream::nProcs());
    allProcParticles[Pstream::myProcNo()] = collectedParticles_;
    Pstream::gatherList(allProcParticles);
    Pstream::scatterList(allProcParticles);
    List<eulerianParticle> allParticles =
        ListListOps::combine<List<eulerianParticle> >
        (
            allProcParticles,
            accessOp<List<eulerianParticle> >()
        );


    // Determine coarse face index (global) and position for each particle
    label nCoarseFaces = globalCoarseFaces_.size();
    List<label> particleCoarseFacei(allParticles.size(), -1);
    List<point> particleCoarseFacePosition(nCoarseFaces, point::min);

    forAll(allParticles, particlei)
    {
        const eulerianParticle& p = allParticles[particlei];
        label globalFaceHiti = p.globalFaceIHit;

        if (globalFaces_.isLocal(globalFaceHiti))
        {
            label localFacei = globalFaces_.toLocal(globalFaceHiti);
            label coarseFacei = fineToCoarseAddr_[localFacei];
            label globalCoarseFacei = globalCoarseFaces_.toGlobal(coarseFacei);

            particleCoarseFacei[particlei] = globalCoarseFacei;
            particleCoarseFacePosition[globalCoarseFacei] =
                coarsePosition_[coarseFacei];
        }
    }
    Pstream::listCombineGather(particleCoarseFacei, maxEqOp<label>());
    Pstream::listCombineGather(particleCoarseFacePosition, maxEqOp<point>());

    // Write the agglomerated particle data to file
    DynamicList<label> processedCoarseFaces;
    if (Pstream::master())
    {
        fileName baseDir(dictBaseFileDir()/name());

        IOdictionary dict
        (
            IOobject
            (
                "particleDistribution",
                obr_.time().timeName(),
                baseDir,
                obr_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            )
        );

        labelListList coarseFaceToParticle =
            invertOneToMany(nCoarseFaces, particleCoarseFacei);

        // Process the allParticles per coarse face
        forAll(coarseFaceToParticle, globalCoarseFacei)
        {
            const List<label>& particleIDs =
                coarseFaceToParticle[globalCoarseFacei];

            const label nParticle = particleIDs.size();

            if (nParticle == 0)
            {
                continue;
            }

            Field<scalar> pd(particleIDs.size());
            scalar sumV = 0;
            vector sumVU = vector::zero;
            scalar startTime = GREAT;
            scalar endTime = -GREAT;
            forAll(particleIDs, i)
            {
                const label particlei = particleIDs[i];
                const eulerianParticle& p = allParticles[particlei];
                scalar pDiameter = cbrt(6*p.V/constant::mathematical::pi);
                pd[i] = pDiameter;
                sumV += p.V;
                sumVU += p.VU;

                startTime = min(startTime, outputTimes_[p.timeIndex]);
                endTime = max(endTime, outputTimes_[p.timeIndex + 1]);
            }

            if (sumV < ROOTVSMALL)
            {
                // Started collecting particle info, but not accumulated any
                // volume yet
                continue;
            }


            distributionModels::binned binnedDiameters
            (
                pd,
                distributionBinWidth_,
                rndGen_
            );

            // Velocity info hard-coded to volume average
            vector Uave = sumVU/sumV;

            dictionary particleDict;
            particleDict.add("startTime", startTime);
            particleDict.add("endTime", endTime);
            particleDict.add("nParticle", nParticle);
            particleDict.add
            (
                "position",
                particleCoarseFacePosition[globalCoarseFacei]
            );
            particleDict.add("volume", sumV);
            particleDict.add("U", Uave);
            particleDict.add
            (
                "binnedDistribution",
                binnedDiameters.writeDict("distribution")
            );
            dict.add
            (
                word("sample" + Foam::name(globalCoarseFacei)),
                particleDict
            );

            processedCoarseFaces.append(globalCoarseFacei);
        }

        dict.regIOobject::write();
    }


    if (resetDistributionOnWrite_)
    {
        // Remove particles from processed coarse faces from collectedParticles_
        Pstream::scatter(processedCoarseFaces);
        labelHashSet processedFaces(processedCoarseFaces);
        DynamicList<eulerianParticle> nonProcessedParticles;
        forAll(collectedParticles_, particlei)
        {
            const eulerianParticle& p = collectedParticles_[particlei];
            label localFacei = globalFaces_.toLocal(p.globalFaceIHit);
            label coarseFacei = fineToCoarseAddr_[localFacei];
            label globalCoarseFacei = globalCoarseFaces_.toGlobal(coarseFacei);
            if (!processedFaces.found(globalCoarseFacei))
            {
                nonProcessedParticles.append(p);
            }
        }
        collectedParticles_.transfer(nonProcessedParticles);
    }
}


// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

Foam::functionObjects::extractEulerianParticleDistribution::
extractEulerianParticleDistribution
(
    const word& name,
    const Time& runTime,
    const dictionary& dict,
    const bool readFields
)
:
    extractEulerianParticleDistribution(name, runTime, dict, false),
    nInjectorLocations_(0),
    resetDistributionOnWrite_(false),
    distributionBinWidth_(0)
    fineToCoarseAddr_(),
    coarseToFineAddr_(),
    coarsePosition_(),
    globalCoarseFaces_(),
    rndGen_(1234, -1)
{
    // We need to cache the collected particles in order to determine the
    // distributions
    cacheCollectedParticles_ = true;

    if (readFields)
    {
        read(dict);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::extractEulerianParticleDistribution::
~extractEulerianParticleDistribution()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


bool Foam::functionObjects::extractEulerianParticleDistribution::read
(
    const dictionary& dict
)
{
    DebugInFunction << endl;

    if (extractEulerianParticles::read(dict))
    {
        dict.lookup("nLocations") >> nInjectorLocations_;
        dict.lookup("distributionBinWidth") >> distributionBinWidth_;
        dict.lookup("resetDistributionOnWrite") >> resetDistributionOnWrite_;

        initialiseBins();

        return true;
    }

    return false;
}


bool Foam::functionObjects::extractEulerianParticleDistribution::execute()
{
    DebugInFunction << endl;

    return extractEulerianParticles::execute();
}


bool Foam::functionObjects::extractEulerianParticleDistribution::write()
{
    DebugInFunction << endl;

    if (extractEulerianParticles::write())
    {
        writeBinnedParticleData();
        return true;
    }

    return false;
}


// ************************************************************************* //

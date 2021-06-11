/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2021 OpenCFD Ltd.
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

#include "PatchParticleHistogram.H"
#include "Pstream.H"
#include "stringListOps.H"
#include "ListOps.H"
#include "ListListOps.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class CloudType>
Foam::label Foam::PatchParticleHistogram<CloudType>::applyToPatch
(
    const label globalPatchi
) const
{
    return patchIDs_.find(globalPatchi);
}


// * * * * * * * * * * * * * protected Member Functions  * * * * * * * * * * //

template<class CloudType>
void Foam::PatchParticleHistogram<CloudType>::write()
{
    forAll(times_, i)
    {
        List<List<scalar>> procTimes(Pstream::nProcs());
        procTimes[Pstream::myProcNo()] = times_[i];
        Pstream::gatherList(procTimes);

        List<List<scalar>> procDiameters(Pstream::nProcs());
        procDiameters[Pstream::myProcNo()] = patchDiameters_[i];
        Pstream::gatherList(procDiameters);

        List<List<scalar>> procParticles(Pstream::nProcs());
        procParticles[Pstream::myProcNo()] = patchParticles_[i];
        Pstream::gatherList(procParticles);

        if (Pstream::master())
        {
            const fvMesh& mesh = this->owner().mesh();

            mkDir(this->writeTimeDir());

            const word& patchName = mesh.boundaryMesh()[patchIDs_[i]].name();

            OFstream patchOutFile
            (
                this->writeTimeDir()/patchName + ".post",
                IOstream::ASCII,
                IOstream::currentVersion,
                mesh.time().writeCompression()
            );

            List<scalar> globalTimes;
            globalTimes = ListListOps::combine<List<scalar>>
            (
                procTimes,
                accessOp<List<scalar>>()
            );

            List<scalar> globalDiameters;
            globalDiameters = ListListOps::combine<List<scalar>>
            (
                procDiameters,
                accessOp<List<scalar>>()
            );

            List<scalar> globalParticles;
            globalParticles = ListListOps::combine<List<scalar>>
            (
                procParticles,
                accessOp<List<scalar>>()
            );

            // Compute histogram
            List<scalar> nParticles(nBins_, Zero);
            forAll(globalDiameters, j)
            {
                const label bini = (globalDiameters[j] - min_)/delta_;
                if (bini >= 0 && bini < nBins_)
                {
                    nParticles[bini] += globalParticles[j];
                    nParticlesCumulative_[i][bini] += globalParticles[j];
                }
            }

            patchOutFile
                << "# nBin=" << nBins_
                << "; min="  << min_
                << "; max="  << max_ << nl
                << "# dEdge1    dEdge2    nParticles    nParticlesCumulative"
                << endl;

            forAll(nParticles, j)
            {
                patchOutFile
                    << binEdges_[j]
                    << ' '
                    << binEdges_[j + 1]
                    << ' '
                    << nParticles[j]
                    << ' '
                    << nParticlesCumulative_[i][j]
                    << nl;
            }
        }

        times_[i].clearStorage();
        patchDiameters_[i].clearStorage();
        patchParticles_[i].clearStorage();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PatchParticleHistogram<CloudType>::PatchParticleHistogram
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    CloudFunctionObject<CloudType>(dict, owner, modelName, typeName),
    nBins_(dict.getCheck<label>("nBins", labelMinMax::ge(1))),
    min_(dict.getScalar("min")),
    max_(dict.getScalar("max")),
    delta_((max_ - min_)/scalar(nBins_)),
    maxStoredParcels_(dict.getScalar("maxStoredParcels")),
    binEdges_(nBins_ + 1),
    patchIDs_(),
    nParticlesCumulative_(),
    times_(),
    patchDiameters_(),
    patchParticles_()
{
    if (min_ >= max_)
    {
        FatalIOErrorInFunction(dict)
            << "Histogram minimum = " << min_
            << ", cannot be larger than histogram maximum = " << max_
            << exit(FatalIOError);
    }

    if (maxStoredParcels_ <= 0)
    {
        FatalIOErrorInFunction(dict)
            << "maxStoredParcels = " << maxStoredParcels_
            << ", cannot be equal to or less than zero"
            << exit(FatalIOError);
    }

    // Compute histogram-bin properties
    binEdges_[0] = min_;
    for (label i = 0; i < nBins_; ++i)
    {
        const scalar next = min_ + (i+1)*delta_;
        binEdges_[i+1] = next;
    }

    // Compute histogram-patch properties
    const wordRes patchMatcher(dict.get<wordRes>("patches"));

    patchIDs_ = patchMatcher.matching(owner.mesh().boundaryMesh().names());

    if (patchIDs_.empty())
    {
        FatalIOErrorInFunction(dict)
            << "No matching patches found: "
            << flatOutput(patchMatcher) << nl
            << exit(FatalIOError);
    }

    nParticlesCumulative_ =
        List<List<scalar>>(patchIDs_.size(), List<scalar>(nBins_, Zero));

    times_.setSize(patchIDs_.size());
    patchDiameters_.setSize(patchIDs_.size());
    patchParticles_.setSize(patchIDs_.size());
}


template<class CloudType>
Foam::PatchParticleHistogram<CloudType>::PatchParticleHistogram
(
    const PatchParticleHistogram<CloudType>& ppm
)
:
    CloudFunctionObject<CloudType>(ppm),
    nBins_(ppm.nBins_),
    min_(ppm.min_),
    max_(ppm.max_),
    delta_(ppm.delta_),
    maxStoredParcels_(ppm.maxStoredParcels_),
    binEdges_(ppm.binEdges_),
    patchIDs_(ppm.patchIDs_),
    nParticlesCumulative_(ppm.nParticlesCumulative_),
    times_(ppm.times_),
    patchDiameters_(ppm.patchDiameters_),
    patchParticles_(ppm.patchParticles_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::PatchParticleHistogram<CloudType>::postPatch
(
    const parcelType& p,
    const polyPatch& pp,
    bool&
)
{
    const label patchi = pp.index();
    const label localPatchi = applyToPatch(patchi);

    if (localPatchi != -1 && times_[localPatchi].size() < maxStoredParcels_)
    {
        times_[localPatchi].append(this->owner().time().value());
        patchDiameters_[localPatchi].append(p.d());
        patchParticles_[localPatchi].append(p.nParticle());
    }
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "PatchPostProcessing.H"
#include "Pstream.H"
#include "stringListOps.H"
#include "ListOps.H"
#include "ListListOps.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class CloudType>
Foam::label Foam::PatchPostProcessing<CloudType>::applyToPatch
(
    const label globalPatchi
) const
{
    return patchIDs_.find(globalPatchi);
}


// * * * * * * * * * * * * * protected Member Functions  * * * * * * * * * * //

template<class CloudType>
void Foam::PatchPostProcessing<CloudType>::write()
{
    forAll(patchData_, i)
    {
        List<List<scalar>> procTimes(Pstream::nProcs());
        procTimes[Pstream::myProcNo()] = times_[i];
        Pstream::gatherList(procTimes);

        List<List<string>> procData(Pstream::nProcs());
        procData[Pstream::myProcNo()] = patchData_[i];
        Pstream::gatherList(procData);

        if (Pstream::master())
        {
            const fvMesh& mesh = this->owner().mesh();

            // Create directory if it doesn't exist
            mkDir(this->writeTimeDir());

            const word& patchName = mesh.boundaryMesh()[patchIDs_[i]].name();

            OFstream patchOutFile
            (
                this->writeTimeDir()/patchName + ".post",
                IOstream::ASCII,
                IOstream::currentVersion,
                mesh.time().writeCompression()
            );

            List<string> globalData;
            globalData = ListListOps::combine<List<string>>
            (
                procData,
                accessOp<List<string>>()
            );

            List<scalar> globalTimes;
            globalTimes = ListListOps::combine<List<scalar>>
            (
                procTimes,
                accessOp<List<scalar>>()
            );

            labelList indices(sortedOrder(globalTimes));

            string header("# Time currentProc " + header_);
            patchOutFile<< header.c_str() << nl;

            forAll(globalTimes, i)
            {
                label dataI = indices[i];

                patchOutFile
                    << globalTimes[dataI] << ' '
                    << globalData[dataI].c_str()
                    << nl;
            }
        }

        patchData_[i].clearStorage();
        times_[i].clearStorage();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PatchPostProcessing<CloudType>::PatchPostProcessing
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    CloudFunctionObject<CloudType>(dict, owner, modelName, typeName),
    maxStoredParcels_(this->coeffDict().getScalar("maxStoredParcels")),
    fields_(),
    patchIDs_(),
    times_(),
    patchData_(),
    header_()
{
    // The "fields" filter is optional
    this->coeffDict().readIfPresent("fields", fields_);

    // The "patches" are required
    const wordRes patchMatcher(this->coeffDict().lookup("patches"));

    patchIDs_ = patchMatcher.matching(owner.mesh().boundaryMesh().names());

    if (patchIDs_.empty())
    {
        WarningInFunction
            << "No matching patches found: "
            << flatOutput(patchMatcher) << nl;
    }

    if (debug)
    {
        Info<< "Post-process fields "
            << flatOutput(fields_) << nl;

        Info<< "On patches (";

        for (const label patchi : patchIDs_)
        {
            Info<< ' ' << owner.mesh().boundaryMesh()[patchi].name();
        }

        Info<< " )" << nl;
    }

    patchData_.setSize(patchIDs_.size());
    times_.setSize(patchIDs_.size());
}


template<class CloudType>
Foam::PatchPostProcessing<CloudType>::PatchPostProcessing
(
    const PatchPostProcessing<CloudType>& ppm
)
:
    CloudFunctionObject<CloudType>(ppm),
    maxStoredParcels_(ppm.maxStoredParcels_),
    fields_(ppm.fields_),
    patchIDs_(ppm.patchIDs_),
    times_(ppm.times_),
    patchData_(ppm.patchData_),
    header_(ppm.header_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::PatchPostProcessing<CloudType>::postPatch
(
    const parcelType& p,
    const polyPatch& pp,
    bool&
)
{
    const label patchi = pp.index();
    const label localPatchi = applyToPatch(patchi);

    if (header_.empty())
    {
        OStringStream data;
        p.writeProperties(data, fields_, " ", true);
        header_ = data.str();
    }

    if (localPatchi != -1 && patchData_[localPatchi].size() < maxStoredParcels_)
    {
        times_[localPatchi].append(this->owner().time().value());

        OStringStream data;
        data<< Pstream::myProcNo();
        p.writeProperties(data, fields_, " ", false);

        patchData_[localPatchi].append(data.str());
    }
}


// ************************************************************************* //

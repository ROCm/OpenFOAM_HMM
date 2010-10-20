/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2010 OpenCFD Ltd.
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

#include "PatchPostProcessing.H"
#include "Pstream.H"
#include "ListListOps.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template <class CloudType>
Foam::label Foam::PatchPostProcessing<CloudType>::applyToPatch
(
    const label globalPatchI
) const
{
    forAll(patchIds_, patchI)
    {
        if (patchIds_[patchI] == globalPatchI)
        {
            return patchI;
        }
    }

    return -1;
}


// * * * * * * * * * * * * * protected Member Functions  * * * * * * * * * * //

template<class CloudType>
void Foam::PatchPostProcessing<CloudType>::write()
{
    forAll(patchData_, patchI)
    {
        List<List<string> > procData(Pstream::nProcs());
        procData[Pstream::myProcNo()] = patchData_[patchI];

        Pstream::gatherList(procData);

        if (Pstream::master())
        {
            fileName outputDir = this->owner().time().path();

            if (Pstream::parRun())
            {
                // Put in undecomposed case (Note: gives problems for
                // distributed data running)
                outputDir =
                    outputDir/".."/"postProcessing"/cloud::prefix/
                    this->owner().name()/this->owner().time().timeName();
            }
            else
            {
                outputDir =
                    outputDir/"postProcessing"/cloud::prefix/
                    this->owner().name()/this->owner().time().timeName();
            }

            // Create directory if it doesn't exist
            mkDir(outputDir);

            OFstream patchOutFile
            (
                outputDir/patchNames_[patchI] + ".post",
                IOstream::ASCII,
                IOstream::currentVersion,
                this->owner().time().writeCompression()
            );

            List<string> globalData;
            globalData = ListListOps::combine<List<string> >
            (
                procData,
                accessOp<List<string> >()
            );
            sort(globalData);

            patchOutFile<< "# Time " + parcelType::propHeader << nl;

            forAll(globalData, i)
            {
                patchOutFile<< globalData[i].c_str() << nl;
            }
        }

        patchData_[patchI].clearStorage();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PatchPostProcessing<CloudType>::PatchPostProcessing
(
    const dictionary& dict,
    CloudType& owner
)
:
    PostProcessingModel<CloudType>(dict, owner, typeName),
    maxStoredParcels_(readLabel(this->coeffDict().lookup("maxStoredParcels"))),
    patchNames_(this->coeffDict().lookup("patches")),
    patchData_(patchNames_.size()),
    patchIds_(patchNames_.size())
{
    const polyBoundaryMesh& bMesh = this->owner().mesh().boundaryMesh();
    forAll(patchNames_, patchI)
    {
        const label id = bMesh.findPatchID(patchNames_[patchI]);
        if (id < 0)
        {
            FatalErrorIn
            (
                "PatchPostProcessing<CloudType>::PatchPostProcessing"
                "("
                    "const dictionary&, "
                    "CloudType& owner"
                ")"
            )<< "Requested patch " << patchNames_[patchI] << " not found" << nl
             << "Available patches are: " << bMesh.names() << nl
             << exit(FatalError);
        }
        patchIds_[patchI] = id;
    }
}


template<class CloudType>
Foam::PatchPostProcessing<CloudType>::PatchPostProcessing
(
    const PatchPostProcessing<CloudType>& ppm
)
:
    PostProcessingModel<CloudType>(ppm),
    maxStoredParcels_(ppm.maxStoredParcels_),
    patchNames_(ppm.patchNames_),
    patchData_(ppm.patchData_),
    patchIds_(ppm.patchIds_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PatchPostProcessing<CloudType>::~PatchPostProcessing()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::PatchPostProcessing<CloudType>::postPatch
(
    const parcelType& p,
    const label patchI
)
{
    label localPatchI = applyToPatch(patchI);
    if (localPatchI >= 0 && patchData_[localPatchI].size() < maxStoredParcels_)
    {
        OStringStream data;
        data<< this->owner().time().timeName() << ' ' << p;
        patchData_[localPatchI].append(data.str());
    }
}


template<class CloudType>
void Foam::PatchPostProcessing<CloudType>::postFace(const parcelType&)
{}


// ************************************************************************* //

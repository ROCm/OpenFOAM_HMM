/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2011 OpenCFD Ltd.
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
#include "stringListOps.H"
#include "ListListOps.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class CloudType>
Foam::label Foam::PatchPostProcessing<CloudType>::applyToPatch
(
    const label globalPatchI
) const
{
    label patchI = 0;
    forAllConstIter(labelHashSet, patchIDs_, iter)
    {
        if (iter.key() == globalPatchI)
        {
            return patchI;
        }

        patchI++;
    }

    return -1;
}


// * * * * * * * * * * * * * protected Member Functions  * * * * * * * * * * //

template<class CloudType>
void Foam::PatchPostProcessing<CloudType>::write()
{
    forAll(patchData_, i)
    {
        List<List<string> > procData(Pstream::nProcs());
        procData[Pstream::myProcNo()] = patchData_[i];

        Pstream::gatherList(procData);

        if (Pstream::master())
        {
            const fvMesh& mesh = this->owner().mesh();

            fileName outputDir = mesh.time().path();

            if (Pstream::parRun())
            {
                // Put in undecomposed case (Note: gives problems for
                // distributed data running)
                outputDir =
                    outputDir/".."/"postProcessing"/cloud::prefix/
                    this->owner().name()/mesh.time().timeName();
            }
            else
            {
                outputDir =
                    outputDir/"postProcessing"/cloud::prefix/
                    this->owner().name()/mesh.time().timeName();
            }

            // Create directory if it doesn't exist
            mkDir(outputDir);

            const word& patchName = mesh.boundaryMesh()[patchIDs_[i]].name();

            OFstream patchOutFile
            (
                outputDir/patchName + ".post",
                IOstream::ASCII,
                IOstream::currentVersion,
                mesh.time().writeCompression()
            );

            List<string> globalData;
            globalData = ListListOps::combine<List<string> >
            (
                procData,
                accessOp<List<string> >()
            );
            sort(globalData);

            patchOutFile<< "# Time " + parcelType::propHeader << nl;

            forAll(globalData, dataI)
            {
                patchOutFile<< globalData[dataI].c_str() << nl;
            }
        }

        patchData_[i].clearStorage();
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
    patchIDs_(),
    patchData_()
{
    const wordList allPatchNames = owner.mesh().boundaryMesh().names();
    wordList patchName(this->coeffDict().lookup("patches"));

    forAllReverse(patchName, i)
    {
        labelList patchIDs = findStrings(patchName[i], allPatchNames);

        if (patchIDs.empty())
        {
            WarningIn
            (
                "Foam::PatchPostProcessing<CloudType>::PatchPostProcessing"
                "("
                    "const dictionary&, "
                    "CloudType& "
                ")"
            )   << "Cannot find any patch names matching " << patchName[i]
                << endl;
        }

        forAll(patchIDs, j)
        {
            patchIDs_.insert(patchIDs[j]);
        }
    }

    patchData_.setSize(patchIDs_.size());
}


template<class CloudType>
Foam::PatchPostProcessing<CloudType>::PatchPostProcessing
(
    const PatchPostProcessing<CloudType>& ppm
)
:
    PostProcessingModel<CloudType>(ppm),
    maxStoredParcels_(ppm.maxStoredParcels_),
    patchIDs_(ppm.patchIDs_),
    patchData_(ppm.patchData_)
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
    const label localPatchI = applyToPatch(patchI);
    if (localPatchI != -1 && patchData_[localPatchI].size() < maxStoredParcels_)
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

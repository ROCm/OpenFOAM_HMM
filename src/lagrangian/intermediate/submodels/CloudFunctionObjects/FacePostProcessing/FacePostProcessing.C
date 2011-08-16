/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "FacePostProcessing.H"
#include "Pstream.H"
#include "ListListOps.H"
#include "surfaceWriter.H"
#include "globalIndex.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class CloudType>
void Foam::FacePostProcessing<CloudType>::applyToFace
(
    const label faceIn,
    label& zoneI,
    label& faceI
) const
{
    const faceZoneMesh& fzm = this->owner().mesh().faceZones();

    forAll(faceZoneIDs_, i)
    {
        const faceZone& fz = fzm[faceZoneIDs_[i]];
        forAll(fz, j)
        {
            if (fz[j] == faceIn)
            {
                zoneI = i;
                faceI = j;
                return;
            }
        }
    }
}


template<class CloudType>
void Foam::FacePostProcessing<CloudType>::makeLogFile
(
    const word& zoneName,
    const label zoneI,
    const label nFaces,
    const scalar totArea
)
{
    // Create the output file if not already created
    if (log_)
    {
        if (debug)
        {
            Info<< "Creating output file." << endl;
        }

        if (Pstream::master())
        {
            const fileName logDir = outputDir_/this->owner().time().timeName();

            // Create directory if does not exist
            mkDir(logDir);

            // Open new file at start up
            outputFilePtr_.set
            (
                zoneI,
                new OFstream(logDir/(type() + '_' + zoneName + ".dat"))
            );

            outputFilePtr_[zoneI]
                << "# Source    : " << type() << nl
                << "# Face zone : " << zoneName << nl
                << "# Faces     : " << nFaces << nl
                << "# Area      : " << totArea << nl
                << "# Time" << tab << "mass" << tab << "massFlux" << endl;
        }
    }
}


template<class CloudType>
void Foam::FacePostProcessing<CloudType>::write()
{
    const fvMesh& mesh = this->owner().mesh();
    const Time& time = mesh.time();
    const faceZoneMesh& fzm = mesh.faceZones();
    const scalar dt = time.deltaTValue();

    totalTime_ += dt;

    const scalar alpha = (totalTime_ - dt)/totalTime_;
    const scalar beta = dt/totalTime_;

    forAll(faceZoneIDs_, zoneI)
    {
        massTotal_[zoneI] += mass_[zoneI];
        massFlux_[zoneI] = alpha*massFlux_[zoneI] + beta*mass_[zoneI]/dt;
    }

    const label procI = Pstream::myProcNo();

    Info<< "particleFaceFlux output:" << nl;

    List<scalarField> zoneMassTotal(mass_.size());
    List<scalarField> zoneMassFlux(massFlux_.size());
    forAll(faceZoneIDs_, zoneI)
    {
        const word& zoneName = fzm[faceZoneIDs_[zoneI]].name();

        scalarListList allProcMass(Pstream::nProcs());
        allProcMass[procI] = massTotal_[zoneI];
        Pstream::gatherList(allProcMass);
        zoneMassTotal[zoneI] =
            ListListOps::combine<scalarList>
            (
                allProcMass, accessOp<scalarList>()
            );
        const scalar sumMassTotal = sum(zoneMassTotal[zoneI]);

        scalarListList allProcMassFlux(Pstream::nProcs());
        allProcMassFlux[procI] = massFlux_[zoneI];
        Pstream::gatherList(allProcMassFlux);
        zoneMassFlux[zoneI] =
            ListListOps::combine<scalarList>
            (
                allProcMassFlux, accessOp<scalarList>()
            );
        const scalar sumMassFlux = sum(zoneMassFlux[zoneI]);

        Info<< "    " << zoneName
            << ": total mass = " << sumMassTotal
            << "; average mass flux = " << sumMassFlux
            << nl;

        if (outputFilePtr_.set(zoneI))
        {
            OFstream& os = outputFilePtr_[zoneI];
            os  << time.timeName() << token::TAB << sumMassTotal << token::TAB
                <<  sumMassFlux<< endl;
        }
    }

    Info<< endl;


    if (surfaceFormat_ != "none")
    {
        forAll(faceZoneIDs_, zoneI)
        {
            const faceZone& fZone = fzm[faceZoneIDs_[zoneI]];

            labelList pointToGlobal;
            labelList uniqueMeshPointLabels;
            autoPtr<globalIndex> globalPointsPtr =
                mesh.globalData().mergePoints
                (
                    fZone().meshPoints(),
                    fZone().meshPointMap(),
                    pointToGlobal,
                    uniqueMeshPointLabels
                );

            pointField uniquePoints(mesh.points(), uniqueMeshPointLabels);
            List<pointField> allProcPoints(Pstream::nProcs());
            allProcPoints[procI] = uniquePoints;
            Pstream::gatherList(allProcPoints);

            faceList faces(fZone().localFaces());
            forAll(faces, i)
            {
                inplaceRenumber(pointToGlobal, faces[i]);
            }
            List<faceList> allProcFaces(Pstream::nProcs());
            allProcFaces[procI] = faces;
            Pstream::gatherList(allProcFaces);

            if (Pstream::master())
            {
                pointField allPoints
                (
                    ListListOps::combine<pointField>
                    (
                        allProcPoints, accessOp<pointField>()
                    )
                );

                faceList allFaces
                (
                    ListListOps::combine<faceList>
                    (
                        allProcFaces, accessOp<faceList>()
                    )
                );

                autoPtr<surfaceWriter> writer
                (
                    surfaceWriter::New(surfaceFormat_)
                );

                writer->write
                (
                    outputDir_/time.timeName(),
                    fZone.name(),
                    allPoints,
                    allFaces,
                    "massTotal",
                    zoneMassTotal[zoneI],
                    false
                );

                writer->write
                (
                    outputDir_/time.timeName(),
                    fZone.name(),
                    allPoints,
                    allFaces,
                    "massFlux",
                    zoneMassFlux[zoneI],
                    false
                );
            }
        }
    }


    if (resetOnWrite_)
    {
        forAll(faceZoneIDs_, zoneI)
        {
            massFlux_[zoneI] = 0.0;
        }
        totalTime_ = 0.0;
    }

    forAll(mass_, zoneI)
    {
        mass_[zoneI] = 0.0;
    }

    // writeProperties();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::FacePostProcessing<CloudType>::FacePostProcessing
(
    const dictionary& dict,
    CloudType& owner
)
:
    CloudFunctionObject<CloudType>(dict, owner, typeName),
    faceZoneIDs_(),
    surfaceFormat_(this->coeffDict().lookup("surfaceFormat")),
    resetOnWrite_(this->coeffDict().lookup("resetOnWrite")),
    totalTime_(0.0),
    mass_(),
    massTotal_(),
    massFlux_(),
    log_(this->coeffDict().lookup("log")),
    outputFilePtr_(),
    outputDir_(owner.mesh().time().path())
{
    wordList faceZoneNames(this->coeffDict().lookup("faceZones"));
    mass_.setSize(faceZoneNames.size());
    massTotal_.setSize(faceZoneNames.size());
    massFlux_.setSize(faceZoneNames.size());

    outputFilePtr_.setSize(faceZoneNames.size());

    if (Pstream::parRun())
    {
        // Put in undecomposed case (Note: gives problems for
        // distributed data running)
        outputDir_ =
            outputDir_/".."/"postProcessing"/cloud::prefix/owner.name();
    }
    else
    {
        outputDir_ =
            outputDir_/"postProcessing"/cloud::prefix/owner.name();
    }

    DynamicList<label> zoneIDs;
    const faceZoneMesh& fzm = owner.mesh().faceZones();
    const surfaceScalarField& magSf = owner.mesh().magSf();
    forAll(faceZoneNames, i)
    {
        const word& zoneName = faceZoneNames[i];
        label zoneI = fzm.findZoneID(zoneName);
        if (zoneI != -1)
        {
            zoneIDs.append(zoneI);
            const faceZone& fz = fzm[zoneI];
            label nFaces = returnReduce(fz.size(), sumOp<label>());
            mass_[i].setSize(nFaces, 0.0);
            massTotal_[i].setSize(nFaces, 0.0);
            massFlux_[i].setSize(nFaces, 0.0);
            Info<< "        " << zoneName << " faces: " << nFaces << nl;

            scalar totArea = 0.0;
            forAll(fz, j)
            {
                totArea += magSf[fz[j]];
            }
            totArea = returnReduce(totArea, sumOp<scalar>());

            makeLogFile(zoneName, i, nFaces, totArea);
        }
    }

    faceZoneIDs_.transfer(zoneIDs);

    // readProperties(); AND initialise mass... fields
}


template<class CloudType>
Foam::FacePostProcessing<CloudType>::FacePostProcessing
(
    const FacePostProcessing<CloudType>& pff
)
:
    CloudFunctionObject<CloudType>(pff),
    faceZoneIDs_(pff.faceZoneIDs_),
    surfaceFormat_(pff.surfaceFormat_),
    resetOnWrite_(pff.resetOnWrite_),
    totalTime_(pff.totalTime_),
    mass_(pff.mass_),
    massTotal_(pff.massTotal_),
    massFlux_(pff.massFlux_),
    log_(pff.log_),
    outputFilePtr_(),
    outputDir_(pff.outputDir_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::FacePostProcessing<CloudType>::~FacePostProcessing()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::FacePostProcessing<CloudType>::postPatch
(
    const parcelType&,
    const label
)
{}


template<class CloudType>
void Foam::FacePostProcessing<CloudType>::postFace(const parcelType& p)
{
    if
    (
        this->owner().solution().output()
     || this->owner().solution().transient()
    )
    {
        label zoneI = -1;
        label faceI = -1;
        applyToFace(p.face(), zoneI, faceI);

        if ((zoneI != -1) && (faceI != -1))
        {
            mass_[zoneI][faceI] += p.mass()*p.nParticle();
        }
    }
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2011 OpenCFD Ltd.
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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class CloudType>
Foam::label Foam::FacePostProcessing<CloudType>::applyToFace
(
    const label faceI
) const
{
    forAll(fZone_, i)
    {
        if (fZone_[i] == faceI)
        {
            return i;
        }
    }

    return -1;
}


template<class CloudType>
void Foam::FacePostProcessing<CloudType>::write()
{
    const fvMesh& mesh = this->owner().mesh();
    const scalar dt = this->owner().time().deltaTValue();

    totalTime_ += dt;

    const scalar alpha = (totalTime_ - dt)/totalTime_;
    const scalar beta = dt/totalTime_;

    massTotal_ += mass_;

    massFlux_ = alpha*massFlux_ + beta*mass_/dt;

    const label procI = Pstream::myProcNo();

    scalarListList allProcMass(Pstream::nProcs());
    allProcMass[procI].setSize(massTotal_.size());
    allProcMass[procI] = massTotal_;
    Pstream::gatherList(allProcMass);
    scalarList allMass
    (
        ListListOps::combine<scalarList>
        (
            allProcMass, accessOp<scalarList>()
        )
    );

    scalarListList allProcMassFlux(Pstream::nProcs());
    allProcMassFlux[procI].setSize(massFlux_.size());
    allProcMassFlux[procI] = massFlux_;
    Pstream::gatherList(allProcMassFlux);
    scalarList allMassFlux
    (
        ListListOps::combine<scalarList>
        (
            allProcMassFlux, accessOp<scalarList>()
        )
    );

    Info<< "particleFaceFlux output:" << nl
        << "    total mass      = " << sum(allMass) << nl
        << "    average mass flux = " << sum(allMassFlux) << nl << endl;


    if (surfaceFormat_ != "none")
    {
        labelList pointToGlobal;
        labelList uniqueMeshPointLabels;
        autoPtr<globalIndex> globalPointsPtr =
            mesh.globalData().mergePoints
            (
                fZone_().meshPoints(),
                fZone_().meshPointMap(),
                pointToGlobal,
                uniqueMeshPointLabels
            );

        pointField uniquePoints(mesh.points(), uniqueMeshPointLabels);
        List<pointField> allProcPoints(Pstream::nProcs());
        allProcPoints[procI].setSize(uniquePoints.size());
        allProcPoints[procI] = uniquePoints;
        Pstream::gatherList(allProcPoints);
        pointField allPoints
        (
            ListListOps::combine<pointField>
            (
                allProcPoints, accessOp<pointField>()
            )
        );

        faceList faces(fZone_().localFaces());
        forAll(faces, i)
        {
            inplaceRenumber(pointToGlobal, faces[i]);
        }
        List<faceList> allProcFaces(Pstream::nProcs());
        allProcFaces[procI].setSize(faces.size());
        allProcFaces[procI] = faces;
        Pstream::gatherList(allProcFaces);
        faceList allFaces
        (
            ListListOps::combine<faceList>
            (
                allProcFaces, accessOp<faceList>()
            )
        );


        if (Pstream::master())
        {
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

            autoPtr<surfaceWriter> writer(surfaceWriter::New(surfaceFormat_));
            writer->write
            (
                outputDir,
                "massTotal",
                allPoints,
                allFaces,
                "massTotal",
                massTotal_,
                false
            );
            writer->write
            (
                outputDir,
                "massFlux",
                allPoints,
                allFaces,
                "massFlux",
                massFlux_,
                false
            );
        }
    }


    if (resetOnWrite_)
    {
        massFlux_ = 0.0;
        totalTime_ = 0.0;
    }

    mass_ = 0.0;

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
    fZone_(owner.mesh().faceZones()[this->coeffDict().lookup("faceZone")]),
    surfaceFormat_(this->coeffDict().lookup("surfaceFormat")),
    resetOnWrite_(this->coeffDict().lookup("resetOnWrite")),
    totalTime_(0.0),
    mass_(),
    massTotal_(),
    massFlux_()
{
    label allFaces = returnReduce(fZone_().size(), sumOp<scalar>());
    Info<< "        Number of faces = " << allFaces << endl;

    mass_.setSize(fZone_.size(), 0.0);

    // readProperties(); AND initialise mass... fields

    massTotal_.setSize(fZone_.size(), 0.0);
    massFlux_.setSize(fZone_.size(), 0.0);
}


template<class CloudType>
Foam::FacePostProcessing<CloudType>::FacePostProcessing
(
    const FacePostProcessing<CloudType>& pff
)
:
    CloudFunctionObject<CloudType>(pff),
    fZone_(pff.fZone_.clone(pff.fZone_.zoneMesh())),
    surfaceFormat_(pff.surfaceFormat_),
    resetOnWrite_(pff.resetOnWrite_),
    totalTime_(pff.totalTime_),
    mass_(pff.mass_),
    massTotal_(pff.massTotal_),
    massFlux_(pff.massFlux_)
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
        const label faceI = applyToFace(p.face());

        if (faceI != -1)
        {
            mass_[faceI] += p.mass()*p.nParticle();
        }
    }
}


// ************************************************************************* //

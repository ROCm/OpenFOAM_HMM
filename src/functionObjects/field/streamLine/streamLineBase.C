/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
    Copyright (C) 2015-2021 OpenCFD Ltd.
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

#include "streamLineBase.H"
#include "fvMesh.H"
#include "ReadFields.H"
#include "sampledSet.H"
#include "globalIndex.H"
#include "mapDistribute.H"
#include "interpolationCellPoint.H"
#include "wallPolyPatch.H"
#include "meshSearchMeshObject.H"
#include "mapPolyMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(streamLineBase, 0);
}
}


const Foam::Enum
<
    Foam::functionObjects::streamLineBase::trackDirType
>
Foam::functionObjects::streamLineBase::trackDirTypeNames
({
    { trackDirType::FORWARD, "forward" },
    { trackDirType::BACKWARD, "backward" },
    { trackDirType::BIDIRECTIONAL, "bidirectional" }
});


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::word&
Foam::functionObjects::streamLineBase::sampledSetAxis() const
{
    if (!sampledSetPtr_)
    {
        sampledSetPoints();
    }

    return sampledSetAxis_;
}


const Foam::sampledSet&
Foam::functionObjects::streamLineBase::sampledSetPoints() const
{
    if (!sampledSetPtr_)
    {
        sampledSetPtr_ = sampledSet::New
        (
            "seedSampleSet",
            mesh_,
            meshSearchMeshObject::New(mesh_),
            dict_.subDict("seedSampleSet")
        );

        sampledSetAxis_ = sampledSetPtr_->axis();
    }

    return *sampledSetPtr_;
}


Foam::autoPtr<Foam::indirectPrimitivePatch>
Foam::functionObjects::streamLineBase::wallPatch() const
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    label nFaces = 0;

    for (const polyPatch& pp : patches)
    {
        //if (!polyPatch::constraintType(pp.type()))
        if (isA<wallPolyPatch>(pp))
        {
            nFaces += pp.size();
        }
    }

    labelList addressing(nFaces);

    nFaces = 0;

    for (const polyPatch& pp : patches)
    {
        //if (!polyPatch::constraintType(pp.type()))
        if (isA<wallPolyPatch>(pp))
        {
            forAll(pp, i)
            {
                addressing[nFaces++] = pp.start()+i;
            }
        }
    }

    return autoPtr<indirectPrimitivePatch>::New
    (
        IndirectList<face>
        (
            mesh_.faces(),
            addressing
        ),
        mesh_.points()
    );
}


void Foam::functionObjects::streamLineBase::initInterpolations
(
    const label nSeeds,
    label& UIndex,
    PtrList<volScalarField>& vsFlds,
    PtrList<interpolation<scalar>>& vsInterp,
    PtrList<volVectorField>& vvFlds,
    PtrList<interpolation<vector>>& vvInterp
)
{
    // Read fields
    label nScalar = 0;
    label nVector = 0;

    for (const word& fieldName : fields_)
    {
        if (foundObject<volScalarField>(fieldName))
        {
            ++nScalar;
        }
        else if (foundObject<volVectorField>(fieldName))
        {
            ++nVector;
        }
        else
        {
            FatalErrorInFunction
                << "Cannot find field " << fieldName << nl
                << "Valid scalar fields are:"
                << flatOutput(mesh_.names(volScalarField::typeName)) << nl
                << "Valid vector fields are:"
                << flatOutput(mesh_.names(volVectorField::typeName))
                << exit(FatalError);
        }
    }
    vsInterp.setSize(nScalar);
    nScalar = 0;
    vvInterp.setSize(nVector);
    nVector = 0;

    for (const word& fieldName : fields_)
    {
        if (foundObject<volScalarField>(fieldName))
        {
            const volScalarField& f = lookupObject<volScalarField>(fieldName);
            vsInterp.set
            (
                nScalar++,
                interpolation<scalar>::New
                (
                    interpolationScheme_,
                    f
                )
            );
        }
        else if (foundObject<volVectorField>(fieldName))
        {
            const volVectorField& f = lookupObject<volVectorField>(fieldName);

            if (f.name() == UName_)
            {
                UIndex = nVector;
            }

            vvInterp.set
            (
                nVector++,
                interpolation<vector>::New
                (
                    interpolationScheme_,
                    f
                )
            );
        }
    }

    // Store the names
    scalarNames_.setSize(vsInterp.size());
    forAll(vsInterp, i)
    {
        scalarNames_[i] = vsInterp[i].psi().name();
    }
    vectorNames_.setSize(vvInterp.size());
    forAll(vvInterp, i)
    {
        vectorNames_[i] = vvInterp[i].psi().name();
    }

    // Check that we know the index of U in the interpolators.

    if (UIndex == -1)
    {
        FatalErrorInFunction
            << "Cannot find field to move particles with : " << UName_ << nl
            << "This field has to be present in the sampled fields " << fields_
            << " and in the objectRegistry."
            << exit(FatalError);
    }

    // Sampled data
    // ~~~~~~~~~~~~

    // Size to maximum expected sizes.
    allTracks_.clear();
    allTracks_.setCapacity(nSeeds);
    allScalars_.setSize(vsInterp.size());
    forAll(allScalars_, i)
    {
        allScalars_[i].clear();
        allScalars_[i].setCapacity(nSeeds);
    }
    allVectors_.setSize(vvInterp.size());
    forAll(allVectors_, i)
    {
        allVectors_[i].clear();
        allVectors_[i].setCapacity(nSeeds);
    }
}


void Foam::functionObjects::streamLineBase::storePoint
(
    const label tracki,

    const scalar w,
    const label lefti,
    const label righti,

    DynamicList<point>& newTrack,
    DynamicList<scalarList>& newScalars,
    DynamicList<vectorList>& newVectors
) const
{
    const label sz = newTrack.size();

    const List<point>& track = allTracks_[tracki];

    newTrack.append((1.0-w)*track[lefti] + w*track[righti]);

    // Scalars
    {
        newScalars.append(scalarList(allScalars_.size()));
        scalarList& newVals = newScalars[sz];

        forAll(allScalars_, scalari)
        {
            const scalarList& trackVals = allScalars_[scalari][tracki];
            newVals[scalari] = (1.0-w)*trackVals[lefti] + w*trackVals[righti];
        }
    }

    // Vectors
    {
        newVectors.append(vectorList(allVectors_.size()));
        vectorList& newVals = newVectors[sz];

        forAll(allVectors_, vectori)
        {
            const vectorList& trackVals = allVectors_[vectori][tracki];
            newVals[vectori] = (1.0-w)*trackVals[lefti] + w*trackVals[righti];
        }
    }
}


// Can split a track into multiple tracks
void Foam::functionObjects::streamLineBase::trimToBox
(
    const treeBoundBox& bb,
    const label tracki,
    PtrList<DynamicList<point>>& newTracks,
    PtrList<DynamicList<scalarList>>& newScalars,
    PtrList<DynamicList<vectorList>>& newVectors
) const
{
    const List<point>& track = allTracks_[tracki];

    if (track.size())
    {
        for
        (
            label segmenti = 1;
            segmenti < track.size();
            segmenti++
        )
        {
            const point& startPt = track[segmenti-1];
            const point& endPt = track[segmenti];

            const vector d(endPt-startPt);
            const scalar magD = mag(d);
            if (magD > ROOTVSMALL)
            {
                if (bb.contains(startPt))
                {
                    // Store 1.0*track[segmenti-1]+0*track[segmenti]
                    storePoint
                    (
                        tracki,

                        0.0,
                        segmenti-1,
                        segmenti,

                        newTracks.last(),
                        newScalars.last(),
                        newVectors.last()
                    );

                    if (!bb.contains(endPt))
                    {
                        point clipPt;
                        if (bb.intersects(endPt, startPt, clipPt))
                        {
                            // End of track. Store point and interpolated
                            // values
                            storePoint
                            (
                                tracki,

                                mag(clipPt-startPt)/magD,
                                segmenti-1,
                                segmenti,

                                newTracks.last(),
                                newScalars.last(),
                                newVectors.last()
                            );

                            newTracks.last().shrink();
                            newScalars.last().shrink();
                            newVectors.last().shrink();
                        }
                    }
                }
                else
                {
                    // startPt outside box. New track. Get starting point

                    point clipPt;
                    if (bb.intersects(startPt, endPt, clipPt))
                    {
                        // New track
                        newTracks.append
                        (
                            new DynamicList<point>(track.size()/10)
                        );
                        newScalars.append
                        (
                            new DynamicList<scalarList>(track.size()/10)
                        );
                        newVectors.append
                        (
                            new DynamicList<vectorList>(track.size()/10)
                        );

                        // Store point and interpolated values
                        storePoint
                        (
                            tracki,

                            mag(clipPt-startPt)/magD,
                            segmenti-1,
                            segmenti,

                            newTracks.last(),
                            newScalars.last(),
                            newVectors.last()
                        );

                        if (!bb.contains(endPt))
                        {
                            bb.intersects
                            (
                                endPt,
                                point(clipPt),
                                clipPt
                            );

                            // Store point and interpolated values
                            storePoint
                            (
                                tracki,

                                mag(clipPt-startPt)/magD,
                                segmenti-1,
                                segmenti,

                                newTracks.last(),
                                newScalars.last(),
                                newVectors.last()
                            );

                            newTracks.last().shrink();
                            newScalars.last().shrink();
                            newVectors.last().shrink();
                        }
                    }
                }
            }
        }

        // Last point
        if (bb.contains(track.last()))
        {
            storePoint
            (
                tracki,

                1.0,
                track.size()-2,
                track.size()-1,

                newTracks.last(),
                newScalars.last(),
                newVectors.last()
            );
        }
    }
}


void Foam::functionObjects::streamLineBase::trimToBox(const treeBoundBox& bb)
{
    // Storage for new tracks. Per track, per sample the coordinate (newTracks)
    // or values for all the sampled fields (newScalars, newVectors)
    PtrList<DynamicList<point>> newTracks;
    PtrList<DynamicList<scalarList>> newScalars;
    PtrList<DynamicList<vectorList>> newVectors;

    forAll(allTracks_, tracki)
    {
        const List<point>& track = allTracks_[tracki];

        if (track.size())
        {
            // New track. Assume it consists of the whole track
            newTracks.append(new DynamicList<point>(track.size()));
            newScalars.append(new DynamicList<scalarList>(track.size()));
            newVectors.append(new DynamicList<vectorList>(track.size()));

            // Trim, split and append to newTracks
            trimToBox(bb, tracki, newTracks, newScalars, newVectors);
        }
    }

    // Transfer newTracks to allTracks_
    allTracks_.setSize(newTracks.size());
    forAll(allTracks_, tracki)
    {
        allTracks_[tracki].transfer(newTracks[tracki]);
    }
    // Replace track scalars
    forAll(allScalars_, scalari)
    {
        DynamicList<scalarList>& fieldVals = allScalars_[scalari];
        fieldVals.setSize(newTracks.size());

        forAll(fieldVals, tracki)
        {
            scalarList& trackVals = allScalars_[scalari][tracki];
            trackVals.setSize(newScalars[tracki].size());
            forAll(trackVals, samplei)
            {
                trackVals[samplei] = newScalars[tracki][samplei][scalari];
            }
        }
    }
    // Replace track vectors
    forAll(allVectors_, vectori)
    {
        DynamicList<vectorList>& fieldVals = allVectors_[vectori];
        fieldVals.setSize(newTracks.size());
        forAll(fieldVals, tracki)
        {
            vectorList& trackVals = allVectors_[vectori][tracki];
            trackVals.setSize(newVectors[tracki].size());
            forAll(trackVals, samplei)
            {
                trackVals[samplei] = newVectors[tracki][samplei][vectori];
            }
        }
    }
}


bool Foam::functionObjects::streamLineBase::writeToFile()
{
    if (Pstream::parRun())
    {
        // Append slave tracks to master ones
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        globalIndex globalTrackIDs(allTracks_.size());

        // Construct a distribution map to pull all to the master.
        labelListList sendMap(Pstream::nProcs());
        labelListList recvMap(Pstream::nProcs());

        if (Pstream::master())
        {
            // Master: receive all. My own first, then consecutive
            // processors.
            label tracki = 0;

            forAll(recvMap, proci)
            {
                labelList& fromProc = recvMap[proci];
                fromProc.setSize(globalTrackIDs.localSize(proci));
                forAll(fromProc, i)
                {
                    fromProc[i] = tracki++;
                }
            }
        }

        labelList& toMaster = sendMap[0];
        toMaster.setSize(globalTrackIDs.localSize());
        forAll(toMaster, i)
        {
            toMaster[i] = i;
        }

        const mapDistribute distMap
        (
            globalTrackIDs.size(),
            std::move(sendMap),
            std::move(recvMap)
        );


        // Distribute the track positions. Note: use scheduled comms
        // to prevent buffering.
        allTracks_.shrink();
        mapDistributeBase::distribute
        (
            Pstream::commsTypes::scheduled,
            distMap.schedule(),
            distMap.constructSize(),
            distMap.subMap(),
            false,
            distMap.constructMap(),
            false,
            allTracks_,
            flipOp()
        );
        allTracks_.setCapacity(allTracks_.size());

        // Distribute the scalars
        forAll(allScalars_, scalari)
        {
            allScalars_[scalari].shrink();
            mapDistributeBase::distribute
            (
                Pstream::commsTypes::scheduled,
                distMap.schedule(),
                distMap.constructSize(),
                distMap.subMap(),
                false,
                distMap.constructMap(),
                false,
                allScalars_[scalari],
                flipOp()
            );
            allScalars_[scalari].setCapacity(allScalars_[scalari].size());
        }
        // Distribute the vectors
        forAll(allVectors_, vectori)
        {
            allVectors_[vectori].shrink();
            mapDistributeBase::distribute
            (
                Pstream::commsTypes::scheduled,
                distMap.schedule(),
                distMap.constructSize(),
                distMap.subMap(),
                false,
                distMap.constructMap(),
                false,
                allVectors_[vectori],
                flipOp()
            );
            allVectors_[vectori].setCapacity(allVectors_[vectori].size());
        }
    }


    // Note: filenames scattered below since used in global call
    fileName scalarVtkFile;
    fileName vectorVtkFile;

    if (Pstream::master())
    {
        if (!bounds_.empty())
        {
            // Clip to bounding box
            trimToBox(treeBoundBox(bounds_));
        }


        label nTracks = 0;
        label n = 0;
        forAll(allTracks_, tracki)
        {
            if (allTracks_[tracki].size())
            {
                nTracks++;
                n += allTracks_[tracki].size();
            }
        }

        Log << "    Tracks:" << nTracks << nl
            << "    Total samples:" << n
            << endl;


        // Massage into form suitable for writers
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // Make output directory

        fileName vtkPath
        (
            time_.globalPath()/functionObject::outputPrefix/"sets"/name()
        );
        if (mesh_.name() != polyMesh::defaultRegion)
        {
            vtkPath = vtkPath/mesh_.name();
        }
        vtkPath = vtkPath/mesh_.time().timeName();

        mkDir(vtkPath);

        // Convert track positions (and compact out empty tracks)

        PtrList<coordSet> tracks(nTracks);
        nTracks = 0;
        labelList oldToNewTrack(allTracks_.size(), -1);

        forAll(allTracks_, tracki)
        {
            if (allTracks_[tracki].size())
            {
                List<point>& points = allTracks_[tracki];
                scalarList dist(points.size());
                dist[0] = 0;
                for (label pointi = 1; pointi < points.size(); ++pointi)
                {
                    dist[pointi] =
                        dist[pointi-1] + mag(points[pointi] - points[pointi-1]);
                }

                tracks.set
                (
                    nTracks,
                    new coordSet
                    (
                        "track" + Foam::name(nTracks),
                        sampledSetAxis(),  // "xyz"
                        std::move(allTracks_[tracki]),
                        std::move(dist)
                    )
                );
                oldToNewTrack[tracki] = nTracks;
                ++nTracks;
            }
        }

        // Convert scalar values

        if (!allScalars_.empty() && !tracks.empty())
        {
            List<List<scalarField>> scalarValues(allScalars_.size());

            forAll(allScalars_, scalari)
            {
                DynamicList<scalarList>& allTrackVals = allScalars_[scalari];
                scalarValues[scalari].setSize(nTracks);

                forAll(allTrackVals, tracki)
                {
                    scalarList& vals = allTrackVals[tracki];
                    if (vals.size())
                    {
                        const label newTracki = oldToNewTrack[tracki];
                        scalarValues[scalari][newTracki].transfer(vals);
                    }
                }
            }

            scalarVtkFile = fileName
            (
                vtkPath
              / scalarFormatterPtr_().getFileName
                (
                    tracks[0],
                    scalarNames_
                )
            );

            Log << "    Writing data to " << scalarVtkFile.path() << endl;

            scalarFormatterPtr_().write
            (
                true,                   // writeTracks
                List<scalarField>(),    // times
                tracks,
                scalarNames_,
                scalarValues,
                OFstream(scalarVtkFile)()
            );
        }

        // Convert vector values

        if (!allVectors_.empty() && !tracks.empty())
        {
            List<List<vectorField>> vectorValues(allVectors_.size());

            forAll(allVectors_, vectori)
            {
                DynamicList<vectorList>& allTrackVals = allVectors_[vectori];
                vectorValues[vectori].setSize(nTracks);

                forAll(allTrackVals, tracki)
                {
                    vectorList& vals = allTrackVals[tracki];
                    if (vals.size())
                    {
                        const label newTracki = oldToNewTrack[tracki];
                        vectorValues[vectori][newTracki].transfer(vals);
                    }
                }
            }

            vectorVtkFile = fileName
            (
                vtkPath
              / vectorFormatterPtr_().getFileName(tracks[0], vectorNames_)
            );

            //Info<< "    Writing vector data to " << vectorVtkFile << endl;

            vectorFormatterPtr_().write
            (
                true,                   // writeTracks
                List<scalarField>(),    // times
                tracks,
                vectorNames_,
                vectorValues,
                OFstream(vectorVtkFile)()
            );
        }
    }


    // File names are generated on the master but setProperty needs to
    // be across all procs
    Pstream::scatter(scalarVtkFile);
    for (const word& fieldName : scalarNames_)
    {
        dictionary propsDict;
        propsDict.add
        (
            "file",
            time_.relativePath(scalarVtkFile, true)
        );
        setProperty(fieldName, propsDict);
    }

    Pstream::scatter(vectorVtkFile);
    for (const word& fieldName : vectorNames_)
    {
        dictionary propsDict;
        propsDict.add
        (
            "file",
            time_.relativePath(vectorVtkFile, true)
        );
        setProperty(fieldName, propsDict);
    }

    return true;
}


void Foam::functionObjects::streamLineBase::resetFieldNames
(
    const word& newUName,
    const wordList& newFieldNames
)
{
    UName_ = newUName;
    fields_ = newFieldNames;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::streamLineBase::streamLineBase
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    dict_(dict),
    fields_()
{}


Foam::functionObjects::streamLineBase::streamLineBase
(
    const word& name,
    const Time& runTime,
    const dictionary& dict,
    const wordList& fieldNames
)
:
    fvMeshFunctionObject(name, runTime, dict),
    dict_(dict),
    fields_(fieldNames)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::streamLineBase::~streamLineBase()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::streamLineBase::read(const dictionary& dict)
{
    if (&dict_ != &dict)
    {
        // Update local copy of dictionary:
        dict_ = dict;
    }

    fvMeshFunctionObject::read(dict);

    Info<< type() << " " << name() << ":" << nl;

    UName_ = dict.getOrDefault<word>("U", "U");

    if (fields_.empty())
    {
        dict.readEntry("fields", fields_);

        if (!fields_.found(UName_))
        {
            FatalIOErrorInFunction(dict)
                << "Velocity field for tracking " << UName_
                << " should be present in the list of fields " << fields_
                << exit(FatalIOError);
        }
    }

    Info<< "    Employing velocity field " << UName_ << endl;

    bool trackForward;
    if (dict.readIfPresent("trackForward", trackForward))
    {
        trackDir_ =
        (
            trackForward
          ? trackDirType::FORWARD
          : trackDirType::BACKWARD
        );

        if (dict.found("direction"))
        {
            FatalIOErrorInFunction(dict)
                << "Cannot specify both 'trackForward' and 'direction'" << nl
                << exit(FatalIOError);
        }
    }
    else
    {
        trackDir_ = trackDirTypeNames.get("direction", dict);
    }
    dict.readEntry("lifeTime", lifeTime_);
    if (lifeTime_ < 1)
    {
        FatalIOErrorInFunction(dict)
            << "Illegal value " << lifeTime_ << " for lifeTime"
            << exit(FatalIOError);
    }


    trackLength_ = VGREAT;
    if (dict.readIfPresent("trackLength", trackLength_))
    {
        Info<< type() << " : fixed track length specified : "
            << trackLength_ << nl << endl;
    }


    bounds_ = boundBox::invertedBox;
    if (dict.readIfPresent("bounds", bounds_) && !bounds_.empty())
    {
        Info<< "    clipping all segments to " << bounds_ << nl << endl;
    }


    interpolationScheme_ = dict.getOrDefault
    (
        "interpolationScheme",
        interpolationCellPoint<scalar>::typeName
    );

    //Info<< "    using interpolation " << interpolationScheme_ << endl;

    cloudName_ = dict.getOrDefault<word>("cloud", type());

    sampledSetPtr_.clear();
    sampledSetAxis_.clear();

    scalarFormatterPtr_ = writer<scalar>::New(dict.get<word>("setFormat"));
    vectorFormatterPtr_ = writer<vector>::New(dict.get<word>("setFormat"));

    return true;
}


bool Foam::functionObjects::streamLineBase::execute()
{
    return true;
}


bool Foam::functionObjects::streamLineBase::write()
{
    Log << type() << " " << name() << " write:" << nl;

    // Do all injection and tracking
    track();

    writeToFile();

    return true;
}


void Foam::functionObjects::streamLineBase::updateMesh(const mapPolyMesh& mpm)
{
    if (&mpm.mesh() == &mesh_)
    {
        read(dict_);
    }
}


void Foam::functionObjects::streamLineBase::movePoints(const polyMesh& mpm)
{
    // Moving mesh affects the search tree
    read(dict_);
}


// ************************************************************************* //

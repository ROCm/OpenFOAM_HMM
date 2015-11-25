/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
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

#include "streamLineBase.H"
#include "fvMesh.H"
#include "ReadFields.H"
#include "sampledSet.H"
#include "globalIndex.H"
#include "mapDistribute.H"
#include "interpolationCellPoint.H"
#include "wallPolyPatch.H"
#include "meshSearchMeshObject.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(streamLineBase, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::autoPtr<Foam::indirectPrimitivePatch>
Foam::streamLineBase::wallPatch() const
{
    const fvMesh& mesh = dynamic_cast<const fvMesh&>(obr_);

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    label nFaces = 0;

    forAll(patches, patchI)
    {
        //if (!polyPatch::constraintType(patches[patchI].type()))
        if (isA<wallPolyPatch>(patches[patchI]))
        {
            nFaces += patches[patchI].size();
        }
    }

    labelList addressing(nFaces);

    nFaces = 0;

    forAll(patches, patchI)
    {
        //if (!polyPatch::constraintType(patches[patchI].type()))
        if (isA<wallPolyPatch>(patches[patchI]))
        {
            const polyPatch& pp = patches[patchI];

            forAll(pp, i)
            {
                addressing[nFaces++] = pp.start()+i;
            }
        }
    }

    return autoPtr<indirectPrimitivePatch>
    (
        new indirectPrimitivePatch
        (
            IndirectList<face>
            (
                mesh.faces(),
                addressing
            ),
            mesh.points()
        )
    );
}


void Foam::streamLineBase::initInterpolations
(
    const label nSeeds,
    label& UIndex,
    PtrList<volScalarField>& vsFlds,
    PtrList<interpolation<scalar> >& vsInterp,
    PtrList<volVectorField>& vvFlds,
    PtrList<interpolation<vector> >& vvInterp
)
{
    const Time& runTime = obr_.time();
    const fvMesh& mesh = dynamic_cast<const fvMesh&>(obr_);

    // Read or lookup fields

    if (loadFromFiles_)
    {
        IOobjectList allObjects(mesh, runTime.timeName());

        IOobjectList objects(2*fields_.size());
        forAll(fields_, i)
        {
            objects.add(*allObjects[fields_[i]]);
        }

        ReadFields(mesh, objects, vsFlds);
        vsInterp.setSize(vsFlds.size());
        forAll(vsFlds, i)
        {
            vsInterp.set
            (
                i,
                interpolation<scalar>::New
                (
                    interpolationScheme_,
                    vsFlds[i]
                )
            );
        }
        ReadFields(mesh, objects, vvFlds);
        vvInterp.setSize(vvFlds.size());
        forAll(vvFlds, i)
        {
            vvInterp.set
            (
                i,
                interpolation<vector>::New
                (
                    interpolationScheme_,
                    vvFlds[i]
                )
            );
        }
    }
    else
    {
        label nScalar = 0;
        label nVector = 0;

        forAll(fields_, i)
        {
            if (mesh.foundObject<volScalarField>(fields_[i]))
            {
                nScalar++;
            }
            else if (mesh.foundObject<volVectorField>(fields_[i]))
            {
                nVector++;
            }
            else
            {
                FatalErrorIn("streamLineBase::track()")
                    << "Cannot find field " << fields_[i] << nl
                    << "Valid scalar fields are:"
                    << mesh.names(volScalarField::typeName) << nl
                    << "Valid vector fields are:"
                    << mesh.names(volVectorField::typeName)
                    << exit(FatalError);
            }
        }
        vsInterp.setSize(nScalar);
        nScalar = 0;
        vvInterp.setSize(nVector);
        nVector = 0;

        forAll(fields_, i)
        {
            if (mesh.foundObject<volScalarField>(fields_[i]))
            {
                const volScalarField& f = mesh.lookupObject<volScalarField>
                (
                    fields_[i]
                );
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
            else if (mesh.foundObject<volVectorField>(fields_[i]))
            {
                const volVectorField& f = mesh.lookupObject<volVectorField>
                (
                    fields_[i]
                );

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
        FatalErrorIn("streamLineBase::track()")
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


void Foam::streamLineBase::storePoint
(
    const label trackI,

    const scalar w,
    const label leftI,
    const label rightI,

    DynamicList<point>& newTrack,
    DynamicList<scalarList>& newScalars,
    DynamicList<vectorList>& newVectors
) const
{
    label sz = newTrack.size();

    const List<point>& track = allTracks_[trackI];

    newTrack.append((1.0-w)*track[leftI] + w*track[rightI]);

    // Scalars
    {
        newScalars.append(scalarList(allScalars_.size()));
        scalarList& newVals = newScalars[sz];

        forAll(allScalars_, scalarI)
        {
            const scalarList& trackVals = allScalars_[scalarI][trackI];
            newVals[scalarI] = (1.0-w)*trackVals[leftI] + w*trackVals[rightI];
        }
    }

    // Vectors
    {
        newVectors.append(vectorList(allVectors_.size()));
        vectorList& newVals = newVectors[sz];

        forAll(allVectors_, vectorI)
        {
            const vectorList& trackVals = allVectors_[vectorI][trackI];
            newVals[vectorI] = (1.0-w)*trackVals[leftI] + w*trackVals[rightI];
        }
    }
}


// Can split a track into multiple tracks
void Foam::streamLineBase::trimToBox
(
    const treeBoundBox& bb,
    const label trackI,
    PtrList<DynamicList<point> >& newTracks,
    PtrList<DynamicList<scalarList> >& newScalars,
    PtrList<DynamicList<vectorList> >& newVectors
) const
{
    const List<point>& track = allTracks_[trackI];
    if (track.size())
    {
        for
        (
            label segmentI = 1;
            segmentI < track.size();
            segmentI++
        )
        {
            const point& startPt = track[segmentI-1];
            const point& endPt = track[segmentI];

            const vector d(endPt-startPt);
            scalar magD = mag(d);
            if (magD > ROOTVSMALL)
            {
                if (bb.contains(startPt))
                {
                    // Store 1.0*track[segmentI-1]+0*track[segmentI]
                    storePoint
                    (
                        trackI,

                        0.0,
                        segmentI-1,
                        segmentI,

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
                                trackI,

                                mag(clipPt-startPt)/magD,
                                segmentI-1,
                                segmentI,

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
                            trackI,

                            mag(clipPt-startPt)/magD,
                            segmentI-1,
                            segmentI,

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
                                trackI,

                                mag(clipPt-startPt)/magD,
                                segmentI-1,
                                segmentI,

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
                trackI,

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


void Foam::streamLineBase::trimToBox(const treeBoundBox& bb)
{
    // Storage for new tracks. Per track, per sample the coordinate (newTracks)
    // or values for all the sampled fields (newScalars, newVectors)
    PtrList<DynamicList<point> > newTracks;
    PtrList<DynamicList<scalarList> > newScalars;
    PtrList<DynamicList<vectorList> > newVectors;

    forAll(allTracks_, trackI)
    {
        const List<point>& track = allTracks_[trackI];

        if (track.size())
        {
            // New track. Assume it consists of the whole track
            newTracks.append(new DynamicList<point>(track.size()));
            newScalars.append(new DynamicList<scalarList>(track.size()));
            newVectors.append(new DynamicList<vectorList>(track.size()));

            // Trim, split and append to newTracks
            trimToBox(bb, trackI, newTracks, newScalars, newVectors);
        }
    }

    // Transfer newTracks to allTracks_
    allTracks_.setSize(newTracks.size());
    forAll(allTracks_, trackI)
    {
        allTracks_[trackI].transfer(newTracks[trackI]);
    }
    // Replace track scalars
    forAll(allScalars_, scalarI)
    {
        DynamicList<scalarList>& fieldVals = allScalars_[scalarI];
        fieldVals.setSize(newTracks.size());

        forAll(fieldVals, trackI)
        {
            scalarList& trackVals = allScalars_[scalarI][trackI];
            trackVals.setSize(newScalars[trackI].size());
            forAll(trackVals, sampleI)
            {
                trackVals[sampleI] = newScalars[trackI][sampleI][scalarI];
            }
        }
    }
    // Replace track vectors
    forAll(allVectors_, vectorI)
    {
        DynamicList<vectorList>& fieldVals = allVectors_[vectorI];
        fieldVals.setSize(newTracks.size());
        forAll(fieldVals, trackI)
        {
            vectorList& trackVals = allVectors_[vectorI][trackI];
            trackVals.setSize(newVectors[trackI].size());
            forAll(trackVals, sampleI)
            {
                trackVals[sampleI] = newVectors[trackI][sampleI][vectorI];
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::streamLineBase::streamLineBase
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    functionObjectState(obr, name),
    dict_(dict),
    obr_(obr),
    loadFromFiles_(loadFromFiles),
    log_(true)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::streamLineBase::~streamLineBase()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::streamLineBase::read(const dictionary& dict)
{
    if (active_)
    {
        log_.readIfPresent("log", dict);

        if (log_) Info<< type() << " " << name_ << ":" << nl;

        dict.lookup("fields") >> fields_;
        if (dict.found("UName"))
        {
            dict.lookup("UName") >> UName_;
        }
        else
        {
            UName_ = "U";
            if (dict.found("U"))
            {
                IOWarningIn("streamLineBase::read(const dictionary&)", dict)
                    << "Using deprecated entry \"U\"."
                    << " Please use \"UName\" instead."
                    << endl;
                dict.lookup("U") >> UName_;
            }
        }

        if (findIndex(fields_, UName_) == -1)
        {
            FatalIOErrorIn("streamLineBase::read(const dictionary&)", dict)
                << "Velocity field for tracking " << UName_
                << " should be present in the list of fields " << fields_
                << exit(FatalIOError);
        }


        dict.lookup("trackForward") >> trackForward_;
        dict.lookup("lifeTime") >> lifeTime_;
        if (lifeTime_ < 1)
        {
            FatalErrorIn(":streamLineBase::read(const dictionary&)")
                << "Illegal value " << lifeTime_ << " for lifeTime"
                << exit(FatalError);
        }


        trackLength_ = VGREAT;
        if (dict.found("trackLength"))
        {
            dict.lookup("trackLength") >> trackLength_;

            if (log_)
            {
                Info<< type() << " : fixed track length specified : "
                    << trackLength_ << nl << endl;
            }
        }


        bounds_ = boundBox::greatBox;
        if (dict.readIfPresent("bounds", bounds_))
        {
            if (log_) Info<< "    clipping all segments to " << bounds_ << nl << endl;
        }


        interpolationScheme_ = dict.lookupOrDefault
        (
            "interpolationScheme",
            interpolationCellPoint<scalar>::typeName
        );

        //if (log_) Info<< "    using interpolation " << interpolationScheme_
        //    << endl;

        cloudName_ = dict.lookupOrDefault<word>("cloudName", type());
        dict.lookup("seedSampleSet") >> seedSet_;

        const fvMesh& mesh = dynamic_cast<const fvMesh&>(obr_);

        const dictionary& coeffsDict = dict.subDict(seedSet_ + "Coeffs");

        sampledSetPtr_ = sampledSet::New
        (
            seedSet_,
            mesh,
            meshSearchMeshObject::New(mesh),
            coeffsDict
        );
        coeffsDict.lookup("axis") >> sampledSetAxis_;

        scalarFormatterPtr_ = writer<scalar>::New(dict.lookup("setFormat"));
        vectorFormatterPtr_ = writer<vector>::New(dict.lookup("setFormat"));
    }
}


void Foam::streamLineBase::execute()
{}


void Foam::streamLineBase::end()
{}


void Foam::streamLineBase::timeSet()
{}


void Foam::streamLineBase::write()
{
    if (active_)
    {
        if (log_) Info<< type() << " " << name_ << " output:" << nl;

        const Time& runTime = obr_.time();
        const fvMesh& mesh = dynamic_cast<const fvMesh&>(obr_);


        // Do all injection and tracking
        track();


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
                label trackI = 0;

                forAll(recvMap, procI)
                {
                    labelList& fromProc = recvMap[procI];
                    fromProc.setSize(globalTrackIDs.localSize(procI));
                    forAll(fromProc, i)
                    {
                        fromProc[i] = trackI++;
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
                sendMap.xfer(),
                recvMap.xfer()
            );


            // Distribute the track positions. Note: use scheduled comms
            // to prevent buffering.
            allTracks_.shrink();
            mapDistributeBase::distribute
            (
                Pstream::scheduled,
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
            forAll(allScalars_, scalarI)
            {
                allScalars_[scalarI].shrink();
                mapDistributeBase::distribute
                (
                    Pstream::scheduled,
                    distMap.schedule(),
                    distMap.constructSize(),
                    distMap.subMap(),
                    false,
                    distMap.constructMap(),
                    false,
                    allScalars_[scalarI],
                    flipOp()
                );
                allScalars_[scalarI].setCapacity(allScalars_[scalarI].size());
            }
            // Distribute the vectors
            forAll(allVectors_, vectorI)
            {
                allVectors_[vectorI].shrink();
                mapDistributeBase::distribute
                (
                    Pstream::scheduled,
                    distMap.schedule(),
                    distMap.constructSize(),
                    distMap.subMap(),
                    false,
                    distMap.constructMap(),
                    false,
                    allVectors_[vectorI],
                    flipOp()
                );
                allVectors_[vectorI].setCapacity(allVectors_[vectorI].size());
            }
        }


        if (Pstream::master())
        {
            if (bounds_ != boundBox::greatBox)
            {
                // Clip to bounding box
                trimToBox(treeBoundBox(bounds_));
            }


            label nTracks = 0;
            label n = 0;
            forAll(allTracks_, trackI)
            {
                if (allTracks_[trackI].size())
                {
                    nTracks++;
                    n += allTracks_[trackI].size();
                }
            }

            if (log_)
            {
                Info<< "    Tracks:" << nTracks << nl
                    << "    Total samples:" << n
                    << endl;
            }


            // Massage into form suitable for writers
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            // Make output directory

            fileName vtkPath
            (
                Pstream::parRun()
              ? runTime.path()/".."/"postProcessing"/"sets"/name()
              : runTime.path()/"postProcessing"/"sets"/name()
            );
            if (mesh.name() != fvMesh::defaultRegion)
            {
                vtkPath = vtkPath/mesh.name();
            }
            vtkPath = vtkPath/mesh.time().timeName();

            mkDir(vtkPath);

            // Convert track positions (and compact out empty tracks)

            PtrList<coordSet> tracks(nTracks);
            nTracks = 0;
            labelList oldToNewTrack(allTracks_.size(), -1);

            forAll(allTracks_, trackI)
            {
                if (allTracks_[trackI].size())
                {
                    tracks.set
                    (
                        nTracks,
                        new coordSet
                        (
                            "track" + Foam::name(nTracks),
                            sampledSetAxis_                 //"xyz"
                        )
                    );
                    oldToNewTrack[trackI] = nTracks;
                    tracks[nTracks].transfer(allTracks_[trackI]);
                    nTracks++;
                }
            }

            // Convert scalar values

            if (allScalars_.size() > 0)
            {
                List<List<scalarField> > scalarValues(allScalars_.size());

                forAll(allScalars_, scalarI)
                {
                    DynamicList<scalarList>& allTrackVals =
                        allScalars_[scalarI];
                    scalarValues[scalarI].setSize(nTracks);

                    forAll(allTrackVals, trackI)
                    {
                        scalarList& vals = allTrackVals[trackI];
                        if (vals.size())
                        {
                            label newTrackI = oldToNewTrack[trackI];
                            scalarValues[scalarI][newTrackI].transfer(vals);
                        }
                    }
                }

                fileName vtkFile
                (
                    vtkPath
                  / scalarFormatterPtr_().getFileName
                    (
                        tracks[0],
                        scalarNames_
                    )
                );

                Info(log_)<< "    Writing data to " << vtkFile.path() << endl;

                scalarFormatterPtr_().write
                (
                    true,           // writeTracks
                    tracks,
                    scalarNames_,
                    scalarValues,
                    OFstream(vtkFile)()
                );

                forAll(scalarNames_, nameI)
                {
                    dictionary propsDict;
                    propsDict.add("file", vtkFile);
                    const word& fieldName = scalarNames_[nameI];
                    setProperty(fieldName, propsDict);
                }
            }

            // Convert vector values

            if (allVectors_.size() > 0)
            {
                List<List<vectorField> > vectorValues(allVectors_.size());

                forAll(allVectors_, vectorI)
                {
                    DynamicList<vectorList>& allTrackVals =
                        allVectors_[vectorI];
                    vectorValues[vectorI].setSize(nTracks);

                    forAll(allTrackVals, trackI)
                    {
                        vectorList& vals = allTrackVals[trackI];
                        if (vals.size())
                        {
                            label newTrackI = oldToNewTrack[trackI];
                            vectorValues[vectorI][newTrackI].transfer(vals);
                        }
                    }
                }

                fileName vtkFile
                (
                    vtkPath
                  / vectorFormatterPtr_().getFileName
                    (
                        tracks[0],
                        vectorNames_
                    )
                );

                //if (log_) Info<< "    Writing vector data to " << vtkFile << endl;

                vectorFormatterPtr_().write
                (
                    true,           // writeTracks
                    tracks,
                    vectorNames_,
                    vectorValues,
                    OFstream(vtkFile)()
                );

                forAll(vectorNames_, nameI)
                {
                    dictionary propsDict;
                    propsDict.add("file", vtkFile);
                    const word& fieldName = vectorNames_[nameI];
                    setProperty(fieldName, propsDict);
                }
            }
        }
    }
}


void Foam::streamLineBase::updateMesh(const mapPolyMesh&)
{
    read(dict_);
}


void Foam::streamLineBase::movePoints(const polyMesh&)
{
    // Moving mesh affects the search tree
    read(dict_);
}


// ************************************************************************* //

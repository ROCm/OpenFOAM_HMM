/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "Pstream.H"
#include "functionObjectList.H"
#include "streamLine.H"
#include "fvMesh.H"
#include "streamLineParticleCloud.H"
#include "ReadFields.H"
#include "meshSearch.H"
#include "sampledSet.H"
#include "globalIndex.H"
#include "mapDistribute.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::streamLine, 0);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::streamLine::track()
{
    const Time& runTime = const_cast<Time&>(obr_.time());
    const fvMesh& mesh = dynamic_cast<const fvMesh&>(obr_);

    IDLList<streamLineParticle> initialParticles;
    streamLineParticleCloud particles
    (
        mesh,
        cloudName_,
        initialParticles
    );

    //Pout<< "Seeding particles." << endl;

    const sampledSet& seedPoints = sampledSetPtr_();
    forAll(seedPoints, i)
    {
        //Pout<< "Seeded particle at " << seedPoints[i]
        //    << " at cell:" << seedPoints.cells()[i]
        //    << endl;

        particles.addParticle
        (
            new streamLineParticle
            (
                particles,
                seedPoints[i],
                seedPoints.cells()[i],
                lifeTime_               // lifetime
            )
        );
    }
    label nSeeds = returnReduce(particles.size(), sumOp<label>());
    Info<< "Seeded " << nSeeds
        << " particles." << endl;


    // Read or lookup fields
    PtrList<volScalarField> vsFlds;
    PtrList<interpolationCellPoint<scalar> > vsInterp;
    PtrList<volVectorField> vvFlds;
    PtrList<interpolationCellPoint<vector> > vvInterp;


    label UIndex = -1;

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
            vsInterp.set(i, new interpolationCellPoint<scalar>(vsFlds[i]));
        }
        ReadFields(mesh, objects, vvFlds);
        vvInterp.setSize(vvFlds.size());
        forAll(vvFlds, i)
        {
            vvInterp.set(i, new interpolationCellPoint<vector>(vvFlds[i]));
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
                FatalErrorIn("streamLine::execute()")
                    << "Cannot find field " << fields_[i] << endl
                    << "Valid scalar fields are:"
                    << mesh.names(volScalarField::typeName) << endl
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
                    new interpolationCellPoint<scalar>(f)
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
                    new interpolationCellPoint<vector>(f)
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
        FatalErrorIn("streamLine::execute()")
            << "Cannot find field to move particles with : " << UName_
            << endl
            << "This field has to be present in the sampled fields "
            << fields_
            << " and in the objectRegistry." << endl
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

    // additional particle info
    streamLineParticle::trackData td
    (
        particles,
        vsInterp,
        vvInterp,
        UIndex,         // index of U in vvInterp
        trackForward_,  // track in +u direction
        allTracks_,
        allScalars_,
        allVectors_
    );

    // Track
    //Pout<< "Tracking particles." << endl;
    particles.move(td);
    //Pout<< "Finished tracking particles." << endl;
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
    dict_(dict),
    name_(name),
    obr_(obr),
    loadFromFiles_(loadFromFiles),
    active_(true)
{
    // Only active if a fvMesh is available
    if (isA<fvMesh>(obr_))
    {
        read(dict_);
    }
    else
    {
        active_ = false;
        WarningIn
        (
            "streamLine::streamLine\n"
            "(\n"
                "const word&,\n"
                "const objectRegistry&,\n"
                "const dictionary&,\n"
                "const bool\n"
            ")"
        )   << "No fvMesh available, deactivating."
            << nl << endl;
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
        //dict_ = dict;
        dict.lookup("fields") >> fields_;
        UName_ = dict.lookupOrDefault<word>("U", "U");
        dict.lookup("trackForward") >> trackForward_;
        dict.lookup("lifeTime") >> lifeTime_;
        cloudName_ = dict.lookupOrDefault<word>("cloudName", "streamLine");
        dict.lookup("seedSampleSet") >> seedSet_;

        const fvMesh& mesh = dynamic_cast<const fvMesh&>(obr_);

        meshSearchPtr_.reset(new meshSearch(mesh, false));

        const dictionary& coeffsDict = dict.subDict(seedSet_ + "Coeffs");
        sampledSetPtr_ = sampledSet::New
        (
            seedSet_,
            mesh,
            meshSearchPtr_(),
            coeffsDict
        );
        coeffsDict.lookup("axis") >> sampledSetAxis_;

        scalarFormatterPtr_ = writer<scalar>::New(dict.lookup("setFormat"));
        vectorFormatterPtr_ = writer<vector>::New(dict.lookup("setFormat"));
    }
}


void Foam::streamLine::execute()
{
//    const Time& runTime = const_cast<Time&>(obr_.time());
//    Pout<< "**streamLine::execute : time:" << runTime.timeName() << endl;
//    Pout<< "**streamLine::execute : time:" << runTime.timeIndex() << endl;
//
//    bool isOutputTime = false;
//
//    const functionObjectList& fobs = runTime.functionObjects();
//
//    forAll(fobs, i)
//    {
//        if (isA<streamLineFunctionObject>(fobs[i]))
//        {
//            const streamLineFunctionObject& fo =
//                dynamic_cast<const streamLineFunctionObject&>(fobs[i]);
//
//            if (fo.name() == name_)
//            {
//                Pout<< "found me:" << i << endl;
//                if (fo.outputControl().output())
//                {
//                    isOutputTime = true;
//                    break;
//                }
//            }
//        }
//    }
//
//
//    if (active_ && isOutputTime)
//    {
//        track();
//    }
}


void Foam::streamLine::end()
{
}


void Foam::streamLine::write()
{
    if (active_)
    {
        const Time& runTime = const_cast<Time&>(obr_.time());
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


            // Distribute the track positions
            mapDistribute::distribute
            (
                Pstream::scheduled,
                distMap.schedule(),
                distMap.constructSize(),
                distMap.subMap(),
                distMap.constructMap(),
                allTracks_
            );

            // Distribute the scalars
            forAll(allScalars_, scalarI)
            {
                mapDistribute::distribute
                (
                    Pstream::scheduled,
                    distMap.schedule(),
                    distMap.constructSize(),
                    distMap.subMap(),
                    distMap.constructMap(),
                    allScalars_[scalarI]
                );
            }
            // Distribute the vectors
            forAll(allVectors_, vectorI)
            {
                mapDistribute::distribute
                (
                    Pstream::scheduled,
                    distMap.schedule(),
                    distMap.constructSize(),
                    distMap.subMap(),
                    distMap.constructMap(),
                    allVectors_[vectorI]
                );
            }
        }


        label n = 0;
        forAll(allTracks_, trackI)
        {
            n += allTracks_[trackI].size();
        }

        Info<< "Tracks:" << allTracks_.size()
            << "  total samples:" << n << endl;


        // Massage into form suitable for writers
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        if (Pstream::master() && allTracks_.size())
        {
            // Make output directory

            fileName vtkPath
            (
                Pstream::parRun()
              ? runTime.path()/".."/"sets"/name()
              : runTime.path()/"sets"/name()
            );
            if (mesh.name() != fvMesh::defaultRegion)
            {
                vtkPath = vtkPath/mesh.name();
            }
            vtkPath = vtkPath/mesh.time().timeName();

            mkDir(vtkPath);


            // Convert track positions

            PtrList<coordSet> tracks(allTracks_.size());
            forAll(allTracks_, trackI)
            {
                tracks.set
                (
                    trackI,
                    new coordSet
                    (
                        "track" + Foam::name(trackI),
                        sampledSetAxis_                 //"xyz"
                    )
                );
                tracks[trackI].transfer(allTracks_[trackI]);
            }



            // Convert scalar values

            if (allScalars_.size() > 0)
            {
                List<List<scalarField> > scalarValues(allScalars_.size());

                forAll(allScalars_, scalarI)
                {
                    DynamicList<scalarList>& allTrackVals =
                        allScalars_[scalarI];
                    scalarValues[scalarI].setSize(allTrackVals.size());

                    forAll(allTrackVals, trackI)
                    {
                        scalarList& trackVals = allTrackVals[trackI];
                        scalarValues[scalarI][trackI].transfer(trackVals);
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

                Info<< "Writing data to " << vtkFile.path() << endl;

                scalarFormatterPtr_().write
                (
                    true,           // writeTracks
                    tracks,
                    scalarNames_,
                    scalarValues,
                    OFstream(vtkFile)()
                );
            }

            // Convert vector values

            if (allVectors_.size() > 0)
            {
                List<List<vectorField> > vectorValues(allVectors_.size());

                forAll(allVectors_, vectorI)
                {
                    DynamicList<vectorList>& allTrackVals =
                        allVectors_[vectorI];
                    vectorValues[vectorI].setSize(allTrackVals.size());

                    forAll(allTrackVals, trackI)
                    {
                        vectorList& trackVals = allTrackVals[trackI];
                        vectorValues[vectorI][trackI].transfer(trackVals);
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

                //Info<< "Writing vector data to " << vtkFile << endl;

                vectorFormatterPtr_().write
                (
                    true,           // writeTracks
                    tracks,
                    vectorNames_,
                    vectorValues,
                    OFstream(vtkFile)()
                );
            }
        }
    }
}


void Foam::streamLine::updateMesh(const mapPolyMesh&)
{
    read(dict_);
}


void Foam::streamLine::movePoints(const pointField&)
{
    // Moving mesh affects the search tree
    read(dict_);
}


//void Foam::streamLine::readUpdate(const polyMesh::readUpdateState state)
//{
//    if (state != UNCHANGED)
//    {
//        read(dict_);
//    }
//}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "wallBoundedStreamLine.H"
#include "fvMesh.H"
#include "wallBoundedStreamLineParticleCloud.H"
#include "sampledSet.H"
#include "faceSet.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(wallBoundedStreamLine, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tetIndices Foam::wallBoundedStreamLine::findNearestTet
(
    const PackedBoolList& isWallPatch,
    const point& seedPt,
    const label cellI
) const
{
    const fvMesh& mesh = dynamic_cast<const fvMesh&>(obr_);

    const cell& cFaces = mesh.cells()[cellI];

    label minFaceI = -1;
    label minTetPtI = -1;
    scalar minDistSqr = sqr(GREAT);

    forAll(cFaces, cFaceI)
    {
        label faceI = cFaces[cFaceI];

        if (isWallPatch[faceI])
        {
            const face& f = mesh.faces()[faceI];
            const label fp0 = mesh.tetBasePtIs()[faceI];
            const point& basePoint = mesh.points()[f[fp0]];

            label fp = f.fcIndex(fp0);
            for (label i = 2; i < f.size(); i++)
            {
                const point& thisPoint = mesh.points()[f[fp]];
                label nextFp = f.fcIndex(fp);
                const point& nextPoint = mesh.points()[f[nextFp]];

                const triPointRef tri(basePoint, thisPoint, nextPoint);

                scalar d2 = magSqr(tri.centre() - seedPt);
                if (d2 < minDistSqr)
                {
                    minDistSqr = d2;
                    minFaceI = faceI;
                    minTetPtI = i-1;
                }
                fp = nextFp;
            }
        }
    }

    // Put particle in tet
    return tetIndices
    (
        cellI,
        minFaceI,
        minTetPtI,
        mesh
    );
}


void Foam::wallBoundedStreamLine::track()
{
    //const Time& runTime = obr_.time();
    const fvMesh& mesh = dynamic_cast<const fvMesh&>(obr_);


    // Determine the 'wall' patches
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // These are the faces that need to be followed

    autoPtr<indirectPrimitivePatch> boundaryPatch(wallPatch());
    PackedBoolList isWallPatch(mesh.nFaces());
    forAll(boundaryPatch().addressing(), i)
    {
        isWallPatch[boundaryPatch().addressing()[i]] = 1;
    }



    // Find nearest wall particle for the seedPoints
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IDLList<wallBoundedStreamLineParticle> initialParticles;
    wallBoundedStreamLineParticleCloud particles
    (
        mesh,
        cloudName_,
        initialParticles
    );

    {
        // Get the seed points
        // ~~~~~~~~~~~~~~~~~~~

        const sampledSet& seedPoints = sampledSetPtr_();


        forAll(seedPoints, i)
        {
            const point& seedPt = seedPoints[i];
            label cellI = seedPoints.cells()[i];

            tetIndices ids(findNearestTet(isWallPatch, seedPt, cellI));

            if (ids.face() != -1 && isWallPatch[ids.face()])
            {
                //Pout<< "Seeding particle :" << nl
                //    << "     seedPt:" << seedPt << nl
                //    << "     face  :" << ids.face() << nl
                //    << "     at    :" << mesh.faceCentres()[ids.face()] << nl
                //    << "     cell  :" << mesh.cellCentres()[ids.cell()] << nl
                //    << endl;

                particles.addParticle
                (
                    new wallBoundedStreamLineParticle
                    (
                        mesh,
                        ids.faceTri(mesh).centre(),
                        ids.cell(),
                        ids.face(),     // tetFace
                        ids.tetPt(),
                        -1,             // not on a mesh edge
                        -1,             // not on a diagonal edge
                        lifeTime_       // lifetime
                    )
                );
            }
            else
            {
                Pout<< type() << " : ignoring seed " << seedPt
                    << " since not in wall cell" << endl;
            }
        }
    }

    label nSeeds = returnReduce(particles.size(), sumOp<label>());

    if (log_) Info<< type() << " : seeded " << nSeeds << " particles" << endl;



    // Read or lookup fields
    PtrList<volScalarField> vsFlds;
    PtrList<interpolation<scalar> > vsInterp;
    PtrList<volVectorField> vvFlds;
    PtrList<interpolation<vector> > vvInterp;
    label UIndex = -1;

    initInterpolations
    (
        nSeeds,
        UIndex,
        vsFlds,
        vsInterp,
        vvFlds,
        vvInterp
    );

    // Additional particle info
    wallBoundedStreamLineParticle::trackingData td
    (
        particles,
        vsInterp,
        vvInterp,
        UIndex,         // index of U in vvInterp
        trackForward_,  // track in +u direction?
        trackLength_,   // fixed track length
        isWallPatch,    // which faces are to follow

        allTracks_,
        allScalars_,
        allVectors_
    );


    // Set very large dt. Note: cannot use GREAT since 1/GREAT is SMALL
    // which is a trigger value for the tracking...
    const scalar trackTime = Foam::sqrt(GREAT);

    // Track
    particles.move(td, trackTime);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallBoundedStreamLine::wallBoundedStreamLine
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    streamLineBase(name, obr, dict, loadFromFiles)
{
    // Check if the available mesh is an fvMesh otherise deactivate
    if (setActive<fvMesh>())
    {
        read(dict_);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::wallBoundedStreamLine::~wallBoundedStreamLine()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::wallBoundedStreamLine::read(const dictionary& dict)
{
    if (active_)
    {
        streamLineBase::read(dict);

        // Make sure that the mesh is trackable
        if (debug)
        {
            const fvMesh& mesh = dynamic_cast<const fvMesh&>(obr_);

            // 1. positive volume decomposition tets
            faceSet faces(mesh, "lowQualityTetFaces", mesh.nFaces()/100+1);
            if
            (
                polyMeshTetDecomposition::checkFaceTets
                (
                    mesh,
                    polyMeshTetDecomposition::minTetQuality,
                    true,
                    &faces
                )
            )
            {
                label nFaces = returnReduce(faces.size(), sumOp<label>());

                WarningInFunction
                    << "Found " << nFaces
                    <<" faces with low quality or negative volume "
                    << "decomposition tets. Writing to faceSet " << faces.name()
                    << endl;
            }

            // 2. all edges on a cell having two faces
            EdgeMap<label> numFacesPerEdge;
            forAll(mesh.cells(), cellI)
            {
                const cell& cFaces = mesh.cells()[cellI];

                numFacesPerEdge.clear();

                forAll(cFaces, cFaceI)
                {
                    label faceI = cFaces[cFaceI];
                    const face& f = mesh.faces()[faceI];
                    forAll(f, fp)
                    {
                        const edge e(f[fp], f.nextLabel(fp));
                        EdgeMap<label>::iterator eFnd =
                            numFacesPerEdge.find(e);
                        if (eFnd != numFacesPerEdge.end())
                        {
                            eFnd()++;
                        }
                        else
                        {
                            numFacesPerEdge.insert(e, 1);
                        }
                    }
                }

                forAllConstIter(EdgeMap<label>, numFacesPerEdge, iter)
                {
                    if (iter() != 2)
                    {
                        FatalErrorInFunction
                            << "problem cell:" << cellI
                            << abort(FatalError);
                    }
                }
            }
        }
    }
}


// ************************************************************************* //

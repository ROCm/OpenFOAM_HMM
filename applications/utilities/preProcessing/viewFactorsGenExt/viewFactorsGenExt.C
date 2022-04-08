/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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

Application
    viewFactorsGenExt

Group
    grpPreProcessingUtilities

Description
    This view factors generation application uses a combined approach of
    double area integral (2AI) and double linear integral (2LI). 2AI is used
    when the two surfaces are 'far' apart and 2LI whenre they are 'close'.
    2LI is integrated along edges using Gaussian quadrature.
    The distance between faces is calculating a ratio between averaged areas
    and the distance between face centres.

    The input from viewFactorsDict are:

        tolGaussQuad              0.1;      // GaussQuad  error
        distTol                   8;        // R/Average(rm)
        alpha                     0.22;     // Use for common edges for 2LI


    For debugging purposes, the following entries can be set in viewFactorsDict:

        writeViewFactorMatrix     true;
        writeFacesAgglomeration   false;
        dumpRays		          false;


        writeViewFactorMatrix   writes the sum of the VF on each face.
        writeFacesAgglomeration writes the agglomeration
        dumpRays                dumps rays


    In order to specify the participants patches in the VF calculation the
    keywaord viewFactorWall should be added to the boundary file.

    floor
    {
        type            wall;
        inGroups        2(wall viewFactorWall);
        nFaces          100;
        startFace       3100;
    }

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "distributedTriSurfaceMesh.H"
#include "meshTools.H"
#include "constants.H"

#include "uindirectPrimitivePatch.H"
#include "DynamicField.H"
#include "unitConversion.H"

#include "scalarMatrices.H"
#include "labelListIOList.H"
#include "scalarListIOList.H"

#include "singleCellFvMesh.H"

#include "IOmapDistribute.H"

#define PBRT_CONSTEXPR constexpr
#define PBRT_THREAD_LOCAL thread_local

#include "pbrt.h"
#include "shape.h"
#include "triangle.h"
#include "geometry.h"
#include "paramset.h"

#include "accelerators/bvh.h"

using namespace Foam;
using namespace Foam::constant;
using namespace Foam::constant::mathematical;
using namespace pbrt;


triSurface triangulate
(
    const polyBoundaryMesh& bMesh,
    const labelHashSet& includePatches,
    const labelListIOList& finalAgglom,
    labelList& triSurfaceToAgglom,
    const globalIndex& globalNumbering,
    const polyBoundaryMesh& coarsePatches
)
{
    const polyMesh& mesh = bMesh.mesh();

    // Storage for surfaceMesh. Size estimate.
    DynamicList<labelledTri> triangles(mesh.nBoundaryFaces());

    label newPatchI = 0;
    label localTriFaceI = 0;

    for (const label patchI : includePatches)
    {
        const polyPatch& patch = bMesh[patchI];
        const pointField& points = patch.points();

        label nTriTotal = 0;

        forAll(patch, patchFaceI)
        {
            const face& f = patch[patchFaceI];

            faceList triFaces(f.nTriangles(points));

            label nTri = 0;

            f.triangles(points, nTri, triFaces);

            forAll(triFaces, triFaceI)
            {
                const face& f = triFaces[triFaceI];

                triangles.append(labelledTri(f[0], f[1], f[2], newPatchI));

                nTriTotal++;

                triSurfaceToAgglom[localTriFaceI++] = globalNumbering.toGlobal
                (
                    Pstream::myProcNo(),
                    finalAgglom[patchI][patchFaceI]
                  + coarsePatches[patchI].start()
                );
            }
        }

        newPatchI++;
    }

    //triSurfaceToAgglom.resize(localTriFaceI-1);

    triangles.shrink();

    triSurface rawSurface(triangles, mesh.points());

    triSurface surface
    (
        rawSurface.localFaces(),
        rawSurface.localPoints()
    );

    // Combine the triSurfaces across all processors
    if (Pstream::parRun())
    {
        List<List<labelledTri>> surfaceProcTris(Pstream::nProcs());
        List<pointField> surfaceProcPoints(Pstream::nProcs());

        surfaceProcTris[Pstream::myProcNo()] = surface;
        surfaceProcPoints[Pstream::myProcNo()] = surface.points();

        Pstream::gatherList(surfaceProcTris);
        Pstream::scatterList(surfaceProcTris);
        Pstream::gatherList(surfaceProcPoints);
        Pstream::scatterList(surfaceProcPoints);

        label nTris = 0;
        forAll(surfaceProcTris, i)
        {
            nTris += surfaceProcTris[i].size();
        }

        List<labelledTri> globalSurfaceTris(nTris);
        label trii = 0;
        label offset = 0;
        forAll(surfaceProcTris, i)
        {
            forAll(surfaceProcTris[i], j)
            {
                globalSurfaceTris[trii] = surfaceProcTris[i][j];
                globalSurfaceTris[trii][0] += offset;
                globalSurfaceTris[trii][1] += offset;
                globalSurfaceTris[trii][2] += offset;
                trii++;
            }
            offset += surfaceProcPoints[i].size();
        }

        label nPoints = 0;
        forAll(surfaceProcPoints, i)
        {
            nPoints += surfaceProcPoints[i].size();
        }

        pointField globalSurfacePoints(nPoints);

        label pointi = 0;
        forAll(surfaceProcPoints, i)
        {
            forAll(surfaceProcPoints[i], j)
            {
                globalSurfacePoints[pointi++] = surfaceProcPoints[i][j];
            }
        }

        surface = triSurface(globalSurfaceTris, globalSurfacePoints);
    }

    // Add patch names to surface
    surface.patches().setSize(newPatchI);

    newPatchI = 0;

    for (const label patchI : includePatches)
    {
        const polyPatch& patch = bMesh[patchI];

        surface.patches()[newPatchI].index() = patchI;
        surface.patches()[newPatchI].name() = patch.name();
        surface.patches()[newPatchI].geometricType() = patch.type();

        newPatchI++;
    }

    return surface;
}


void writeRays
(
    const fileName& fName,
    const pointField& compactCf,
    const pointField& myFc,
    const labelListList& visibleFaceFaces
)
{
    OFstream str(fName);
    label vertI = 0;

    Pout<< "Dumping rays to " << str.name() << endl;

    forAll(myFc, faceI)
    {
        const labelList visFaces = visibleFaceFaces[faceI];
        forAll(visFaces, faceRemote)
        {
            label compactI = visFaces[faceRemote];
            const point& remoteFc = compactCf[compactI];

            meshTools::writeOBJ(str, myFc[faceI]);
            vertI++;
            meshTools::writeOBJ(str, remoteFc);
            vertI++;
            str << "l " << vertI-1 << ' ' << vertI << nl;
        }
    }
    str.flush();
}


scalar calculateViewFactorFij2AI
(
    const vector& i,
    const vector& j,
    const vector& dAi,
    const vector& dAj
)
{
    vector r = i - j;
    scalar rMag = mag(r);

    if (rMag > SMALL)
    {
        scalar dAiMag = mag(dAi);
        scalar dAjMag = mag(dAj);

        vector ni = dAi/dAiMag;
        vector nj = dAj/dAjMag;
        scalar cosThetaJ = mag(nj & r)/rMag;
        scalar cosThetaI = mag(ni & r)/rMag;

        return
        (
            (cosThetaI*cosThetaJ*dAjMag*dAiMag)
           /(sqr(rMag)*constant::mathematical::pi)
        );
    }
    else
    {
        return 0;
    }
}


void insertMatrixElements
(
    const globalIndex& globalNumbering,
    const label fromProcI,
    const labelListList& globalFaceFaces,
    const scalarListList& viewFactors,
    scalarSquareMatrix& matrix
)
{
    forAll(viewFactors, faceI)
    {
        const scalarList& vf = viewFactors[faceI];
        const labelList& globalFaces = globalFaceFaces[faceI];

        label globalI = globalNumbering.toGlobal(fromProcI, faceI);
        forAll(globalFaces, i)
        {
            matrix[globalI][globalFaces[i]] = vf[i];
        }
    }
}

scalar GaussQuad
(

    const scalarList& w,
    const scalarList& p,
    const scalar& magSi,
    const scalar& magSj,
    const vector& di,
    const vector& dj,
    const vector& ci,
    const vector& cj,
    const scalar cosij,
    const scalar alpha,
    label gi
)
{
    scalar dIntFij = 0;
    if (gi == 0)
    {
        vector r(ci - cj);
        if (mag(r) < SMALL)
        {
            r = (alpha*magSi)*di;
        }
        dIntFij = max(cosij*Foam::log(r&r)*magSi*magSj, 0);
    }
    else
    {
        List<vector> pi(w.size());
        forAll (pi, i)
        {
            pi[i] = ci + p[i]*(magSi/2)*di;
        }


        List<vector> pj(w.size());
        forAll (pj, i)
        {
            pj[i] = cj + p[i]*(magSj/2)*dj;
        }

        forAll (w, i)
        {
            forAll (w, j)
            {
                vector r(pi[i] - pj[j]);

                if (mag(r) < SMALL)
                {
                    r = (alpha*magSi)*di;
                     dIntFij +=
                        cosij*w[i]*w[j]*Foam::log(r&r);
                }
                else
                {
                    dIntFij +=
                        cosij*w[i]*w[j]*Foam::log(r&r);
                }

            }
        }

        dIntFij *= (magSi/2) * (magSj/2);

    }
    return dIntFij;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Calculate view factors from face agglomeration array."
        " The finalAgglom generated by faceAgglomerate utility."
    );

    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createNamedMesh.H"

    // Read view factor dictionary
    IOdictionary viewFactorDict
    (
       IOobject
       (
            "viewFactorsDict",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
       )
    );

    const word viewFactorWall("viewFactorWall");

    const bool writeViewFactors =
        viewFactorDict.getOrDefault("writeViewFactorMatrix", false);

    const bool dumpRays =
        viewFactorDict.getOrDefault("dumpRays", false);

    const label debug = viewFactorDict.getOrDefault<label>("debug", 0);

    const scalar tolGaussQuad =
        viewFactorDict.getOrDefault<scalar>("tolGaussQuad", 0.01);

    const scalar distTol =
         viewFactorDict.getOrDefault<scalar>("distTol", 8);

    const scalar alpha =
         viewFactorDict.getOrDefault<scalar>("alpha", 0.21);

    // Read agglomeration map
    labelListIOList finalAgglom
    (
        IOobject
        (
            "finalAgglom",
            mesh.facesInstance(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    // Create the coarse mesh  using agglomeration
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (debug)
    {
        Pout << "\nCreating single cell mesh..." << endl;
    }

    singleCellFvMesh coarseMesh
    (
        IOobject
        (
            "coarse:" + mesh.name(),
            runTime.timeName(),
            runTime,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        finalAgglom
    );

    if (debug)
    {
        Pout << "\nCreated single cell mesh..." << endl;
    }


    // Calculate total number of fine and coarse faces
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    label nCoarseFaces = 0;      //total number of coarse faces
    label nFineFaces = 0;        //total number of fine faces

    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    const polyBoundaryMesh& coarsePatches = coarseMesh.boundaryMesh();

    labelList viewFactorsPatches(patches.indices(viewFactorWall));
    for (const label patchi : viewFactorsPatches)
    {
        nCoarseFaces += coarsePatches[patchi].size();
        nFineFaces += patches[patchi].size();
    }

    // total number of coarse faces
    label totalNCoarseFaces = nCoarseFaces;

    reduce(totalNCoarseFaces, sumOp<label>());

    if (Pstream::master())
    {
        Info << "\nTotal number of coarse faces: "<< totalNCoarseFaces << endl;
    }

    if (Pstream::master() && debug)
    {
        Pout << "\nView factor patches included in the calculation : "
             << viewFactorsPatches << endl;
    }

    // Collect local Cf and Sf on coarse mesh
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    DynamicList<point> localCoarseCf(nCoarseFaces);
    DynamicList<point> localCoarseSf(nCoarseFaces);
    DynamicList<label> localAgg(nCoarseFaces);
    labelHashSet includePatches;

    for (const label patchID : viewFactorsPatches)
    {
        const polyPatch& pp = patches[patchID];
        const labelList& agglom = finalAgglom[patchID];

        includePatches.insert(patchID);

        if (agglom.size() > 0)
        {
            label nAgglom = max(agglom)+1;
            labelListList coarseToFine(invertOneToMany(nAgglom, agglom));
            const labelList& coarsePatchFace =
                coarseMesh.patchFaceMap()[patchID];

            const pointField& coarseCf =
                coarseMesh.Cf().boundaryField()[patchID];
            const pointField& coarseSf =
                coarseMesh.Sf().boundaryField()[patchID];

            forAll(coarseCf, faceI)
            {
                point cf = coarseCf[faceI];

                const label coarseFaceI = coarsePatchFace[faceI];
                const labelList& fineFaces = coarseToFine[coarseFaceI];
                const label agglomI =
                    agglom[fineFaces[0]] + coarsePatches[patchID].start();

                // Construct single face
                uindirectPrimitivePatch upp
                (
                    UIndirectList<face>(pp, fineFaces),
                    pp.points()
                );

                List<point> availablePoints
                (
                    upp.faceCentres().size()
                  + upp.localPoints().size()
                );

                SubList<point>
                (
                    availablePoints,
                    upp.faceCentres().size()
                ) = upp.faceCentres();

                SubList<point>
                (
                    availablePoints,
                    upp.localPoints().size(),
                    upp.faceCentres().size()
                ) = upp.localPoints();

                point cfo = cf;
                scalar dist = GREAT;
                forAll(availablePoints, iPoint)
                {
                    point cfFine = availablePoints[iPoint];
                    if (mag(cfFine-cfo) < dist)
                    {
                        dist = mag(cfFine-cfo);
                        cf = cfFine;
                    }
                }

                point sf = coarseSf[faceI];
                localCoarseCf.append(cf);
                localCoarseSf.append(sf);
                localAgg.append(agglomI);

            }
        }
    }

    // Distribute local coarse Cf and Sf for shooting rays
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    List<pointField> remoteCoarseCf(Pstream::nProcs());
    List<pointField> remoteCoarseSf(Pstream::nProcs());
    List<labelField> remoteCoarseAgg(Pstream::nProcs());

    remoteCoarseCf[Pstream::myProcNo()] = localCoarseCf;
    remoteCoarseSf[Pstream::myProcNo()] = localCoarseSf;
    remoteCoarseAgg[Pstream::myProcNo()] = localAgg;

    Pstream::gatherList(remoteCoarseCf);
    Pstream::scatterList(remoteCoarseCf);
    Pstream::gatherList(remoteCoarseSf);
    Pstream::scatterList(remoteCoarseSf);
    Pstream::gatherList(remoteCoarseAgg);
    Pstream::scatterList(remoteCoarseAgg);


    globalIndex globalNumbering(nCoarseFaces);

    // Set up searching engine for obstacles
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #include "searchingEngineExt.H"

    // Determine rays between coarse face centres
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DynamicList<label> rayStartFace(nCoarseFaces + 0.01*nCoarseFaces);

    DynamicList<label> rayEndFace(rayStartFace.size());

    // Return rayStartFace in local index and rayEndFace in global index
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     #include "shootRaysExt.H"

    // Calculate number of visible faces from local index
    labelList nVisibleFaceFaces(nCoarseFaces, Zero);

    forAll(rayStartFace, i)
    {
        nVisibleFaceFaces[rayStartFace[i]]++;
    }

    labelListList visibleFaceFaces(nCoarseFaces);

    label nViewFactors = 0;
    forAll(nVisibleFaceFaces, faceI)
    {
        visibleFaceFaces[faceI].setSize(nVisibleFaceFaces[faceI]);
        nViewFactors += nVisibleFaceFaces[faceI];
    }

    // - Construct compact numbering
    // - return map from remote to compact indices
    //   (per processor (!= myProcNo) a map from remote index to compact index)
    // - construct distribute map
    // - renumber rayEndFace into compact addressing

    List<Map<label>> compactMap(Pstream::nProcs());

    mapDistribute map(globalNumbering, rayEndFace, compactMap);

    // visibleFaceFaces has:
    //    (local face, local viewed face) = compact viewed face
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    nVisibleFaceFaces = 0;
    forAll(rayStartFace, i)
    {
        label faceI = rayStartFace[i];
        label compactI = rayEndFace[i];
        visibleFaceFaces[faceI][nVisibleFaceFaces[faceI]++] = compactI;
    }

    // Construct data in compact addressing
    // (2AA) need coarse (Ai), fine Sf (dAi) and fine Cf(r) to calculate Fij
    // (2LI) need edges (li)
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    pointField compactCoarseCf(map.constructSize(), Zero);
    pointField compactCoarseSf(map.constructSize(), Zero);
    List<List<point>> compactFineSf(map.constructSize());
    List<List<point>> compactFineCf(map.constructSize());

    DynamicList<List<point>> compactPoints(map.constructSize());

    DynamicList<label> compactPatchId(map.constructSize());

    // Insert my coarse local values
    SubList<point>(compactCoarseSf, nCoarseFaces) = localCoarseSf;
    SubList<point>(compactCoarseCf, nCoarseFaces) = localCoarseCf;

    const faceList& faces = mesh.faces();

    // Insert my fine local values
    label compactI = 0;
    forAll(viewFactorsPatches, i)
    {
        label patchID = viewFactorsPatches[i];

        const labelList& agglom = finalAgglom[patchID];
        if (agglom.size() > 0)
        {
            label nAgglom = max(agglom)+1;
            labelListList coarseToFine(invertOneToMany(nAgglom, agglom));
            const labelList& coarsePatchFace =
                coarseMesh.patchFaceMap()[patchID];

            const polyPatch& pp = patches[patchID];

            forAll(coarseToFine, coarseI)
            {
                compactPatchId.append(patchID);
                List<point>& fineCf = compactFineCf[compactI];
                List<point>& fineSf = compactFineSf[compactI];

                label startFace = pp.start();

                const vectorField locPoints
                (
                    mesh.points(),
                    faces[coarseI + startFace]
                );

                const label coarseFaceI = coarsePatchFace[coarseI];
                const labelList& fineFaces = coarseToFine[coarseFaceI];

                fineCf.setSize(fineFaces.size());
                fineSf.setSize(fineFaces.size());

                compactPoints.append(locPoints);

                fineCf = UIndirectList<point>
                (
                    mesh.Cf().boundaryField()[patchID],
                    coarseToFine[coarseFaceI]
                );
                fineSf = UIndirectList<point>
                (
                    mesh.Sf().boundaryField()[patchID],
                    coarseToFine[coarseFaceI]
                );

                compactI++;
            }
        }
    }

    if (Pstream::master())
    {
        Info<< "map distribute..."  << endl;
    }

    // Do all swapping
    map.distribute(compactCoarseSf);
    map.distribute(compactCoarseCf);
    map.distribute(compactFineCf);
    map.distribute(compactFineSf);

    map.distribute(compactPoints);

    map.distribute(compactPatchId);

    // Plot all rays between visible faces.
    if (dumpRays)
    {
        writeRays
        (
            runTime.path()/"allVisibleFaces.obj",
            compactCoarseCf,
            remoteCoarseCf[Pstream::myProcNo()],
            visibleFaceFaces
        );
    }


    // Fill local view factor matrix
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    scalarListIOList F2LI
    (
        IOobject
        (
            "F",
            mesh.facesInstance(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        nCoarseFaces
    );

    label totalPatches = coarsePatches.size();
    reduce(totalPatches, maxOp<label>());

    // Matrix sum in j(Fij) for each i (if enclosure sum = 1)
    scalarSquareMatrix sumViewFactorPatch
    (
        totalPatches,
        0.0
    );

    scalarList patchArea(totalPatches, Zero);

    if (Pstream::master())
    {
        Info<< "\nCalculating view factors..." << endl;
    }

    List<scalarList> GaussPoints(5);
    GaussPoints[0].setSize(1);
    GaussPoints[0] = 0;

    GaussPoints[1].setSize(2);
    GaussPoints[1][0] =  1/std::sqrt(3);
    GaussPoints[1][1] = -1/std::sqrt(3);

    GaussPoints[2].setSize(3);
    GaussPoints[2][0] =  0;
    GaussPoints[2][1] =  0.774597;
    GaussPoints[2][2] = -0.774597;

    GaussPoints[3].setSize(4);
    GaussPoints[3][0] =  0.339981;
    GaussPoints[3][1] = -0.339981;
    GaussPoints[3][2] = 0.861136;
    GaussPoints[3][3] = -0.861136;

    GaussPoints[4].setSize(5);
    GaussPoints[4][0] =  0;
    GaussPoints[4][1] = 0.538469;
    GaussPoints[4][2] = -0.538469;
    GaussPoints[4][3] = 0.90618;
    GaussPoints[4][4] = -0.90618;


    List<scalarList> GaussWeights(5);
    GaussWeights[0].setSize(1);
    GaussWeights[0] = 2;

    GaussWeights[1].setSize(2);
    GaussWeights[1][0] =  1;
    GaussWeights[1][1] =  1;

    GaussWeights[2].setSize(3);
    GaussWeights[2][0] =  0.888889;
    GaussWeights[2][1] =  0.555556;
    GaussWeights[2][2] =  0.555556;

    GaussWeights[3].setSize(4);
    GaussWeights[3][0] =  0.652145;
    GaussWeights[3][1] =  0.652145;
    GaussWeights[3][2] =  0.347855;
    GaussWeights[3][3] =  0.347855;

    GaussWeights[4].setSize(5);
    GaussWeights[4][0] =  0.568889;
    GaussWeights[4][1] =  0.478629;
    GaussWeights[4][2] =  0.478629;
    GaussWeights[4][3] =  0.236927;
    GaussWeights[4][4] =  0.236927;

    label maxQuadOrder = 5;

    if (mesh.nSolutionD() == 3)
    {
        forAll(localCoarseSf, coarseFaceI)
        {
            const List<point>& localFineSf = compactFineSf[coarseFaceI];
            const vector Ai = sum(localFineSf);
            const List<point>& localFineCf = compactFineCf[coarseFaceI];
            const label fromPatchId = compactPatchId[coarseFaceI];

            const List<point>& lPoints = compactPoints[coarseFaceI];

            patchArea[fromPatchId] += mag(Ai);

            const labelList& visCoarseFaces = visibleFaceFaces[coarseFaceI];

            forAll(visCoarseFaces, visCoarseFaceI)
            {
                //F2AI[coarseFaceI].setSize(visCoarseFaces.size());
                F2LI[coarseFaceI].setSize(visCoarseFaces.size());
                label compactJ = visCoarseFaces[visCoarseFaceI];
                const List<point>& remoteFineSj = compactFineSf[compactJ];
                const List<point>& remoteFineCj = compactFineCf[compactJ];

                const List<point>& rPointsCj = compactPoints[compactJ];

                const label toPatchId = compactPatchId[compactJ];

                bool far(false);
                // Relative distance
                forAll(localFineSf, i)
                {
                    const scalar dAi =
                        Foam::sqrt
                        (
                            mag(localFineSf[i])/constant::mathematical::pi
                        );
                    const vector& dCi = localFineCf[i];

                    forAll(remoteFineSj, j)
                    {
                        const scalar dAj =
                            Foam::sqrt
                            (
                                mag(remoteFineSj[j])/constant::mathematical::pi
                            );
                        const vector& dCj = remoteFineCj[j];

                        const scalar dist = mag(dCi - dCj)/((dAi + dAj)/2);

                        if (dist > distTol)
                        {
                            far = true;
                        }
                    }
                }

                if (far)
                {
                    // 2AI method
                    scalar F2AIij = 0;

                    forAll(localFineSf, i)
                    {
                        const vector& dAi = localFineSf[i];
                        const vector& dCi = localFineCf[i];

                        forAll(remoteFineSj, j)
                        {
                            const vector& dAj = remoteFineSj[j];
                            const vector& dCj = remoteFineCj[j];

                            scalar dIntFij = calculateViewFactorFij2AI
                            (
                                dCi,
                                dCj,
                                dAi,
                                dAj
                            );

                            F2AIij += dIntFij;
                        }
                    }
                    F2LI[coarseFaceI][visCoarseFaceI] = F2AIij/mag(Ai);
                }
                else
                {
                    // 2LI method
                    label nLocal = lPoints.size();
                    label nRemote = rPointsCj.size();

                    // Using sub-divisions (quadrature)
                    scalar oldEToeInt = 0;
                    for (label gi=0; gi < maxQuadOrder; gi++)
                    {
                        scalar F2LIij = 0;
                        for(label i=0; i<nLocal; i++)
                        {
                            vector si;
                            vector ci;

                            vector sj;
                            vector cj;

                            if (i == 0)
                            {
                                si = lPoints[i] - lPoints[nLocal-1];
                                ci = (lPoints[i] + lPoints[nLocal-1])/2;
                            }
                            else
                            {
                                si = lPoints[i] - lPoints[i-1];
                                ci = (lPoints[i] + lPoints[i-1])/2;
                            }

                            for(label j=0; j<nRemote; j++)
                            {
                                if (j == 0)
                                {
                                    sj = rPointsCj[j]-rPointsCj[nRemote-1];
                                    cj = (rPointsCj[j]+rPointsCj[nRemote-1])/2;
                                }
                                else
                                {
                                    sj = rPointsCj[j] - rPointsCj[j-1];
                                    cj = (rPointsCj[j] + rPointsCj[j-1])/2;
                                }


                                scalar magSi = mag(si);
                                scalar magSj = mag(sj);
                                scalar cosij = (si & sj)/(magSi * magSj);

                                vector di = si/magSi;
                                vector dj = sj/magSj;

                                label quadOrder = gi;
                                const vector r(ci - cj);
                                // Common edges use n = 0
                                if (mag(r) < SMALL)
                                {
                                    quadOrder = 0;
                                }

                                scalar dIntFij =
                                    GaussQuad
                                    (
                                        GaussWeights[gi],
                                        GaussPoints[gi],
                                        magSi,
                                        magSj,
                                        di,
                                        dj,
                                        ci,
                                        cj,
                                        cosij,
                                        alpha,
                                        quadOrder
                                    );

                                F2LIij += dIntFij;
                            }
                        }

                        scalar err = (F2LIij-oldEToeInt)/F2LIij;

                        if
                        (
                            (mag(err) < tolGaussQuad && gi > 0)
                         || gi == maxQuadOrder-1
                        )
                        {
                            //DebugVar(gi)
                            F2LI[coarseFaceI][visCoarseFaceI] =
                                F2LIij/mag(Ai)/4/constant::mathematical::pi;
                            break;
                        }
                        else
                        {
                            oldEToeInt = F2LIij;
                        }
                    }
                }

                sumViewFactorPatch[fromPatchId][toPatchId] +=
                    F2LI[coarseFaceI][visCoarseFaceI]*mag(Ai);
            }
        }
    }
    else
    {
         FatalErrorInFunction
            << " View factors are not available in 2D "
            << exit(FatalError);
    }



    // Write view factors matrix in listlist form
    F2LI.write();

    reduce(sumViewFactorPatch, sumOp<scalarSquareMatrix>());
    reduce(patchArea, sumOp<scalarList>());

    if (Pstream::master() && debug)
    {
        forAll(viewFactorsPatches, i)
        {
            label patchI =  viewFactorsPatches[i];
            for (label j=i; j<viewFactorsPatches.size(); j++)
            {
                label patchJ =  viewFactorsPatches[j];

                Info << "F" << patchI << patchJ << ": "
                     << sumViewFactorPatch[patchI][patchJ]/patchArea[patchI]
                     << endl;
            }
        }
    }

    if (writeViewFactors)
    {
        if (Pstream::master())
        {
            Info << "Writing view factor matrix..." << endl;
        }

        volScalarField viewFactorField
        (
            IOobject
            (
                "viewFactorField",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar(dimless, Zero)
        );

        label compactI = 0;

        volScalarField::Boundary& vfbf = viewFactorField.boundaryFieldRef();
        forAll(viewFactorsPatches, i)
        {
            label patchID = viewFactorsPatches[i];
            const labelList& agglom = finalAgglom[patchID];
            if (agglom.size() > 0)
            {
                label nAgglom = max(agglom)+1;
                labelListList coarseToFine(invertOneToMany(nAgglom, agglom));
                const labelList& coarsePatchFace =
                    coarseMesh.patchFaceMap()[patchID];

                forAll(coarseToFine, coarseI)
                {
                    scalar FiSum = sum(F2LI[compactI]);

                    const label coarseFaceID = coarsePatchFace[coarseI];
                    const labelList& fineFaces = coarseToFine[coarseFaceID];
                    forAll(fineFaces, fineId)
                    {
                        const label faceID = fineFaces[fineId];
                        vfbf[patchID][faceID] = FiSum;
                    }
                    compactI++;
                }
            }
        }
        viewFactorField.write();
    }


    // Invert compactMap (from processor+localface to compact) to go
    // from compact to processor+localface (expressed as a globalIndex)
    // globalIndex globalCoarFaceNum(coarseMesh.nFaces());
    labelList compactToGlobal(map.constructSize());

    // Local indices first (note: are not in compactMap)
    for (label i = 0; i < globalNumbering.localSize(); i++)
    {
        compactToGlobal[i] = globalNumbering.toGlobal(i);
    }


    forAll(compactMap, procI)
    {
        const Map<label>& localToCompactMap = compactMap[procI];

        forAllConstIters(localToCompactMap, iter)
        {
            compactToGlobal[*iter] = globalNumbering.toGlobal
            (
                procI,
                iter.key()
            );
        }
    }


    labelListList globalFaceFaces(visibleFaceFaces.size());

    // Create globalFaceFaces needed to insert view factors
    // in F to the global matrix Fmatrix
    forAll(globalFaceFaces, faceI)
    {
        globalFaceFaces[faceI] = renumber
        (
            compactToGlobal,
            visibleFaceFaces[faceI]
        );
    }

    labelListIOList IOglobalFaceFaces
    (
        IOobject
        (
            "globalFaceFaces",
            mesh.facesInstance(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        std::move(globalFaceFaces)
    );
    IOglobalFaceFaces.write();


    IOmapDistribute IOmapDist
    (
        IOobject
        (
            "mapDist",
            mesh.facesInstance(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        std::move(map)
    );

    IOmapDist.write();

    Info<< "End\n" << endl;
    return 0;
}


// ************************************************************************* //

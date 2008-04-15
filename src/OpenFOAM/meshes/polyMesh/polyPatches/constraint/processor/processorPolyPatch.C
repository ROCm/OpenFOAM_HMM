/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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

#include "processorPolyPatch.H"
#include "polyBoundaryMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "dictionary.H"
#include "SubField.H"
#include "demandDrivenData.H"
#include "matchPoints.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(processorPolyPatch, 0);
    addToRunTimeSelectionTable(polyPatch, processorPolyPatch, dictionary);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::processorPolyPatch::processorPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const int myProcNo,
    const int neighbProcNo
)
:
    coupledPolyPatch(name, size, start, index, bm),
    myProcNo_(myProcNo),
    neighbProcNo_(neighbProcNo),
    neighbFaceCentres_(),
    neighbFaceAreas_(),
    neighbFaceCellCentres_(),
    neighbPointsPtr_(NULL),
    neighbEdgesPtr_(NULL)
{}


Foam::processorPolyPatch::processorPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm
)
:
    coupledPolyPatch(name, dict, index, bm),
    myProcNo_(readLabel(dict.lookup("myProcNo"))),
    neighbProcNo_(readLabel(dict.lookup("neighbProcNo"))),
    neighbFaceCentres_(),
    neighbFaceAreas_(),
    neighbFaceCellCentres_(),
    neighbPointsPtr_(NULL),
    neighbEdgesPtr_(NULL)
{}


Foam::processorPolyPatch::processorPolyPatch
(
    const processorPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    coupledPolyPatch(pp, bm),
    myProcNo_(pp.myProcNo_),
    neighbProcNo_(pp.neighbProcNo_),
    neighbFaceCentres_(),
    neighbFaceAreas_(),
    neighbFaceCellCentres_(),
    neighbPointsPtr_(NULL),
    neighbEdgesPtr_(NULL)
{}


Foam::processorPolyPatch::processorPolyPatch
(
    const processorPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart
)
:
    coupledPolyPatch(pp, bm, index, newSize, newStart),
    myProcNo_(pp.myProcNo_),
    neighbProcNo_(pp.neighbProcNo_),
    neighbFaceCentres_(),
    neighbFaceAreas_(),
    neighbFaceCellCentres_(),
    neighbPointsPtr_(NULL),
    neighbEdgesPtr_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::processorPolyPatch::~processorPolyPatch()
{
    deleteDemandDrivenData(neighbPointsPtr_);
    deleteDemandDrivenData(neighbEdgesPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::processorPolyPatch::initGeometry()
{
    if (Pstream::parRun())
    {
        OPstream toNeighbProc
        (
            Pstream::blocking,
            neighbProcNo(),
          + 3*(sizeof(label) + size()*sizeof(vector))
        );

        toNeighbProc
            << faceCentres()
            << faceAreas()
            << faceCellCentres();
    }
}


void Foam::processorPolyPatch::calcGeometry()
{
    if (Pstream::parRun())
    {
        {
            IPstream fromNeighbProc
            (
                Pstream::blocking,
                neighbProcNo(),
                3*(sizeof(label) + size()*sizeof(vector))
            );
            fromNeighbProc
                >> neighbFaceCentres_
                >> neighbFaceAreas_
                >> neighbFaceCellCentres_;
        }

        scalarField magSf = mag(faceAreas());

        forAll(magSf, facei)
        {
            scalar nmagSf = mag(neighbFaceAreas_[facei]);
            scalar avSf = (magSf[facei] + nmagSf)/2.0;

            if (avSf > VSMALL && mag(magSf[facei] - nmagSf)/avSf > 1e-4)
            {
                FatalErrorIn
                (
                    "processorPolyPatch::calcGeometry()"
                )   << "face " << facei << " area does not match neighbour by "
                    << 100*mag(magSf[facei] - nmagSf)/avSf
                    << "% -- possible face ordering problem." << endl
                    << "patch:" << name() << " mesh face:" << start()+facei
                    << " face centre:" << faceCentres()[facei]
                    << " my area:" << magSf[facei]
                    << " neighbour area:" << nmagSf
                    << exit(FatalError);
            }
        }

        calcTransformTensors
        (
            faceCentres(),
            neighbFaceCentres_,
            faceNormals(),
            neighbFaceAreas_/(mag(neighbFaceAreas_)+VSMALL)
        );
    }
}


void Foam::processorPolyPatch::initMovePoints(const pointField& p)
{
    polyPatch::movePoints(p);
    processorPolyPatch::initGeometry();
}


void Foam::processorPolyPatch::movePoints(const pointField&)
{
    processorPolyPatch::calcGeometry();
}


void Foam::processorPolyPatch::initUpdateMesh()
{
    polyPatch::initUpdateMesh();

    deleteDemandDrivenData(neighbPointsPtr_);
    deleteDemandDrivenData(neighbEdgesPtr_);

    if (Pstream::parRun())
    {
        // Express all points as patch face and index in face.
        labelList patchFace(nPoints());
        labelList indexInFace(nPoints());

        for (label patchPointI = 0; patchPointI < nPoints(); patchPointI++)
        {
            label faceI = pointFaces()[patchPointI][0];

            patchFace[patchPointI] = faceI;

            const face& f = localFaces()[faceI];

            indexInFace[patchPointI] = findIndex(f, patchPointI);
        }

        OPstream toNeighbProc
        (
            Pstream::blocking,
            neighbProcNo(),
            3*sizeof(label)
          + 2*nPoints()*sizeof(label)
          + nEdges()*sizeof(edge)
        );

        toNeighbProc
            << patchFace
            << indexInFace
            << edges();
    }
}


void Foam::processorPolyPatch::updateMesh()
{
    // For completeness
    polyPatch::updateMesh();

    if (Pstream::parRun())
    {
        labelList nbrPatchFace(nPoints());
        labelList nbrIndexInFace(nPoints());
        edgeList nbrEdges(nEdges());

        {
            // Note cannot predict exact size since edgeList not (yet) sent as
            // binary entity but as List of edges.
            IPstream fromNeighbProc(Pstream::blocking, neighbProcNo());

            fromNeighbProc
                >> nbrPatchFace
                >> nbrIndexInFace
                >> nbrEdges;
        }

        if (nbrPatchFace.size() == nPoints() && nbrEdges.size() == nEdges())
        {
            // Convert neighbour edges and indices into face back into
            // my edges and points.
            neighbPointsPtr_ = new labelList(nPoints());
            labelList& neighbPoints = *neighbPointsPtr_;

            // Inverse of neighbPoints so from neighbour point to current point.
            labelList nbrToThis(nPoints(), -1);

            forAll(nbrPatchFace, nbrPointI)
            {
                // Find face and index in face on this side.
                const face& f = localFaces()[nbrPatchFace[nbrPointI]];
                label index = (f.size() - nbrIndexInFace[nbrPointI]) % f.size();
                label patchPointI = f[index];

                neighbPoints[patchPointI] = nbrPointI;
                nbrToThis[nbrPointI] = patchPointI;
            }

            // Convert edges.
            neighbEdgesPtr_ = new labelList(nEdges());
            labelList& neighbEdges = *neighbEdgesPtr_;

            forAll(nbrEdges, nbrEdgeI)
            {
                const edge& nbrEdge = nbrEdges[nbrEdgeI];

                // Get edge in local point numbering
                edge e(nbrToThis[nbrEdge[0]], nbrToThis[nbrEdge[1]]);

                // Find the edge.
                const labelList& pEdges = pointEdges()[e[0]];

                label edgeI = -1;

                forAll(pEdges, i)
                {
                    if (edges()[pEdges[i]] == e)
                    {
                        edgeI = pEdges[i];
                        break;
                    }
                }

                if (edgeI == -1)
                {
                    if (debug)
                    {
                        WarningIn("processorPolyPatch::updateMesh()")
                            << "Patch:" << name()
                            << "Cannot find patch edge with vertices " << e
                            << " coords:"
                            << localPoints()[e[0]]<< localPoints()[e[1]]
                            << " on patch " << name() << nl
                            << "Can only find edges "
                            << IndirectList<edge>(edges(), pEdges)()
                            << " connected to first vertex" << nl
                            << "Either your mesh is incorrect or this patch"
                            << " was constructed from part of a cyclic patch."
                            << nl
                            << "Not calculating edge neighbour addressing."
                            << endl;
                    }
                    deleteDemandDrivenData(neighbEdgesPtr_);
                    break;
                }

                neighbEdges[edgeI] = nbrEdgeI;
            }
        }
        else
        {
            // Differing number of points or edges. Probably patch includes
            // part of a cyclic.
            neighbPointsPtr_ = NULL;
            neighbEdgesPtr_ = NULL;

            if (debug)
            {
                {
                    fileName nm(name()+"_faces.obj");
                    Pout<< "processorPolyPatch::order : Writing my " << size()
                        << " faces to OBJ file " << nm << endl;
                    writeOBJ(nm, *this, points());
                }

                WarningIn("processorPolyPatch::updateMesh()")
                    << "Patch:" << name()
                    << "my nPoints:" << nPoints()
                    << " edges:" << nEdges()
                    << " nbrPatchFace:" << nbrPatchFace.size()
                    << " nbrEdges:" << nbrEdges.size()
                    << endl;
            }
        }
    }
}


const Foam::labelList& Foam::processorPolyPatch::neighbPoints() const
{
    if (!neighbPointsPtr_)
    {
        // Was probably created from cyclic patch and hence the
        // number of edges or points might differ on both
        // sides of the processor patch since one side might have
        // it merged with another bit of geometry

        FatalErrorIn("processorPolyPatch::neighbPoints() const")
            << "No extended addressing calculated for patch " << name()
            << nl
            << "This can happen if the number of points or edges on both"
            << " sides of the two coupled patches differ." << nl
            << "This happens if the processorPatch was constructed from"
            << " part of a cyclic patch."
            << abort(FatalError);
    }
    return *neighbPointsPtr_;
}


const Foam::labelList& Foam::processorPolyPatch::neighbEdges() const
{
    if (!neighbEdgesPtr_)
    {
        // Was probably created from cyclic patch and hence the
        // number of edges or points might differ on both
        // sides of the processor patch since one side might have
        // it merged with another bit of geometry

        FatalErrorIn("processorPolyPatch::neighbEdges() const")
            << "No extended addressing calculated for patch " << name()
            << nl
            << "This can happen if the number of points or edges on both"
            << " sides of the two coupled patches differ." << nl
            << "This happens if the processorPatch was constructed from"
            << " part of a cyclic patch."
            << abort(FatalError);
    }
    return *neighbEdgesPtr_;
}


void Foam::processorPolyPatch::initOrder(const primitivePatch& pp) const
{
    if (!Pstream::parRun())
    {
        return;
    }

    if (debug)
    {
        fileName nm(name()+"_faces.obj");
        Pout<< "processorPolyPatch::order : Writing my " << pp.size()
            << " faces to OBJ file " << nm << endl;
        writeOBJ(nm, pp, pp.points());

        // Calculate my face centres
        pointField ctrs(calcFaceCentres(pp, pp.points()));

        OFstream localStr(name() + "_localFaceCentres.obj");
        Pout<< "processorPolyPatch::order : "
            << "Dumping " << ctrs.size()
            << " local faceCentres to " << localStr.name() << endl;

        forAll(ctrs, faceI)
        {
            writeOBJ(localStr, ctrs[faceI]);
        }
    }

    const bool isMaster = Pstream::myProcNo() < neighbProcNo();

    // Check (on old patch!) for weirdness.
    if (separated() || !parallel())
    {
        WarningIn
        (
            "processorPolyPatch::initOrder(const primitivePatch&) const"
        )   << "in patch:" << name() << " : "
            << "using geometric matching on this processor patch might fail"
            << " since it has 'separated' faces or is not 'parallel'"
            << endl;
    }


    if (isMaster)
    {
        pointField ctrs(calcFaceCentres(pp, pp.points()));

        pointField anchorPoints(getAnchorPoints(pp, pp.points()));

        // Now send all info over to the neighbour
        OPstream toNeighbour(Pstream::blocking, neighbProcNo());
        toNeighbour << ctrs << anchorPoints;
    }
}


// Return new ordering. Ordering is -faceMap: for every face index
// the new face -rotation:for every new face the clockwise shift
// of the original face. Return false if nothing changes (faceMap
// is identity, rotation is 0)
bool Foam::processorPolyPatch::order
(
    const primitivePatch& pp,
    labelList& faceMap,
    labelList& rotation
) const
{
    if (!Pstream::parRun())
    {
        return false;
    }

    faceMap.setSize(pp.size());
    faceMap = -1;

    rotation.setSize(pp.size());
    rotation = 0;

    const bool isMaster = Pstream::myProcNo() < neighbProcNo();

    if (isMaster)
    {
        // Do nothing (i.e. identical mapping, zero rotation).
        // See comment at top.
        forAll(faceMap, patchFaceI)
        {
            faceMap[patchFaceI] = patchFaceI;
        }

        return false;
    }
    else
    {
        vectorField masterCtrs;
        vectorField masterAnchorPoints;

        // Receive data from neighbour
        {
            IPstream fromNeighbour(Pstream::blocking, neighbProcNo());
            fromNeighbour >> masterCtrs >> masterAnchorPoints;
        }

        // Calculate my face centres
        pointField ctrs(calcFaceCentres(pp, pp.points()));

        // Calculate typical distance from face centre
        scalarField tols(calcFaceTol(pp, pp.points(), ctrs));

        if (debug || masterCtrs.size() != pp.size())
        {
            {
                OFstream nbrStr(name() + "_nbrFaceCentres.obj");
                Pout<< "processorPolyPatch::order : "
                    << "Dumping neighbour faceCentres to " << nbrStr.name()
                    << endl;
                forAll(masterCtrs, faceI)
                {
                    writeOBJ(nbrStr, masterCtrs[faceI]);
                }
            }

            if (masterCtrs.size() != pp.size())
            {
                FatalErrorIn
                (
                    "processorPolyPatch::order(const primitivePatch&"
                    ", labelList&, labelList&) const"
                )   << "in patch:" << name() << " : "
                    << "Local size of patch is " << pp.size() << " (faces)."
                    << endl
                    << "Received from neighbour " << masterCtrs.size()
                    << " faceCentres!"
                    << abort(FatalError);
            }
        }

        // Geometric match of face centre vectors
        bool matchedAll = matchPoints(ctrs, masterCtrs, tols, true, faceMap);

        if (debug)
        {
            OFstream ccStr(name() + "_faceCentresConnections.obj");

            Pout<< "processorPolyPatch::order : "
                << "Dumping newly found match as lines between"
                << " corresponding face centres to OBJ file " << ccStr.name()
                << endl;

            label vertI = 0;

            forAll(ctrs, faceI)
            {
                label masterFaceI = faceMap[faceI];

                if (masterFaceI != -1)
                {
                    const point& c0 = masterCtrs[masterFaceI];
                    const point& c1 = ctrs[faceI];
                    writeOBJ(ccStr, c0, c1, vertI);
                }
            }
        }

        if (!matchedAll)
        {
            SeriousErrorIn
            (
                "processorPolyPatch::order(const primitivePatch&"
                ", labelList&, labelList&) const"
            )   << "in patch:" << name() << " : "
                << "Cannot match vectors to faces on both sides of patch"
                << endl
                << "masterCtrs[0]:" << masterCtrs[0] << endl
                << "ctrs[0]:" << ctrs[0] << endl
                << "Please check your topology changes or maybe you have"
                << " multiple separated (from cyclics) processor patches"
                << endl
                << "Continuing with incorrect face ordering from now on!"
                << endl;

            return false;
        }

        // Set rotation.
        forAll(faceMap, oldFaceI)
        {
            // The face f will be at newFaceI (after morphing) and we want its
            // anchorPoint (= f[0]) to align with the anchorpoint for the
            // corresponding face on the other side.

            label newFaceI = faceMap[oldFaceI];

            const point& wantedAnchor = masterAnchorPoints[newFaceI];

            rotation[newFaceI] = getRotation
            (
                pp.points(),
                pp[oldFaceI],
                wantedAnchor,
                tols[oldFaceI]
            );

            if (rotation[newFaceI] == -1)
            {
                SeriousErrorIn
                (
                    "processorPolyPatch::order(const primitivePatch&"
                    ", labelList&, labelList&) const"
                )   << "in patch:" << name()
                    << " : "
                    << "Cannot find point on face " << pp[oldFaceI]
                    << " with vertices:"
                    << IndirectList<point>(pp.points(), pp[oldFaceI])()
                    << " that matches point " << wantedAnchor
                    << " when matching the halves of processor patch " << name()
                    << "Continuing with incorrect face ordering from now on!"
                    << endl;

                return false;
            }
        }

        forAll(faceMap, faceI)
        {
            if (faceMap[faceI] != faceI || rotation[faceI] != 0)
            {
                return true;
            }
        }

        return false;
    }
}


void Foam::processorPolyPatch::write(Ostream& os) const
{
    polyPatch::write(os);
    os.writeKeyword("myProcNo") << myProcNo_
        << token::END_STATEMENT << nl;
    os.writeKeyword("neighbProcNo") << neighbProcNo_
        << token::END_STATEMENT << nl;
}


// ************************************************************************* //

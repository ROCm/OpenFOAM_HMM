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

#include "processorPolyPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "dictionary.H"
#include "SubField.H"
#include "demandDrivenData.H"
#include "matchPoints.H"
#include "OFstream.H"
#include "polyMesh.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(processorPolyPatch, 0);
    addToRunTimeSelectionTable(polyPatch, processorPolyPatch, dictionary);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::processorPolyPatch::checkSubPatches() const
{
    if (starts_.size() != patchIDs_.size()+1)
    {
        FatalErrorIn("processorPolyPatch::checkPatches() const")
            << "starts should have one more element than patchIDs." << endl
            << "starts:" << starts_ << " patchIDs:" << patchIDs_
            << exit(FatalError);
    }
    if (starts_[0] != 0)
    {
        FatalErrorIn("processorPolyPatch::checkPatches() const")
            << "starts[0] should be 0." << endl
            << "starts:" << starts_ << " patchIDs:" << patchIDs_
            << exit(FatalError);
    }
    if (starts_[starts_.size()-1] != size())
    {
        FatalErrorIn("processorPolyPatch::checkPatches() const")
            << "Last element in starts should be the size." << endl
            << "starts:" << starts_ << " size:" << size()
            << exit(FatalError);
    }
    // If there are any faces from internal they should be first.
    for (label i = 1; i < patchIDs_.size(); i++)
    {
        if (patchIDs_[i] == -1)
        {
            FatalErrorIn("processorPolyPatch::checkPatches() const")
            << "Faces originating from internal faces (subPatch is -1)"
            << " should always be first in the patch."
            << " The current list of subpatches " << patchIDs_
            << exit(FatalError);
        }
    }
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
    const int neighbProcNo,
    const labelList& patchIDs,
    const labelList& starts
)
:
    coupledPolyPatch(name, size, start, index, bm),
    myProcNo_(myProcNo),
    neighbProcNo_(neighbProcNo),
    patchIDs_(patchIDs),
    starts_(starts),
    neighbFaceCentres_(),
    neighbFaceAreas_(),
    neighbFaceCellCentres_()
//    neighbPointsPtr_(NULL),
//    neighbEdgesPtr_(NULL)
{
    checkSubPatches();
}


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
    patchIDs_(dict.lookup("patchIDs")),
    starts_(dict.lookup("starts")),
    neighbFaceCentres_(),
    neighbFaceAreas_(),
    neighbFaceCellCentres_()
//    neighbPointsPtr_(NULL),
//    neighbEdgesPtr_(NULL)
{
    checkSubPatches();
}


Foam::processorPolyPatch::processorPolyPatch
(
    const processorPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    coupledPolyPatch(pp, bm),
    myProcNo_(pp.myProcNo_),
    neighbProcNo_(pp.neighbProcNo_),
    patchIDs_(pp.patchIDs_),
    starts_(pp.starts_),
    neighbFaceCentres_(),
    neighbFaceAreas_(),
    neighbFaceCellCentres_()
//    neighbPointsPtr_(NULL),
//    neighbEdgesPtr_(NULL)
{}


Foam::processorPolyPatch::processorPolyPatch
(
    const processorPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart,
    const labelList& patchIDs,
    const labelList& starts
)
:
    coupledPolyPatch(pp, bm, index, newSize, newStart),
    myProcNo_(pp.myProcNo_),
    neighbProcNo_(pp.neighbProcNo_),
    patchIDs_(patchIDs),
    starts_(starts),
    neighbFaceCentres_(),
    neighbFaceAreas_(),
    neighbFaceCellCentres_()
//    neighbPointsPtr_(NULL),
//    neighbEdgesPtr_(NULL)
{
    checkSubPatches();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::processorPolyPatch::~processorPolyPatch()
{
//    deleteDemandDrivenData(neighbPointsPtr_);
//    deleteDemandDrivenData(neighbEdgesPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::processorPolyPatch::transformPosition(pointField& l) const
{
    if (l.size() != size())
    {
        FatalErrorIn("processorPolyPatch::transformPosition(pointField&) const")
            << "Size of field " << l.size() << " differs from patch size "
            << size() << abort(FatalError);
    }

    forAll(patchIDs_, subI)
    {
        label patchI = patchIDs_[subI];

        if (patchI != -1)
        {
            // Get field on patch
            SubField<point> subFld(subSlice(l, subI));

            refCast<const coupledPolyPatch>
            (
                boundaryMesh()[patchI]
            ).transformPosition(reinterpret_cast<pointField&>(subFld));
        }
    }
}


void Foam::processorPolyPatch::initGeometry()
{
Pout<< "**processorPolyPatch::initGeometry()" << endl;
    if (Pstream::parRun())
    {
        pointField fc(size());
        vectorField fa(size());
        pointField cc(size());

        forAll(patchIDs_, i)
        {
            label patchI = patchIDs_[i];

            autoPtr<primitivePatch> subPp = subPatch(i);

            SubField<point> subFc(subSlice(fc, i));
            SubField<vector> subFa(subSlice(fa, i));
            SubField<point> subCc(subSlice(cc, i));

            subFc.assign(subSlice(faceCentres(), i));
            subFa.assign(subSlice(faceAreas(), i));
            subCc.assign(subSlice(faceCellCentres()(), i));

            if (patchI != -1)
            {
                coupledPolyPatch& pp = const_cast<coupledPolyPatch&>
                (
                    refCast<const coupledPolyPatch>
                    (
                        boundaryMesh()[patchI]
                    )
                );

                Pout<< name() << " calling initGeometry on " << pp.name()
                    << endl;

                pp.initGeometry(subPp, subFc, subFa, subCc);

                Pout<< name() << " from patchI:" << patchI
                    << " calculated fc:" << subFc
                    << " fa:" << subFa << " subCC:" << subCc << endl;
                Pout<< name() << " fc:" << fc << endl;
            }
        }

        Pout<< name() << " fc:" << fc << " fa:" << fa << " cc:" << cc << endl;


        OPstream toNeighbProc
        (
            Pstream::blocking,
            neighbProcNo(),
            3*(sizeof(label) + size()*sizeof(vector) + sizeof(scalar))
        );

        toNeighbProc << fc << fa << cc;
    }
}


void Foam::processorPolyPatch::calcGeometry()
{
Pout<< "processorPolyPatch::calcGeometry() for " << name() << endl;
    if (Pstream::parRun())
    {
        {
            IPstream fromNeighbProc
            (
                Pstream::blocking,
                neighbProcNo(),
                3*(sizeof(label) + size()*sizeof(vector) + sizeof(scalar))
            );
            fromNeighbProc
                >> neighbFaceCentres_
                >> neighbFaceAreas_
                >> neighbFaceCellCentres_;
        }

Pout<< "processorPolyPatch::calcGeometry() : received data for "
    << neighbFaceCentres_.size() << " faces." << endl;

        forAll(patchIDs_, i)
        {
            label patchI = patchIDs_[i];

            if (patchI == -1)
            {
                // Anything needs doing for ex-internal faces?
            }
            else
            {
                coupledPolyPatch& pp = const_cast<coupledPolyPatch&>
                (
                    refCast<const coupledPolyPatch>
                    (
                        boundaryMesh()[patchI]
                    )
                );

Pout<< "processorPolyPatch::calcGeometry() : referring to " << pp.name()
    << " for faces size:" << starts_[i+1]-starts_[i]
    << " start:" << starts_[i]
    << endl;

                pp.calcGeometry
                (
                    subPatch(i),
                    subSlice(faceCentres(), i),
                    subSlice(faceAreas(), i),
                    subSlice(faceCellCentres()(), i),
                    subSlice(neighbFaceCentres_, i),
                    subSlice(neighbFaceAreas_, i),
                    subSlice(neighbFaceCellCentres_, i)
                );
            }
        }
    }
Pout<< "**neighbFaceCentres_:" << neighbFaceCentres_ << endl;
Pout<< "**neighbFaceAreas_:" << neighbFaceAreas_ << endl;
Pout<< "**neighbFaceCellCentres_:" << neighbFaceCellCentres_ << endl;

}


void Foam::processorPolyPatch::initMovePoints(const pointField& p)
{
    // Update local polyPatch quantities
    polyPatch::movePoints(p);
    // Recalculate geometry
    processorPolyPatch::initGeometry();
}


void Foam::processorPolyPatch::movePoints(const pointField&)
{
    processorPolyPatch::calcGeometry();
}


void Foam::processorPolyPatch::initUpdateMesh()
{
    polyPatch::initUpdateMesh();

//    deleteDemandDrivenData(neighbPointsPtr_);
//    deleteDemandDrivenData(neighbEdgesPtr_);
//
//    if (Pstream::parRun())
//    {
//        // Express all points as patch face and index in face.
//        labelList pointFace(nPoints());
//        labelList pointIndex(nPoints());
//
//        for (label patchPointI = 0; patchPointI < nPoints(); patchPointI++)
//        {
//            label faceI = pointFaces()[patchPointI][0];
//
//            pointFace[patchPointI] = faceI;
//
//            const face& f = localFaces()[faceI];
//
//            pointIndex[patchPointI] = findIndex(f, patchPointI);
//        }
//
//        // Express all edges as patch face and index in face.
//        labelList edgeFace(nEdges());
//        labelList edgeIndex(nEdges());
//
//        for (label patchEdgeI = 0; patchEdgeI < nEdges(); patchEdgeI++)
//        {
//            label faceI = edgeFaces()[patchEdgeI][0];
//
//            edgeFace[patchEdgeI] = faceI;
//
//            const labelList& fEdges = faceEdges()[faceI];
//
//            edgeIndex[patchEdgeI] = findIndex(fEdges, patchEdgeI);
//        }
//
//        OPstream toNeighbProc
//        (
//            Pstream::blocking,
//            neighbProcNo(),
//            8*sizeof(label)             // four headers of labelList
//          + 2*nPoints()*sizeof(label)   // two point-based labellists
//          + 2*nEdges()*sizeof(label)    // two edge-based labelLists
//        );
//
//        toNeighbProc
//            << pointFace
//            << pointIndex
//            << edgeFace
//            << edgeIndex;
//    }
}


void Foam::processorPolyPatch::updateMesh()
{
    // For completeness
    polyPatch::updateMesh();

//    if (Pstream::parRun())
//    {
//        labelList nbrPointFace;
//        labelList nbrPointIndex;
//        labelList nbrEdgeFace;
//        labelList nbrEdgeIndex;
//
//        {
//            // Note cannot predict exact size since opposite nPoints might
//            // be different from one over here.
//            IPstream fromNeighbProc(Pstream::blocking, neighbProcNo());
//
//            fromNeighbProc
//                >> nbrPointFace
//                >> nbrPointIndex
//                >> nbrEdgeFace
//                >> nbrEdgeIndex;
//        }
//
//        // Convert neighbour faces and indices into face back into
//        // my edges and points.
//
//        // Convert points.
//        // ~~~~~~~~~~~~~~~
//
//        neighbPointsPtr_ = new labelList(nPoints(), -1);
//        labelList& neighbPoints = *neighbPointsPtr_;
//
//        forAll(nbrPointFace, nbrPointI)
//        {
//            // Find face and index in face on this side.
//            const face& f = localFaces()[nbrPointFace[nbrPointI]];
//            label index = (f.size() - nbrPointIndex[nbrPointI]) % f.size();
//            label patchPointI = f[index];
//
//            if (neighbPoints[patchPointI] == -1)
//            {
//                // First reference of point
//                neighbPoints[patchPointI] = nbrPointI;
//            }
//            else if (neighbPoints[patchPointI] >= 0)
//            {
//                // Point already visited. Mark as duplicate.
//                neighbPoints[patchPointI] = -2;
//            }
//        }
//
//        // Reset all duplicate entries to -1.
//        forAll(neighbPoints, patchPointI)
//        {
//            if (neighbPoints[patchPointI] == -2)
//            {
//                neighbPoints[patchPointI] = -1;
//            }
//        }
//
//        // Convert edges.
//        // ~~~~~~~~~~~~~~
//
//        neighbEdgesPtr_ = new labelList(nEdges(), -1);
//        labelList& neighbEdges = *neighbEdgesPtr_;
//
//        forAll(nbrEdgeFace, nbrEdgeI)
//        {
//            // Find face and index in face on this side.
//            const labelList& f = faceEdges()[nbrEdgeFace[nbrEdgeI]];
//            label index = (f.size() - nbrEdgeIndex[nbrEdgeI] - 1) % f.size();
//            label patchEdgeI = f[index];
//
//            if (neighbEdges[patchEdgeI] == -1)
//            {
//                // First reference of edge
//                neighbEdges[patchEdgeI] = nbrEdgeI;
//            }
//            else if (neighbEdges[patchEdgeI] >= 0)
//            {
//                // Edge already visited. Mark as duplicate.
//                neighbEdges[patchEdgeI] = -2;
//            }
//        }
//
//        // Reset all duplicate entries to -1.
//        forAll(neighbEdges, patchEdgeI)
//        {
//            if (neighbEdges[patchEdgeI] == -2)
//            {
//                neighbEdges[patchEdgeI] = -1;
//            }
//        }
//
//        // Remove any addressing used for shared points/edges calculation
//        primitivePatch::clearOut();
//    }
}


//const Foam::labelList& Foam::processorPolyPatch::neighbPoints() const
//{
//    if (!neighbPointsPtr_)
//    {
//        FatalErrorIn("processorPolyPatch::neighbPoints() const")
//            << "No extended addressing calculated for patch " << name()
//            << abort(FatalError);
//    }
//    return *neighbPointsPtr_;
//}
//
//
//const Foam::labelList& Foam::processorPolyPatch::neighbEdges() const
//{
//    if (!neighbEdgesPtr_)
//    {
//        FatalErrorIn("processorPolyPatch::neighbEdges() const")
//            << "No extended addressing calculated for patch " << name()
//            << abort(FatalError);
//    }
//    return *neighbEdgesPtr_;
//}


const Foam::labelListList& Foam::processorPolyPatch::subMeshPoints() const
{
    if (!subMeshPointsPtr_.valid())
    {
        subMeshPointsPtr_.reset(new labelListList(patchIDs_.size()));
        labelListList& meshPoints = subMeshPointsPtr_();

        forAll(patchIDs_, subI)
        {
            meshPoints[subI] = subPatch(subI)().meshPoints();
        }
    }
    return subMeshPointsPtr_();
}


const Foam::labelListList& Foam::processorPolyPatch::reverseSubMeshPoints()
 const
{
    if (!reverseSubMeshPointsPtr_.valid())
    {
        reverseSubMeshPointsPtr_.reset(new labelListList(patchIDs_.size()));
        labelListList& meshPoints = reverseSubMeshPointsPtr_();

        forAll(patchIDs_, subI)
        {
            label subStart = starts()[subI];
            label subSize = starts()[subI+1]-subStart;

            faceList reverseFaces(subSize);
            forAll(reverseFaces, i)
            {
                reverseFaces[i] = operator[](subStart+i).reverseFace();
            }
            meshPoints[subI] = primitivePatch
            (
                faceSubList
                (
                    reverseFaces,
                    reverseFaces.size()
                ),
                points()
            ).meshPoints();
        }
    }
    return reverseSubMeshPointsPtr_();
}


Foam::label Foam::processorPolyPatch::whichSubPatch(const label facei) const
{
    label index = findSortedIndex(starts_, facei);

    if (index < 0 || index >= starts_.size())
    {
        FatalErrorIn("processorPolyPatch::whichSubPatch(const label) const")
            << "Illegal local face index " << facei << " for patch " << name()
            << endl << "Face index should be between 0 and "
            << starts_[starts_.size()-1]-1 << abort(FatalError);
    }
    return patchIDs_[index];
}


void Foam::processorPolyPatch::initOrder(const primitivePatch& pp) const
{
    if (!Pstream::parRun())
    {
        return;
    }

    if (debug)
    {
        fileName nm
        (
            boundaryMesh().mesh().time().path()
           /name()+"_faces.obj"
        );
        Pout<< "processorPolyPatch::order : Writing my " << pp.size()
            << " faces to OBJ file " << nm << endl;
        writeOBJ(nm, pp, pp.points());

        // Calculate my face centres
        pointField ctrs(calcFaceCentres(pp, pp.points()));

        OFstream localStr
        (
            boundaryMesh().mesh().time().path()
           /name() + "_localFaceCentres.obj"
        );
        Pout<< "processorPolyPatch::order : "
            << "Dumping " << ctrs.size()
            << " local faceCentres to " << localStr.name() << endl;

        forAll(ctrs, faceI)
        {
            writeOBJ(localStr, ctrs[faceI]);
        }
    }

    if (owner())
    {
        pointField ctrs(calcFaceCentres(pp, pp.points()));

        pointField anchors(getAnchorPoints(pp, pp.points()));

        // Now send all info over to the neighbour
        OPstream toNeighbour(Pstream::blocking, neighbProcNo());
        toNeighbour << ctrs << anchors;
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
    // Note: we only get the faces that originate from internal faces.

    if (!Pstream::parRun())
    {
        return false;
    }

    faceMap.setSize(pp.size());
    faceMap = -1;

    rotation.setSize(pp.size());
    rotation = 0;

    if (owner())
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
        vectorField masterAnchors;

        // Receive data from neighbour
        {
            IPstream fromNeighbour(Pstream::blocking, neighbProcNo());
            fromNeighbour >> masterCtrs >> masterAnchors;
        }

        // Calculate my face centres
        pointField ctrs(calcFaceCentres(pp, pp.points()));

        // Calculate typical distance from face centre
        scalarField tols(calcFaceTol(pp, pp.points(), ctrs));

        if (debug || masterCtrs.size() != pp.size())
        {
            {
                OFstream nbrStr
                (
                    boundaryMesh().mesh().time().path()
                   /name() + "_nbrFaceCentres.obj"
                );
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
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // 1. Try existing ordering and transformation
        bool matchedAll = matchPoints(ctrs, masterCtrs, tols, true, faceMap);

        if (!matchedAll || debug)
        {
            // Dump faces
            fileName str
            (
                boundaryMesh().mesh().time().path()
               /name()/name()+"_faces.obj"
            );
            Pout<< "processorPolyPatch::order :"
                << " Writing faces to OBJ file " << str.name() << endl;
            writeOBJ(str, pp, pp.points());

            OFstream ccStr
            (
                boundaryMesh().mesh().time().path()
               /name() + "_faceCentresConnections.obj"
            );

            Pout<< "processorPolyPatch::order :"
                << " Dumping newly found match as lines between"
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
                << "    masterCtrs[0]:" << masterCtrs[0] << endl
                << "    ctrs[0]:" << ctrs[0] << endl
                << "    Please check your topology changes or maybe you have"
                << " multiple separated (from cyclics) processor patches"
                << endl
                << "    Continuing with incorrect face ordering from now on!"
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

            const point& wantedAnchor = masterAnchors[newFaceI];

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
                )   << "in patch " << name()
                    << " : "
                    << "Cannot find point on face " << pp[oldFaceI]
                    << " with vertices "
                    << UIndirectList<point>(pp.points(), pp[oldFaceI])()
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
    os.writeKeyword("patchIDs") << patchIDs_
        << token::END_STATEMENT << nl;
    os.writeKeyword("starts") << starts_
        << token::END_STATEMENT << nl;
}


// ************************************************************************* //

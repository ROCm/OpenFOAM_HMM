/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd.
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

#include "wallBoundedParticle.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const std::size_t Foam::wallBoundedParticle::sizeofFields_
(
    sizeof(wallBoundedParticle) - sizeof(particle)
);


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::tetIndices Foam::wallBoundedParticle::currentTetIndices() const
{
    // Replacement for particle::currentTetIndices that avoids error
    // upon invalid tetBasePtIs

    const faceList& pFaces = mesh_.faces();
    const labelList& pOwner = mesh_.faceOwner();

    const Foam::face& f = pFaces[tetFacei_];

    bool own = (pOwner[tetFacei_] == celli_);

    label faceBasePtI = mesh_.tetBasePtIs()[tetFacei_];
    if (faceBasePtI == -1)
    {
        //WarningInFunction
        //    << "No base point for face " << tetFacei_ << ", " << f
        //    << ", produces a decomposition that has a minimum "
        //    << "volume greater than tolerance."
        //    << endl;
        faceBasePtI = 0;
    }

    label facePtI = (tetPti_ + faceBasePtI) % f.size();
    label otherFacePtI = f.fcIndex(facePtI);

    label facePtAI;
    label facePtBI;
    if (own)
    {
        facePtAI = facePtI;
        facePtBI = otherFacePtI;
    }
    else
    {
        facePtAI = otherFacePtI;
        facePtBI = facePtI;
    }

    return tetIndices
    (
        celli_,
        tetFacei_,
        faceBasePtI,
        facePtAI,
        facePtBI,
        tetPti_
    );
}


Foam::edge Foam::wallBoundedParticle::currentEdge() const
{
    if ((meshEdgeStart_ != -1) == (diagEdge_ != -1))
    {
        FatalErrorInFunction
            << "Particle:"
            << info()
            << "cannot both be on a mesh edge and a face-diagonal edge."
            << " meshEdgeStart_:" << meshEdgeStart_
            << " diagEdge_:" << diagEdge_
            << abort(FatalError);
    }

    const Foam::face& f = mesh_.faces()[tetFace()];

    if (meshEdgeStart_ != -1)
    {
        return edge(f[meshEdgeStart_], f.nextLabel(meshEdgeStart_));
    }
    else
    {
        label faceBasePti = mesh_.tetBasePtIs()[tetFace()];
        if (faceBasePti == -1)
        {
            //FatalErrorInFunction
            //WarningInFunction
            //    << "No base point for face " << tetFace() << ", " << f
            //    << ", produces a decomposition that has a minimum "
            //    << "volume greater than tolerance."
            //    //<< abort(FatalError);
            //    << endl;
            faceBasePti = 0;
        }

        label diagPti = (faceBasePti+diagEdge_)%f.size();

        return edge(f[faceBasePti], f[diagPti]);
    }
}


void Foam::wallBoundedParticle::crossEdgeConnectedFace
(
    const label& celli,
    label& tetFacei,
    label& tetPti,
    const edge& e
)
{
    const faceList& pFaces = mesh_.faces();
    const cellList& pCells = mesh_.cells();

    const Foam::face& f = pFaces[tetFacei];

    const Foam::cell& thisCell = pCells[celli];

    forAll(thisCell, cFI)
    {
        // Loop over all other faces of this cell and
        // find the one that shares this edge

        label fI = thisCell[cFI];

        if (tetFacei == fI)
        {
            continue;
        }

        const Foam::face& otherFace = pFaces[fI];

        label edDir = otherFace.edgeDirection(e);

        if (edDir == 0)
        {
            continue;
        }
        else if (f == pFaces[fI])
        {
            // This is a necessary condition if using duplicate baffles
            // (so coincident faces). We need to make sure we don't cross into
            // the face with the same vertices since we might enter a tracking
            // loop where it never exits. This test should be cheap
            // for most meshes so can be left in for 'normal' meshes.
            continue;
        }
        else
        {
            //Found edge on other face
            tetFacei = fI;

            label eIndex = -1;

            if (edDir == 1)
            {
                // Edge is in the forward circulation of this face, so
                // work with the start point of the edge
                eIndex = findIndex(otherFace, e.start());
            }
            else
            {
                // edDir == -1, so the edge is in the reverse
                // circulation of this face, so work with the end
                // point of the edge
                eIndex = findIndex(otherFace, e.end());
            }

            label tetBasePtI = mesh_.tetBasePtIs()[fI];

            if (tetBasePtI == -1)
            {
                //FatalErrorInFunction
                //WarningInFunction
                //    << "No base point for face " << fI << ", " << f
                //    << ", produces a decomposition that has a minimum "
                //    << "volume greater than tolerance."
                //    //<< abort(FatalError);
                //    << endl;
                tetBasePtI = 0;
            }

            // Find eIndex relative to the base point on new face
            eIndex -= tetBasePtI;

            if (neg(eIndex))
            {
                eIndex = (eIndex + otherFace.size()) % otherFace.size();
            }

            if (eIndex == 0)
            {
                // The point is the base point, so this is first tet
                // in the face circulation
                tetPti = 1;
            }
            else if (eIndex == otherFace.size() - 1)
            {
                // The point is the last before the base point, so
                // this is the last tet in the face circulation
                tetPti = otherFace.size() - 2;
            }
            else
            {
                tetPti = eIndex;
            }

            break;
        }
    }
}


void Foam::wallBoundedParticle::crossEdgeConnectedFace(const edge& meshEdge)
{
    // Update tetFace, tetPt
    crossEdgeConnectedFace(cell(), tetFace(), tetPt(), meshEdge);

    // Update face to be same as tracking one
    face() = tetFace();

    // And adapt meshEdgeStart_.
    const Foam::face& f = mesh_.faces()[tetFace()];
    label fp = findIndex(f, meshEdge[0]);

    if (f.nextLabel(fp) == meshEdge[1])
    {
        meshEdgeStart_ = fp;
    }
    else
    {
        label fpMin1 = f.rcIndex(fp);

        if (f[fpMin1] == meshEdge[1])
        {
            meshEdgeStart_ = fpMin1;
        }
        else
        {
            FatalErrorInFunction
                << "Problem :"
                << " particle:"
                << info()
                << "face:" << tetFace()
                << " verts:" << f
                << " meshEdge:" << meshEdge
                << abort(FatalError);
        }
    }

    diagEdge_ = -1;

    // Check that still on same mesh edge
    const edge eNew(f[meshEdgeStart_], f.nextLabel(meshEdgeStart_));
    if (eNew != meshEdge)
    {
        FatalErrorInFunction
            << "Problem" << abort(FatalError);
    }
}


void Foam::wallBoundedParticle::crossDiagonalEdge()
{
    if (diagEdge_ == -1)
    {
        FatalErrorInFunction
            << "Particle:"
            << info()
            << "not on a diagonal edge" << abort(FatalError);
    }
    if (meshEdgeStart_ != -1)
    {
        FatalErrorInFunction
            << "Particle:"
            << info()
            << "meshEdgeStart_:" << meshEdgeStart_ << abort(FatalError);
    }

    const Foam::face& f = mesh_.faces()[tetFace()];

    // tetPti starts from 1, goes up to f.size()-2

    if (tetPt() == diagEdge_)
    {
        tetPt() = f.rcIndex(tetPt());
    }
    else
    {
        label nextTetPt = f.fcIndex(tetPt());
        if (diagEdge_ == nextTetPt)
        {
            tetPt() = nextTetPt;
        }
        else
        {
            FatalErrorInFunction
                << "Particle:"
                << info()
                << "tetPt:" << tetPt()
                << " diagEdge:" << diagEdge_ << abort(FatalError);
        }
    }

    meshEdgeStart_ = -1;
}


Foam::scalar Foam::wallBoundedParticle::trackFaceTri
(
    const vector& n,
    const vector& endPosition,
    label& minEdgei
)
{
    // Track p from position to endPosition
    const triFace tri(currentTetIndices().faceTriIs(mesh_));

    // Check which edge intersects the trajectory.
    // Project trajectory onto triangle
    minEdgei = -1;
    scalar minS = 1;        // end position

    edge currentE(-1, -1);
    if (meshEdgeStart_ != -1 || diagEdge_ != -1)
    {
        currentE = currentEdge();
    }

    // Determine path along line position+s*d to see where intersections are.
    forAll(tri, i)
    {
        label j = tri.fcIndex(i);

        const point& pt0 = mesh_.points()[tri[i]];
        const point& pt1 = mesh_.points()[tri[j]];

        if (edge(tri[i], tri[j]) == currentE)
        {
            // Do not check particle is on
            continue;
        }

        // Outwards pointing normal
        vector edgeNormal = (pt1-pt0)^n;

        edgeNormal /= mag(edgeNormal)+VSMALL;

        // Determine whether position and end point on either side of edge.
        scalar sEnd = (endPosition-pt0)&edgeNormal;
        if (sEnd >= 0)
        {
            // endPos is outside triangle. position() should always be
            // inside.
            scalar sStart = (position()-pt0)&edgeNormal;
            if (mag(sEnd-sStart) > VSMALL)
            {
                scalar s = sStart/(sStart-sEnd);

                if (s >= 0 && s < minS)
                {
                    minS = s;
                    minEdgei = i;
                }
            }
        }
    }

    if (minEdgei != -1)
    {
        position() += minS*(endPosition-position());
    }
    else
    {
        // Did not hit any edge so tracked to the end position. Set position
        // without any calculation to avoid truncation errors.
        position() = endPosition;
        minS = 1.0;
    }

    // Project position() back onto plane of triangle
    const point& triPt = mesh_.points()[tri[0]];
    position() -= ((position()-triPt)&n)*n;

    return minS;
}


bool Foam::wallBoundedParticle::isTriAlongTrack
(
    const vector& n,
    const point& endPosition
) const
{
    const triFace triVerts(currentTetIndices().faceTriIs(mesh_));
    const edge currentE = currentEdge();

    if
    (
        currentE[0] == currentE[1]
     || findIndex(triVerts, currentE[0]) == -1
     || findIndex(triVerts, currentE[1]) == -1
    )
    {
        FatalErrorInFunction
            << "Edge " << currentE << " not on triangle " << triVerts
            << info()
            << abort(FatalError);
    }


    const vector dir = endPosition-position();

    forAll(triVerts, i)
    {
        label j = triVerts.fcIndex(i);
        const point& pt0 = mesh_.points()[triVerts[i]];
        const point& pt1 = mesh_.points()[triVerts[j]];

        if (edge(triVerts[i], triVerts[j]) == currentE)
        {
            vector edgeNormal = (pt1-pt0)^n;
            return (dir&edgeNormal) < 0;
        }
    }

    FatalErrorInFunction
        << "Problem" << abort(FatalError);

    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallBoundedParticle::wallBoundedParticle
(
    const polyMesh& mesh,
    const vector& position,
    const label celli,
    const label tetFacei,
    const label tetPti,
    const label meshEdgeStart,
    const label diagEdge
)
:
    particle(mesh, position, celli, tetFacei, tetPti),
    meshEdgeStart_(meshEdgeStart),
    diagEdge_(diagEdge)
{}


Foam::wallBoundedParticle::wallBoundedParticle
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields,
    bool newFormat
)
:
    particle(mesh, is, readFields, newFormat)
{
    if (readFields)
    {
        if (is.format() == IOstream::ASCII)
        {
            is  >> meshEdgeStart_ >> diagEdge_;
        }
        else
        {
            is.read(reinterpret_cast<char*>(&meshEdgeStart_), sizeofFields_);
        }
    }

    is.check(FUNCTION_NAME);
}


Foam::wallBoundedParticle::wallBoundedParticle
(
    const wallBoundedParticle& p
)
:
    particle(p),
    meshEdgeStart_(p.meshEdgeStart_),
    diagEdge_(p.diagEdge_)
{}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const wallBoundedParticle& p
)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const particle&>(p)
            << token::SPACE << p.meshEdgeStart_
            << token::SPACE << p.diagEdge_;
    }
    else
    {
        os  << static_cast<const particle&>(p);
        os.write
        (
            reinterpret_cast<const char*>(&p.meshEdgeStart_),
            wallBoundedParticle::sizeofFields_
        );
    }

    return os;
}


Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const InfoProxy<wallBoundedParticle>& ip
)
{
    const wallBoundedParticle& p = ip.t_;

    tetPointRef tpr(p.currentTetIndices().tet(p.mesh()));

    os  << "    " << static_cast<const particle&>(p) << nl
        << "    tet:" << nl;
    os  << "    ";
    meshTools::writeOBJ(os, tpr.a());
    os  << "    ";
    meshTools::writeOBJ(os, tpr.b());
    os  << "    ";
    meshTools::writeOBJ(os, tpr.c());
    os  << "    ";
    meshTools::writeOBJ(os, tpr.d());
    os  << "    l 1 2" << nl
        << "    l 1 3" << nl
        << "    l 1 4" << nl
        << "    l 2 3" << nl
        << "    l 2 4" << nl
        << "    l 3 4" << nl;
    os  << "    ";
    meshTools::writeOBJ(os, p.position());

    return os;
}


// ************************************************************************* //

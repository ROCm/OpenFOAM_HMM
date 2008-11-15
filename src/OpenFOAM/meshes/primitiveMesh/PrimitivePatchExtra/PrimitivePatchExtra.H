/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
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

Class
    Foam::PrimitivePatchExtra

Description
    PrimitivePatch with some extra functionality.

SourceFiles
    PrimitivePatchExtra.C

\*---------------------------------------------------------------------------*/

#ifndef PrimitivePatchExtra_H
#define PrimitivePatchExtra_H

#include "PrimitivePatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class PrimitivePatchExtra Declaration
\*---------------------------------------------------------------------------*/

template
<
    class Face,
    template<class> class FaceList,
    class PointField,
    class PointType=point
>
class PrimitivePatchExtra
:
    public PrimitivePatch<Face, FaceList, PointField, PointType>
{

public:

    // Public typedefs

        typedef Face FaceType;
        typedef FaceList<Face> FaceListType;
        typedef PointField PointFieldType;

private:

    // Private typedefs
    typedef PrimitivePatch<Face, FaceList, PointField, PointType> TemplateType;

    // Private data

    // Demand driven private data

        //- Edge-face addressing (sorted)
        mutable labelListList* sortedEdgeFacesPtr_;

        //- Label of face that 'owns' edge
        //  i.e. e.vec() is righthanded walk along face
        mutable labelList* edgeOwnerPtr_;


    // Private Member Functions

        //- Calculate sorted edgeFaces
        void calcSortedEdgeFaces() const;

        //- Calculate owner
        void calcEdgeOwner() const;

public:

    // Constructors

        //- Construct from components
        PrimitivePatchExtra
        (
            const FaceList<Face>& faces,
            const Field<PointType>& points
        );

        //- Construct as copy
        PrimitivePatchExtra
        (
            const PrimitivePatchExtra<Face, FaceList, PointField, PointType>&
        );


    // Destructor

        virtual ~PrimitivePatchExtra();

        void clearOut();

        void clearTopology();


    // Member Functions

    // Access functions for demand driven data

    // Topological data; no mesh required.

        //- Return edge-face addressing sorted
        //  (for edges with more than 2 faces) according to the
        //  angle around the edge.
        //  Orientation is anticlockwise looking from
        //  edge.vec(localPoints())
        const labelListList& sortedEdgeFaces() const;

        //- If 2 face neighbours: label of face where ordering of edge
        //  is consistent with righthand walk.
        //  If 1 neighbour: label of only face.
        //  If >2 neighbours: undetermined.
        const labelList& edgeOwner() const;


    // Addressing into mesh


    // Other patch operations


    // Check

        //- Check multiply-connected edges.
        void checkEdges(const bool verbose) const;

        //- Check orientation (normals) and normals of neighbouring faces
        boolList checkOrientation(const bool verbose) const;


    // Edit

        //- Fill faceZone with currentZone for every face reachable
        //  from faceI without crossing edge marked in borderEdge.
        //  Note: faceZone has to be sized nFaces before calling.
        void markZone
        (
            const UList<bool>& borderEdge,
            const label faceI,
            const label currentZone,
            labelList& faceZone
        ) const;

        //- (size and) fills faceZone with zone of face.
        //  Zone is area reachable by edge crossing without crossing borderEdge
        //  (bool for every edge in surface). Returns number of zones.
        label markZones
        (
            const UList<bool>& borderEdge,
            labelList& faceZone
        ) const;

        //- Determine the mapping for a sub-mesh.
        //  Only include faces for which bool List entry is true
        //  @param[out] pointMap mapping new to old localPoints
        //  @param[out] faceMap  mapping new to old faces
        void subsetMap
        (
            const UList<bool>& include,
            labelList& pointMap,
            labelList& faceMap
        ) const;


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "PrimitivePatchExtra.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    Foam::localPointRegion

Description
    Takes mesh with 'baffles' (= boundary faces sharing points).
    Determines for selected points on boundary faces the 'point region' it is
    connected to. Each region can be visited by a cell-face-cell walk.
    Used in duplicating points after splitting baffles.

    Regions are not consecutive per processor. They will be -1..nRegions_.

    Note: coupled boundaries (cyclics, parallel) not fully tested.

SourceFiles
    localPointRegion.C

\*---------------------------------------------------------------------------*/

#ifndef localPointRegion_H
#define localPointRegion_H

#include "typeInfo.H"
#include "Map.H"
#include "labelList.H"
#include "labelHashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class primitiveMesh;
class polyMesh;
class face;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                           Class localPointRegion Declaration
\*---------------------------------------------------------------------------*/

class localPointRegion
{
    // Private data

        //- (global) number of regions
        label nRegions_;

        //- From mesh face to local face
        Map<label> localFaces_;

        //- From local face addressing to mesh faces.
        labelList meshFaces_;

        //- Per local face the region
        labelList faceRegion_;

        //- Per local point (i.e. index in pointsToBeDuplicated) the list
        //  of regions
        labelListList pointRegions_;


    // Private Member Functions

        //- Set face (and its connected points) to global region nRegions_.
        void setFaceRegion
        (
            const polyMesh&,
            const Map<label>&,
            const label faceI,
            labelListList& pointFaceRegion // in pointFaces addressing.
        );

        //- Set point (and its connected faces) to global region nRegions_
        void setPointRegion
        (
            const polyMesh&,
            const Map<label>&,
            const label pointI,
            const label localPointI,
            const label faceI,
            labelListList& pointFaceRegion
        );

        //- Gets local face regions (local to each point) and merges them into
        //  global face region (always possible since faces cannot cross
        //  a baffle.
        void mergeLocalRegions
        (
            const polyMesh&,
            const labelList& pointsToBeDuplicated,
            labelListList& pointFaceRegion
        );

        //- Change region
        void changeRegion
        (
            const polyMesh&,
            const label faceI,
            const label oldRegionI,
            const label newRegionI,
            labelHashSet& changedFaces
        );

        //- Check if two faces are equal. If forward = false checks f1 in
        //  reverse order.
        static bool isDuplicate
        (
            const face& f0,
            const face& f1,
            const bool forward
        );

public:

    //- Runtime type information
    ClassName("localPointRegion");


    // Constructors

        //- Construct from mesh and points that should be separated.
        localPointRegion
        (
            const polyMesh& mesh,
            const labelList& pointsToBeDuplicated
        );


    // Member Functions

        // Static Member Functions

            //- Helper routine to find baffles (two boundary faces using the
            //  same points but in reverse order)
            //  Gets list of (boundary!) faces to check. Returns labelList
            //  of same size as the input list
            //  with -1 or index of other face in the input list.
            //  Does not handle duplicate faces on both sides of processor patch
            static labelList findDuplicateFaces
            (
                const primitiveMesh&,
                const labelList&
            );

            //- Are all the cells using pointI reachable through (processor
            //  local) cell-face-cell walk.
            static bool isSingleCellRegion
            (
                const primitiveMesh& mesh,
                const label pointI
            );

            //- See walkPointConnectedFaces. Temporarily exposed so as to
            //  use outside since par problem with localPointRegion (which is
            //  probably fixed now with plusEqOp v.s. sumOp in gather/scatter)
            static label walkCellFaceCell
            (
                const primitiveMesh&,
                const label startCellI,
                const label startPointI,
                const label regionI,
                labelList& regionPerFace
            );


        // Access

            //- Global number of regions
            label nRegions() const
            {
                return nRegions_;
            }

            //- From local face addressing to mesh faces.
            const labelList& meshFaces() const
            {   
                return meshFaces_;
            }

            //- From mesh face to local face
            const Map<label>& localFaces() const
            {   
                return localFaces_;
            }

            //- Per local face the global region. Regions are not consecutive
            //  per processor. They will be -1..nRegions_.
            const labelList& faceRegion() const
            {
                return faceRegion_;
            }

            //- Per pointsToBeDuplicated the list of regions
            const labelListList& pointRegions() const
            {
                return pointRegions_;
            }

        // Edit

            //- Force recalculation of locally stored data on topological change
            void updateMesh(const mapPolyMesh&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

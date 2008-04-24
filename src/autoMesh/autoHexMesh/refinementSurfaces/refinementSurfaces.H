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

Class
    refinementSurfaces

Description
    Container for triSurfaces used for surface-driven refinement.
    These contain all the data about the level of refinement needed per
    surface.

SourceFiles
    refinementSurfaces.C

\*---------------------------------------------------------------------------*/

#ifndef refinementSurfaces_H
#define refinementSurfaces_H

#include "triSurfaceMeshes.H"
#include "PackedList.H"
#include "triSurfaceGeoMesh.H"
#include "triSurfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class searchableSurface;

/*---------------------------------------------------------------------------*\
                           Class refinementSurfaces Declaration
\*---------------------------------------------------------------------------*/

class refinementSurfaces
:
    public triSurfaceMeshes
{
    // Private data

        //- Name (word)
        wordList names_;

        //- Per surface whether is closed
        boolList closed_;

        //- Per 'interface' surface : name of faceZone to put faces into
        wordList faceZoneNames_;

        //- Per 'interface' surface : name of cellZone to put cells into
        wordList cellZoneNames_;

        //- Per 'interface' surface : (only used if surface is closed)
        //  whether to zone cells inside or outside surface.
        boolList zoneInside_;


        //- From local region number to global region number
        labelList regionOffset_;

        //- From global region number to refinement level
        labelList minLevel_;

        //- From global region number to refinement level
        labelList maxLevel_;

        //- Per surface refinement level adapted for shells.
        PtrList<triSurfaceLabelField> minLevelFields_;


    // Private Member Functions

        static fileNameList extractFileNames(const PtrList<dictionary>&);

        //- Disallow default bitwise copy construct
        refinementSurfaces(const refinementSurfaces&);

        //- Disallow default bitwise assignment
        void operator=(const refinementSurfaces&);


public:

    // Constructors

        //- Construct from directory (io.instance()) and dictionaries
        refinementSurfaces
        (
            const IOobject& io,
            const PtrList<dictionary>&
        );


    // Member Functions

        // Access

            //- Names of surfaces
            const wordList& names() const
            {
                return names_;
            }

            //- Per surface whether is closed
            const boolList& closed() const
            {
                return closed_;
            }

            //- Per 'interface' surface : name of faceZone to put faces into
            const wordList& faceZoneNames() const
            {
                return faceZoneNames_;
            }

            //- Per 'interface' surface : name of cellZone to put cells into
            const wordList& cellZoneNames() const
            {
                return cellZoneNames_;
            }

            //- Per 'interface' surface : if closed: zonify cells inside surface
            const boolList& zoneInside() const
            {
                return zoneInside_;
            }

            //- Get indices of named surfaces (surfaces with cellZoneName)
            labelList getNamedSurfaces() const;


            //- From local region number to global region number
            const labelList& regionOffset() const
            {
                return regionOffset_;
            }

            //- From global region number to refinement level
            const labelList& minLevel() const
            {
                return minLevel_;
            }

            //- From global region number to refinement level
            const labelList& maxLevel() const
            {
                return maxLevel_;
            }


        // Helper

            //- From surface and region on surface to global region
            label globalRegion(const label surfI, const label regionI) const
            {
                return regionOffset_[surfI]+regionI;
            }

            //- From triangle on surface to global region
            label triangleRegion(const label surfI, const label triI) const
            {
                const triSurface& s = operator[](surfI);

                return globalRegion(surfI, s[triI].region());
            }

            //- Min level for surface and region on surface
            label minLevel(const label surfI, const label regionI) const
            {
                return minLevel_[globalRegion(surfI, regionI)];
            }

            //- Max level for surface and region on surface
            label maxLevel(const label surfI, const label regionI) const
            {
                return maxLevel_[globalRegion(surfI, regionI)];
            }

            label nRegions() const
            {
                return minLevel_.size();
            }

            //- Minlevel updated for refinement shells
            const triSurfaceLabelField& minLevelField(const label surfI) const
            {
                return minLevelFields_[surfI];
            }

            //- Calculate minLevelFields
            void setMinLevelFields
            (
                const PtrList<searchableSurface>& shells,
                const labelList& shellLevels,
                const boolList& shellRefineInside
            );

            //- Helper: is surface closed?
            static bool isSurfaceClosed(const triSurface&);

            //- Helper: orient (closed only) surfaces so keepPoint is outside.
            static void orientSurface
            (
                const point& keepPoint,
                triSurfaceMesh& surface
            );

            //- Helper: count number of triangles per region
            static labelList countRegions(const triSurface&);


        // Searching

            //- Find intersection of edge. Return -1 or first surface
            //  with higher (than currentLevel) minlevel.
            //  Return surface number and hit info.
            label findHigherIntersection
            (
                const point& start,
                const point& end,
                const label currentLevel,   // current cell refinement level
                pointIndexHit&
            ) const;

            //- Find intersection with max level. Return -1 or the surface
            //  with the highest maxLevel above currentLevel
            label findHighestIntersection
            (
                const point& start,
                const point& end,
                const label currentLevel,   // current cell refinement level
                pointIndexHit&
            ) const;

            //- Detect if a point is 'inside' (depending on zoneInside flag) a
            //  zoneable surface. Returns -1 if not, returns first surface it
            //  is.
            label insideZone(const labelList& surfaces, const point& pt) const;


            //- Mark for all points whether it is inside any closed surface
            //  Return number of inside points.
            label markInsidePoints(const pointField&, PackedList<1>& isInside)
             const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

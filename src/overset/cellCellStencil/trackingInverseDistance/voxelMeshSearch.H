/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
     \\/     M anipulation  |
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

Class
    Foam::voxelMeshSearch

Description
    Fast, non-parallel searching in mesh without use of octree.

SourceFiles
    voxelMeshSearch.C

\*---------------------------------------------------------------------------*/

#ifndef voxelMeshSearch_H
#define voxelMeshSearch_H

#include "boundBox.H"
#include "labelVector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
class polyMesh;
class OBJstream;

/*---------------------------------------------------------------------------*\
                        Class voxelMeshSearch Declaration
\*---------------------------------------------------------------------------*/

class voxelMeshSearch
{
    // Private data

        const polyMesh& mesh_;

        //- Local mesh bounding box
        boundBox localBb_;

        //- Number of points in all directions (for local mesh only)
        labelVector nDivs_;

        //- Voxel to seed cell
        labelList seedCell_;


   // Private Member Functions

        //- Find nearest cell (on same processor patch as seed face)
        label searchProcPatch(const label seedFacei, const point&) const;

        //- Find the face on the cell that gets intersected
        label findIntersectedFace(const label celli, const point&) const;


public:

    // Declare name of the class and its debug switch
    ClassName("voxelMeshSearch");


    // Constructors

        //- Construct from mesh; voxels estimated from local number of cells
        voxelMeshSearch(const polyMesh&, const bool doUpdate = true);

        //- Construct from mesh and voxel discretisation
        voxelMeshSearch
        (
            const polyMesh&,
            const labelVector&,
            const bool doUpdate = true
        );


    // Member functions

        //- Number of voxels for local mesh
        const labelVector& nDivs() const
        {
            return nDivs_;
        }

        //- Update lookup tables for geometry changes
        bool update();

        //- Find a cell
        label findCell(const point&) const;

        //- Find cells. Returns number of cells found
        //label findCells(const UList<point>&, labelList&) const;


        //Voxel helper functions

            //- Coordinate to voxel indices
            static labelVector index3
            (
                const boundBox& bb,
                const labelVector& nDivs,
                const point& p
            );

            //- Coordinate to combined voxel index. If clip makes sure
            //  components are all inside. If not clip returns -1 if outside bb.
            static label index
            (
                const boundBox& bb,
                const labelVector& nDivs,
                const point& p,
                const bool clip
            );

            //- Voxel index to voxel centre
            static point centre
            (
                const boundBox& bb,
                const labelVector& nDivs,
                const labelVector& voxel
            );

            //- Fill voxels indicated by bounding box
            template<class Container, class Type>
            static void fill
            (
                Container& elems,
                const boundBox& bb,
                const labelVector& nDivs,
                const boundBox& subBb,
                const Type val
            );

            //- Check if any voxel inside bounding box is set to val or
            //  not set to val (isNot = true)
            template<class Container, class Type>
            static bool overlaps
            (
                const boundBox& bb,
                const labelVector& nDivs,
                const boundBox& subBb,
                const Container& elems,
                const Type val,
                const bool isNot = false
            );

            //- Debug: write points for every set element
            template<class Container, class Type>
            static void write
            (
                OBJstream&,
                const boundBox& bb,
                const labelVector& nDivs,
                const Container& elems,
                const Type val,
                const bool isNot = false
            );

            //- Debug: write all edges
            static void writeGrid
            (
                OBJstream&,
                const boundBox&,
                const labelVector&
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "voxelMeshSearchTemplates.C"
#endif

#endif

// ************************************************************************* //

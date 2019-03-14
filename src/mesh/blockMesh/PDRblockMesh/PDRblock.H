/*---------------------------------------------------------------------------* \
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::PDRblock

Description
    A single block x-y-z rectilinear mesh addressable as i,j,k with
    simplified creation. Some of the input is similar to blockMeshDict,
    but since this specialization is for a single-block that is aligned
    with the x-y-z directions, it provides a different means of specifying
    the mesh.

    Dictionary controls
    \table
        Property    | Description                          | Required | Default
        x           | X-direction grid specification       | yes |
        y           | Y-direction grid specification       | yes |
        z           | Z-direction grid specification       | yes |
        scale       | Point scaling                        | no  | 1.0
        boundary    | Boundary patches                     | yes |
        defaultPatch | Default patch specification         | no  |
    \endtable

    Grid coordinate controls
    \table
        Property| Description                          | Required | Default
        points  | Locations defining the mesh segment  | yes |
        nCells  | Divisions per mesh segment           | yes |
        ratios  | Expansion ratio (end/start) per segment | no | uniform
    \endtable

    The expansion ratios are defined as per blockMesh and represent the ratio
    of end-size / start-size for the section. A negative value is trapped
    and treated as its reciprocal.

SourceFiles
    PDRblockI.H
    PDRblock.C
    PDRblockCreate.C

\*---------------------------------------------------------------------------*/

#ifndef PDRblock_H
#define PDRblock_H

#include "ijkMesh.H"
#include "boundBox.H"
#include "pointField.H"
#include "faceList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class IOobject;
class polyMesh;

/*---------------------------------------------------------------------------*\
                           Class PDRblock Declaration
\*---------------------------------------------------------------------------*/

class PDRblock
:
    public ijkMesh
{

public:

    // Public Classes

        //- Grid locations in an axis direction.
        //  The number of points is one larger than the number of elements
        //  it represents
        struct location
        :
            public scalarList
        {
            //- The number of cells in this direction.
            inline label nCells() const;

            //- The number of points in this direction.
            inline label nPoints() const;

            //- The number of divisions (nCells) in this direction.
            inline label size() const;

            //- Check that element index is within valid range.
            inline void checkIndex(const label i) const;

            //- Cell size at element position.
            inline scalar width(const label i) const;

            //- Cell centre at element position.
            inline scalar C(const label i) const;
        };


private:

    // Private Classes

        //- Extracted patch settings
        struct boundaryEntry
        {
            //- The patch name
            word name_;

            //- The patch type
            word type_;

            //- The patch size
            label size_;

            //- The associated block face ids [0,5]
            labelList faces_;
        };


    // Private Data

        //- Verbosity
        bool verbose_;

        //- The grid points in all (i,j,k / x,y,z) directions.
        Vector<location> grid_;

        //- The mesh bounding box
        boundBox bounds_;

        //- The boundary patch information
        PtrList<boundaryEntry> patches_;


    // Private Member Functions

        //- Read and define grid points in given direction
        void readGridControl
        (
            const direction cmpt,
            const dictionary& dict,
            const scalar scaleFactor = -1
        );

        //- Read "boundary" information
        void readBoundary(const dictionary& dict);


        //- Populate point field for the block
        void createPoints(pointField& pts) const;

        //- Add internal faces to lists.
        //  Lists must be properly sized!
        //  \return the number of faces added
        label addInternalFaces
        (
            faceList::iterator& faceIter,
            labelList::iterator& ownIter,
            labelList::iterator& neiIter
        ) const;


        //- Add boundary faces for the shape face to lists
        //  Lists must be properly sized!
        //  \return the number of faces added
        label addBoundaryFaces
        (
            const direction shapeFacei,
            faceList::iterator& faceIter,
            labelList::iterator& ownIter
        ) const;


public:

    // Constructors

        //- Construct zero-size
        inline PDRblock();

        //- Construct from dictionary
        explicit PDRblock(const dictionary& dict, bool verbose=false);


    // Member Functions

        //- Read dictionary
        bool read(const dictionary& dict);


    // Access

        //- The grid point locations in the i,j,k (x,y,z) directions.
        inline const Vector<location>& grid() const;


    // Mesh information

        //- The mesh bounding box
        inline const boundBox& bounds() const;

        //- Cell size in x-direction at i position.
        inline scalar dx(const label i) const;

        //- Cell size in x-direction at i position.
        inline scalar dx(const labelVector& ijk) const;

        //- Cell size in y-direction at j position.
        inline scalar dy(const label j) const;

        //- Cell size in y-direction at j position.
        inline scalar dy(const labelVector& ijk) const;

        //- Cell size in z-direction at k position.
        inline scalar dz(const label k) const;

        //- Cell size in z-direction at k position.
        inline scalar dz(const labelVector& ijk) const;

        //- Cell dimensions at i,j,k position.
        inline vector span(const label i, const label j, const label k) const;

        //- Cell dimensions at i,j,k position.
        inline vector span(const labelVector& ijk) const;

        //- Cell centre at i,j,k position.
        inline point C(const label i, const label j, const label k) const;

        //- Cell centre at i,j,k position.
        inline point C(const labelVector& ijk) const;

        //- Cell volume at i,j,k position.
        inline scalar V(const label i, const label j, const label k) const;

        //- Cell volume at i,j,k position.
        inline scalar V(const labelVector& ijk) const;

        //- Characteristic cell size at i,j,k position.
        //  This is the cubic root of the volume
        inline scalar width(const label i, const label j, const label k) const;

        //- Characteristic cell size at i,j,k position.
        //  This is the cubic root of the volume
        inline scalar width(const labelVector& ijk) const;


    // Searching

        //- Return i,j,k index for cell enclosing this location
        //  The value (-1,-1,-1) is returned for out-of-bounds (not found).
        labelVector findCell(const point& pt) const;


    // Mesh generation

        //- Create polyMesh for grid definition and patch information
        autoPtr<polyMesh> mesh(const IOobject& io) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "PDRblockI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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
    Foam::cellSizeControlSurfaces

Description

SourceFiles
    cellSizeControlSurfacesI.H
    cellSizeControlSurfaces.C

\*---------------------------------------------------------------------------*/

#ifndef cellSizeControlSurfaces_H
#define cellSizeControlSurfaces_H

#include "searchableSurfaces.H"
#include "searchableSurfacesQueries.H"
#include "conformationSurfaces.H"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template < class GT, class Vb = CGAL::Triangulation_vertex_base_3<GT> >
class vertexWithInfo
:
    public Vb
{
public:

    Foam::scalar value_;

    Foam::tensor alignment_;

    typedef typename Vb::Vertex_handle  Vertex_handle;
    typedef typename Vb::Cell_handle    Cell_handle;
    typedef typename Vb::Point          Point;

    template < class TDS2 >
    struct Rebind_TDS
    {
        typedef typename Vb::template Rebind_TDS<TDS2>::Other   Vb2;
        typedef vertexWithInfo<GT, Vb2>                         Other;
    };

    vertexWithInfo()
    :
        Vb(),
        value_(0.0),
        alignment_(Foam::tensor::zero)
    {}

    explicit vertexWithInfo(const Point& p)
    :
        Vb(p),
        value_(0.0),
        alignment_(Foam::tensor::zero)
    {}

    vertexWithInfo(const Point& p, const Foam::scalar& value)
    :
        Vb(p),
        value_(value),
        alignment_(Foam::tensor::zero)
    {}

    vertexWithInfo
    (
        const Point& p,
        const Foam::scalar& value,
        const Foam::tensor& alignment
    )
    :
        Vb(p),
        value_(value),
        alignment_(alignment)
    {}


    const Foam::scalar& value() const
    {
        return value_;
    }

    void value(const Foam::scalar& value)
    {
        value_ = value;
    }

    const Foam::tensor& alignment() const
    {
        return alignment_;
    }

    void alignment(const Foam::tensor& alignment)
    {
        alignment_ = alignment;
    }
};


namespace Foam
{

// Forward declaration of classes
class conformalVoronoiMesh;

class cellSizeFunction;


/*---------------------------------------------------------------------------*\
                   Class cellSizeControlSurfaces Declaration
\*---------------------------------------------------------------------------*/

class cellSizeControlSurfaces
{
    typedef CGAL::Exact_predicates_inexact_constructions_kernel             K;
    typedef CGAL::Triangulation_data_structure_3<vertexWithInfo<K> >        Tds;

    typedef CGAL::Delaunay_triangulation_3<K, Tds, CGAL::Fast_location>
        Delaunay;

    typedef Delaunay::Cell_handle    Cell_handle;
    typedef Delaunay::Vertex_handle  Vertex_handle;
    typedef Delaunay::Locate_type    Locate_type;
    typedef Delaunay::Point          Point;

    Delaunay T_;

    mutable Cell_handle oldCellHandle_;


    // Private data

        //- Reference to the searchableSurfaces object holding all geometry data
        const searchableSurfaces& allGeometry_;

        const conformationSurfaces& geometryToConformTo_;

        //- Indices of surfaces in allGeometry that are to be conformed to
        labelList surfaces_;

        //- A list of all of the cellSizeFunction objects
        PtrList<cellSizeFunction> cellSizeFunctions_;

        //- Regions where no cell size function is specified use defaultCellSize
        scalar defaultCellSize_;

        //- Assigning a priority to all requests for cell sizes, the highest
        //  overrules
        label defaultPriority_;


    // Private Member Functions

        //- Evaluate the cell size functions, returning a bool stating if a
        //  function was found or not
        bool evalCellSizeFunctions
        (
            const point& pt,
            scalar& minSize
        ) const;

        bool checkCoplanarTet(Cell_handle c, const scalar tol) const;

        bool checkClosePoints(Cell_handle c, const scalar tol) const;

        label refineTriangulation(const scalar factor);

        void writeRefinementTriangulation();

        //- Disallow default bitwise copy construct
        cellSizeControlSurfaces(const cellSizeControlSurfaces&);

        //- Disallow default bitwise assignment
        void operator=(const cellSizeControlSurfaces&);


public:

    //- Runtime type information
    ClassName("cellSizeControlSurfaces");


    // Constructors

        //- Construct from dictionary and references to conformalVoronoiMesh and
        //  searchableSurfaces
        cellSizeControlSurfaces
        (
            const searchableSurfaces& allGeometry,
            const conformationSurfaces& geometryToConformTo,
            const dictionary& motionControlDict
        );

    //- Destructor
    ~cellSizeControlSurfaces();


    // Member Functions

        // Access

            //- Return reference to the searchableSurfaces object containing all
            //  of the geometry
            inline const searchableSurfaces& geometry() const;

            //- Return the surface indices
            inline const labelList& surfaces() const;

            //- In regions where no cell size function is specified
            // use defaultCellSize
            inline scalar defaultCellSize() const;

        // Query

            //- Return the cell size at the given location
            scalar cellSize(const point& pt) const;

            //- Return the cell size at the given locations
            scalarField cellSize(const pointField& pts) const;

        // Edit

            //- Set the cell sizes for each point
            void setCellSizes(const pointField& pts);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "cellSizeControlSurfacesI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

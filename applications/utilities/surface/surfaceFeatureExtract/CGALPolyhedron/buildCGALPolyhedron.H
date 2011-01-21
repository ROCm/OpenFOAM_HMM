/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2011 OpenCFD Ltd.
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
    Foam::buildCGALPolyhedron

Description
    Convert a triSurface into a CGAL Polyhedron

SourceFiles
    buildCGALPolyhedron.C

\*---------------------------------------------------------------------------*/

#ifndef buildCGALPolyhedron_H
#define buildCGALPolyhedron_H

#include "triSurface.H"
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Polyhedron_3.h>

typedef CGAL::Simple_cartesian<double>               Kernel;
typedef CGAL::Polyhedron_3<Kernel>                   Polyhedron;
typedef Polyhedron::HalfedgeDS                       HalfedgeDS;
typedef Polyhedron::Vertex                           Vertex;
typedef Polyhedron::Vertex_iterator                  Vertex_iterator;
typedef Kernel::Point_3                              Point_3;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class buildCGALPolyhedron Declaration
\*---------------------------------------------------------------------------*/

class buildCGALPolyhedron
:
    public CGAL::Modifier_base<HalfedgeDS>
{
    // Private data

        //- Reference to triSurface to convert
        const Foam::triSurface& surf_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        buildCGALPolyhedron(const buildCGALPolyhedron&);

        //- Disallow default bitwise assignment
        void operator=(const buildCGALPolyhedron&);


public:

    // Constructors

        //- Construct with reference to triSurface
        buildCGALPolyhedron(const triSurface& surf);


    //- Destructor
    ~buildCGALPolyhedron();


    // Member Operators

        //- operator() of this `modifier' called by delegate function of
        //  Polyhedron
        void operator()(HalfedgeDS& hds);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

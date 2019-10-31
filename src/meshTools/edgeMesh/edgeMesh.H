/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2014 OpenFOAM Foundation
    Copyright (C) 2017-2018 OpenCFD Ltd.
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
    Foam::edgeMesh

Description
    Points connected by edges.

    Can be read from fileName based on extension. Uses ::New factory method
    to select the reader and transfer the result.

SourceFiles
    edgeMeshI.H
    edgeMesh.C
    edgeMeshIO.C
    edgeMeshNew.C

\*---------------------------------------------------------------------------*/

#ifndef edgeMesh_H
#define edgeMesh_H

#include "pointField.H"
#include "edgeList.H"
#include "edgeMeshFormatsCore.H"
#include "runTimeSelectionTables.H"
#include "memberFunctionSelectionTables.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations

class edgeMesh;
class Istream;
class Ostream;

Istream& operator>>(Istream& is, edgeMesh& em);
Ostream& operator<<(Ostream& os, const edgeMesh& em);


/*---------------------------------------------------------------------------*\
                           Class edgeMesh Declaration
\*---------------------------------------------------------------------------*/

class edgeMesh
:
    public fileFormats::edgeMeshFormatsCore
{
    // Private data

        //- Vertices of the edges
        pointField points_;

        //- The edges defining the boundary
        edgeList edges_;

        //- From point to edges
        mutable autoPtr<labelListList> pointEdgesPtr_;


    // Private Member Functions

        //- Calculate point-edge addressing (inverse of edges)
        void calcPointEdges() const;


protected:

    // Protected Member Functions

        //- Non-const access to global points
        inline pointField& storedPoints();

        //- Non-const access to the edges
        inline edgeList& storedEdges();


public:

        //- Runtime type information
        TypeName("edgeMesh");


    // Static

        //- Can we read this file format?
        static bool canRead(const fileName& name, bool verbose=false);

        //- Can we read this file format?
        static bool canReadType(const word& ext, bool verbose=false);

        //- Can we write this file format type?
        static bool canWriteType(const word& ext, bool verbose=false);

        static wordHashSet readTypes();
        static wordHashSet writeTypes();


    // Constructors

        //- Construct null
        inline edgeMesh();

        //- Copy construct
        inline edgeMesh(const edgeMesh& em);

        //- Move construct
        inline edgeMesh(edgeMesh&& em);

        //- Copy construct from components
        inline edgeMesh(const pointField& points, const edgeList& edges);

        //- Move construct from components
        inline edgeMesh(pointField&& pointLst, edgeList&& edgeLst);

        //- Construct from file name (uses extension to determine type)
        edgeMesh(const fileName& name);

        //- Construct from file name (uses extension to determine type)
        edgeMesh(const fileName& name, const word& ext);


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            edgeMesh,
            fileExtension,
            (
                const fileName& name
            ),
            (name)
        );


    // Selectors

        //- Select constructed from filename (explicit extension)
        static autoPtr<edgeMesh> New
        (
            const fileName& name,
            const word& ext
        );

        //- Select constructed from filename (implicit extension)
        static autoPtr<edgeMesh> New(const fileName& name);


    //- Destructor
    virtual ~edgeMesh() = default;


    // Member Function Selectors

        declareMemberFunctionSelectionTable
        (
            void,
            edgeMesh,
            write,
            fileExtension,
            (
                const fileName& name,
                const edgeMesh& mesh
            ),
            (name, mesh)
        );

        //- Write to file
        static void write(const fileName& name, const edgeMesh& mesh);


    // Member Functions

        //- Transfer the contents of the argument and annul the argument
        void transfer(edgeMesh& mesh);


    // Read

        //- Read from file. Chooses reader based on explicit extension
        bool read(const fileName& name, const word& ext);

        //- Read from file. Chooses reader based on detected extension
        virtual bool read(const fileName& name);


    // Access

        //- Return points
        inline const pointField& points() const;

        //- Return edges
        inline const edgeList& edges() const;

        //- Return edges
        inline const labelListList& pointEdges() const;

        //- Find connected regions. Set region number per edge.
        //  Returns number of regions.
        label regions(labelList& edgeRegion) const;


    // Edit

        //- Clear all storage
        virtual void clear();

        //- Scale points. A non-positive factor is ignored
        virtual void scalePoints(const scalar scaleFactor);

        //- Geometric merge points (points within mergeDist) prior to
        //  automatically calling mergeEdges().
        virtual void mergePoints(const scalar mergeDist);

        //- Merge duplicate edges and eliminate unused points.
        virtual void mergeEdges();


    // Write

        virtual void writeStats(Ostream&) const;

        //- Generic write routine. Chooses writer based on extension.
        virtual void write(const fileName& name) const
        {
            write(name, *this);
        }


    // Member Operators

        //- Copy assignment
        inline void operator=(const edgeMesh& rhs);

        //- Move assignment
        inline void operator=(edgeMesh&& rhs);


    // Ostream Operator

        friend Ostream& operator<<(Ostream& os, const edgeMesh& em);
        friend Istream& operator>>(Istream& is, edgeMesh& em);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "edgeMeshI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

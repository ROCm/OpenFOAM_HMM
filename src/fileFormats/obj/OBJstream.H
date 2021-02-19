/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2015 OpenFOAM Foundation
    Copyright (C) 2017-2021 OpenCFD Ltd.
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
    Foam::OBJstream

Description
    OFstream that keeps track of vertices

SourceFiles
    OBJstream.C

\*---------------------------------------------------------------------------*/

#ifndef OBJstream_H
#define OBJstream_H

#include "OFstream.H"
#include "point.H"
#include "edge.H"
#include "face.H"
#include "triPointRef.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class treeBoundBox;

/*---------------------------------------------------------------------------*\
                          Class OBJstream Declaration
\*---------------------------------------------------------------------------*/

class OBJstream
:
    public OFstream
{
    // Private Data

        bool startOfLine_;

        label nVertices_;


    // Private Member Functions

        void writeAndCheck(const char c);

public:

    //- Declare type-name (with debug switch)
    ClassName("OBJstream");


    // Constructors

        //- Construct from pathname
        explicit OBJstream
        (
            const fileName& pathname,
            IOstreamOption streamOpt = IOstreamOption()
        );

        //- Construct from pathname
        OBJstream
        (
            const fileName& pathname,
            IOstreamOption::streamFormat fmt,
            IOstreamOption::versionNumber ver = IOstreamOption::currentVersion,
            IOstreamOption::compressionType cmp = IOstreamOption::UNCOMPRESSED
        )
        :
            OBJstream(pathname, IOstreamOption(fmt, ver, cmp))
        {}


    //- Destructor
    ~OBJstream() = default;


    // Member Functions

        // Access

            //- Return the number of vertices written
            label nVertices() const
            {
                return nVertices_;
            }


        // Ostream implementation

            //- Inherit write from Ostream
            using Ostream::write;

            //- Write character
            virtual Ostream& write(const char c);

            //- Write character string
            virtual Ostream& write(const char* str);

            //- Write word
            virtual Ostream& write(const word& str);

            //- Write string
            virtual Ostream& write(const string& str);

            //- Write std::string surrounded by quotes.
            //  Optional write without quotes.
            virtual Ostream& writeQuoted
            (
                const std::string& str,
                const bool quoted=true
            );


        // Direct write functionality

            //- Write point
            Ostream& write(const point& pt);

            //- Write point and vector normal ('vn')
            Ostream& write(const point& pt, const vector& n);

            //- Write edge as points and line
            Ostream& write(const edge& e, const UList<point>& points);

            //- Write line
            Ostream& write(const linePointRef& ln);

            //- Write line with points and vector normals ('vn')
            Ostream& write
            (
                const linePointRef& ln,
                const vector& n0,
                const vector& n1
            );

            //- Write triangle as points and lines/filled-polygon
            Ostream& write(const triPointRef& f, const bool lines = true);

            //- Write face as points and lines/filled-polygon
            Ostream& write
            (
                const face& f,
                const UList<point>& points,
                const bool lines = true
            );

            //- Write patch faces as points and lines/filled-polygon
            Ostream& write
            (
                const UList<face>& faces,
                const pointField& points,
                const bool lines = true
            );

            //- Write edges as points and lines.
            //  Optionally eliminate unused points.
            Ostream& write
            (
                const UList<edge>& edges,
                const UList<point>& points,
                const bool compact = false
            );

            //- Write tree-bounding box as lines/filled-polygons
            Ostream& write
            (
                const treeBoundBox& bb,
                const bool lines = true
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

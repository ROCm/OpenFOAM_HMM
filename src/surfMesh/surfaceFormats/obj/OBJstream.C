/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2016 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "OBJstream.H"
#include "primitivePatch.H"
#include "treeBoundBox.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(OBJstream, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::OBJstream::writeAndCheck(const char c)
{
    if (c == '\n')
    {
        startOfLine_ = true;
    }
    else if (startOfLine_)
    {
        startOfLine_ = false;
        if (c == 'v')
        {
            ++nVertices_;
        }
    }

    OFstream::write(c);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::OBJstream::OBJstream
(
    const fileName& pathname,
    streamFormat format,
    versionNumber version,
    compressionType compression
)
:
    OFstream(pathname, format, version, compression),
    startOfLine_(true),
    nVertices_(0)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::Ostream& Foam::OBJstream::write(const char c)
{
    writeAndCheck(c);
    return *this;
}


Foam::Ostream& Foam::OBJstream::write(const char* str)
{
    for (const char* iter = str; *iter; ++iter)
    {
        writeAndCheck(*iter);
    }
    return *this;
}


Foam::Ostream& Foam::OBJstream::write(const word& str)
{
    return writeQuoted(str, false);
}


Foam::Ostream& Foam::OBJstream::write(const string& str)
{
    return writeQuoted(str, true);
}


Foam::Ostream& Foam::OBJstream::writeQuoted
(
    const std::string& str,
    const bool quoted
)
{
    if (!quoted)
    {
        // Output unquoted, only advance line number on newline
        for (auto iter = str.cbegin(); iter != str.cend(); ++iter)
        {
            writeAndCheck(*iter);
        }
        return *this;
    }


    OFstream::write(token::BEGIN_STRING);

    unsigned backslash = 0;
    for (auto iter = str.cbegin(); iter != str.cend(); ++iter)
    {
        const char c = *iter;

        if (c == '\\')
        {
            ++backslash;
            continue; // only output after escaped character is known
        }
        else if (c == token::NL)
        {
            ++lineNumber_;
            ++backslash;    // backslash escape for newline
        }
        else if (c == token::END_STRING)
        {
            ++backslash;    // backslash escape for quote
        }

        // output all pending backslashes
        while (backslash)
        {
            OFstream::write('\\');
            --backslash;
        }

        writeAndCheck(c);
    }

    // silently drop any trailing backslashes
    // they would otherwise appear like an escaped end-quote
    OFstream::write(token::END_STRING);

    return *this;
}


Foam::Ostream& Foam::OBJstream::write(const point& pt)
{
    write("v ") << pt.x() << ' ' << pt.y() << ' ' << pt.z() << nl;
    return *this;
}


Foam::Ostream& Foam::OBJstream::write(const point& pt, const vector& n)
{
    write(pt);
    OFstream::write("vn ") << n.x() << ' ' << n.y() << ' ' << n.z() << nl;
    return *this;
}


Foam::Ostream& Foam::OBJstream::write(const edge& e, const UList<point>& points)
{
    write(points[e[0]]);
    write(points[e[1]]);
    write("l ") << nVertices_-1 << ' ' << nVertices_ << nl;
    return *this;
}


Foam::Ostream& Foam::OBJstream::write(const linePointRef& ln)
{
    write(ln.start());
    write(ln.end());
    write("l ") << nVertices_-1 << ' ' << nVertices_ << nl;
    return *this;
}


Foam::Ostream& Foam::OBJstream::write
(
    const linePointRef& ln,
    const vector& n0,
    const vector& n1
)
{
    write(ln.start(), n0);
    write(ln.end(), n1);
    write("l ") << nVertices_-1 << ' ' << nVertices_ << nl;
    return *this;
}


Foam::Ostream& Foam::OBJstream::write
(
    const triPointRef& f,
    const bool lines
)
{
    const label start = nVertices_+1;  // 1-offset for obj included here
    write(f.a());
    write(f.b());
    write(f.c());
    if (lines)
    {
        write('l');
        for (int i = 0; i < 3; i++)
        {
            write(' ') << i+start;
        }
        write(' ') << start << '\n';
    }
    else
    {
        write('f');
        for (int i = 0; i < 3; i++)
        {
            write(' ') << i+start;
        }
        write('\n');
    }
    return *this;
}


Foam::Ostream& Foam::OBJstream::write
(
    const face& f,
    const UList<point>& points,
    const bool lines
)
{
    const label start = nVertices_+1;  // 1-offset for obj included here
    forAll(f, i)
    {
        write(points[f[i]]);
    }
    if (lines)
    {
        write('l');
        forAll(f, i)
        {
            write(' ') << i+start;
        }
        write(' ') << start << '\n';
    }
    else
    {
        write('f');
        forAll(f, i)
        {
            write(' ') << i+start;
        }
        write('\n');
    }
    return *this;
}


Foam::Ostream& Foam::OBJstream::write
(
    const UList<face>& faces,
    const pointField& points,
    const bool lines
)
{
    SubList<face> allFcs(faces, faces.size());

    primitivePatch pp(allFcs, points);

    const pointField& localPoints = pp.localPoints();
    const faceList& localFaces = pp.localFaces();

    const label start = nVertices_+1;  // 1-offset for obj included here

    forAll(localPoints, i)
    {
        write(localPoints[i]);
    }

    if (lines)
    {
        const edgeList& edges = pp.edges();
        forAll(edges, edgeI)
        {
            const edge& e = edges[edgeI];

            write("l ") << e[0]+start << ' ' << e[1]+start << nl;
        }
    }
    else
    {
        forAll(localFaces, facei)
        {
            const face& f = localFaces[facei];
            write('f');
            forAll(f, i)
            {
                write(' ') << f[i]+start;
            }
            write('\n');
        }
    }
    return *this;
}


Foam::Ostream& Foam::OBJstream::write
(
    const UList<edge>& edges,
    const UList<point>& points,
    const bool compact
)
{
    if (compact)
    {
        // Code similar to PrimitivePatch::calcMeshData()
        // Unsorted version

        label objPointId = nVertices_+1;  // 1-offset for obj included here

        Map<label> markedPoints(2*edges.size());
        forAll(edges, edgei)
        {
            const edge& e = edges[edgei];

            if (markedPoints.insert(e[0], objPointId))
            {
                write(points[e[0]]);
                ++objPointId;
            }
            if (markedPoints.insert(e[1], objPointId))
            {
                write(points[e[1]]);
                ++objPointId;
            }
        }

        forAll(edges, edgei)
        {
            const edge& e = edges[edgei];

            write("l ")
                << markedPoints[e[0]] << ' '
                << markedPoints[e[1]] << nl;
        }
    }
    else
    {
        const label start = nVertices_+1;  // 1-offset for obj included here

        forAll(points, i)
        {
            write(points[i]);
        }

        forAll(edges, edgei)
        {
            const edge& e = edges[edgei];

            write("l ")
                << e[0]+start << ' ' << e[1]+start << nl;
        }
    }

    return *this;
}


Foam::Ostream& Foam::OBJstream::write
(
    const treeBoundBox& bb,
    const bool lines
)
{
    const label start = nVertices_+1;  // 1-offset for obj included here

    pointField points(bb.points());
    forAll(points, i)
    {
        write(points[i]);
    }

    if (lines)
    {
        forAll(treeBoundBox::edges, edgei)
        {
            const edge& e = treeBoundBox::edges[edgei];

            write("l ") << e[0]+start << ' ' << e[1]+start << nl;
        }
    }
    else
    {
        forAll(treeBoundBox::faces, facei)
        {
            const face& f = treeBoundBox::faces[facei];

            write('f');
            forAll(f, i)
            {
                write(' ') << f[i]+start;
            }
            write('\n');
        }
    }

    return *this;
}


// ************************************************************************* //

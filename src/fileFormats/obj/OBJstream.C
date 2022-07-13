/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
    Copyright (C) 2020-2022 OpenCFD Ltd.
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
    defineTypeName(OBJstream);
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
    IOstreamOption streamOpt
)
:
    OFstream(pathname, streamOpt),
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

    // Output with surrounding quotes and backslash escaping
    OFstream::write(static_cast<char>(token::DQUOTE));

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
        else if (c == token::DQUOTE)
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
    OFstream::write(static_cast<char>(token::DQUOTE));

    return *this;
}


Foam::Ostream& Foam::OBJstream::write(const point& p)
{
    write('v') << ' ' << p.x() << ' ' << p.y() << ' ' << p.z() << nl;
    return *this;
}


Foam::Ostream& Foam::OBJstream::write(const point& p, const vector& n)
{
    write(p);
    OFstream::write("vn ") << n.x() << ' ' << n.y() << ' ' << n.z() << nl;
    return *this;
}


Foam::Ostream& Foam::OBJstream::write(const UList<point>& points)
{
    for (const point& p : points)
    {
        write('v') << ' ' << p.x() << ' ' << p.y() << ' ' << p.z() << nl;
    }
    return *this;
}


Foam::Ostream& Foam::OBJstream::write(const edge& e, const UList<point>& points)
{
    write(points[e.first()]);
    write(points[e.second()]);
    write('l') << ' ' << nVertices_-1 << ' ' << nVertices_ << nl;
    return *this;
}


Foam::Ostream& Foam::OBJstream::write(const linePointRef& ln)
{
    write(ln.first());
    write(ln.second());
    write('l') << ' ' << nVertices_-1 << ' ' << nVertices_ << nl;
    return *this;
}


Foam::Ostream& Foam::OBJstream::write
(
    const linePointRef& ln,
    const vector& n0,
    const vector& n1
)
{
    write(ln.first(), n0);
    write(ln.second(), n1);
    write('l') << ' ' << nVertices_-1 << ' ' << nVertices_ << nl;
    return *this;
}


Foam::Ostream& Foam::OBJstream::writeLine
(
    const point& p0,
    const point& p1
)
{
    write(p0);
    write(p1);
    write('l') << ' ' << nVertices_-1 << ' ' << nVertices_ << nl;
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
        for (int i = 0; i < 3; ++i)
        {
            write(' ') << i+start;
        }
        write(' ') << start << '\n';
    }
    else
    {
        write('f');
        for (int i = 0; i < 3; ++i)
        {
            write(' ') << i+start;
        }
        write('\n');
    }
    return *this;
}


Foam::Ostream& Foam::OBJstream::writeFace
(
    const UList<point>& points,
    const bool lines
)
{
    const label start = nVertices_+1;  // 1-offset for obj included here

    write(points);

    if (lines)
    {
        write('l');
        forAll(points, i)
        {
            write(' ') << i+start;
        }
        write(' ') << start << '\n';
    }
    else
    {
        write('f');
        forAll(points, i)
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

    for (const label fp : f)
    {
        write(points[fp]);
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
    primitivePatch pp(SubList<face>(faces), points);

    const label start = nVertices_+1;  // 1-offset for obj included here

    write(pp.localPoints());

    if (lines)
    {
        for (const edge& e : pp.edges())
        {
            write('l') << ' '
                << e.first()+start << ' '
                << e.second()+start << nl;
        }
    }
    else
    {
        for (const face& f : pp.localFaces())
        {
            write('f');
            for (const label fp : f)
            {
                write(' ') << fp+start;
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

        for (const edge& e : edges)
        {
            if (markedPoints.insert(e.first(), objPointId))
            {
                write(points[e.first()]);
                ++objPointId;
            }
            if (markedPoints.insert(e.second(), objPointId))
            {
                write(points[e.second()]);
                ++objPointId;
            }
        }

        for (const edge& e : edges)
        {
            write('l') << ' '
                << markedPoints[e.first()] << ' '
                << markedPoints[e.second()] << nl;
        }
    }
    else
    {
        const label start = nVertices_+1;  // 1-offset for obj included here

        write(points);

        for (const edge& e : edges)
        {
            write('l') << ' '
                << e.first()+start << ' '
                << e.second()+start << nl;
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

    write(bb.points());

    if (lines)
    {
        for (const edge& e : treeBoundBox::edges)
        {
            write('l') << ' '
                << e.first()+start << ' '
                << e.second()+start << nl;
        }
    }
    else
    {
        for (const face& f : treeBoundBox::faces)
        {
            write('f');
            for (const label fp : f)
            {
                write(' ') << fp+start;
            }
            write('\n');
        }
    }

    return *this;
}


// ************************************************************************* //

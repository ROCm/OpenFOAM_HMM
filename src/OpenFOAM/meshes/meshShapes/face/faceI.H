/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd.
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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Edge to the right of face vertex i
inline Foam::label Foam::face::right(const label i) const
{
    return i;
}


// Edge to the left of face vertex i
inline Foam::label Foam::face::left(const label i) const
{
    return rcIndex(i);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::face::face()
{}


inline Foam::face::face(const label sz)
:
    labelList(sz, -1)
{}


inline Foam::face::face(const labelUList& lst)
:
    labelList(lst)
{}


template<unsigned Size>
inline Foam::face::face(const FixedList<label, Size>& lst)
:
    labelList(lst)
{}


inline Foam::face::face(std::initializer_list<label> lst)
:
    labelList(lst)
{}


inline Foam::face::face(const Xfer<labelList>& lst)
:
    labelList(lst)
{}


inline Foam::face::face(labelList&& lst)
:
    labelList(std::move(lst))
{}


inline Foam::face::face(Istream& is)
{
    is >> *this;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::pointField Foam::face::points
(
    const UList<point>& meshPoints
) const
{
    // There are as many points as there are labels for them
    pointField p(size());

    // For each point in list, set it to the point in 'pnts' addressed
    // by 'labs'
    label i = 0;
    for (const label pointi : *this)
    {
        p[i++] = meshPoints[pointi];
    }

    // Return list
    return p;
}


inline Foam::scalar Foam::face::mag(const UList<point>& p) const
{
    return ::Foam::mag(normal(p));
}


inline Foam::label Foam::face::nEdges() const
{
    // for a closed polygon a number of edges is the same as number of points
    return size();
}


inline Foam::edge Foam::face::faceEdge(const label n) const
{
    return edge(operator[](n), operator[](fcIndex(n)));
}


inline bool Foam::face::found(const label globalIndex) const
{
    return which(globalIndex) != -1;
}


inline Foam::label Foam::face::nextLabel(const label i) const
{
    return this->fcValue(i);
}


inline Foam::label Foam::face::prevLabel(const label i) const
{
    return this->rcValue(i);
}


inline Foam::label Foam::face::nTriangles() const
{
    return size() - 2;
}


// * * * * * * * * * * * * * * * Global Operators   * * * * * * * * * * * * * //

inline bool Foam::operator==(const face& a, const face& b)
{
    return face::compare(a,b) != 0;
}

inline bool Foam::operator!=(const face& a, const face& b)
{
    return face::compare(a,b) == 0;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

inline Foam::Istream& Foam::operator>>(Istream& is, face& f)
{
    if (is.version() == IOstream::originalVersion)
    {
        // Read starting (
        is.readBegin("face");

        // Read the 'name' token for the face
        token t(is);

        // Read labels
        is >> static_cast<labelList&>(f);

        // Read end)
        is.readEnd("face");
    }
    else
    {
        is >> static_cast<labelList&>(f);
    }

    is.check(FUNCTION_NAME);
    return is;
}

// ************************************************************************* //

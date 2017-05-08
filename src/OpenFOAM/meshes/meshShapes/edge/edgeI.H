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

#include "IOstreams.H"
#include "Swap.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

inline int Foam::edge::compare(const edge& a, const edge& b)
{
    if (a[0] == b[0] && a[1] == b[1])
    {
        return 1;
    }
    else if (a[0] == b[1] && a[1] == b[0])
    {
        return -1;
    }
    else
    {
        return 0;
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class InputIter>
inline Foam::label Foam::edge::insertMultiple
(
    const InputIter begIter,
    const InputIter endIter
)
{
    // Available slots.
    // Don't use count() since it has special treatment for duplicates
    const int maxChange = (start() < 0 ? 1 : 0) + (end() < 0 ? 1 : 0);

    int changed = 0;
    if (maxChange)
    {
        for (InputIter iter = begIter; iter != endIter; ++iter)
        {
            if (insert(*iter))
            {
                if (++changed >= maxChange)
                {
                    break;
                }
            }
        }
    }

    return changed;
}


template<class InputIter>
inline Foam::label Foam::edge::eraseMultiple
(
    const InputIter begIter,
    const InputIter endIter
)
{
    // Occupied slots.
    // Don't use count() since it has special treatment for duplicates
    const int maxChange = (start() >= 0 ? 1 : 0) + (end() >= 0 ? 1 : 0);

    int changed = 0;
    if (maxChange)
    {
        for (InputIter iter = begIter; iter != endIter; ++iter)
        {
            changed += erase(*iter);
            if (changed >= maxChange)
            {
                break;
            }
        }
    }

    return changed;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::edge::edge()
:
    FixedList<label, 2>(-1)
{}


inline Foam::edge::edge(const label from, const label to)
{
    start() = from;
    end()   = to;
}


inline Foam::edge::edge(const label from, const label to, const bool doSort)
{
    if (doSort && from > to)
    {
        start() = to;
        end()   = from;
    }
    else
    {
        start() = from;
        end()   = to;
    }
}


inline Foam::edge::edge(const labelPair& pair)
{
    start() = pair.first();
    end()   = pair.second();
}


inline Foam::edge::edge(const FixedList<label, 2>& lst)
{
    start() = lst[0];
    end()   = lst[1];
}


inline Foam::edge::edge(const FixedList<label, 2>& lst, const bool doSort)
{
    if (doSort && lst[0] > lst[1])
    {
        start() = lst[1];
        end()   = lst[0];
    }
    else
    {
        start() = lst[0];
        end()   = lst[1];
    }
}


inline Foam::edge::edge(Istream& is)
:
    FixedList<label, 2>(is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::label Foam::edge::start() const
{
    return operator[](0);
}

inline Foam::label& Foam::edge::start()
{
    return operator[](0);
}


inline Foam::label Foam::edge::end() const
{
    return operator[](1);
}

inline Foam::label& Foam::edge::end()
{
    return operator[](1);
}


inline Foam::label Foam::edge::minVertex() const
{
    return (start() < end() ? start() : end());
}


inline Foam::label Foam::edge::maxVertex() const
{
    return (start() > end() ? start() : end());
}


inline bool Foam::edge::found(const label pointLabel) const
{
    // -1: always false
    return (pointLabel >= 0 && (pointLabel == start() || pointLabel == end()));
}


inline Foam::label Foam::edge::which(const label pointLabel) const
{
    // -1: always false
    if (pointLabel >= 0)
    {
        if (pointLabel == start())
        {
            return 0;
        }
        if (pointLabel == end())
        {
            return 1;
        }
    }
    return -1;
}


inline bool Foam::edge::connects(const edge& other) const
{
    return (other.found(start()) || other.found(end()));
}


inline Foam::label Foam::edge::commonVertex(const edge& other) const
{
    if (other.found(start()))
    {
        return start();
    }
    else if (other.found(end()))
    {
        return end();
    }
    else
    {
        // No shared vertex.
        return -1;
    }
}


inline Foam::label Foam::edge::otherVertex(const label index) const
{
    if (index == start())
    {
        return end();
    }
    else if (index == end())
    {
        return start();
    }
    else
    {
        // The given vertex is not on the edge in the first place.
        return -1;
    }
}


inline Foam::label Foam::edge::collapse()
{
    // Cannot resize FixedList, so mark duplicates with '-1'
    // (the lower vertex is retained)
    // catch any '-1' (eg, if called multiple times)

    label n = 2;
    if (start() == end() || end() < 0)
    {
        end() = -1;
        --n;
    }
    if (start() < 0)
    {
        --n;
    }

    return n;
}


inline bool Foam::edge::sorted() const
{
    return (start() < end());
}


inline void Foam::edge::sort()
{
    if (start() > end())
    {
        flip();
    }
}


inline void Foam::edge::flip()
{
    Swap(operator[](0), operator[](1));
}


inline Foam::edge Foam::edge::reverseEdge() const
{
    return edge(end(), start());
}


inline void Foam::edge::clear()
{
    start() = -1;
    end()   = -1;
}


inline Foam::label Foam::edge::count() const
{
    label n = 2;
    if (start() == end() || end() < 0)
    {
        --n;
    }
    if (start() < 0)
    {
        --n;
    }

    return n;
}


inline bool Foam::edge::empty() const
{
    return (start() < 0 && end() < 0);
}


inline bool Foam::edge::insert(const label index)
{
    if (index < 0)
    {
        // Cannot insert invalid point labels (use direct assignment for that)
        return false;
    }

    if (start() < 0)
    {
        // Store at [0], if not duplicate of [1]
        if (index != end())
        {
            start() = index;
            return true;
        }
    }
    else if (end() < 0)
    {
        // Store at [1], if not duplicate of [0]
        if (index != start())
        {
            end() = index;
            return true;
        }
    }

    return false;
}


inline Foam::label Foam::edge::insert(const UList<label>& lst)
{
    return insertMultiple(lst.begin(), lst.end());
}


template<unsigned AnySize>
inline Foam::label Foam::edge::insert(const FixedList<label, AnySize>& lst)
{
    return insertMultiple(lst.begin(), lst.end());
}


inline Foam::label Foam::edge::insert(std::initializer_list<label> lst)
{
    return insertMultiple(lst.begin(), lst.end());
}


inline Foam::label Foam::edge::erase(const label index)
{
    if (index < 0)
    {
        // Can never remove invalid point labels!
        return 0;
    }

    label n = 0;
    if (index == start())
    {
        start() = -1;
        ++n;
    }

    // Automatically handle duplicates, which should not have been there anyhow
    if (index == end())
    {
        end() = -1;
        ++n;
    }

    return n;
}


inline Foam::label Foam::edge::erase(const UList<label>& lst)
{
    return eraseMultiple(lst.begin(), lst.end());
}


template<unsigned AnySize>
inline Foam::label Foam::edge::erase(const FixedList<label, AnySize>& lst)
{
    return eraseMultiple(lst.begin(), lst.end());
}


inline Foam::label Foam::edge::erase(std::initializer_list<label> lst)
{
    return eraseMultiple(lst.begin(), lst.end());
}


// Geometric

inline Foam::point Foam::edge::centre(const UList<point>& pts) const
{
    #ifdef FULLDEBUG
    if (start() < 0 || end() < 0)
    {
        FatalErrorInFunction
            << "negative point index on edge " << *this
            << abort(FatalError);
    }
    #endif

    return 0.5*(pts[start()] + pts[end()]);
}


inline Foam::vector Foam::edge::vec(const UList<point>& pts) const
{
    #ifdef FULLDEBUG
    if (start() < 0 || end() < 0)
    {
        FatalErrorInFunction
            << "negative point index on edge " << *this
            << abort(FatalError);
    }
    #endif

    return pts[end()] - pts[start()];
}


inline Foam::vector Foam::edge::unitVec(const UList<point>& pts) const
{
    #ifdef FULLDEBUG
    if (start() < 0 || end() < 0)
    {
        FatalErrorInFunction
            << "negative point index on edge " << *this
            << abort(FatalError);
    }
    #endif

    Foam::vector v = pts[end()] - pts[start()];
    v /= ::Foam::mag(v) + VSMALL;

    return v;
}


inline Foam::scalar Foam::edge::mag(const UList<point>& pts) const
{
    return ::Foam::mag(vec(pts));
}


inline Foam::linePointRef Foam::edge::line(const UList<point>& pts) const
{
    #ifdef FULLDEBUG
    if (start() < 0 || end() < 0)
    {
        FatalErrorInFunction
            << "negative point index on edge " << *this
            << abort(FatalError);
    }
    #endif

    return linePointRef(pts[start()], pts[end()]);
}


// * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * * //

inline bool Foam::operator==(const edge& a, const edge& b)
{
    return edge::compare(a,b) != 0;
}


inline bool Foam::operator!=(const edge& a, const edge& b)
{
    return edge::compare(a,b) == 0;
}


// ************************************************************************* //

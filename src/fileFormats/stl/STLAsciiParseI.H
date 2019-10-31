/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018 OpenCFD Ltd.
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

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

inline void Foam::Detail::STLAsciiParse::beginSolid(word solidName)
{
    if (solidName.empty())
    {
        solidName = "solid";  // Could also use solid0, solid1, ...
    }

    auto iter = nameLookup_.cfind(solidName);
    if (iter.found())
    {
        if (groupId_ != iter.val())
        {
            sorted_ = false;  // Group appeared out of order
            groupId_ = iter.val();
        }
    }
    else
    {
        groupId_ = sizes_.size();
        if (nameLookup_.insert(solidName, groupId_))
        {
            names_.append(solidName);
            sizes_.append(0);
        }
        else
        {
            FatalErrorInFunction<< "Duplicate solid-name: " << solidName
                << exit(FatalError);
        }
    }
}


inline void Foam::Detail::STLAsciiParse::beginFacet()
{
    nFacetPoints_ = 0;
    nVertexCmpt_ = 0;
}


inline void Foam::Detail::STLAsciiParse::resetVertex()
{
    nVertexCmpt_ = 0;
}


inline bool Foam::Detail::STLAsciiParse::addVertexComponent(float val)
{
    currVertex_[nVertexCmpt_] = val;

    if (++nVertexCmpt_ == 3)
    {
        points_.append(currVertex_);
        nVertexCmpt_ = 0;
        ++nFacetPoints_;
    }

    return !nVertexCmpt_;
}


inline bool Foam::Detail::STLAsciiParse::addVertexComponent(const char* text)
{
    //-> safer, but slower: readFloat(text, currVertex_[nVertexCmpt_]);
    currVertex_[nVertexCmpt_] = ::atof(text);

    if (++nVertexCmpt_ == 3)
    {
        points_.append(currVertex_);
        nVertexCmpt_ = 0;
        ++nFacetPoints_;
    }

    return !nVertexCmpt_;
}


inline void Foam::Detail::STLAsciiParse::endFacet()
{
    if (nFacetPoints_ == 3)
    {
        facets_.append(groupId_);
        sizes_[groupId_]++;
    }
    else
    {
        if (nFacetPoints_ > 3)
        {
            nFacetPoints_ -= 3;
        }

        if (nFacetPoints_)
        {
            points_.resize(points_.size() - nFacetPoints_);
        }
    }

    nFacetPoints_ = 0;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::Detail::STLAsciiParse::STLAsciiParse(const label approxNpoints)
:
    sorted_(true),
    groupId_(-1),
    lineNum_(1),
    nFacetPoints_(0),
    nVertexCmpt_(0),
    points_(approxNpoints),
    facets_(approxNpoints/2)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline void Foam::Detail::STLAsciiParse::clear()
{
    sorted_ = true;
    groupId_ = -1;
    lineNum_ = 0;

    nFacetPoints_ = 0;
    nVertexCmpt_ = 0;

    points_.clear();
    facets_.clear();
    names_.clear();
    sizes_.clear();
    nameLookup_.clear();
}


inline bool Foam::Detail::STLAsciiParse::sorted() const
{
    return sorted_;
}


inline Foam::DynamicList<Foam::STLpoint>& Foam::Detail::STLAsciiParse::points()
{
    return points_;
}


inline Foam::DynamicList<Foam::label>& Foam::Detail::STLAsciiParse::facets()
{
    return facets_;
}


inline Foam::DynamicList<Foam::word>& Foam::Detail::STLAsciiParse::names()
{
    return names_;
}


inline Foam::DynamicList<Foam::label>& Foam::Detail::STLAsciiParse::sizes()
{
    return sizes_;
}


// ************************************************************************* //

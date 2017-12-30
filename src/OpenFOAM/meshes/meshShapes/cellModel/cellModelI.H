/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::word& Foam::cellModel::name() const
{
    return name_;
}


inline Foam::label Foam::cellModel::index() const
{
    return index_;
}


inline Foam::label Foam::cellModel::nPoints() const
{
    return nPoints_;
}


inline Foam::label Foam::cellModel::nEdges() const
{
    return edges_.size();
}


inline Foam::label Foam::cellModel::nFaces() const
{
    return faces_.size();
}


inline const Foam::edgeList& Foam::cellModel::modelEdges() const
{
    return edges_;
}


inline const Foam::faceList& Foam::cellModel::modelFaces() const
{
    return faces_;
}


//  Return the faces of a cellModel by untangling the geometry
//  supplied in terms of the face labels
inline Foam::edgeList Foam::cellModel::edges
(
    const labelUList& pointLabels
) const
{
    edgeList e(edges_.size());

    // Translate model labels into global labels
    forAll(edges_, edgei)
    {
         e[edgei] =
             edge
             (
                 pointLabels[edges_[edgei].start()],
                 pointLabels[edges_[edgei].end()]
             );
    }

    return e;
}


//  Return the faces of a cellModel by untangling the geometry
//  supplied in terms of the face labels
inline Foam::faceList Foam::cellModel::faces
(
    const labelUList& pointLabels
) const
{
    faceList f(faces_.size());

    // Translate model labels into global labels
    forAll(faces_, facei)
    {
         const labelList& curModelLabels = faces_[facei];

         face& curFace = f[facei];

         curFace.setSize(curModelLabels.size());

         forAll(curModelLabels, labeli)
         {
             curFace[labeli] = pointLabels[curModelLabels[labeli]];
         }
    }

    return f;
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

inline bool Foam::operator==(const cellModel& lhs, const cellModel& rhs)
{
    return (&lhs == &rhs);
}


inline bool Foam::operator!=(const cellModel& lhs, const cellModel& rhs)
{
    return (&lhs != &rhs);
}


// ************************************************************************* //

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
    pointFeatureEdgesTypes

Description
    struct for holding information on the types of feature edges attached to
    feature points

\*---------------------------------------------------------------------------*/

#ifndef pointFeatureEdgesTypes_H
#define pointFeatureEdgesTypes_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class pointFeatureEdgesTypes Declaration
\*---------------------------------------------------------------------------*/

//- Hold the types of feature edges attached to the point.
class pointFeatureEdgesTypes
:
    public HashTable<label, extendedFeatureEdgeMesh::edgeStatus>
{
    label pointLabel_;

public:

    pointFeatureEdgesTypes(const label pointLabel)
    :
        HashTable<label, extendedFeatureEdgeMesh::edgeStatus>(),
        pointLabel_(pointLabel)
    {}


    friend Ostream& operator<<(Ostream& os, const pointFeatureEdgesTypes& p)
    {
        os  << "Point = " << p.pointLabel_ << endl;

        for
        (
            HashTable<label, extendedFeatureEdgeMesh::edgeStatus>
                ::const_iterator iter = p.cbegin();
            iter != p.cend();
            ++iter
        )
        {
            os  << "    "
                << extendedFeatureEdgeMesh::edgeStatusNames_[iter.key()]
                << " = "
                << iter()
                << endl;
        }

        return os;
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

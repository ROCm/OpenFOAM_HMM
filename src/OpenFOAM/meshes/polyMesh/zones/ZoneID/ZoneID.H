/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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
    Foam::ZoneID

Description
    A class that holds the data needed to identify a zone in a dynamic mesh.

    The zone is identified by name.
    Its index in the zoneMesh is updated if the mesh has changed.

\*---------------------------------------------------------------------------*/

#ifndef ZoneID_H
#define ZoneID_H

#include "keyType.H"
#include "labelList.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class ZoneType, class MeshType> class ZoneMesh;

// Forward declaration of friend functions and operators

template<class ZoneType> class ZoneID;

template<class ZoneType>
Ostream& operator<<(Ostream&, const ZoneID<ZoneType>&);

/*---------------------------------------------------------------------------*\
                           Class ZoneID Declaration
\*---------------------------------------------------------------------------*/

template<class ZoneType>
class ZoneID
{
    // Private data

        //- Zone name
        keyType key_;

        //- Zone indices
        labelList indices_;


public:

    // Constructors

        //- Construct from name
        ZoneID(const keyType& key, const ZoneMesh<ZoneType, polyMesh>& zm)
        :
            key_(key),
            indices_(zm.findIndices(key_))
        {}

        //- Construct from Istream
        ZoneID(Istream& is, const ZoneMesh<ZoneType, polyMesh>& zm)
        :
            key_(is),
            indices_(zm.findIndices(key_))
        {}


    // Destructor - default


    // Member Functions

        // Access

            //- Return name
            const keyType& name() const
            {
                return key_;
            }

            //- Return indices of matching zones
            const labelList& indices() const
            {
                return indices_;
            }

            //- Return index of first matching zone
            label index() const
            {
                return indices_.empty() ? -1 : indices_[0];
            }

            //- Has the zone been found
            bool active() const
            {
                return !indices_.empty();
            }


        // Edit

            //- Update
            void update(const ZoneMesh<ZoneType, polyMesh>& zm)
            {
                indices_ = zm.findIndices(key_);
            }


    // IOstream Operators

        friend Ostream& operator<< <ZoneType>
        (Ostream&, const ZoneID<ZoneType>&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class ZoneType>
Ostream& operator<<(Ostream& os, const ZoneID<ZoneType>& p)
{
    os  << token::BEGIN_LIST
        << p.name() << token::SPACE << p.index()
        << token::END_LIST;

    // Check state of Ostream
    os.check("Ostream& operator<<(Ostream&, const ZoneID<ZoneType>&)");

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

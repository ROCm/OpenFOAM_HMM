/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2018 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2011-2014 OpenFOAM Foundation
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
    Foam::volumeType

Description
    An enumeration wrapper for classification of a location as being
    inside/outside of a volume.

SourceFiles
    volumeType.C

\*---------------------------------------------------------------------------*/

#ifndef volumeType_H
#define volumeType_H

#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class dictionary;
class volumeType;
Istream& operator>>(Istream& is, volumeType& vt);
Ostream& operator<<(Ostream& os, const volumeType& vt);


/*---------------------------------------------------------------------------*\
                         Class volumeType Declaration
\*---------------------------------------------------------------------------*/

class volumeType
{
public:

    //- Volume classification types
    enum type : char
    {
        UNKNOWN = 0,    //!< Unknown state
        INSIDE  = 0x1,  //!< A location inside the volume
        OUTSIDE = 0x2,  //!< A location outside the volume
        MIXED   = 0x3   //!< A location that is partly inside and outside
    };

    // Static data

    //- Names for the classification enumeration
    static const Enum<volumeType::type> names;


private:

    // Private data

        //- Volume type
        type t_;


public:

    // Constructors

        //- Construct null as UNKNOWN state
        volumeType()
        :
            t_(UNKNOWN)
        {}

        //- Construct from enumeration
        volumeType(type t)
        :
            t_(t)
        {}

        //- Construct as lookupOrDefault by name from dictionary
        volumeType(const word& key, const dictionary& dict, const type deflt);

        //- Construct from integer
        explicit volumeType(const int t)
        :
            t_(static_cast<volumeType::type>(t))
        {}


    // Member Functions

        //- Return the enumeration
        operator type() const
        {
            return t_;
        }

        //- The string representation of the volume type enumeration
        const word& str() const;


    // IOstream operators

        friend Istream& operator>>(Istream& is, volumeType& vt);
        friend Ostream& operator<<(Ostream& os, const volumeType& vt);
};


//- Data associated with volumeType type are contiguous
template<>
inline bool contiguous<volumeType>() {return true;}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

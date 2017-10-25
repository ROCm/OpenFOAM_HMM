/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

Class
    Foam::Map

Description
    A HashTable to objects of type \<T\> with a label key.

See also
    PtrMap

\*---------------------------------------------------------------------------*/

#ifndef Map_H
#define Map_H

#include "HashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class Map Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class Map
:
    public HashTable<T, label, Hash<label>>
{
public:

    //- The template instance used for this Map
    typedef Map<T> this_type;

    //- The template instance used for the parent HashTable
    typedef HashTable<T, label, Hash<label>> parent_type;

    using iterator = typename parent_type::iterator;
    using const_iterator = typename parent_type::const_iterator;


    // Constructors

        //- Construct given initial size
        Map(const label size = 128)
        :
            parent_type(size)
        {}

        //- Construct from Istream
        Map(Istream& is)
        :
            parent_type(is)
        {}

        //- Copy construct
        Map(const this_type& map)
        :
            parent_type(map)
        {}

        //- Move construct
        Map(this_type&& map)
        :
            parent_type(std::move(map))
        {}

        //- Construct by transferring the parameter contents
        Map(const Xfer<this_type>& map)
        :
            parent_type(map)
        {}

        //- Construct by transferring the parameter contents
        Map(const Xfer<HashTable<T, label, Hash<label>>>& map)
        :
            parent_type(map)
        {}

        //- Construct from an initializer list
        Map(std::initializer_list<std::pair<label, T>> map)
        :
            parent_type(map)
        {}


    // Member Operators

        using parent_type::operator=;

        //- Copy assignment
        void operator=(const this_type& rhs)
        {
            parent_type::operator=(rhs);
        }

        //- Move assignment
        void operator=(this_type&& rhs)
        {
            parent_type::operator=(std::move(rhs));
        }


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

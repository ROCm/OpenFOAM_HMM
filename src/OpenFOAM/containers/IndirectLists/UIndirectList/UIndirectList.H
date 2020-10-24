/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2019 OpenCFD Ltd.
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
    Foam::UIndirectList

Description
    A List with indirect addressing.

    Like IndirectList but does not store addressing.

    Note the const_cast of the list values. This is so we can use it both
    on const and non-const lists. Alternative would be to have a const_
    variant etc.

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef UIndirectList_H
#define UIndirectList_H

#include "List.H"
#include "IndirectListBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class T> class UIndirectList;

// Common list types
typedef UIndirectList<bool> boolUIndList;       //!< UIndirectList of bools
typedef UIndirectList<label> labelUIndList;     //!< UIndirectList of labels

/*---------------------------------------------------------------------------*\
                        Class UIndirectList Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class UIndirectList
:
    public IndirectListBase<T, labelUList>
{
public:

    // Constructors

        //- Shallow copy values and addressing
        UIndirectList(const UList<T>& values, const labelUList& addr)
        :
            IndirectListBase<T, labelUList>(values, addr)
        {}

        //- Copy construct (shallow copy values and addressing)
        UIndirectList(const UIndirectList<T>& list)
        :
            UIndirectList<T>(list.values(), list.addressing())
        {}


    // Member Operators

        //- Use standard assignment operations
        using IndirectListBase<T, labelUList>::operator=;

        //- Deep copy values, Fatal if list sizes are not identical
        void operator=(const UIndirectList<T>& rhs)
        {
            this->copyList(rhs);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

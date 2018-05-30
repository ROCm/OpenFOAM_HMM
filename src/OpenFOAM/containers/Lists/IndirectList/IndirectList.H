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
    Foam::IndirectList

Description
    A List with indirect addressing.

See also
    Foam::UIndirectList for a version without any allocation for the
    addressing.

SourceFiles
    IndirectListI.H

\*---------------------------------------------------------------------------*/

#ifndef IndirectList_H
#define IndirectList_H

#include "UIndirectList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class IndirectListAddressing Declaration
\*---------------------------------------------------------------------------*/

//- A helper class for storing addresses.
class IndirectListAddressing
{
    // Private data

        //- Storage for the list addressing
        List<label> addressing_;


    // Private Member Functions

        //- No copy construct
        IndirectListAddressing(const IndirectListAddressing&) = delete;

        //- No copy assignment
        void operator=(const IndirectListAddressing&) = delete;


protected:

    // Constructors

        //- Copy construct from addressing array
        inline explicit IndirectListAddressing(const labelUList& addr);

        //- Move construct from addressing array
        inline explicit IndirectListAddressing(List<label>&& addr);


    // Member Functions

        //- Return the list addressing
        inline const List<label>& addressing() const;

        //- Copy reset addressing
        inline void resetAddressing(const labelUList& addr);

        //- Move reset addressing
        inline void resetAddressing(List<label>&& addr);

};


/*---------------------------------------------------------------------------*\
                        Class IndirectList Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class IndirectList
:
    private IndirectListAddressing,
    public  UIndirectList<T>
{
    // Private Member Functions

        //- No copy assignment
        void operator=(const IndirectList<T>&) = delete;

        //- No copy assignment from UIndirectList
        void operator=(const UIndirectList<T>&) = delete;


public:

    // Constructors

        //- Copy construct addressing, shallow copy values list reference.
        inline IndirectList(const UList<T>& values, const labelUList& addr);

        //- Move construct addressing, shallow copy values list reference.
        inline IndirectList(const UList<T>& values, List<label>&& addr);

        //- Copy construct addressing, shallow copy values list reference.
        inline IndirectList(const IndirectList<T>& list);

        //- Copy construct addressing, shallow copy values list reference
        //- from UIndirectList
        inline explicit IndirectList(const UIndirectList<T>& list);


    // Member Functions

        //- Return the list addressing
        using UIndirectList<T>::addressing;

        //- Reset the list addressing
        using IndirectListAddressing::resetAddressing;


    // Member Operators

        //- Assignment operator
        using UIndirectList<T>::operator=;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "IndirectListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

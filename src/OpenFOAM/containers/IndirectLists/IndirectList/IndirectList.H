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
    Foam::IndirectList

Description
    A List with indirect addressing.

See also
    Foam::UIndirectList for a version without addressing allocation.

SourceFiles
    IndirectListI.H

\*---------------------------------------------------------------------------*/

#ifndef IndirectList_H
#define IndirectList_H

#include "IndirectListAddressing.H"
#include "UIndirectList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class IndirectList Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class IndirectList
:
    private IndirectListAddressing<labelList>,
    public UIndirectList<T>
{
public:

    // Constructors

        //- Copy construct addressing, shallow copy values reference
        inline IndirectList(const UList<T>& values, const labelUList& addr);

        //- Move construct addressing, shallow copy values reference
        inline IndirectList(const UList<T>& values, labelList&& addr);

        //- Zero-sized addressing, shallow copy values reference
        inline IndirectList(const UList<T>& values, const Foam::zero);

        //- Copy construct addressing, shallow copy values reference
        inline IndirectList(const IndirectList<T>& list);

        //- Move construct addressing, shallow copy values reference
        inline IndirectList(IndirectList<T>&& list);

        //- Copy construct addressing, shallow copy values reference
        inline explicit IndirectList(const UIndirectList<T>& list);


    // Member Functions

        //- The list addressing
        using IndirectListAddressing::addressing;


    // Member Operators

        //- Assignment operator
        using UIndirectList<T>::operator=;

        //- Deep copy values, Fatal if list sizes are not identical
        void operator=(const IndirectList<T>& rhs)
        {
            this->copyList(rhs);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "IndirectListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

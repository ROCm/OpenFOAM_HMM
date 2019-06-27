/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018-2019 OpenCFD Ltd.
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
    Foam::Detail::PtrListDetail

Description
    A rudimentary list of pointers used for PtrList, UPtrList, etc.
    This class is considered implementation detail and should not normally
    be used other than by OpenFOAM container classes.

SourceFiles
    PtrListDetailI.H
    PtrListDetail.C
    PtrListDetailIO.C

\*---------------------------------------------------------------------------*/

#ifndef PtrListDetail_H
#define PtrListDetail_H

#include "List.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class Ostream;

namespace Detail
{

/*---------------------------------------------------------------------------*\
                    Class Detail::PtrListDetail Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class PtrListDetail
:
    public List<T*>
{
public:

    // Constructors

        //- Construct null
        inline constexpr PtrListDetail() noexcept;

        //- Construct with specified size, each element initialized to nullptr
        inline PtrListDetail(const label len);

        //- Copy construct (shallow copies addresses)
        inline PtrListDetail(const PtrListDetail<T>& list);

        //- Move construct
        inline PtrListDetail(PtrListDetail<T>&& list);

        //- Copy or move (reuse) construct as specified
        inline PtrListDetail(PtrListDetail<T>& list, bool reuse);


    // Member Functions

        //- Return the count of non-nullptr entries
        label count() const;

        //- Locate the first null entry, -1 if there are not any
        label findNull() const;

        //- Delete the allocated entries, but retain the list size.
        void free();

        //- Make a copy by cloning each of the list pointers.
        template<class... Args>
        PtrListDetail<T> clone(Args&&... args) const;

        //- Reset size of list.
        //  New entries are initialized to nullptr.
        inline void resize(const label newLen);

        //- Override size to be inconsistent with allocated storage.
        //  Use with care
        inline void setAddressableSize(const label n);

        //- Write output, optionally silently trimming nullptrs
        Ostream& write(Ostream& os, const bool trimNull=false) const;


    // Member Operators

        //- Copy assignment (shallow copies addresses)
        inline void operator=(const PtrListDetail<T>& list);

        //- Move assignment
        inline void operator=(PtrListDetail<T>&& list);


    // Housekeeping

        // Just use resize().
        void setSize(const label) = delete;
        void setSize(const label, const T&) = delete;
        void setSize(const label, const T*) = delete;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Detail
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "PtrListDetailI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "PtrListDetail.C"
    #include "PtrListDetailIO.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //

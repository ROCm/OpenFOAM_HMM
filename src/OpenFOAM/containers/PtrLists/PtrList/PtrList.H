/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2019 OpenCFD Ltd.
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
    Foam::PtrList

Description
    A list of pointers to objects of type \<T\>, with allocation/deallocation
    management of the pointers.
    The operator[] returns a reference to the object, not the pointer.

See Also
    Foam::UPtrList
    Foam::PtrDynList

SourceFiles
    PtrListI.H
    PtrList.C
    PtrListIO.C

\*---------------------------------------------------------------------------*/

#ifndef PtrList_H
#define PtrList_H

#include "UPtrList.H"
#include "SLPtrListFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations

template<class T> class autoPtr;
template<class T> class tmp;
template<class T> class PtrList;

template<class T> Istream& operator>>(Istream& is, PtrList<T>& list);


/*---------------------------------------------------------------------------*\
                           Class PtrList Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class PtrList
:
    public UPtrList<T>
{
protected:

    // Protected Member Functions

        //- Read from Istream using Istream constructor class
        template<class INew>
        void readIstream(Istream& is, const INew& inew);

        //- Delete the allocated entries, but retain the list size.
        inline void free();

public:

    // Constructors

        //- Construct null
        inline constexpr PtrList() noexcept;

        //- Construct with specified size, each element initialized to nullptr
        inline explicit PtrList(const label len);

        //- Copy construct using 'clone()' method on each element
        inline PtrList(const PtrList<T>& list);

        //- Move construct
        inline PtrList(PtrList<T>&& list);

        //- Take ownerskip of pointers in the list, set old pointers to null.
        inline explicit PtrList(UList<T*>& list);

        //- Copy construct using 'clone()' method on each element
        template<class CloneArg>
        inline PtrList(const PtrList<T>& list, const CloneArg& cloneArgs);

        //- Construct as copy or re-use as specified
        PtrList(PtrList<T>& list, bool reuse);

        //- Copy construct using 'clone()' on each element of SLPtrList\<T\>
        explicit PtrList(const SLPtrList<T>& list);

        //- Construct from Istream using given Istream constructor class
        template<class INew>
        PtrList(Istream& is, const INew& inew);

        //- Construct from Istream using default Istream constructor class
        PtrList(Istream& is);


    //- Destructor
    ~PtrList();


    // Member Functions

        //- Make a copy by cloning each of the list elements.
        template<class... Args>
        PtrList<T> clone(Args&&... args) const;

        //- Clear the PtrList. Delete allocated entries and set size to zero.
        inline void clear();

        //- Adjust size of PtrList.
        //  New entries are initialized to nullptr, removed entries are deleted
        void resize(const label newLen);

        //- Same as resize()
        inline void setSize(const label newLen);

        //- Append an element to the end of the list
        inline void append(T* ptr);

        //- Move append an element to the end of the list
        inline void append(autoPtr<T>& aptr);

        //- Move append an element to the end of the list
        inline void append(autoPtr<T>&& aptr);

        //- Move or clone append a tmp to the end of the list
        inline void append(const tmp<T>& tptr);

        //- Transfer into this list and annul the argument list
        inline void transfer(PtrList<T>& list);

        //- Return const pointer to element (if set) or nullptr.
        //  The return value can be tested as a bool.
        inline const T* set(const label i) const;

        //- Set element to given pointer and return old element (can be null)
        //  No-op if the new pointer value is identical to the current content.
        inline autoPtr<T> set(const label i, T* ptr);

        //- Set element to given autoPtr and return old element
        inline autoPtr<T> set(const label i, autoPtr<T>& aptr);

        //- Set element to given autoPtr and return old element
        inline autoPtr<T> set(const label i, autoPtr<T>&& aptr);

        //- Set element to given tmp and return old element
        inline autoPtr<T> set(const label i, const tmp<T>& tptr);

        //- Release ownership of the pointer at the given position.
        //  Out of bounds addressing is a no-op and returns nullptr.
        inline autoPtr<T> release(const label i);


    // Member Operators

        //- Copy assignment.
        //  For existing list entries, values are copied from the list.
        //  For new list entries, pointers are cloned from the list.
        void operator=(const PtrList<T>& list);

        //- Move assignment
        inline void operator=(PtrList<T>&& list);


    // IOstream operator

        //- Read from Istream, discarding contents of existing list
        friend Istream& operator>> <T>(Istream& is, PtrList<T>& list);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "PtrListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "PtrList.C"
    #include "PtrListIO.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

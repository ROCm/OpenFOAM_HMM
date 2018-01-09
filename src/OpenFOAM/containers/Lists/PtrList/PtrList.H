/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2018 OpenCFD Ltd.
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
    Foam::PtrList

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

// Forward declarations

template<class T> class autoPtr;
template<class T> class tmp;
template<class T> class PtrList;

template<class T> Istream& operator>>(Istream& is, PtrList<T>& lst);


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

        //- Read from Istream using given Istream constructor class
        template<class INew>
        void read(Istream& is, const INew& inew);


public:

    // Constructors

        //- Construct null
        PtrList() = default;

        //- Construct with specified size, each element initialized to nullptr
        explicit inline PtrList(const label nElem);

        //- Copy construct using 'clone()' method on each element
        PtrList(const PtrList<T>& lst);

        //- Move construct
        inline PtrList(PtrList<T>&& lst);

        //- Copy construct with additional argument for 'clone()'
        template<class CloneArg>
        PtrList(const PtrList<T>& lst, const CloneArg& cloneArg);

        //- Construct by transferring the parameter contents
        PtrList(const Xfer<PtrList<T>>& lst);

        //- Construct as copy or re-use as specified
        PtrList(PtrList<T>& lst, bool reuse);

        //- Construct as copy of SLPtrList<T>
        explicit PtrList(const SLPtrList<T>& lst);

        //- Construct from Istream using given Istream constructor class
        template<class INew>
        PtrList(Istream& is, const INew& inew);

        //- Construct from Istream using default Istream constructor class
        PtrList(Istream& is);


    //- Destructor
    ~PtrList();


    // Member Functions

        //- Clear the PtrList. Set size to zero and delete allocated entries
        void clear();

        //- Reset size of PtrList.
        //  New entries are initialized to nullptr, removed entries are deleted
        void setSize(const label newLen);

        //- Reset size of PtrList.
        //  New entries are initialized to nullptr, removed entries are deleted
        inline void resize(const label newLen);

        //- Append an element at the end of the list
        using UPtrList<T>::append;

        //- Append an element at the end of the list
        inline void append(const autoPtr<T>& aptr);

        //- Append an element at the end of the list
        inline void append(const tmp<T>& tptr);

        //- Transfer into this list and annul the argument list
        inline void transfer(PtrList<T>& lst);

        //- Transfer contents to the Xfer container
        inline Xfer<PtrList<T>> xfer();

        //- Return true if element is set (ie, not a nullptr)
        inline bool set(const label i) const;

        //- Set element to given pointer and return old element (can be null)
        inline autoPtr<T> set(const label i, T* ptr);

        //- Set element to given autoPtr and return old element
        inline autoPtr<T> set(const label i, const autoPtr<T>& aptr);

        //- Set element to given tmp and return old element
        inline autoPtr<T> set(const label i, const tmp<T>& tptr);


    // Member operators

        //- Copy assignment.
        //  For existing list entries, values are copied from the list.
        //  For new list entries, pointers are cloned from the list.
        void operator=(const PtrList<T>& lst);

        //- Move assignment
        inline void operator=(PtrList<T>&& lst);


    // IOstream operator

        //- Read PtrList from Istream, discarding contents of existing PtrList
        friend Istream& operator>> <T>(Istream& is, PtrList<T>& lst);
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

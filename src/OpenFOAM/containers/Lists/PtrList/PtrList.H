/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::PtrList

Description
    A templated 1D list of pointers to objects of type \<T\>, where the
    size of the array is known and used for subscript bounds checking, etc.

    The element operator [] returns a reference to the object rather than a
    pointer.

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

        //- Null Constructor
        PtrList();

        //- Construct with size specified
        explicit PtrList(const label len);

        //- Copy constructor, using 'clone()' method on each element
        PtrList(const PtrList<T>& a);

        //- Copy constructor with additional argument for clone
        template<class CloneArg>
        PtrList(const PtrList<T>& a, const CloneArg& cloneArg);

        //- Construct by transferring the parameter contents
        PtrList(const Xfer<PtrList<T>>& lst);

        //- Construct as copy or re-use as specified
        PtrList(PtrList<T>& a, bool reuse);

        //- Construct as copy of SLPtrList<T>
        explicit PtrList(const SLPtrList<T>& sll);

        //- Construct from Istream using given Istream constructor class
        template<class INew>
        PtrList(Istream& is, const INew& inew);

        //- Construct from Istream using default Istream constructor class
        PtrList(Istream& is);


    //- Destructor
    ~PtrList();


    // Member functions

        // Edit

            //- Reset size of PtrList. If extending the PtrList, new entries are
            //  set to nullptr. If truncating the PtrList, removed entries are
            //  deleted
            void setSize(const label);

            //- Alias for setSize(const label)
            inline void resize(const label);

            //- Clear the PtrList, i.e. set size to zero deleting all the
            //  allocated entries
            void clear();

            //- Append an element at the end of the list
            inline void append(T*);
            inline void append(const autoPtr<T>&);
            inline void append(const tmp<T>&);

            //- Transfer the contents of the argument PtrList into this PtrList
            //  and annul the argument list
            void transfer(PtrList<T>&);

            //- Transfer contents to the Xfer container
            inline Xfer<PtrList<T>> xfer();

            //- Is element set
            inline bool set(const label) const;

            //- Set element to given T* and return old element (can be nullptr)
            inline autoPtr<T> set(const label, T*);

            //- Set element to given autoPtr<T> and return old element
            inline autoPtr<T> set(const label, const autoPtr<T>&);

            //- Set element to given tmp<T> and return old element
            inline autoPtr<T> set(const label, const tmp<T>&);

            //- Reorders elements. Ordering does not have to be done in
            //  ascending or descending order. Reordering has to be unique.
            //  (is shuffle)
            void reorder(const labelUList&);


    // Member operators

        //- Assignment
        void operator=(const PtrList<T>&);


    // IOstream operator

        //- Read PtrList from Istream, discarding contents of existing PtrList
        friend Istream& operator>> <T>(Istream&, PtrList<T>&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "PtrListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "PtrList.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

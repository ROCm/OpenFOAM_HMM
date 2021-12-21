/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2021 OpenCFD Ltd.
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
    Foam::PtrDynList

Description
    A dynamically resizable PtrList with allocation management.

See Also
    Foam::UPtrList
    Foam::PtrList

SourceFiles
    PtrDynListI.H

\*---------------------------------------------------------------------------*/

#ifndef PtrDynList_H
#define PtrDynList_H

#include "PtrList.H"
#include <type_traits>
#include <memory>
#include <utility>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class T, int SizeMin> class PtrDynList;

/*---------------------------------------------------------------------------*\
                         Class PtrDynList Declaration
\*---------------------------------------------------------------------------*/

template<class T, int SizeMin=64>
class PtrDynList
:
    public PtrList<T>
{
    static_assert(SizeMin > 0, "Invalid min size parameter");

    // Private Data

        //- The capacity (allocated size) of the list.
        label capacity_;


public:

    // Constructors

        //- Default construct
        inline constexpr PtrDynList() noexcept;

        //- Construct with given capacity.
        inline explicit PtrDynList(const label len);

        //- Copy construct using 'clone()' method on each element
        inline PtrDynList(const PtrDynList<T, SizeMin>& list);

        //- Move construct
        inline PtrDynList(PtrDynList<T, SizeMin>&& list);

        //- Take ownership of pointers in the list, set old pointers to null.
        inline explicit PtrDynList(UList<T*>& list);


    //- Destructor
    ~PtrDynList() = default;


    // Member Functions

    // Access

        //- Size of the underlying storage.
        inline label capacity() const noexcept;

        //- Return const pointer to element (can be nullptr),
        //  with bounds checking.
        inline const T* get(const label i) const;

        //- Return const pointer to element (if set) or nullptr,
        //  with bounds checking.
        //  The return value can be tested as a bool.
        const T* set(const label i) const { return this->get(i); }


    // Sizing

        //- Reserve allocation space for at least this size.
        inline void reserve(const label len);

        //- Alter the addressed list size.
        inline void resize(const label newLen);

        //- Same as resize()
        void setSize(const label n) { this->resize(n); }

        //- Clear the addressed list, i.e. set the size to zero.
        //  Allocated size does not change
        inline void clear();

        //- Clear the list and delete storage.
        inline void clearStorage();

        //- Expand the addressable size to fit the allocated capacity.
        //  Returns the previous addressable size.
        inline label expandStorage() noexcept;

        //- Shrink the allocated space to the number of elements used.
        inline void shrink();


    // Edit

        //- Delete the allocated entries, but retain the list size.
        using PtrList<T>::free;

        //- Squeeze out intermediate nullptr entries in the list of pointers
        //- and adjust the addressable size accordingly.
        //  \return the number of non-null entries
        inline label squeezeNull();

        //- Swap content, independent of sizing parameter
        template<int AnySizeMin>
        inline void swap(PtrDynList<T, AnySizeMin>& other);

        //- Construct and append an element to the end of the list
        template<class... Args>
        inline void emplace_append(Args&&... args);

        //- Append an element to the end of the list
        inline void append(T* ptr);

        //- Move append an element to the end of the list
        inline void append(autoPtr<T>& ptr);

        //- Move append an element to the end of the list
        inline void append(autoPtr<T>&& ptr);

        //- Move append an element to the end of the list
        inline void append(std::unique_ptr<T>&& ptr);

        //- Move or clone append a tmp to the end of the list
        inline void append(const refPtr<T>& ptr);

        //- Move or clone append a tmp to the end of the list
        inline void append(const tmp<T>& ptr);

        //- Move append another list to the end of this list.
        inline void append(PtrList<T>&& other);

        //- Move append another list to the end of this list.
        template<int AnySizeMin>
        inline void append(PtrDynList<T, AnySizeMin>&& other);

        //- Remove and return the top element
        inline autoPtr<T> remove();

        //- Construct and set an element
        template<class... Args>
        inline autoPtr<T> emplace(const label i, Args&&... args);

        //- Set element to given pointer and return old element (can be null)
        inline autoPtr<T> set(const label i, T* ptr);

        //- Set element to given autoPtr and return old element
        inline autoPtr<T> set(const label i, autoPtr<T>& ptr);

        //- Set element to given autoPtr and return old element
        inline autoPtr<T> set(const label i, autoPtr<T>&& ptr);

        //- Set element to given pointer and return old element
        inline autoPtr<T> set(const label i, std::unique_ptr<T>&& ptr);

        //- Set element to given refPtr and return old element
        inline autoPtr<T> set(const label i, const refPtr<T>& ptr);

        //- Set element to given tmp and return old element
        inline autoPtr<T> set(const label i, const tmp<T>& ptr);

        //- Reorder elements. Reordering must be unique (ie, shuffle).
        inline void reorder(const labelUList& oldToNew);


    // Member Operators

        //- Copy (clone) assignment
        inline void operator=(const PtrList<T>& list);

        //- Copy (clone) assignment
        inline void operator=(const PtrDynList<T, SizeMin>& list);

        //- Copy (clone) assignment with different sizing parameters
        template<int AnySizeMin>
        inline void operator=(const PtrDynList<T, AnySizeMin>& list);

        //- Move assignment
        inline void operator=(PtrList<T>&& list);

        //- Move assignment
        inline void operator=(PtrDynList<T, SizeMin>&& list);

        //- Move assignment with different sizing parameters
        template<int AnySizeMin>
        inline void operator=(PtrDynList<T, AnySizeMin>&& list);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "PtrDynListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

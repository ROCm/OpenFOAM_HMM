/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2018-2020 OpenCFD Ltd.
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
    Foam::CompactIOList

Description
    A List of objects of type \<T\> with automated input and output using
    a compact storage. Behaves like IOList except when binary output in
    case it writes a CompactListList.

    Useful for lists of small sublists e.g. faceList, cellList.

SourceFiles
    CompactIOList.C

\*---------------------------------------------------------------------------*/

#ifndef CompactIOList_H
#define CompactIOList_H

#include "IOList.H"
#include "regIOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class Istream;
class Ostream;
template<class T, class BaseType> class CompactIOList;

template<class T, class BaseType> Istream& operator>>
(
    Istream&,
    CompactIOList<T, BaseType>&
);
template<class T, class BaseType> Ostream& operator<<
(
    Ostream&,
    const CompactIOList<T, BaseType>&
);

/*---------------------------------------------------------------------------*\
                           Class CompactIOList Declaration
\*---------------------------------------------------------------------------*/

template<class T, class BaseType>
class CompactIOList
:
    public regIOobject,
    public List<T>
{
    // Private Member Functions

        //- Read according to header type
        void readFromStream();

        //- Has too many elements in it?
        bool overflows() const;

public:

    //- Runtime type information
    TypeName("CompactList");


    // Constructors

        //- Default copy construct
        CompactIOList(const CompactIOList&) = default;

        //- Construct from IOobject
        explicit CompactIOList(const IOobject& io);

        //- Construct from IOobject and default length of CompactIOList
        CompactIOList(const IOobject& io, const label len);

        //- Construct from IOobject and List content
        CompactIOList(const IOobject& io, const UList<T>& content);

        //- Construct by transferring the List content
        CompactIOList(const IOobject& io, List<T>&& content);

    // Destructor

        virtual ~CompactIOList() = default;


    // Member Functions

        //- Write using stream options. Checks for overflow in binary
        virtual bool writeObject
        (
            IOstreamOption streamOpt,
            const bool valid
        ) const;

        virtual bool writeData(Ostream&) const;


    // Member Operators

        //- Copy assignment of entries
        void operator=(const CompactIOList<T, BaseType>& rhs);

        //- Copy or move assignment of entries
        using List<T>::operator=;


    // IOstream operators

        //- Read List from Istream, discarding contents of existing List.
        friend Istream& operator>> <T, BaseType>
        (
            Istream&,
            CompactIOList<T, BaseType>&
        );

        // Write List to Ostream.
        friend Ostream& operator<< <T, BaseType>
        (
            Ostream&,
            const CompactIOList<T, BaseType>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "CompactIOList.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

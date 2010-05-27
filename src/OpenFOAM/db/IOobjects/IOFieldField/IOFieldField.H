/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2010-2010 OpenCFD Ltd.
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
    Foam::IOFieldField

Description
    A Field of objects of type \<T\> with automated input and output.

SourceFiles
    IOFieldField.C

\*---------------------------------------------------------------------------*/

#ifndef IOFieldField_H
#define IOFieldField_H

#include "IOField.H"
#include "regIOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class Istream;
class Ostream;

// Forward declaration of friend functions and operators
template<class T, class BaseType> class IOFieldField;

template<class T, class BaseType> Istream& operator>>
(
    Istream&,
    IOFieldField<T, BaseType>&
);
template<class T, class BaseType> Ostream& operator<<
(
    Ostream&,
    const IOFieldField<T, BaseType>&
);

/*---------------------------------------------------------------------------*\
                        Class IOFieldField Declaration
\*---------------------------------------------------------------------------*/

template<class T, class BaseType>
class IOFieldField
:
    public regIOobject,
    public Field<T>
{
    // Private Member Functions

        //- Read according to header type
        void readFromStream();

public:

    //- Runtime type information
    TypeName("FieldField");


    // Constructors

        //- Construct from IOobject
        IOFieldField(const IOobject&);

        //- Construct from IOobject and size of IOFieldField
        IOFieldField(const IOobject&, const label);

        //- Construct from IOobject and a Field
        IOFieldField(const IOobject&, const Field<T>&);

        //- Construct by transferring the Field contents
        IOFieldField(const IOobject&, const Xfer<Field<T> >&);


    // Destructor

        virtual ~IOFieldField();


    // Member functions

        virtual bool writeObject
        (
            IOstream::streamFormat,
            IOstream::versionNumber,
            IOstream::compressionType
        ) const;

        virtual bool writeData(Ostream&) const;


    // Member operators

        void operator=(const IOFieldField<T, BaseType>&);

        void operator=(const Field<T>&);


    // IOstream operators

        //- Read Field from Istream, discarding contents of existing Field.
        friend Istream& operator>> <T, BaseType>
        (
            Istream&,
            IOFieldField<T, BaseType>&
        );

        // Write Field to Ostream.
        friend Ostream& operator<< <T, BaseType>
        (
            Ostream&,
            const IOFieldField<T, BaseType>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "IOFieldField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

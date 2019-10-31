/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2018 OpenCFD Ltd.
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
    Foam::IOField

Description
    A primitive field of type \<T\> with automated input and output.

SourceFiles
    IOField.C

\*---------------------------------------------------------------------------*/

#ifndef IOField_H
#define IOField_H

#include "regIOobject.H"
#include "Field.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class IOField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class IOField
:
    public regIOobject,
    public Field<Type>
{
public:

    TypeName("Field");


    // Constructors

        //- Default copy construct
        IOField(const IOField&) = default;

        //- Construct from IOobject
        explicit IOField(const IOobject& io);

        //- Construct from IOobject; does local processor require reading?
        IOField(const IOobject& io, const bool valid);

        //- Construct from IOobject and size (does not set values)
        IOField(const IOobject& io, const label size);

        //- Construct from IOobject and a List/Field content
        IOField(const IOobject& io, const UList<Type>& content);

        //- Construct by transferring the Field content
        IOField(const IOobject& io, Field<Type>&& content);

        //- Construct by copying/moving tmp content
        IOField(const IOobject& io, const tmp<Field<Type>>& tfld);


    //- Destructor
    virtual ~IOField() = default;


    // Member Functions

        bool writeData(Ostream& os) const;


    // Member Operators

        //- Copy assignment of entries
        void operator=(const IOField<Type>& rhs);

        //- Copy or move assignment of entries
        using Field<Type>::operator=;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "IOField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

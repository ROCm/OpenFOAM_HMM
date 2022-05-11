/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2022 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "IOField.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
bool Foam::IOField<Type>::readContents()
{
    if
    (
        (
            readOpt() == IOobject::MUST_READ
         || readOpt() == IOobject::MUST_READ_IF_MODIFIED
        )
     || (readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    )
    {
        readStream(typeName) >> *this;
        close();
        return true;
    }

    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::IOField<Type>::IOField(const IOobject& io)
:
    regIOobject(io)
{
    // Check for MUST_READ_IF_MODIFIED
    warnNoRereading<IOField<Type>>();

    readContents();
}


template<class Type>
Foam::IOField<Type>::IOField(const IOobject& io, const bool valid)
:
    regIOobject(io)
{
    // Check for MUST_READ_IF_MODIFIED
    warnNoRereading<IOField<Type>>();

    if
    (
        io.readOpt() == IOobject::MUST_READ
     || io.readOpt() == IOobject::MUST_READ_IF_MODIFIED
    )
    {
        Istream& is = readStream(typeName, valid);

        if (valid)
        {
            is >> *this;
        }
        close();
    }
    else if (io.readOpt() == IOobject::READ_IF_PRESENT)
    {
        bool haveFile = headerOk();

        Istream& is = readStream(typeName, haveFile && valid);

        if (valid && haveFile)
        {
            is >> *this;
        }
        close();
    }
}


template<class Type>
Foam::IOField<Type>::IOField(const IOobject& io, Foam::zero)
:
    regIOobject(io)
{
    // Check for MUST_READ_IF_MODIFIED
    warnNoRereading<IOField<Type>>();

    readContents();
}


template<class Type>
Foam::IOField<Type>::IOField(const IOobject& io, const label len)
:
    regIOobject(io)
{
    // Check for MUST_READ_IF_MODIFIED
    warnNoRereading<IOField<Type>>();

    if (!readContents())
    {
        Field<Type>::resize(len);
    }
}


template<class Type>
Foam::IOField<Type>::IOField(const IOobject& io, const UList<Type>& content)
:
    regIOobject(io)
{
    // Check for MUST_READ_IF_MODIFIED
    warnNoRereading<IOField<Type>>();

    if (!readContents())
    {
        Field<Type>::operator=(content);
    }
}


template<class Type>
Foam::IOField<Type>::IOField(const IOobject& io, Field<Type>&& content)
:
    regIOobject(io)
{
    // Check for MUST_READ_IF_MODIFIED
    warnNoRereading<IOField<Type>>();

    Field<Type>::transfer(content);

    readContents();
}


template<class Type>
Foam::IOField<Type>::IOField(const IOobject& io, const tmp<Field<Type>>& tfld)
:
    regIOobject(io)
{
    const bool reuse = tfld.movable();

    if (reuse)
    {
        Field<Type>::transfer(tfld.ref());
    }

    if (!readContents() && !reuse)
    {
        Field<Type>::operator=(tfld());
    }

    tfld.clear();
}


template<class Type>
Foam::IOFieldRef<Type>::IOFieldRef
(
    const IOobject& io,
    const Field<Type>& content
)
:
    regIOobject(io),
    contentRef_(content)  // cref
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::IOField<Type>::writeData(Ostream& os) const
{
    os << static_cast<const Field<Type>&>(*this);
    return os.good();
}


template<class Type>
bool Foam::IOFieldRef<Type>::writeData(Ostream& os) const
{
    os << contentRef_.cref();
    return os.good();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::IOField<Type>::operator=(const IOField<Type>& rhs)
{
    Field<Type>::operator=(rhs);
}


// ************************************************************************* //

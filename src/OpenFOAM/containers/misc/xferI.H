/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T>
Foam::xfer<T>::xfer()
:
    ref_()
{}


template<class T>
Foam::xfer<T>::xfer(T& t)
:
    ref_()
{
    transfer(t);
}


template<class T>
Foam::xfer<T>::xfer(xfer<T>& t)
:
    ref_()
{
    transfer(t);
}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class T>
Foam::xfer<T>::~xfer()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
T& Foam::xfer<T>::ref()
{
    return ref_;
}


template<class T>
void Foam::xfer<T>::transfer(T& t)
{
    ref_.transfer(t);
}


template<class T>
void Foam::xfer<T>::transfer(xfer<T>& t)
{
    ref_.transfer(t.ref_);
}


template<class T>
void Foam::xfer<T>::yield(T& t)
{
    t.transfer(ref_);
}


template<class T>
void Foam::xfer<T>::yield(xfer<T>& t)
{
    t.ref_.transfer(ref_);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


template<class T>
void Foam::xfer<T>::operator=(T& t)
{
    transfer(t);
}


template<class T>
void Foam::xfer<T>::operator=(xfer<T>& t)
{
    if (this == &t)
    {
        FatalErrorIn
        (
            "Foam::xfer<T>::operator="
            "(const Foam::xfer<T>&)"
        )   << "Attempted assignment to self"
            << abort(FatalError);
    }

    transfer(t);
}


template<class T>
Foam::xfer<T>& Foam::xfer<T>::operator()() const
{
    return const_cast<xfer<T>&>(*this);
}



// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class T>
Foam::Istream& Foam::operator>>(Istream& is, xfer<T>& t)
{
    is >> t.ref_;
    return is;
}


template<class T>
Foam::Ostream& Foam::operator<<(Ostream& os, const xfer<T>& t)
{
    os << t.ref_;
    return os;
}

// ************************************************************************* //

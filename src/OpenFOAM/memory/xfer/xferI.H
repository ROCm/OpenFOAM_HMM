/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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
inline Foam::xfer<T>::xfer(T* p)
:
    ptr_(p ? p : new T)
{}


template<class T>
inline Foam::xfer<T>::xfer(T& t, bool allowTransfer)
:
    ptr_(new T)
{
    if (allowTransfer)
    {
        ptr_->transfer(t);
    }
    else
    {
        ptr_->operator=(t);
    }
}


template<class T>
inline Foam::xfer<T>::xfer(const T& t)
:
    ptr_(new T)
{
    ptr_->operator=(t);
}


template<class T>
inline Foam::xfer<T>::xfer(const xfer<T>& t)
:
    ptr_(new T)
{
    ptr_->transfer(*(t.ptr_));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class T>
inline Foam::xfer<T>::~xfer()
{
    delete ptr_;
    ptr_ = 0;
}


// * * * * * * * * * * * * *  Member Functions * * * * * * * * * * * * * * * //

template<class T>
inline const Foam::xfer<T>& Foam::xfer<T>::null()
{
    xfer<T>* nullPtr = reinterpret_cast<xfer<T>*>(0);
    return *nullPtr;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
inline void Foam::xfer<T>::operator=(T& t)
{
    ptr_->transfer(t);
}


template<class T>
inline void Foam::xfer<T>::operator=(const xfer<T>& t)
{
    // silently ignore attempted copy to self
    if (this != &t)
    {
        ptr_->transfer(*(t.ptr_));
    }
}


template<class T>
inline T& Foam::xfer<T>::operator()() const
{
    return *ptr_;
}


template<class T>
inline T* Foam::xfer<T>::operator->() const
{
    return ptr_;
}


// * * * * * * * * * * * * *  Helper Functions * * * * * * * * * * * * * * * //


template<class T>
inline Foam::xfer<T> Foam::xferCopy(const T& t)
{
    return Foam::xfer<T>(t);
}

template<class T>
inline Foam::xfer<T> Foam::xferMove(T& t)
{
    return Foam::xfer<T>(t, true);
}


template<class T>
inline Foam::xfer<T> Foam::xferTmp(Foam::tmp<T>& tt)
{
    return Foam::xfer<T>(tt(), tt.isTmp());
}


template<class To, class From>
inline Foam::xfer<To> Foam::xferCopyTo(const From& t)
{
    Foam::xfer<To> xf;
    xf() = t;
    return xf;
}


template<class To, class From>
inline Foam::xfer<To> Foam::xferMoveTo(From& t)
{
    Foam::xfer<To> xf;
    xf().transfer(t);
    return xf;
}


// ************************************************************************* //

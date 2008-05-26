/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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

#include "timeSeries.H"
#include "Istream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<typename T>
Foam::timeSeries<T>::timeSeries(const bounds bound)
:
    List<Tuple2<scalar, T> >(),
    bounding_(bound)
{}


template<typename T>
Foam::timeSeries<T>::timeSeries(const word& bound)
:
    List<Tuple2<scalar, T> >(),
    bounding_(timeSeries::WARN)
{
    bounding(bound);
}


template<typename T>
Foam::timeSeries<T>::timeSeries(Istream& is, const bounds bound)
:
    List<Tuple2<scalar, T> >(is),
    bounding_(bound)
{}


template<typename T>
Foam::timeSeries<T>::timeSeries(Istream& is, const word& bound)
:
    List<Tuple2<scalar, T> >(is),
    bounding_(timeSeries::WARN)
{
    bounding(bound);
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<typename T>
Foam::timeSeries<T>::~timeSeries()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<typename T>
Foam::word Foam::timeSeries<T>::bounding() const
{
    word enumName("warn");

    switch (bounding_)
    {
        case timeSeries::ERROR:
            enumName = "error";
            break;

        case timeSeries::WARN:
            enumName = "warn";
            break;

        case timeSeries::CLAMP:
            enumName = "clamp";
            break;

        case timeSeries::REPEAT:
            enumName = "repeat";
            break;
    }

    return enumName;
}


template<typename T>
void Foam::timeSeries<T>::bounding(const word& bound)
{
    if (bound == "error")
    {
        bounding_ = timeSeries::ERROR;
    }
    else if (bound == "warn")
    {
        bounding_ = timeSeries::WARN;
    }
    else if (bound == "clamp")
    {
        bounding_ = timeSeries::CLAMP;
    }
    else if (bound == "repeat")
    {
        bounding_ = timeSeries::REPEAT;
    }
    else
    {
        WarningIn("Foam::timeSeries<T>::boundingEnum(const word&)")
            << "bad bounding specifier " << bound << " using 'warn'" << endl;

        bounding_ = timeSeries::WARN;
    }
}


template<typename T>
void Foam::timeSeries<T>::check() const
{
    label n = size();
    scalar prevTime = List<Tuple2<scalar, T> >::operator[](0).first();

    for (label i = 1; i < n; ++i)
    {
        const scalar currTime = List<Tuple2<scalar, T> >::operator[](i).first();

        // avoid duplicate times (divide-by-zero error)
        if (currTime <= prevTime)
        {
            FatalErrorIn
            (
                "Foam::timeSeries<T>::checkOrder() const"
            )   << "out-of-order time: "
                << currTime << " at index " << i << nl
                << exit(FatalError);
        }
        prevTime = currTime;
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


template<typename T>
const Foam::Tuple2<Foam::scalar, T>&
Foam::timeSeries<T>::operator[](const label i) const
{
    label ii = i;
    label n  = size();

    if (n <= 1)
    {
        ii = 0;
    }
    else if (ii < 0)
    {
        switch (bounding_)
        {
            case timeSeries::ERROR:
                FatalErrorIn
                (
                    "Foam::timeSeries<T>::operator[](const label) const"
                )   << "index (" << ii << ") underflow" << nl
                    << exit(FatalError);
                break;

            case timeSeries::WARN:
                WarningIn
                (
                    "Foam::timeSeries<T>::operator[](const label) const"
                )   << "index (" << ii << ") underflow" << nl
                    << "    Continuing with the first entry"
                    << endl;
                // fall-through to 'CLAMP'

            case timeSeries::CLAMP:
                ii = 0;
                break;

            case timeSeries::REPEAT:
                while (ii < 0)
                {
                    ii += n;
                }
                break;
        }
    }
    else if (ii >= n)
    {
        switch (bounding_)
        {
            case timeSeries::ERROR:
                FatalErrorIn
                (
                    "Foam::timeSeries<T>::operator[](const label) const"
                )   << "index (" << ii << ") overflow" << nl
                    << exit(FatalError);
                break;

            case timeSeries::WARN:
                WarningIn
                (
                    "Foam::timeSeries<T>::operator[](const label) const"
                )   << "index (" << ii << ") overflow" << nl
                    << "    Continuing with the last entry"
                    << endl;
                // fall-through to 'CLAMP'

            case timeSeries::CLAMP:
                ii = n - 1;
                break;

            case timeSeries::REPEAT:
                while (ii >= n)
                {
                    ii -= n;
                }
                break;
        }
    }

    return List<Tuple2<scalar, T> >::operator[](ii);
}


template<typename T>
T Foam::timeSeries<T>::operator()(const scalar timeValue) const
{
    label n = size();

    if (n <= 1)
    {
        return List<Tuple2<scalar, T> >::operator[](0).second();
    }

    scalar minTime = List<Tuple2<scalar, T> >::operator[](0).first();
    scalar maxTime = List<Tuple2<scalar, T> >::operator[](n-1).first();
    scalar lookupTime = timeValue;

    if (lookupTime < minTime)
    {
        switch (bounding_)
        {
            case timeSeries::ERROR:
                FatalErrorIn
                (
                    "Foam::timeSeries<T>::operator[](const scalar) const"
                )   << "time (" << lookupTime << ") underflow" << nl
                    << exit(FatalError);
                break;

            case timeSeries::WARN:
                WarningIn
                (
                    "Foam::timeSeries<T>::operator[](const scalar) const"
                )   << "time (" << lookupTime << ") underflow" << nl
                    << "    Continuing with the first entry"
                    << endl;
                // fall-through to 'CLAMP'

            case timeSeries::CLAMP:
                return List<Tuple2<scalar, T> >::operator[](0).second();
                break;

            case timeSeries::REPEAT:
                // adjust lookupTime to >= 0
                while (lookupTime < 0)
                {
                    lookupTime += maxTime;
                }
                break;
        }
    }
    else if (lookupTime >= maxTime)
    {
        switch (bounding_)
        {
            case timeSeries::ERROR:
                FatalErrorIn
                (
                    "Foam::timeSeries<T>::operator[](const label) const"
                )   << "time (" << lookupTime << ") overflow" << nl
                    << exit(FatalError);
                break;

            case timeSeries::WARN:
                WarningIn
                (
                    "Foam::timeSeries<T>::operator[](const label) const"
                )   << "time (" << lookupTime << ") overflow" << nl
                    << "    Continuing with the last entry"
                    << endl;
                // fall-through to 'CLAMP'

            case timeSeries::CLAMP:
                return List<Tuple2<scalar, T> >::operator[](n-1).second();
                break;

            case timeSeries::REPEAT:
                // adjust lookupTime <= maxTime
                while (lookupTime > maxTime)
                {
                    lookupTime -= maxTime;
                }
                break;
        }
    }

    label lo = 0;
    label hi = 0;

    // look for the correct range
    for (label i = 0; i < n; ++i)
    {
        if (lookupTime >= List<Tuple2<scalar, T> >::operator[](i).first())
        {
            lo = hi = i;
        }
        else
        {
            hi = i;
            break;
        }
    }

    if (lo == hi)
    {
        // we are at the end of the table - or there is only a single entry
        return List<Tuple2<scalar, T> >::operator[](hi).second();
    }
    else if (hi == 0)
    {
        // this treatment should should only occur under these condition:
        //  -> the 'REPEAT' treatment
        //  -> (0 <= time <= minTime)
        //  -> minTime > 0
        // Use the value at maxTime as the value for time=0
        lo = n - 1;

        return
        (
            List<Tuple2<scalar, T> >::operator[](lo).second()
          +
            (
                List<Tuple2<scalar, T> >::operator[](hi).second()
              - List<Tuple2<scalar, T> >::operator[](lo).second()
            )
          * (lookupTime / minTime)
        );
    }
    else
    {
        // normal interpolation
        return
        (
            List<Tuple2<scalar, T> >::operator[](lo).second()
          +
            (
                List<Tuple2<scalar, T> >::operator[](hi).second()
              - List<Tuple2<scalar, T> >::operator[](lo).second()
            )
          *
            (
                lookupTime
              - List<Tuple2<scalar, T> >::operator[](lo).first()
            )
          /
            (
                List<Tuple2<scalar, T> >::operator[](hi).first()
              - List<Tuple2<scalar, T> >::operator[](lo).first()
            )
        );
    }
}


// ************************************************************************* //

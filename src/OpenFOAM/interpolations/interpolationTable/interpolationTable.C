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

#include "interpolationTable.H"
#include "IFstream.H"
#include "objectRegistry.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::interpolationTable<Type>::interpolationTable()
:
    List<Tuple2<scalar, Type> >(),
    dict_(dictionary::null),
    boundAction_(interpolationTable::WARN),
    fileName_("undefined_fileName")
{}


template<class Type>
Foam::interpolationTable<Type>::interpolationTable
(
    const objectRegistry& obr,
    const dictionary& dict
)
:
    List<Tuple2<scalar, Type> >(),
    dict_(dict),
    boundAction_(wordToBoundAction(dict.lookup("boundAction"))),
    fileName_(dict.lookup("fileName"))
{
    fileName_.expand();

    // Correct for relative path
    if (fileName_[0] != '/')
    {
        fileName_ = obr.db().path()/fileName_;
    }

    // Read data from file
    IFstream(fileName_)() >> *this;

    // Check that the data is okay
    check();

    if (this->size() == 0)
    {
        FatalErrorIn
        (
            "Foam::interpolationTable<Type>::interpolationTable"
            "(const dictionary& dict)"
        )   << "table is empty" << nl
            << exit(FatalError);
    }
}


template<class Type>
Foam::interpolationTable<Type>::interpolationTable
(
     const interpolationTable& interpTable
)
:
    List<Tuple2<scalar, Type> >(interpTable),
    dict_(interpTable.dict_),
    boundAction_(interpTable.boundAction_),
    fileName_(interpTable.fileName_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::interpolationTable<Type>::~interpolationTable()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::word Foam::interpolationTable<Type>::boundActionToWord
(
     const boundActions& bound
) const
{
    word enumName("warn");

    switch (bound)
    {
        case interpolationTable::ERROR:
        {
            enumName = "error";
            break;
        }
        case interpolationTable::WARN:
        {
            enumName = "warn";
            break;
        }
        case interpolationTable::CLAMP:
        {
            enumName = "clamp";
            break;
        }
        case interpolationTable::REPEAT:
        {
            enumName = "repeat";
            break;
        }
    }

    return enumName;
}


template<class Type>
typename Foam::interpolationTable<Type>::boundActions
Foam::interpolationTable<Type>::wordToBoundAction
(
    const word& bound
) const
{
    if (bound == "error")
    {
        return interpolationTable::ERROR;
    }
    else if (bound == "warn")
    {
        return interpolationTable::WARN;
    }
    else if (bound == "clamp")
    {
        return interpolationTable::CLAMP;
    }
    else if (bound == "repeat")
    {
        return interpolationTable::REPEAT;
    }
    else
    {
        WarningIn
        (
            "Foam::interpolationTable<Type>::wordToBoundAction(const word&)"
        )   << "bad bounding specifier " << bound << " using 'warn'" << endl;

        return interpolationTable::WARN;
    }
}


template<class Type>
void Foam::interpolationTable<Type>::check() const
{
    label n = size();
    scalar prevValue = List<Tuple2<scalar, Type> >::operator[](0).first();

    for (label i=1; i<n; ++i)
    {
        const scalar currValue =
            List<Tuple2<scalar, Type> >::operator[](i).first();

        // avoid duplicate values (divide-by-zero error)
        if (currValue <= prevValue)
        {
            FatalErrorIn
            (
                "Foam::interpolationTable<Type>::checkOrder() const"
            )   << "out-of-order value: "
                << currValue << " at index " << i << nl
                << exit(FatalError);
        }
        prevValue = currValue;
    }
}


template<class Type>
typename Foam::interpolationTable<Type>::boundActions
Foam::interpolationTable<Type>::boundAction
(
    const boundActions& bound
)
{
    boundActions prev = boundAction_;
    boundAction_ = bound;
    return prev;
}


template<class Type>
void Foam::interpolationTable<Type>::write(Ostream& os) const
{
    os.writeKeyword("fileName")
        << fileName_ << token::END_STATEMENT << nl;
    os.writeKeyword("boundAction")
        << boundActionToWord(boundAction_) << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
const Foam::Tuple2<Foam::scalar, Type>&
Foam::interpolationTable<Type>::operator[](const label i) const
{
    label ii = i;
    label n  = size();

    if (n <= 1)
    {
        ii = 0;
    }
    else if (ii < 0)
    {
        switch (boundAction_)
        {
            case interpolationTable::ERROR:
            {
                FatalErrorIn
                (
                    "Foam::interpolationTable<Type>::operator[]"
                    "(const label) const"
                )   << "index (" << ii << ") underflow" << nl
                    << exit(FatalError);
                break;
            }
            case interpolationTable::WARN:
            {
                WarningIn
                (
                    "Foam::interpolationTable<Type>::operator[]"
                    "(const label) const"
                )   << "index (" << ii << ") underflow" << nl
                    << "    Continuing with the first entry"
                    << endl;
                // fall-through to 'CLAMP'
            }
            case interpolationTable::CLAMP:
            {
                ii = 0;
                break;
            }
            case interpolationTable::REPEAT:
            {
                while (ii < 0)
                {
                    ii += n;
                }
                break;
            }
        }
    }
    else if (ii >= n)
    {
        switch (boundAction_)
        {
            case interpolationTable::ERROR:
            {
                FatalErrorIn
                (
                    "Foam::interpolationTable<Type>::operator[]"
                    "(const label) const"
                )   << "index (" << ii << ") overflow" << nl
                    << exit(FatalError);
                break;
            }
            case interpolationTable::WARN:
            {
                WarningIn
                (
                    "Foam::interpolationTable<Type>::operator[]"
                    "(const label) const"
                )   << "index (" << ii << ") overflow" << nl
                    << "    Continuing with the last entry"
                    << endl;
                // fall-through to 'CLAMP'
            }
            case interpolationTable::CLAMP:
            {
                ii = n - 1;
                break;
            }
            case interpolationTable::REPEAT:
            {
                while (ii >= n)
                {
                    ii -= n;
                }
                break;
            }
        }
    }

    return List<Tuple2<scalar, Type> >::operator[](ii);
}


template<class Type>
Type Foam::interpolationTable<Type>::operator()(const scalar value) const
{
    label n = size();

    if (n <= 1)
    {
        return List<Tuple2<scalar, Type> >::operator[](0).second();
    }

    scalar minLimit = List<Tuple2<scalar, Type> >::operator[](0).first();
    scalar maxLimit = List<Tuple2<scalar, Type> >::operator[](n-1).first();
    scalar lookupValue = value;

    if (lookupValue < minLimit)
    {
        switch (boundAction_)
        {
            case interpolationTable::ERROR:
            {
                FatalErrorIn
                (
                    "Foam::interpolationTable<Type>::operator[]"
                    "(const scalar) const"
                )   << "value (" << lookupValue << ") underflow" << nl
                    << exit(FatalError);
                break;
            }
            case interpolationTable::WARN:
            {
                WarningIn
                (
                    "Foam::interpolationTable<Type>::operator[]"
                    "(const scalar) const"
                )   << "value (" << lookupValue << ") underflow" << nl
                    << "    Continuing with the first entry"
                    << endl;
                // fall-through to 'CLAMP'
            }
            case interpolationTable::CLAMP:
            {
                return List<Tuple2<scalar, Type> >::operator[](0).second();
                break;
            }
            case interpolationTable::REPEAT:
            {
                // adjust lookupValue to >= 0
                while (lookupValue < 0)
                {
                    lookupValue += maxLimit;
                }
                break;
            }
        }
    }
    else if (lookupValue >= maxLimit)
    {
        switch (boundAction_)
        {
            case interpolationTable::ERROR:
            {
                FatalErrorIn
                (
                    "Foam::interpolationTable<Type>::operator[]"
                    "(const label) const"
                )   << "value (" << lookupValue << ") overflow" << nl
                    << exit(FatalError);
                break;
            }
            case interpolationTable::WARN:
            {
                WarningIn
                (
                    "Foam::interpolationTable<Type>::operator[]"
                    "(const label) const"
                )   << "value (" << lookupValue << ") overflow" << nl
                    << "    Continuing with the last entry"
                    << endl;
                // fall-through to 'CLAMP'
            }
            case interpolationTable::CLAMP:
            {
                return List<Tuple2<scalar, Type> >::operator[](n-1).second();
                break;
            }
            case interpolationTable::REPEAT:
            {
                // adjust lookupValue <= maxLimit
                while (lookupValue > maxLimit)
                {
                    lookupValue -= maxLimit;
                }
                break;
            }
        }
    }

    label lo = 0;
    label hi = 0;

    // look for the correct range
    for (label i = 0; i < n; ++i)
    {
        if (lookupValue >= List<Tuple2<scalar, Type> >::operator[](i).first())
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
        return List<Tuple2<scalar, Type> >::operator[](hi).second();
    }
    else if (hi == 0)
    {
        // this treatment should should only occur under these conditions:
        //  -> the 'REPEAT' treatment
        //  -> (0 <= value <= minLimit)
        //  -> minLimit > 0
        // Use the value at maxLimit as the value for value=0
        lo = n - 1;

        return
        (
            List<Tuple2<scalar, Type> >::operator[](lo).second()
          + (
                List<Tuple2<scalar, Type> >::operator[](hi).second()
              - List<Tuple2<scalar, Type> >::operator[](lo).second()
            )
           *(lookupValue / minLimit)
        );
    }
    else
    {
        // normal interpolation
        return
        (
            List<Tuple2<scalar, Type> >::operator[](lo).second()
          + (
                List<Tuple2<scalar, Type> >::operator[](hi).second()
              - List<Tuple2<scalar, Type> >::operator[](lo).second()
            )
           *(
                lookupValue
              - List<Tuple2<scalar, Type> >::operator[](lo).first()
            )
           /(
                List<Tuple2<scalar, Type> >::operator[](hi).first()
              - List<Tuple2<scalar, Type> >::operator[](lo).first()
            )
        );
    }
}


// ************************************************************************* //

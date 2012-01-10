/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "IFstream.H"
#include "openFoamTableReader.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class Type>
void Foam::interpolation2DTable<Type>::readTable()
{
    fileName fName(fileName_);
    fName.expand();

    // Read data from file
    reader_()(fName, *this);
    //IFstream(fName)() >> *this;

    if (this->empty())
    {
        FatalErrorIn
        (
            "Foam::interpolation2DTable<Type>::readTable()"
        )   << "table read from " << fName << " is empty" << nl
            << exit(FatalError);
    }

    // Check that the data are okay
    check();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::interpolation2DTable<Type>::interpolation2DTable()
:
    List<Tuple2<scalar, List<Tuple2<scalar, Type> > > >(),
    boundsHandling_(interpolation2DTable::WARN),
    fileName_("fileNameIsUndefined"),
    reader_(NULL)
{}


template<class Type>
Foam::interpolation2DTable<Type>::interpolation2DTable
(
    const List<Tuple2<scalar, List<Tuple2<scalar, Type> > > >& values,
    const boundsHandling bounds,
    const fileName& fName
)
:
    List<Tuple2<scalar, List<Tuple2<scalar, Type> > > >(values),
    boundsHandling_(bounds),
    fileName_(fName),
    reader_(NULL)
{}


template<class Type>
Foam::interpolation2DTable<Type>::interpolation2DTable(const fileName& fName)
:
    List<Tuple2<scalar, List<Tuple2<scalar, Type> > > >(),
    boundsHandling_(interpolation2DTable::WARN),
    fileName_(fName),
    reader_(new openFoamTableReader<Type>(dictionary()))
{
    readTable();
}


template<class Type>
Foam::interpolation2DTable<Type>::interpolation2DTable(const dictionary& dict)
:
    List<Tuple2<scalar, List<Tuple2<scalar, Type> > > >(),
    boundsHandling_(wordToBoundsHandling(dict.lookup("outOfBounds"))),
    fileName_(dict.lookup("fileName")),
    reader_(tableReader<Type>::New(dict))
{
    readTable();
}


template<class Type>
Foam::interpolation2DTable<Type>::interpolation2DTable
(
     const interpolation2DTable& interpTable
)
:
    List<Tuple2<scalar, List<Tuple2<scalar, Type> > > >(interpTable),
    boundsHandling_(interpTable.boundsHandling_),
    fileName_(interpTable.fileName_),
    reader_(interpTable.reader_)    // note: steals reader. Used in write().
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Type Foam::interpolation2DTable<Type>::interpolateValue
(
    const List<Tuple2<scalar, Type> >& data,
    const scalar lookupValue
) const
{
    label n = data.size();

    scalar minLimit = data[0].first();
    scalar maxLimit = data[n-1].first();

    if (lookupValue < minLimit)
    {
        switch (boundsHandling_)
        {
            case interpolation2DTable::ERROR:
            {
                FatalErrorIn
                (
                    "Foam::interpolation2DTable<Type>::interpolateValue("
                    "List<Tuple2<scalar, Type> > data,"
                    "const scalar lookupValue)"
                )   << "value (" << lookupValue << ") underflow" << nl
                    << exit(FatalError);
                break;
            }
            case interpolation2DTable::WARN:
            {
                WarningIn
                (
                    "Foam::interpolation2DTable<Type>::interpolateValue("
                    "List<Tuple2<scalar, Type> > data,"
                    "const scalar lookupValue)"
                )   << "value (" << lookupValue << ") underflow" << nl
                    << "    Continuing with the first entry"
                    << endl;
                // fall-through to 'CLAMP'
            }
            case interpolation2DTable::CLAMP:
            {
                return data[0].second();
                break;
            }
        }
    }
    else if (lookupValue >= maxLimit)
    {
        switch (boundsHandling_)
        {
            case interpolation2DTable::ERROR:
            {
                FatalErrorIn
                (
                    "Foam::interpolation2DTable<Type>::interpolateValue("
                    "List<Tuple2<scalar, Type> > data,"
                    "const scalar lookupValue)"
                )   << "value (" << lookupValue << ") overflow" << nl
                    << exit(FatalError);
                break;
            }
            case interpolation2DTable::WARN:
            {
                WarningIn
                (
                    "Foam::interpolation2DTable<Type>::interpolateValue("
                    "List<Tuple2<scalar, Type> > data,"
                    "const scalar lookupValue)"
                )   << "value (" << lookupValue << ") overflow" << nl
                    << "    Continuing with the last entry"
                    << endl;
                // fall-through to 'CLAMP'
            }
            case interpolation2DTable::CLAMP:
            {
                return data[n-1].second();
                break;
            }
        }
    }

    // look for the correct range in X
    label lo = 0;
    label hi = 0;

    for (label i = 0; i < n; ++i)
    {
        if (lookupValue >= data[i].first())
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
        return data[lo].second();
    }
    else
    {
        // normal interpolation
        return
        (
            data[lo].second()
          + (
                data[hi].second()
              - data[lo].second()
            )
           *(
                lookupValue
              - data[lo].first()
            )
           /(
                data[hi].first()
              - data[lo].first()
            )
        );
    }

}


template<class Type>
Foam::word Foam::interpolation2DTable<Type>::boundsHandlingToWord
(
     const boundsHandling& bound
) const
{
    word enumName("warn");

    switch (bound)
    {
        case interpolation2DTable::ERROR:
        {
            enumName = "error";
            break;
        }
        case interpolation2DTable::WARN:
        {
            enumName = "warn";
            break;
        }
        case interpolation2DTable::CLAMP:
        {
            enumName = "clamp";
            break;
        }
    }

    return enumName;
}


template<class Type>
typename Foam::interpolation2DTable<Type>::boundsHandling
Foam::interpolation2DTable<Type>::wordToBoundsHandling
(
    const word& bound
) const
{
    if (bound == "error")
    {
        return interpolation2DTable::ERROR;
    }
    else if (bound == "warn")
    {
        return interpolation2DTable::WARN;
    }
    else if (bound == "clamp")
    {
        return interpolation2DTable::CLAMP;
    }
    else
    {
        WarningIn
        (
            "Foam::interpolation2DTable<Type>::wordToBoundsHandling(const word&)"
        )   << "bad outOfBounds specifier " << bound << " using 'warn'" << endl;

        return interpolation2DTable::WARN;
    }
}


template<class Type>
typename Foam::interpolation2DTable<Type>::boundsHandling
Foam::interpolation2DTable<Type>::outOfBounds
(
    const boundsHandling& bound
)
{
    boundsHandling prev = boundsHandling_;
    boundsHandling_ = bound;
    return prev;
}


template<class Type>
void Foam::interpolation2DTable<Type>::check() const
{
    label n = this->size();
    typedef  List<Tuple2<scalar, List<Tuple2<scalar, Type> > > > matrix;
    scalar prevValue = matrix::operator[](0).first();

    for (label i=1; i<n; ++i)
    {
        const scalar currValue = matrix::operator[](i).first();

        // avoid duplicate values (divide-by-zero error)
        if (currValue <= prevValue)
        {
            FatalErrorIn
            (
                "Foam::interpolation2DTable<Type>::check() const"
            )   << "out-of-order value: "
                << currValue << " at index " << i << nl
                << exit(FatalError);
        }
        prevValue = currValue;
    }
}


template<class Type>
void Foam::interpolation2DTable<Type>::write(Ostream& os) const
{
    os.writeKeyword("fileName")
        << fileName_ << token::END_STATEMENT << nl;
    os.writeKeyword("outOfBounds")
        << boundsHandlingToWord(boundsHandling_) << token::END_STATEMENT << nl;

    if (reader_.valid())
    {
        reader_->write(os);
    }
    //*this >> os;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
const Foam::List<Foam::Tuple2<Foam::scalar, Type> >&
Foam::interpolation2DTable<Type>::operator[](const label i) const
{
    label ii = i;
    label n  = this->size();

    if (n <= 1)
    {
        ii = 0;
    }
    else if (ii < 0)
    {
        switch (boundsHandling_)
        {
            case interpolation2DTable::ERROR:
            {
                FatalErrorIn
                (
                    "Foam::interpolation2DTable<Type>::operator[]"
                    "(const label) const"
                )   << "index (" << ii << ") underflow" << nl
                    << exit(FatalError);
                break;
            }
            case interpolation2DTable::WARN:
            {
                WarningIn
                (
                    "Foam::interpolation2DTable<Type>::operator[]"
                    "(const label) const"
                )   << "index (" << ii << ") underflow" << nl
                    << "    Continuing with the first entry"
                    << endl;
                // fall-through to 'CLAMP'
            }
            case interpolation2DTable::CLAMP:
            {
                ii = 0;
                break;
            }
        }
    }
    else if (ii >= n)
    {
        switch (boundsHandling_)
        {
            case interpolation2DTable::ERROR:
            {
                FatalErrorIn
                (
                    "Foam::interpolation2DTable<Type>::operator[]"
                    "(const label) const"
                )   << "index (" << ii << ") overflow" << nl
                    << exit(FatalError);
                break;
            }
            case interpolation2DTable::WARN:
            {
                WarningIn
                (
                    "Foam::interpolation2DTable<Type>::operator[]"
                    "(const label) const"
                )   << "index (" << ii << ") overflow" << nl
                    << "    Continuing with the last entry"
                    << endl;
                // fall-through to 'CLAMP'
            }
            case interpolation2DTable::CLAMP:
            {
                ii = n - 1;
                break;
            }
        }
    }

    return List<Tuple2<scalar, List<Tuple2<scalar, Type> > > >::operator[](ii);
}


template<class Type>
Type Foam::interpolation2DTable<Type>::operator()
(
    const scalar valueX,
    const scalar valueY
) const
{
    typedef  List<Tuple2<scalar, List<Tuple2<scalar, Type> > > > matrix;
    label nX = this->size();

    if (nX <= 1)
    {
        const List<Tuple2<scalar, Type> >& dataY =
            matrix::operator[](0).second();

        return interpolateValue(dataY, valueY);
    }


    scalar minLimit = matrix::operator[](0).first();
    scalar maxLimit = matrix::operator[](nX-1).first();
    scalar lookupValue = valueX;

    if (lookupValue < minLimit)
    {
        switch (boundsHandling_)
        {
            case interpolation2DTable::ERROR:
            {
                FatalErrorIn
                (
                    "Foam::interpolation2DTable<Type>::operator[]"
                    "(const scalar, const scalar) const"
                )   << "value (" << lookupValue << ") underflow" << nl
                    << exit(FatalError);
                break;
            }
            case interpolation2DTable::WARN:
            {
                WarningIn
                (
                    "Foam::interpolation2DTable<Type>::operator[]"
                    "(const scalar, const scalar) const"
                )   << "value (" << lookupValue << ") underflow" << nl
                    << "    Continuing with the first entry"
                    << endl;
                // fall-through to 'CLAMP'
            }
            case interpolation2DTable::CLAMP:
            {
                return interpolateValue
                (
                    matrix::operator[](0).second(), valueY
                );
                break;
            }
        }
    }
    else if (lookupValue >= maxLimit)
    {
        switch (boundsHandling_)
        {
            case interpolation2DTable::ERROR:
            {
                FatalErrorIn
                (
                    "Foam::interpolation2DTable<Type>::operator[]"
                    "(const label, const scalar) const"
                )   << "value (" << lookupValue << ") overflow" << nl
                    << exit(FatalError);
                break;
            }
            case interpolation2DTable::WARN:
            {
                WarningIn
                (
                    "Foam::interpolation2DTable<Type>::operator[]"
                    "(const label, const scalar) const"
                )   << "value (" << lookupValue << ") overflow" << nl
                    << "    Continuing with the last entry"
                    << endl;
                // fall-through to 'CLAMP'
            }
            case interpolation2DTable::CLAMP:
            {
                return interpolateValue
                (
                    matrix::operator[](nX-1).second(), valueY
                );
                break;
            }
        }
    }

    label loX = 0;
    label hiX = 0;

    // look for the correct range in X
    for (label i = 0; i < nX; ++i)
    {
        if (lookupValue >= matrix::operator[](i).first())
        {
            loX = hiX = i;
        }
        else
        {
            hiX = i;
            break;
        }
    }

    // look for the correct range in y
    lookupValue = valueY;
    label loY1 = 0;
    label hiY1 = 0;

    label nY = matrix::operator[](loX).second().size();

    minLimit = matrix::operator[](loX).second()[0].first();
    maxLimit = matrix::operator[](loX).second()[nY-1].first();

    if (lookupValue < minLimit)
    {
        switch (boundsHandling_)
        {
            case interpolation2DTable::ERROR:
            {
                FatalErrorIn
                (
                    "Foam::interpolation2DTable<Type>::operator[]"
                    "(const scalar, const scalar) const"
                )   << "value (" << lookupValue << ") underflow" << nl
                    << exit(FatalError);
                break;
            }
            case interpolation2DTable::WARN:
            {
                WarningIn
                (
                    "Foam::interpolation2DTable<Type>::operator[]"
                    "(const scalar, const scalar) const"
                )   << "value (" << lookupValue << ") underflow" << nl
                    << "    Continuing with the first entry"
                    << endl;
                // fall-through to 'CLAMP'
            }
            case interpolation2DTable::CLAMP:
            {
                hiY1 = 0;
                hiY1 = 1;
                break;
            }
        }
    }
    else if (lookupValue >= maxLimit)
    {
        switch (boundsHandling_)
        {
            case interpolation2DTable::ERROR:
            {
                FatalErrorIn
                (
                    "Foam::interpolation2DTable<Type>::operator[]"
                    "(const scalar, const scalar) const"
                )   << "value (" << lookupValue << ") overflow" << nl
                    << exit(FatalError);
                break;
            }
            case interpolation2DTable::WARN:
            {
                WarningIn
                (
                    "Foam::interpolation2DTable<Type>::operator[]"
                    "(const scalar, const scalar) const"
                )   << "value (" << lookupValue << ") overflow" << nl
                    << "    Continuing with the last entry"
                    << endl;
                // fall-through to 'CLAMP'
            }
            case interpolation2DTable::CLAMP:
            {
                hiY1 = nY-1;
                hiY1 = nY;
                break;
            }
        }
    }
    else
    {
        // Finds the lo and hi of Y on the lowest x
        for (label i = 0; i < nY; ++i)
        {
            if
            (
                lookupValue >= matrix::operator[](loX).second()[i].first()
            )
            {
                loY1 = hiY1 = i;
            }
            else
            {
                hiY1 = i;
                break;
            }
        }
    }

    label loY2 = 0;
    label hiY2 = 0;

    nY = matrix::operator[](hiX).second().size();

    minLimit = matrix::operator[](loX).second()[0].first();
    maxLimit = matrix::operator[](loX).second()[nY-1].first();

    if (lookupValue < minLimit)
    {
        switch (boundsHandling_)
        {
            case interpolation2DTable::ERROR:
            {
                FatalErrorIn
                (
                    "Foam::interpolation2DTable<Type>::operator[]"
                    "(const scalar, const scalar) const"
                )   << "value (" << lookupValue << ") underflow" << nl
                    << exit(FatalError);
                break;
            }
            case interpolation2DTable::WARN:
            {
                WarningIn
                (
                    "Foam::interpolation2DTable<Type>::operator[]"
                    "(const scalar, const scalar) const"
                )   << "value (" << lookupValue << ") underflow" << nl
                    << "    Continuing with the first entry"
                    << endl;
                // fall-through to 'CLAMP'
            }
            case interpolation2DTable::CLAMP:
            {
                loY2 = 0;
                loY2 = 1;
                break;
            }
        }
    }
    else if (lookupValue >= maxLimit)
    {
        switch (boundsHandling_)
        {
            case interpolation2DTable::ERROR:
            {
                FatalErrorIn
                (
                    "Foam::interpolation2DTable<Type>::operator[]"
                    "(const scalar, const scalar) const"
                )   << "value (" << lookupValue << ") overflow" << nl
                    << exit(FatalError);
                break;
            }
            case interpolation2DTable::WARN:
            {
                WarningIn
                (
                    "Foam::interpolation2DTable<Type>::operator[]"
                    "(const scalar, const scalar) const"
                )   << "value (" << lookupValue << ") overflow" << nl
                    << "    Continuing with the last entry"
                    << endl;
                // fall-through to 'CLAMP'
            }
            case interpolation2DTable::CLAMP:
            {
                loY2 = nY-1;
                loY2 = nY;
                break;
            }
        }
    }
    else
    {
        // Finds the lo and hi of Y on the high x
        for (label i = 0; i < nY; ++i)
        {
            if
            (
                lookupValue >= matrix::operator[](hiX).second()[i].first()
            )
            {
                loY2 = hiY2 = i;
            }
            else
            {
                hiY2 = i;
                break;
            }
        }
    }

    if (loX == hiX)
    {
        // we are at the end of the table - or there is only a single entry
        return (interpolateValue(matrix::operator[](hiX).second(), valueY));
    }
    else
    {
        Type loXData = matrix::operator[](loX).second()[loY1].second();
        Type hiXData = matrix::operator[](hiX).second()[loY1].second();

        Type hiYData = matrix::operator[](loX).second()[hiY1].second();

        Type refValue = matrix::operator[](loX).second()[loY1].second();

        // normal interpolation on x
        refValue +=
            (
                hiXData
              - loXData
            )
            *(
                valueX
              - matrix::operator[](loX).first()
             )
            /(
                matrix::operator[](hiX).first()
              - matrix::operator[](loX).first()
            );

        // normal interpolation on y
        refValue +=
            (
                hiYData
              - loXData
            )
            *(
                valueY
              - matrix::operator[](loX).second()[loY1].first()
            )
            /(
                matrix::operator[](loX).second()[hiY1].first()
              - matrix::operator[](loX).second()[loY1].first()
            );

        return refValue;
    }
}


// ************************************************************************* //

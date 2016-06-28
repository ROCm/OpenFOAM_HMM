/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2016 OpenCFD Ltd.
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

#include <iomanip>
#include <sstream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::ensightSurfaceReader::readSkip
(
    IFstream& is,
    const label nSkip,
    Type& value
) const
{
    skip(nSkip, is);

    is  >> value;
}


template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::ensightSurfaceReader::readField
(
    const label timeIndex,
    const label fieldIndex
) const
{
    if (debug)
    {
        InfoInFunction<< endl;
    }

    const word& fieldName(fieldNames_[fieldIndex]);
    const label fileIndex = timeStartIndex_ + timeIndex*timeIncrement_;

    fileName fieldFileName(fieldFileNames_[fieldIndex]);

    std::ostringstream oss;
    oss << std::setfill('0') << std::setw(4) << fileIndex;
    const word indexStr = oss.str();
    fieldFileName.replace("****", indexStr);

    IFstream is(baseDir_/fieldFileName);

    if (!is.good())
    {
        FatalErrorInFunction
            << "Cannot read file " << is.name()
            << " for field " << fieldName
            << exit(FatalError);
    }

    // Check that data type is as expected
    word primitiveType(is);

    if (debug)
    {
        Info<< "primitiveType: " << primitiveType << endl;
    }

    if (primitiveType != pTraits<Type>::typeName)
    {
        FatalIOErrorInFunction(is)
            << "Expected " << pTraits<Type>::typeName << "values "
            << "but found type " << primitiveType
            << exit(FatalIOError);
    }

    tmp<Field<Type> > tField(new Field<Type>());

    label n;
    if (surfPtr_.valid())
    {
        n = surfPtr_->size();
    }
    else
    {
        n = 1000;
    }

    Type value;
    word wValue;
    label iValue;

    // Read header info: part index, e.g. part 1
    is  >> wValue >> iValue;

    // Read data file
    // - Assume that file contains a mix of words and numbers, and that all
    //   numbers relate to face values, e.g. header comprises of words and
    //   element types are also words, e.g. tria3, quad4, nsided
    DynamicList<Type> values(n);
    while (is.good())
    {
        token t(is);

        if (is.eof())
        {
            break;
        }

        if (t.isWord())
        {
            wValue = t.wordToken();
        }
        else
        {
            is.putBack(t);
            is  >> value;
            values.append(value);
        }
    }

    tField().transfer(values);

    return tField;
}


// ************************************************************************* //

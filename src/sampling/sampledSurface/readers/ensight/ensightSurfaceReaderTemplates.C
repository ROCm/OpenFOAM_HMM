/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2017 OpenCFD Ltd.
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
void Foam::ensightSurfaceReader::readFromLine
(
    const label nSkip,
    IStringStream& is,
    Type& value
) const
{
    skip(nSkip, is);

    is  >> value;
}


template<class Type>
void Foam::ensightSurfaceReader::readFromLine
(
    const label nSkip,
    const string& buffer,
    Type& value
) const
{
    IStringStream is(buffer);

    readFromLine(nSkip, is, value);
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::ensightSurfaceReader::readField
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
    label nMask = 0;
    for (size_t chari = 0; chari < fieldFileName.size(); ++chari)
    {
        if (fieldFileName[chari] == '*')
        {
            nMask++;
        }
    }

    const std::string maskStr(nMask, '*');
    oss << std::setfill('0') << std::setw(nMask) << fileIndex;
    const word indexStr = oss.str();
    fieldFileName.replace(maskStr, indexStr);


    ensightReadFile is(baseDir_/fieldFileName, streamFormat_);

    if (!is.good())
    {
        FatalErrorInFunction
            << "Cannot read file " << is.name()
            << " for field " << fieldName
            << exit(FatalError);
    }

    // Check that data type is as expected
    string primitiveType;
    is.read(primitiveType);


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

    scalar value;
    string strValue;
    label iValue;

    // Read header info: part index, e.g. part 1
    is.read(strValue);
    is.read(iValue);

    // Allocate storage for data as a list per component
    List<DynamicList<scalar>> values(pTraits<Type>::nComponents);
    label n = surfPtr_->size();
    forAll(values, cmptI)
    {
        values[cmptI].setCapacity(n);
    }

    // Read data file using schema generated while reading the surface
    forAll(schema_, i)
    {
        if (debug)
        {
            const string& faceType = schema_[i].first();
            Info<< "Reading face type " << faceType << " data" << endl;
        }

        const label nFace = schema_[i].second();

        if (nFace != 0)
        {
            is.read(strValue);

            for
            (
                direction cmptI=0;
                cmptI < pTraits<Type>::nComponents;
                ++cmptI
            )
            {
                for (label faceI = 0; faceI < nFace; ++faceI)
                {
                    is.read(value);
                    values[cmptI].append(value);
                }
            }
        }
    }

    tmp<Field<Type>> tField(new Field<Type>(n, pTraits<Type>::zero));
    Field<Type>& field = tField.ref();

    for
    (
        direction cmptI=0;
        cmptI < pTraits<Type>::nComponents;
        ++cmptI
    )
    {
        field.replace(cmptI, values[cmptI]);
        values[cmptI].clear();
    }

    return tField;
}


// ************************************************************************* //

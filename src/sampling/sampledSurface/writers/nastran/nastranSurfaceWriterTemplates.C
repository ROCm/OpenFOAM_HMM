/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015-2017 OpenCFD Ltd.
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

#include "OFstream.H"
#include "IOmanip.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::Ostream& Foam::nastranSurfaceWriter::writeValue
(
    Ostream& os,
    const Type& value
) const
{
    switch (writeFormat_)
    {
        case fieldFormat::SHORT :
        {
            os  << setw(8) << value;
            break;
        }

        case fieldFormat::LONG :
        {
            os  << setw(16) << value;
            break;
        }

        case fieldFormat::FREE :
        {
            os  << value;
            break;
        }
    }

    return os;
}


template<class Type>
Foam::Ostream& Foam::nastranSurfaceWriter::writeFaceValue
(
    Ostream& os,
    const loadFormat format,
    const Type& value,
    const label EID
) const
{
    // Fixed short/long formats supporting PLOAD2 and PLOAD4:

    // PLOAD2:
    // 1 descriptor : PLOAD2
    // 2 SID        : load set ID
    // 3 data value : load value - MUST be singular
    // 4 EID        : element ID

    // PLOAD4:
    // 1 descriptor : PLOAD4
    // 2 SID        : load set ID
    // 3 EID        : element ID
    // 4 onwards    : load values

    label SID = 1;

    Type scaledValue = scale_*value;

    // Write keyword
    writeKeyword(os, loadFormatNames_[format])  << separator_;

    // Write load set ID
    os.setf(std::ios_base::right);

    writeValue(os, SID) << separator_;

    switch (format)
    {
        case loadFormat::PLOAD2 :
        {
            if (pTraits<Type>::nComponents == 1)
            {
                writeValue(os, scaledValue) << separator_;
            }
            else
            {
                WarningInFunction
                    << loadFormatNames_[format] << " requires scalar values "
                    << "and cannot be used for higher rank values"
                    << endl;

                writeValue(os, scalar(0)) << separator_;
            }

            writeValue(os, EID);
            break;
        }

        case loadFormat::PLOAD4 :
        {
            writeValue(os, EID);

            for (direction d = 0; d < pTraits<Type>::nComponents; ++d)
            {
                os  << separator_;
                writeValue(os, component(scaledValue, d));
            }
            break;
        }
    }

    os.unsetf(std::ios_base::right);

    os << nl;

    return os;
}


template<class Type>
Foam::fileName Foam::nastranSurfaceWriter::writeTemplate
(
    const fileName& outputDir,
    const fileName& surfaceName,
    const meshedSurf& surf,
    const word& fieldName,
    const Field<Type>& values,
    const bool isNodeValues,
    const bool verbose
) const
{
    if (!fieldMap_.found(fieldName))
    {
        FatalErrorInFunction
            << "No mapping found between field " << fieldName
            << " and corresponding Nastran field.  Available types are:"
            << fieldMap_
            << exit(FatalError);

        return fileName::null;
    }

    const loadFormat& format(fieldMap_[fieldName]);

    if (!isDir(outputDir/fieldName))
    {
        mkDir(outputDir/fieldName);
    }

    // const scalar timeValue = Foam::name(this->mesh().time().timeValue());
    const scalar timeValue = 0.0;

    OFstream os(outputDir/fieldName/surfaceName + ".nas");
    fileFormats::NASCore::setPrecision(os, writeFormat_);

    if (verbose)
    {
        Info<< "Writing nastran file to " << os.name() << endl;
    }

    os  << "TITLE=OpenFOAM " << surfaceName.c_str()
        << " " << fieldName << " data" << nl
        << "$" << nl
        << "TIME " << timeValue << nl
        << "$" << nl
        << "BEGIN BULK" << nl;

    List<DynamicList<face>> decomposedFaces;
    writeGeometry(os, surf, decomposedFaces);

    os  << "$" << nl
        << "$ Field data" << nl
        << "$" << nl;

    label elemId = 0;

    if (isNodeValues)
    {
        for (const DynamicList<face>& dFaces : decomposedFaces)
        {
            for (const face& f : dFaces)
            {
                Type v = Zero;

                for (const label verti : f)
                {
                    v += values[verti];
                }
                v /= f.size();

                writeFaceValue(os, format, v, ++elemId);
            }
        }
    }
    else
    {
        for (const DynamicList<face>& dFaces : decomposedFaces)
        {
            forAll(dFaces, facei)
            {
                writeFaceValue(os, format, values[facei], ++elemId);
            }
        }
    }

    writeFooter(os, surf)
        << "ENDDATA" << endl;

    return os.name();
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015-2016 OpenCFD Ltd.
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
void Foam::nastranSurfaceWriter::writeValue
(
    const Type& value,
    Ostream& os
) const
{
    switch (writeFormat_)
    {
        case wfShort:
        {
            os  << setw(8) << value;
            break;
        }
        case wfLong:
        {
            os  << setw(16) << value;
            break;
        }
        case wfFree:
        {
            os  << value;
            break;
        }
    }
}


template<class Type>
void Foam::nastranSurfaceWriter::writeFaceValue
(
    const dataFormat& format,
    const Type& value,
    const label EID,
    Ostream& os
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

    // Write Keyword
    writeKeyword(dataFormatNames_[format], os);

    os  << separator_;

    // Write load set ID
    os.setf(ios_base::right);
    writeValue(SID, os);

    os  << separator_;

    switch (format)
    {
        case dfPLOAD2:
        {
            if (pTraits<Type>::nComponents == 1)
            {
                writeValue(scaledValue, os);
            }
            else
            {
                WarningInFunction
                    << dataFormatNames_[format] << " requires scalar values "
                    << "and cannot be used for higher rank values"
                    << endl;

                writeValue(scalar(0), os);
            }

            os  << separator_;
            writeValue(EID, os);
            break;
        }

        case dfPLOAD4:
        {
            writeValue(EID, os);

            for (direction dirI = 0; dirI < pTraits<Type>::nComponents; dirI++)
            {
                os  << separator_;
                writeValue(component(scaledValue, dirI), os);
            }
            break;
        }

        default:
        {
            FatalErrorInFunction
                << "Unhandled enumeration " << dataFormatNames_[format]
                << exit(FatalError);
        }
    }

    os.unsetf(ios_base::right);

    os << nl;
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

    const dataFormat& format(fieldMap_[fieldName]);

    if (!isDir(outputDir/fieldName))
    {
        mkDir(outputDir/fieldName);
    }

    // const scalar timeValue = Foam::name(this->mesh().time().timeValue());
    const scalar timeValue = 0.0;

    OFstream os(outputDir/fieldName/surfaceName + ".nas");
    formatOS(os);

    if (verbose)
    {
        Info<< "Writing nastran file to " << os.name() << endl;
    }

    os  << "TITLE=OpenFOAM " << surfaceName.c_str() << " " << fieldName
        << " data" << nl
        << "$" << nl
        << "TIME " << timeValue << nl
        << "$" << nl
        << "BEGIN BULK" << nl;

    List<DynamicList<face>> decomposedFaces;
    writeGeometry(surf, decomposedFaces, os);

    os  << "$" << nl
        << "$ Field data" << nl
        << "$" << nl;

    if (isNodeValues)
    {
        label n = 0;

        forAll(decomposedFaces, i)
        {
            const DynamicList<face>& dFaces = decomposedFaces[i];
            forAll(dFaces, facei)
            {
                Type v = Zero;
                const face& f = dFaces[facei];

                forAll(f, fptI)
                {
                    v += values[f[fptI]];
                }
                v /= f.size();

                writeFaceValue(format, v, ++n, os);
            }
        }
    }
    else
    {
        label n = 0;

        forAll(decomposedFaces, i)
        {
            const DynamicList<face>& dFaces = decomposedFaces[i];

            forAll(dFaces, facei)
            {
                writeFaceValue(format, values[facei], ++n, os);
            }
        }
    }

    writeFooter(os, surf);

    os  << "ENDDATA" << endl;

    return os.name();
}


// ************************************************************************* //

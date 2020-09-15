/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
    Copyright (C) 2015-2020 OpenCFD Ltd.
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
#include "ListOps.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::Ostream& Foam::surfaceWriters::nastranWriter::writeValue
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
Foam::Ostream& Foam::surfaceWriters::nastranWriter::writeFaceValue
(
    Ostream& os,
    const loadFormat format,
    const Type& value,
    const label elemId
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

    const label setId = 1;

    // Write keyword
    writeKeyword(os, fileFormats::NASCore::loadFormatNames[format])
        << separator_;

    // Write load set ID
    os.setf(std::ios_base::right);

    writeValue(os, setId) << separator_;

    switch (format)
    {
        case loadFormat::PLOAD2 :
        {
            if (pTraits<Type>::nComponents == 1)
            {
                writeValue(os, value);
            }
            else
            {
                writeValue(os, mag(value));
            }

            os << separator_;
            writeValue(os, elemId);
            break;
        }

        case loadFormat::PLOAD4 :
        {
            writeValue(os, elemId);

            for (direction d = 0; d < pTraits<Type>::nComponents; ++d)
            {
                os  << separator_;
                writeValue(os, component(value, d));
            }
            break;
        }
    }

    os.unsetf(std::ios_base::right);

    os << nl;

    return os;
}


template<class Type>
Foam::fileName Foam::surfaceWriters::nastranWriter::writeTemplate
(
    const word& fieldName,
    const Field<Type>& localValues
)
{
    // Geometry changed since last output? Capture now before any merging.
    /// const bool geomChanged = (!upToDate_);

    // Separate geometry, when commonGeometry = true
    if (!wroteGeom_ && commonGeometry_)
    {
        write();
    }

    checkOpen();

    const loadFormat format
    (
        fieldMap_.lookup
        (
            fieldName,
            // Default format
            (
                pTraits<Type>::nComponents == 1
              ? loadFormat::PLOAD2
              : loadFormat::PLOAD4
            )
        )
    );

    if
    (
        !std::is_integral<Type>::value  // Handle 'Ids' etc silently
     && !fieldMap_.empty()
     && !fieldMap_.found(fieldName)
    )
    {
        WarningInFunction
            << "No mapping found between field " << fieldName
            << " and corresponding Nastran field.  Available types:"
            << fieldMap_ << nl;
    }

    // Emit any common warnings
    if (format == loadFormat::PLOAD2 && pTraits<Type>::nComponents != 1)
    {
        WarningInFunction
            << fileFormats::NASCore::loadFormatNames[format]
            << " cannot be used for higher rank values"
            << " - reverting to mag()" << endl;
    }


    // Common geometry
    // Field:  rootdir/<TIME>/<field>_surfaceName.bdf

    // Embedded geometry
    // Field:  rootdir/<TIME>/<field>/surfaceName.bdf

    fileName outputFile = outputPath_.path();
    if (useTimeDir() && !timeName().empty())
    {
        // Splice in time-directory
        outputFile /= timeName();
    }

    fileName geomFileName;
    if (commonGeometry_)
    {
        // Common geometry
        geomFileName = outputPath_.name().ext("nas");

        // Append <field>_surfaceName.bdf
        outputFile /= fieldName + '_' + outputPath_.name();
    }
    else
    {
        // Embedded geometry
        // Use sub-directory
        outputFile /= fieldName / outputPath_.name();
    }
    outputFile.ext("bdf");


    // Output scaling for the variable, but not for integer types.
    // could also solve with clever templating

    const scalar varScale =
    (
        std::is_integral<Type>::value
      ? scalar(1)
      : fieldScale_.getOrDefault<scalar>(fieldName, 1)
    );

    if (verbose_)
    {
        Info<< "Writing field " << fieldName;
        if (!equal(varScale, 1))
        {
            Info<< " (scaling " << varScale << ')';
        }
        Info<< " to " << outputFile << endl;
    }


    // Implicit geometry merge()
    tmp<Field<Type>> tfield = mergeField(localValues) * varScale;

    const meshedSurf& surf = surface();

    if (Pstream::master() || !parallel_)
    {
        const auto& values = tfield();

        if (!isDir(outputFile.path()))
        {
            mkDir(outputFile.path());
        }

        const scalar timeValue(0);

        // Additional bookkeeping for decomposing non tri/quad
        labelList decompOffsets;
        DynamicList<face> decompFaces;


        OFstream os(outputFile);
        fileFormats::NASCore::setPrecision(os, writeFormat_);

        os  << "TITLE=OpenFOAM " << outputFile.name()
            << token::SPACE << fieldName << " data" << nl;

        if (useTimeDir() && !timeName().empty())
        {
            os  << '$' << nl
                << "$ TIME " << timeName() << nl;
        }

        os  << "TIME " << timeValue << nl
            << nl
            << "BEGIN BULK" << nl;

        if (commonGeometry_)
        {
            os  << "INCLUDE '" << geomFileName.c_str() << "'" << nl;

            // Geometry already written (or suppressed)
            // - still need decomposition information
            fileFormats::NASCore::faceDecomposition
            (
                surf.points(),
                surf.faces(),
                decompOffsets,
                decompFaces
            );
        }
        else
        {
            // Write geometry
            writeGeometry(os, surf, decompOffsets, decompFaces);
        }

        // Write field
        os  << '$' << nl
            << "$ Field data" << nl
            << '$' << nl;


        // Regular (undecomposed) faces
        const faceList& faces = surf.faces();
        const labelUList& elemIds = surf.faceIds();

        // Possible to use faceIds?
        const bool useOrigFaceIds =
        (
            elemIds.size() == faces.size()
         && !ListOps::found(elemIds, lessOp1<label>(0))
         && decompFaces.empty()
        );


        label elemId = 0;

        if (this->isPointData())
        {
            forAll(faces, facei)
            {
                if (useOrigFaceIds)
                {
                    elemId = elemIds[facei];
                }

                const label beginElemId = elemId;

                // Any face decomposition
                for
                (
                    label decompi = decompOffsets[facei];
                    decompi < decompOffsets[facei+1];
                    ++decompi
                )
                {
                    const face& f = decompFaces[decompi];

                    Type v = Zero;
                    for (const label verti : f)
                    {
                        v += values[verti];
                    }
                    v /= f.size();

                    writeFaceValue(os, format, v, ++elemId);
                }


                // Face not decomposed
                if (beginElemId == elemId)
                {
                    const face& f = faces[facei];

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
            auto valIter = values.cbegin();

            forAll(faces, facei)
            {
                if (useOrigFaceIds)
                {
                    elemId = elemIds[facei];
                }

                const Type v(*valIter);
                ++valIter;

                label nValues =
                    max
                    (
                        label(1),
                        (decompOffsets[facei+1] - decompOffsets[facei])
                    );

                while (nValues--)
                {
                    writeFaceValue(os, format, v, ++elemId);
                }
            }
        }

        os  << "ENDDATA" << endl;
    }

    wroteGeom_ = true;
    return outputFile;
}


// ************************************************************************* //

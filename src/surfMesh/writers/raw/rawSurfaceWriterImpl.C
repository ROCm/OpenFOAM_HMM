/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2011, 2015-2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2011-2014 OpenFOAM Foundation
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
#include "OSspecific.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{
    // Emit each component
    template<class Type>
    static inline void writeData(Ostream& os, const Type& val)
    {
        for (direction i=0; i < pTraits<Type>::nComponents; ++i)
        {
            os  << ' ' << component(val, i);
        }
        os  << nl;
    }

    template<class Type>
    static inline void writeHeader(Ostream& os, const word& fieldName) {}

    template<>
    void writeHeader<label>(Ostream& os, const word& fieldName)
    {
        os  << "# x  y  z"
            << "  " << fieldName << nl;
    }

    template<>
    void writeHeader<scalar>(Ostream& os, const word& fieldName)
    {
        os  << "# x  y  z"
            << "  " << fieldName << nl;
    }

    template<>
    void writeHeader<vector>(Ostream& os, const word& fieldName)
    {
        os  << "# x  y  z"
            << "  " << fieldName << "_x"
            << "  " << fieldName << "_y"
            << "  " << fieldName << "_z"
            << nl;
    }

    template<>
    void writeHeader<sphericalTensor>(Ostream& os, const word& fieldName)
    {
        os  << "# x  y  z"
            << "  " << fieldName << "_ii" << nl;
    }

    template<>
    void writeHeader<symmTensor>(Ostream& os, const word& fieldName)
    {
        // This may not be quite right (not sure what people might need)
        os  << "#  xx  xy  xz  yy  yz  zz";
        for (direction i=0; i < pTraits<symmTensor>::nComponents; ++i)
        {
            os  << "  " << fieldName << '_' << int(i);
        }
        os  << nl;
    }

    template<>
    void writeHeader<tensor>(Ostream& os, const word& fieldName)
    {
        // This may not be quite right (not sure what people might need)
        os  << "#  xx  xy  xz  yx  yy  yz  zx  zy  zz";
        for (direction i=0; i < pTraits<tensor>::nComponents; ++i)
        {
            os  << "  " << fieldName << '_' << int(i);
        }
        os  << nl;
    }

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::fileName Foam::surfaceWriters::rawWriter::writeTemplate
(
    const word& fieldName,
    const Field<Type>& localValues
)
{
    checkOpen();

    // Field:  rootdir/<TIME>/<field>_surfaceName.raw

    fileName outputFile = outputPath_.path();
    if (useTimeDir() && !timeName().empty())
    {
        // Splice in time-directory
        outputFile /= timeName();
    }

    // Append <field>_surfaceName.raw
    outputFile /= fieldName + '_' + outputPath_.name();
    outputFile.ext("raw");

    if (verbose_)
    {
        Info<< "Writing field " << fieldName << " to " << outputFile << endl;
    }


    // geometry merge() implicit
    tmp<Field<Type>> tfield = mergeField(localValues);

    const meshedSurf& surf = surface();

    if (Pstream::master() || !parallel_)
    {
        const auto& values = tfield();
        const pointField& points = surf.points();
        const faceList& faces = surf.faces();

        if (!isDir(outputFile.path()))
        {
            mkDir(outputFile.path());
        }

        OFstream os
        (
            outputFile,
            IOstream::ASCII,
            IOstream::currentVersion,
            writeCompression_
        );

        // Header
        {
            os  << "# " << fieldName;
            if (this->isPointData())
            {
                os  << "  POINT_DATA ";
            }
            else
            {
                os << "  FACE_DATA ";
            }
            os << values.size() << nl;

            // #  x  y  z  field
            writeHeader<Type>(os, fieldName);
        }


        if (this->isPointData())
        {
            // Node values
            forAll(values, elemi)
            {
                writePoint(os, points[elemi]);
                writeData(os, values[elemi]);
            }
        }
        else
        {
            // Face values
            forAll(values, elemi)
            {
                writePoint(os, faces[elemi].centre(points));
                writeData(os,  values[elemi]);
            }
        }
    }

    wroteGeom_ = true;
    return outputFile;
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015-2018 OpenCFD Ltd.
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
Foam::fileName Foam::rawSurfaceWriter::writeTemplate
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
    // field:  rootdir/time/<field>_surfaceName.raw

    const pointField& points = surf.points();
    const faceList&    faces = surf.faces();

    if (!isDir(outputDir))
    {
        mkDir(outputDir);
    }

    OFstream os(outputDir/fieldName + '_' + surfaceName + ".raw");

    if (verbose)
    {
        Info<< "Writing field " << fieldName << " to " << os.name() << endl;
    }

    // Header
    os  << "# " << fieldName;
    if (isNodeValues)
    {
        os  << "  POINT_DATA " << values.size() << nl;
    }
    else
    {
        os  << "  FACE_DATA " << values.size() << nl;
    }

    // Header
    // #  x  y  z  field
    writeHeader<Type>(os, fieldName);

    if (isNodeValues)
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
            writeData(os, values[elemi]);
        }
    }

    return os.name();
}


// ************************************************************************* //

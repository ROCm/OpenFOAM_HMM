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

#include "foamVtkOutputOptions.H"
#include "OSspecific.H"
#include <fstream>

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{
    template<class Type>
    static inline void beginFloatField
    (
        vtk::formatter& format,
        const word& fieldName,
        const Field<Type>& values
    )
    {
        // Use 'double' instead of 'float' for ASCII output (issue #891)
        if (format.opts().ascii())
        {
            format.os()
                << fieldName << ' '
                << int(pTraits<Type>::nComponents) << ' '
                << values.size() << " double" << nl;
        }
        else
        {
            format.os()
                << fieldName << ' '
                << int(pTraits<Type>::nComponents) << ' '
                << values.size() << " float" << nl;
        }
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// // Unspecialized (unknown) data type - map as zeros
// void Foam::vtkSurfaceWriter::writeZeros
// (
//     vtk::formatter& format,
//     const label count
// )
// {
//     const float val(0);
//
//     for (label i=0; i < count; ++i)
//     {
//         format.write(val);
//     }
//     format.flush();
// }


template<class Type>
void Foam::vtkSurfaceWriter::writeField
(
    vtk::formatter& format,
    const Field<Type>& values
)
{
    vtk::writeList(format, values);
    format.flush();
}


template<class Type>
Foam::fileName Foam::vtkSurfaceWriter::writeTemplate
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
    // field:  rootdir/time/<field>_surfaceName.{vtk|vtp}

    fileName outputFile(outputDir/fieldName + '_' + surfaceName + ".vtk");

    if (!isDir(outputDir))
    {
        mkDir(outputDir);
    }
    if (verbose)
    {
        Info<< "Writing field " << fieldName << " to "
            << outputFile << endl;
    }

    // As vtk::outputOptions
    vtk::outputOptions opts(static_cast<vtk::formatType>(fmtType_));
    opts.legacy(true);
    opts.precision(precision_);

    std::ofstream os(outputFile);

    auto format = opts.newFormatter(os);

    writeGeometry(*format, surf);

    if (isNodeValues)
    {
        format->os()
            << "POINT_DATA "
            << values.size() << nl
            << "FIELD attributes 1" << nl;
    }
    else
    {
        format->os()
            << "CELL_DATA "
            << values.size() << nl
            << "FIELD attributes 1" << nl;
    }

    beginFloatField(*format, fieldName, values);
    writeField(*format, values);

    return outputFile;
}


// ************************************************************************* //

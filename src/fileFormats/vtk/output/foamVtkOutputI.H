/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * Specializations * * * * * * * * * * * * * * * //

namespace Foam
{
namespace foamVtkOutput
{

//- Template specialization for label
template<>
inline void write<label>(foamVtkOutput::formatter& fmt, const label& val)
{
    fmt.write(val);
}


//- Template specialization for float
template<>
inline void write<float>(foamVtkOutput::formatter& fmt, const float& val)
{
    fmt.write(val);
}


//- Template specialization for double
template<>
inline void write<double>(foamVtkOutput::formatter& fmt, const double& val)
{
    fmt.write(val);
}

} // End namespace foamVtkOutput
} // End namespace Foam


// * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * * //

inline std::ostream& Foam::foamVtkOutput::legacy::fileHeader
(
    std::ostream& os,
    const std::string& title,
    const bool binary
)
{
    os  << "# vtk DataFile Version 2.0" << nl
        << title << nl
        << (binary ? "BINARY" : "ASCII") << nl;

    return os;
}


inline void Foam::foamVtkOutput::legacy::fileHeader
(
    foamVtkOutput::formatter& fmt,
    const std::string& title,
    const vtkFileTag& contentTypeTag
)
{
    fileHeader(fmt, title, contentNames[contentTypeTag]);
}


inline void Foam::foamVtkOutput::legacy::beginPoints
(
    std::ostream& os,
    const label nPoints
)
{
    os  << "POINTS " << nPoints << " float" << nl;
}


inline void Foam::foamVtkOutput::legacy::beginPolys
(
    std::ostream& os,
    const label nPolys,
    const label nConnectivity
)
{
    os  << "POLYGONS " << nPolys << ' ' << (nPolys + nConnectivity) << nl;
}


inline void Foam::foamVtkOutput::legacy::dataHeader
(
    std::ostream& os,
    const vtkFileTag& dataTypeTag,
    const label nEntries,
    const label nFields
)
{
    os  << dataTypeNames[dataTypeTag] << ' ' << nEntries << nl
        << "FIELD attributes " << nFields << nl;
}


inline void Foam::foamVtkOutput::legacy::floatField
(
    std::ostream& os,
    const word& fieldName,
    const label nCmpt,
    const label nEntries
)
{
    os  << fieldName << ' ' << nCmpt << ' ' << nEntries << " float" << nl;
}


inline void Foam::foamVtkOutput::legacy::intField
(
    std::ostream& os,
    const word& fieldName,
    const label nCmpt,
    const label nEntries
)
{
    os  << fieldName << ' ' << nCmpt << ' ' << nEntries << " int" << nl;
}


// ************************************************************************* //

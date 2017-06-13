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
namespace vtk
{

//- Template specialization for label
template<>
inline void write<label>(vtk::formatter& fmt, const label& val)
{
    fmt.write(val);
}


//- Template specialization for float
template<>
inline void write<float>(vtk::formatter& fmt, const float& val)
{
    fmt.write(val);
}


//- Template specialization for double
template<>
inline void write<double>(vtk::formatter& fmt, const double& val)
{
    fmt.write(val);
}

} // End namespace vtk
} // End namespace Foam


// * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * * //

inline std::ostream& Foam::vtk::legacy::fileHeader
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


inline void Foam::vtk::legacy::fileHeader
(
    vtk::formatter& fmt,
    const std::string& title,
    const vtk::fileTag& contentTypeTag
)
{
    fileHeader(fmt, title, contentNames[contentTypeTag]);
}


inline void Foam::vtk::legacy::beginPoints
(
    std::ostream& os,
    const label nPoints
)
{
    os  << "POINTS " << nPoints << " float" << nl;
}


inline void Foam::vtk::legacy::beginPolys
(
    std::ostream& os,
    const label nPolys,
    const label nConnectivity
)
{
    os  << "POLYGONS " << nPolys << ' ' << (nPolys + nConnectivity) << nl;
}


inline void Foam::vtk::legacy::dataHeader
(
    std::ostream& os,
    const vtk::fileTag& dataTypeTag,
    const label nEntries,
    const label nFields
)
{
    os  << dataTypeNames[dataTypeTag] << ' ' << nEntries << nl
        << "FIELD attributes " << nFields << nl;
}


inline void Foam::vtk::legacy::floatField
(
    std::ostream& os,
    const word& fieldName,
    const int nCmpt,
    const label nEntries
)
{
    os  << fieldName << ' ' << nCmpt << ' ' << nEntries << " float" << nl;
}


inline void Foam::vtk::legacy::intField
(
    std::ostream& os,
    const word& fieldName,
    const int nCmpt,
    const label nEntries
)
{
    os  << fieldName << ' ' << nCmpt << ' ' << nEntries << " int" << nl;
}


// ************************************************************************* //

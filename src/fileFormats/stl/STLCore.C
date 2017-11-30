/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016-2017 OpenCFD Ltd.
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

#include "STLCore.H"
#include "gzstream.h"
#include "OSspecific.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

//! \cond fileScope

// The number of bytes in the STL binary header
static const unsigned STLHeaderSize = 80;

// Check if "SOLID" or "solid" appears as the first non-space content.
// Assume that any leading space is less than 75 chars or so, otherwise
// it is really bad input.
static bool startsWithSolid(const char header[STLHeaderSize])
{
    unsigned pos = 0;
    while (std::isspace(header[pos]) && pos < STLHeaderSize)
    {
        ++pos;
    }

    return
    (
        pos < (STLHeaderSize-5)  // At least 5 chars remaining
     && std::toupper(header[pos+0]) == 'S'
     && std::toupper(header[pos+1]) == 'O'
     && std::toupper(header[pos+2]) == 'L'
     && std::toupper(header[pos+3]) == 'I'
     && std::toupper(header[pos+4]) == 'D'
    );
}
//! \endcond


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fileFormats::STLCore::isBinaryName
(
    const fileName& filename,
    const STLFormat format
)
{
    return
    (
        format == STLFormat::UNKNOWN
      ? (filename.ext() == "stlb")
      : format == STLFormat::BINARY
    );
}


// Check binary by getting the header and number of facets
// this seems to work better than the old token-based method
// - using wordToken can cause an abort if non-word (binary) content
//   is detected ... this is not exactly what we want.
// - some programs (eg, PROSTAR) have 'solid' as the first word in
//   the binary header. This is just wrong and not our fault.
int Foam::fileFormats::STLCore::detectBinaryHeader
(
    const fileName& filename
)
{
    bool compressed = false;
    autoPtr<std::istream> streamPtr
    (
        new std::ifstream(filename, std::ios::binary)
    );

    // If the file is compressed, decompress it before further checking.
    if (!streamPtr->good() && isFile(filename + ".gz", false))
    {
        compressed = true;
        streamPtr.reset(new igzstream((filename + ".gz").c_str()));
    }
    std::istream& is = streamPtr();

    if (!is.good())
    {
        FatalErrorInFunction
            << "Cannot read file " << filename
            << " or file " << filename + ".gz"
            << exit(FatalError);
    }

    // Read the STL header
    char header[STLHeaderSize];
    is.read(header, STLHeaderSize);

    // If the stream is bad, it can't be a binary STL
    if (!is.good() || startsWithSolid(header))
    {
        return 0;
    }


    // Read the number of triangles in the STL file
    // (note: read as signed so we can check whether >2^31)
    int32_t nTris;
    is.read(reinterpret_cast<char*>(&nTris), sizeof(int32_t));

    // Check that stream is OK and number of triangles is positive,
    // if not this may be an ASCII file
    //
    // Also compare the file size with that expected from the number of tris
    // If the comparison is still not sensible then it may be an ASCII file
    if (!is || nTris < 0)
    {
        return 0;
    }
    else if (!compressed)
    {
        const off_t dataFileSize = Foam::fileSize(filename);

        if
        (
            nTris < int(dataFileSize - STLHeaderSize)/50
         || nTris > int(dataFileSize - STLHeaderSize)/25
        )
        {
            return 0;
        }
    }

    // looks like it might be BINARY, return number of triangles
    return nTris;
}


Foam::autoPtr<std::istream>
Foam::fileFormats::STLCore::readBinaryHeader
(
    const fileName& filename,
    label& nTrisEstimated
)
{
    bool bad = false;
    bool compressed = false;
    nTrisEstimated = 0;

    autoPtr<std::istream> streamPtr
    (
        new std::ifstream(filename, std::ios::binary)
    );

    // If the file is compressed, decompress it before reading.
    if (!streamPtr->good() && isFile(filename + ".gz", false))
    {
        compressed = true;
        streamPtr.reset(new igzstream((filename + ".gz").c_str()));
    }
    std::istream& is = streamPtr();

    if (!is.good())
    {
        streamPtr.clear();

        FatalErrorInFunction
            << "Cannot read file " << filename
            << " or file " << filename + ".gz"
            << exit(FatalError);
    }


    // Read the STL header
    char header[STLHeaderSize];
    is.read(header, STLHeaderSize);

    // Check that stream is OK, if not this may be an ASCII file
    if (!is.good()) // could check again: startsWithSolid(header)
    {
        streamPtr.clear();

        FatalErrorInFunction
            << "problem reading header, perhaps file is not binary "
            << exit(FatalError);
    }

    // Read the number of triangles in the STl file
    // (note: read as int so we can check whether >2^31)
    int32_t nTris;
    is.read(reinterpret_cast<char*>(&nTris), sizeof(int32_t));

    // Check that stream is OK and number of triangles is positive,
    // if not this maybe an ASCII file
    //
    // Also compare the file size with that expected from the number of tris
    // If the comparison is not sensible then it may be an ASCII file
    if (!is || nTris < 0)
    {
        bad = true;
    }
    else if (!compressed)
    {
        const off_t dataFileSize = Foam::fileSize(filename);

        if
        (
            nTris < int(dataFileSize - STLHeaderSize)/50
         || nTris > int(dataFileSize - STLHeaderSize)/25
        )
        {
            bad = true;
        }
    }

    if (bad)
    {
        streamPtr.clear();

        FatalErrorInFunction
            << "problem reading number of triangles, perhaps file is not binary"
            << exit(FatalError);
    }

    nTrisEstimated = nTris;
    return streamPtr;
}


void Foam::fileFormats::STLCore::writeBinaryHeader
(
    ostream& os,
    uint32_t nTris
)
{
    // STL header with extra information about nTris
    char header[STLHeaderSize];
    sprintf(header, "STL binary file %u facets", nTris);

    // avoid trailing junk
    for (size_t i = strlen(header); i < STLHeaderSize; ++i)
    {
        header[i] = 0;
    }

    os.write(header, STLHeaderSize);
    os.write(reinterpret_cast<char*>(&nTris), sizeof(uint32_t));
}


// ************************************************************************* //

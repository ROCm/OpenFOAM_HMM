/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "ensightSurfaceReader.H"
#include "stringOps.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ensightSurfaceReader, 0);
    addToRunTimeSelectionTable(surfaceReader, ensightSurfaceReader, fileName);
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Read and discard specified number of elements
template<class Type>
static inline void discard(label n, ensightReadFile& is)
{
    Type val;

    while (n > 0)
    {
        is.read(val);
        --n;
    }
}

} // End namespace Foam


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::ensightSurfaceReader::skip(const label n, Istream& is) const
{
    label i = 0;
    token tok;
    while (is.good() && (i < n))
    {
        is >> tok;
        ++i;

        DebugInfo
            << "Skipping token " << tok << nl;
    }

    if (i != n)
    {
        WarningInFunction
            << "Requested to skip " << n << " tokens, but stream exited after "
            << i << " tokens. Last token read: " << tok
            << nl;
    }
}


void Foam::ensightSurfaceReader::readLine(IFstream& is, string& line) const
{
    do
    {
        is.getLine(line);

        // Trim out any '#' comments
        const auto pos = line.find('#');
        if (pos != std::string::npos)
        {
            line.erase(pos);
        }
        stringOps::inplaceTrimRight(line);
    }
    while (line.empty() && is.good());
}


void Foam::ensightSurfaceReader::debugSection
(
    const word& expected,
    IFstream& is
) const
{
    string actual;
    readLine(is, actual);

    if (expected != actual)
    {
        FatalIOErrorInFunction(is)
            << "Expected section header '" << expected
            << "' but read " << actual << nl
            << exit(FatalIOError);
    }

    DebugInfo
        << "Read section header: " << expected << nl;
}


Foam::fileName Foam::ensightSurfaceReader::replaceMask
(
    const fileName& fName,
    const label timeIndex
)
{
    fileName result(fName);

    const auto nMask = stringOps::count(fName, '*');

    // If there are any '*' chars, they are assumed to be contiguous
    // Eg, data/******/geometry

    if (nMask)
    {
        std::ostringstream oss;
        oss << std::setfill('0') << std::setw(nMask) << timeIndex;

        const std::string maskStr(nMask, '*');
        const std::string indexStr = oss.str();
        result.replace(maskStr, indexStr);
    }

    return result;
}


Foam::Pair<Foam::ensightSurfaceReader::idTypes>
Foam::ensightSurfaceReader::readGeometryHeader(ensightReadFile& is) const
{
    // Binary flag string if applicable
    is.readBinaryHeader();

    string buffer;

    Pair<idTypes> idHandling(idTypes::NONE, idTypes::NONE);

    // Ensight Geometry File
    is.read(buffer);
    DebugInfo<< "buffer [" << buffer.length() << "] " << buffer << nl;

    // Description - 1
    is.read(buffer);
    DebugInfo<< "buffer [" << buffer.length() << "] " << buffer << nl;

    // "node id (off|assign|given|ignore)" - "given" is not actually supported
    is.read(buffer);
    DebugInfo<< "buffer [" << buffer.length() << "] " << buffer << nl;

    if (buffer.find("ignore") != std::string::npos)
    {
        idHandling.first() = idTypes::IGNORE;
    }
    else if (buffer.find("given") != std::string::npos)
    {
        idHandling.first() = idTypes::GIVEN;
    }

    // "element id (off|assign|given|ignore)"
    is.read(buffer);
    DebugInfo<< "buffer [" << buffer.length() << "] " << buffer << nl;

    if (buffer.find("ignore") != std::string::npos)
    {
        idHandling.second() = idTypes::IGNORE;
    }
    else if (buffer.find("given") != std::string::npos)
    {
        idHandling.second() = idTypes::GIVEN;
    }


    // "part" - but could also be an optional "extents"
    is.read(buffer);
    DebugInfo<< "buffer [" << buffer.length() << "] " << buffer << nl;

    if (buffer.find("extents") != std::string::npos)
    {
        // Optional extents - read and discard 6 floats
        // (xmin, xmax, ymin, ymax, zmin, zmax)

        discard<scalar>(6, is);

        // Part
        is.read(buffer);
        DebugInfo<< "buffer [" << buffer.length() << "] " << buffer << nl;
    }

    // The part number
    label ivalue;
    is.read(ivalue);
    DebugInfo<< "ivalue: " << ivalue << nl;

    // Part description / name
    is.read(buffer);
    DebugInfo<< "buffer [" << buffer.length() << "] " << buffer << nl;

    // "coordinates"
    is.read(buffer);
    DebugInfo<< "buffer [" << buffer.length() << "] " << buffer << nl;

    return idHandling;
}


void Foam::ensightSurfaceReader::readCase(IFstream& is)
{
    DebugInFunction << endl;

    if (!is.good())
    {
        FatalErrorInFunction
            << "Cannot read file " << is.name()
            << exit(FatalError);
    }

    string buffer;

    // Read the file

    debugSection("FORMAT", is);
    readLine(is, buffer);  // type: ensight gold

    debugSection("GEOMETRY", is);
    readLine(is, buffer);

    // GEOMETRY with any of these
    //     model: 1 xxx.0000.mesh
    //     model: xxx.0000.mesh
    //     model: data/directory/geometry
    //
    // - use the last entry
    meshFileName_ = stringOps::splitSpace(buffer).last().str();

    DebugInfo << "mesh file:" << meshFileName_ << endl;

    debugSection("VARIABLE", is);

    // Read the field description
    DynamicList<word> fieldNames(16);
    DynamicList<string> fieldFileNames(16);

    while (is.good())
    {
        readLine(is, buffer);

        if (buffer == "TIME")
        {
            break;
        }

        // Read the field name and associated file name. Eg,
        //     scalar per element: 1  p  data/********/p

        const auto parsed = stringOps::splitSpace(buffer);

        if (buffer.find(':') == string::npos || parsed.size() < 4)
        {
            WarningInFunction
                << "Error reading field file name. Current buffer: "
                << buffer << endl;
            continue;
        }
        else if (debug)
        {
            Info<< "variable line: " << parsed.size();
            for (const auto& s : parsed)
            {
                Info<< " " << s.str();
            }
            Info<< nl;
        }

        fieldNames.append(parsed[parsed.size()-2].str());
        fieldFileNames.append(parsed.last().str());
    }
    fieldNames_.transfer(fieldNames);
    fieldFileNames_.transfer(fieldFileNames);

    DebugInfo
        << "fieldNames: " << fieldNames_ << nl
        << "fieldFileNames: " << fieldFileNames_ << nl;

    // Start reading time information
    readLine(is, buffer);                       // time set: <int>

    readLine(is, buffer);
    readFromLine(3, buffer, nTimeSteps_);       // number of steps: <int>
    readLine(is, buffer);
    readFromLine(3, buffer, timeStartIndex_);   // filename start number: <int>
    readLine(is, buffer);
    readFromLine(2, buffer, timeIncrement_);    // filename increment: <int>

    DebugInfo
        << "nTimeSteps: " << nTimeSteps_ << nl
        << "timeStartIndex: " << timeStartIndex_ << nl
        << "timeIncrement: " << timeIncrement_ << nl;

    // Read the time values
    readLine(is, buffer); // time values:
    timeValues_.setSize(nTimeSteps_);
    for (label i = 0; i < nTimeSteps_; ++i)
    {
        scalar t(readScalar(is));

        timeValues_[i].value() = t;
        // TODO: use character representation of t directly instead of
        // regenerating from scalar value
        timeValues_[i].name() = Foam::name(t);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ensightSurfaceReader::ensightSurfaceReader(const fileName& fName)
:
    surfaceReader(fName),
    streamFormat_(IOstream::ASCII),
    baseDir_(fName.path()),
    meshFileName_(),
    fieldNames_(),
    fieldFileNames_(),
    nTimeSteps_(0),
    timeStartIndex_(0),
    timeIncrement_(1),
    timeValues_(),
    surfPtr_(nullptr)
{
    IFstream is(fName);
    readCase(is);
}


// * * * * * * * * * * * * * Public Member Functions   * * * * * * * * * * * //

const Foam::meshedSurface& Foam::ensightSurfaceReader::geometry
(
    const label timeIndex
)
{
    DebugInFunction << endl;

    if (!surfPtr_)
    {
        fileName meshInstance(replaceMask(meshFileName_, timeIndex));
        IFstream isBinary(baseDir_/meshInstance, IOstream::BINARY);

        if (!isBinary.good())
        {
            FatalErrorInFunction
                << "Cannot read file " << isBinary.name()
                << exit(FatalError);
        }

        streamFormat_ = IOstream::BINARY;
        {
            istream& iss = isBinary.stdStream();

            // Binary string is *exactly* 80 characters
            string buf(size_t(80), '\0');
            iss.read(&buf[0], 80);

            if (!iss)
            {
                // Truncated?
                buf.erase(iss.gcount());
            }

            // Truncate at the first embedded '\0'
            const auto endp = buf.find('\0');
            if (endp != std::string::npos)
            {
                buf.erase(endp);
            }

            // Contains "C Binary" ?
            if
            (
                (buf.find("binary") == std::string::npos)
             && (buf.find("Binary") == std::string::npos)
            )
            {
                streamFormat_ = IOstream::ASCII;
            }
        }

        if (debug)
        {
            Info<< "stream format: ";
            if (streamFormat_ == IOstream::ASCII)
            {
                Info<< "ascii" << endl;
            }
            else
            {
                Info<< "binary" << endl;
            }
        }


        ensightReadFile is(baseDir_/meshInstance, streamFormat_);

        DebugInfo
            << "File: " << is.name() << nl;

        Pair<idTypes> idHandling = readGeometryHeader(is);

        label nPoints;
        is.read(nPoints);

        DebugInfo
            << "nPoints: " << nPoints << nl;

        if (idHandling.first() == idTypes::GIVEN)
        {
            WarningInFunction
                << "Treating node id 'given' as being 'ignore'" << nl
                << "If something fails, this could be the reason" << nl
                << endl;

            idHandling.first() = idTypes::IGNORE;
        }

        if (idHandling.first() == idTypes::IGNORE)
        {
            DebugInfo
                << "Ignore " << nPoints << " node ids" << nl;

            // Read and discard labels
            discard<label>(nPoints, is);
        }

        pointField points(nPoints);
        for (direction cmpt = 0; cmpt < vector::nComponents; ++cmpt)
        {
            for (point& pt : points)
            {
                is.read(pt[cmpt]);
            }
        }


        // Read faces - may be a mix of tris, quads and polys
        DynamicList<face> faces(ceil(nPoints/3));
        DynamicList<Tuple2<string, label>> schema(faces.size());
        string faceType;
        while (is.good()) // (is.peek() != EOF)
        {
            is.read(faceType);

            if (!is.good())
            {
                break;
            }

            label nFace = 0;

            if (faceType == "tria3")
            {
                is.read(nFace);

                DebugInfo
                    << "faceType <" << faceType.c_str() << "> count: "
                    << nFace << nl;

                if
                (
                    idHandling.second() == idTypes::IGNORE
                 || idHandling.second() == idTypes::GIVEN
                )
                {
                    DebugInfo
                        << "Ignore " << nFace << " element ids" << nl;

                    // Read and discard labels
                    discard<label>(nFace, is);
                }

                face f(3);
                for (label facei = 0; facei < nFace; ++facei)
                {
                    for (label& fp : f)
                    {
                        is.read(fp);
                    }

                    faces.append(f);
                }
            }
            else if (faceType == "quad4")
            {
                is.read(nFace);

                DebugInfo
                    << "faceType <" << faceType.c_str() << "> count: "
                    << nFace << nl;

                if
                (
                    idHandling.second() == idTypes::IGNORE
                 || idHandling.second() == idTypes::GIVEN
                )
                {
                    DebugInfo
                        << "Ignore " << nFace << " element ids" << nl;

                    // Read and discard labels
                    discard<label>(nFace, is);
                }

                face f(4);
                for (label facei = 0; facei < nFace; ++facei)
                {
                    for (label& fp : f)
                    {
                        is.read(fp);
                    }

                    faces.append(f);
                }
            }
            else if (faceType == "nsided")
            {
                is.read(nFace);

                DebugInfo
                    << "faceType <" << faceType.c_str() << "> count: "
                    << nFace << nl;

                if
                (
                    idHandling.second() == idTypes::IGNORE
                 || idHandling.second() == idTypes::GIVEN
                )
                {
                    DebugInfo
                        << "Ignore " << nFace << " element ids" << nl;

                    // Read and discard labels
                    discard<label>(nFace, is);
                }

                labelList np(nFace);
                for (label facei = 0; facei < nFace; ++facei)
                {
                    is.read(np[facei]);
                }
                for (label facei = 0; facei < nFace; ++facei)
                {
                    face f(np[facei]);
                    for (label& fp : f)
                    {
                        is.read(fp);
                    }

                    faces.append(f);
                }
            }
            else
            {
                if (debug)
                {
                    WarningInFunction
                        << "Unknown face type: <" << faceType.c_str()
                        << ">. Stopping read and continuing with current "
                        << "elements only" << endl;
                }
                break;
            }
            schema.append(Tuple2<string, label>(faceType, nFace));
        }

        schema_.transfer(schema);

        DebugInfo
            << "read nFaces: " << faces.size() << nl
            << "file schema: " << schema_ << nl;

        // Convert from 1-based Ensight addressing to 0-based OF addressing
        for (face& f : faces)
        {
            for (label& pointi : f)
            {
                --pointi;
            }
        }

        surfPtr_.reset(new meshedSurface(std::move(points), std::move(faces)));
    }

    return *surfPtr_;
}


Foam::instantList Foam::ensightSurfaceReader::times() const
{
    return timeValues_;
}


Foam::wordList Foam::ensightSurfaceReader::fieldNames
(
    const label timeIndex
) const
{
    return fieldNames_;
}


Foam::tmp<Foam::Field<Foam::scalar>> Foam::ensightSurfaceReader::field
(
    const label timeIndex,
    const label fieldIndex,
    const scalar& refValue
) const
{
    return readField<scalar>(timeIndex, fieldIndex);
}


Foam::tmp<Foam::Field<Foam::vector>> Foam::ensightSurfaceReader::field
(
    const label timeIndex,
    const label fieldIndex,
    const vector& refValue
) const
{
    return readField<vector>(timeIndex, fieldIndex);
}


Foam::tmp<Foam::Field<Foam::sphericalTensor>>
Foam::ensightSurfaceReader::field
(
    const label timeIndex,
    const label fieldIndex,
    const sphericalTensor& refValue
) const
{
    return readField<sphericalTensor>(timeIndex, fieldIndex);
}


Foam::tmp<Foam::Field<Foam::symmTensor>> Foam::ensightSurfaceReader::field
(
    const label timeIndex,
    const label fieldIndex,
    const symmTensor& refValue
) const
{
    return readField<symmTensor>(timeIndex, fieldIndex);
}


Foam::tmp<Foam::Field<Foam::tensor>> Foam::ensightSurfaceReader::field
(
    const label timeIndex,
    const label fieldIndex,
    const tensor& refValue
) const
{
    return readField<tensor>(timeIndex, fieldIndex);
}


// ************************************************************************* //

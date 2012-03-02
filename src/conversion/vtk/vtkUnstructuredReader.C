/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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

#include "vtkUnstructuredReader.H"
#include "labelIOField.H"
#include "scalarIOField.H"
#include "stringIOList.H"
#include "cellModeller.H"
#include "vectorIOField.H"

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

defineTypeNameAndDebug(Foam::vtkUnstructuredReader, 0);

template<>
const char*
Foam::NamedEnum<Foam::vtkUnstructuredReader::vtkDataType, 3>::names[] =
{
    "int",
    "float",
    "string"
};
const Foam::NamedEnum<Foam::vtkUnstructuredReader::vtkDataType, 3>
Foam::vtkUnstructuredReader::vtkDataTypeNames;


template<>
const char*
Foam::NamedEnum<Foam::vtkUnstructuredReader::vtkDataSetType, 3>::names[] =
{
    "FIELD",
    "SCALARS",
    "VECTORS"
};
const Foam::NamedEnum<Foam::vtkUnstructuredReader::vtkDataSetType, 3>
Foam::vtkUnstructuredReader::vtkDataSetTypeNames;


template<>
const char*
Foam::NamedEnum<Foam::vtkUnstructuredReader::parseMode, 5>::names[] =
{
    "NOMODE",
    "UNSTRUCTURED_GRID",
    "POLYDATA",
    "CELL_DATA",
    "POINT_DATA"
};
const Foam::NamedEnum<Foam::vtkUnstructuredReader::parseMode, 5>
Foam::vtkUnstructuredReader::parseModeNames;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class T>
void Foam::vtkUnstructuredReader::readBlock
(
    Istream& inFile,
    const label n,
    List<T>& lst
) const
{
    lst.setSize(n);
    forAll(lst, i)
    {
        inFile >> lst[i];
    }
}


void Foam::vtkUnstructuredReader::warnUnhandledType
(
    Istream& inFile,
    const label type,
    labelHashSet& warningGiven
) const
{
    if (warningGiven.insert(type))
    {
        IOWarningIn("vtkUnstructuredReader::warnUnhandledType(..)", inFile)
            << "Skipping unknown cell type " << type << endl;
    }
}


// Split cellTypes into cells, faces and lines
void Foam::vtkUnstructuredReader::extractCells
(
    Istream& inFile,
    const labelList& cellTypes,
    const labelList& cellVertData
)
{
    const cellModel& hex = *(cellModeller::lookup("hex"));
    const cellModel& prism = *(cellModeller::lookup("prism"));
    const cellModel& pyr = *(cellModeller::lookup("pyr"));
    const cellModel& tet = *(cellModeller::lookup("tet"));

    labelList tetPoints(4);
    labelList pyrPoints(5);
    labelList prismPoints(6);
    labelList hexPoints(8);


    cells_.setSize(cellTypes.size());
    faces_.setSize(cellTypes.size());
    lines_.setSize(cellTypes.size());

    label dataIndex = 0;
    label cellI = 0;
    label faceI = 0;
    label lineI = 0;


    // To mark whether unhandled type has been visited.
    labelHashSet warningGiven;

    forAll(cellTypes, i)
    {
        switch (cellTypes[i])
        {
            case VTK_VERTEX:
            {
                warnUnhandledType(inFile, cellTypes[i], warningGiven);
                label nRead = cellVertData[dataIndex++];
                if (nRead != 1)
                {
                    FatalIOErrorIn
                    (
                        "vtkUnstructuredReader::extractCells(..)",
                        inFile
                    )   << "Expected size 1 for VTK_VERTEX but found "
                        << nRead << exit(FatalIOError);
                }
                dataIndex += nRead;
            }
            break;

            case VTK_POLY_VERTEX:
            {
                warnUnhandledType(inFile, cellTypes[i], warningGiven);
                label nRead = cellVertData[dataIndex++];
                dataIndex += nRead;
            }
            break;

            case VTK_LINE:
            {
                //warnUnhandledType(inFile, cellTypes[i], warningGiven);
                label nRead = cellVertData[dataIndex++];
                if (nRead != 2)
                {
                    FatalIOErrorIn
                    (
                        "vtkUnstructuredReader::extractCells(..)",
                        inFile
                    )   << "Expected size 2 for VTK_LINE but found "
                        << nRead << exit(FatalIOError);
                }
                labelList& segment = lines_[lineI++];
                segment.setSize(2);
                segment[0] = cellVertData[dataIndex++];
                segment[1] = cellVertData[dataIndex++];
            }
            break;

            case VTK_POLY_LINE:
            {
                //warnUnhandledType(inFile, cellTypes[i], warningGiven);
                label nRead = cellVertData[dataIndex++];
                labelList& segment = lines_[lineI++];
                segment.setSize(nRead);
                forAll(segment, i)
                {
                    segment[i] = cellVertData[dataIndex++];
                }
            }
            break;

            case VTK_TRIANGLE:
            {
                face& f = faces_[faceI++];
                f.setSize(3);
                label nRead = cellVertData[dataIndex++];
                if (nRead != 3)
                {
                    FatalIOErrorIn
                    (
                        "vtkUnstructuredReader::extractCells(..)",
                        inFile
                    )   << "Expected size 3 for VTK_TRIANGLE but found "
                        << nRead << exit(FatalIOError);
                }
                f[0] = cellVertData[dataIndex++];
                f[1] = cellVertData[dataIndex++];
                f[2] = cellVertData[dataIndex++];
            }
            break;

            case VTK_QUAD:
            {
                face& f = faces_[faceI++];
                f.setSize(4);
                label nRead = cellVertData[dataIndex++];
                if (nRead != 4)
                {
                    FatalIOErrorIn
                    (
                        "vtkUnstructuredReader::extractCells(..)",
                        inFile
                    )   << "Expected size 4 for VTK_QUAD but found "
                        << nRead << exit(FatalIOError);
                }
                f[0] = cellVertData[dataIndex++];
                f[1] = cellVertData[dataIndex++];
                f[2] = cellVertData[dataIndex++];
                f[3] = cellVertData[dataIndex++];
            }
            break;

            case VTK_POLYGON:
            {
                face& f = faces_[faceI++];
                label nRead = cellVertData[dataIndex++];
                f.setSize(nRead);
                forAll(f, fp)
                {
                    f[fp] = cellVertData[dataIndex++];
                }
            }
            break;

            case VTK_TETRA:
            {
                label nRead = cellVertData[dataIndex++];
                if (nRead != 4)
                {
                    FatalIOErrorIn
                    (
                        "vtkUnstructuredReader::extractCells(..)",
                        inFile
                    )   << "Expected size 4 for VTK_TETRA but found "
                        << nRead << exit(FatalIOError);
                }
                tetPoints[0] = cellVertData[dataIndex++];
                tetPoints[1] = cellVertData[dataIndex++];
                tetPoints[2] = cellVertData[dataIndex++];
                tetPoints[3] = cellVertData[dataIndex++];
                cells_[cellI++] = cellShape(tet, tetPoints, true);
            }
            break;

            case VTK_PYRAMID:
            {
                label nRead = cellVertData[dataIndex++];
                if (nRead != 5)
                {
                    FatalIOErrorIn
                    (
                        "vtkUnstructuredReader::extractCells(..)",
                        inFile
                    )   << "Expected size 5 for VTK_PYRAMID but found "
                        << nRead << exit(FatalIOError);
                }
                pyrPoints[0] = cellVertData[dataIndex++];
                pyrPoints[1] = cellVertData[dataIndex++];
                pyrPoints[2] = cellVertData[dataIndex++];
                pyrPoints[3] = cellVertData[dataIndex++];
                pyrPoints[4] = cellVertData[dataIndex++];
                cells_[cellI++] = cellShape(pyr, pyrPoints, true);
            }
            break;

            case VTK_WEDGE:
            {
                label nRead = cellVertData[dataIndex++];
                if (nRead != 6)
                {
                    FatalIOErrorIn
                    (
                        "vtkUnstructuredReader::extractCells(..)",
                        inFile
                    )   << "Expected size 6 for VTK_WEDGE but found "
                        << nRead << exit(FatalIOError);
                }
                prismPoints[0] = cellVertData[dataIndex++];
                prismPoints[1] = cellVertData[dataIndex++];
                prismPoints[2] = cellVertData[dataIndex++];
                prismPoints[3] = cellVertData[dataIndex++];
                prismPoints[4] = cellVertData[dataIndex++];
                prismPoints[5] = cellVertData[dataIndex++];
                cells_[cellI++] = cellShape(prism, prismPoints, true);
            }
            break;

            case VTK_HEXAHEDRON:
            {
                label nRead = cellVertData[dataIndex++];
                if (nRead != 8)
                {
                    FatalIOErrorIn
                    (
                        "vtkUnstructuredReader::extractCells(..)",
                        inFile
                    )   << "Expected size 8 for VTK_HEXAHEDRON but found "
                        << nRead << exit(FatalIOError);
                }
                hexPoints[0] = cellVertData[dataIndex++];
                hexPoints[1] = cellVertData[dataIndex++];
                hexPoints[2] = cellVertData[dataIndex++];
                hexPoints[3] = cellVertData[dataIndex++];
                hexPoints[4] = cellVertData[dataIndex++];
                hexPoints[5] = cellVertData[dataIndex++];
                hexPoints[6] = cellVertData[dataIndex++];
                hexPoints[7] = cellVertData[dataIndex++];
                cells_[cellI++] = cellShape(hex, hexPoints, true);
            }
            break;

            default:
                warnUnhandledType(inFile, cellTypes[i], warningGiven);
                label nRead = cellVertData[dataIndex++];
                dataIndex += nRead;
        }
    }

    if (debug)
    {
        Info<< "Read " << cellI << " cells;" << faceI << " faces." << endl;
    }
    cells_.setSize(cellI);
    faces_.setSize(faceI);
    lines_.setSize(lineI);
}


// Read single field and stores it on the objectRegistry.
void Foam::vtkUnstructuredReader::readField
(
    ISstream& inFile,
    objectRegistry& obj,
    const word& arrayName,
    const word& dataType,
    const label size
) const
{
    switch (vtkDataTypeNames[dataType])
    {
        case VTK_INT:
        {
            autoPtr<labelIOField> fieldVals
            (
                new labelIOField
                (
                    IOobject
                    (
                        arrayName,
                        "",
                        obj
                    ),
                    size
                )
            );
            readBlock(inFile, fieldVals().size(), fieldVals());
            regIOobject::store(fieldVals);
        }
        break;

        case VTK_FLOAT:
        {
            autoPtr<scalarIOField> fieldVals
            (
                new scalarIOField
                (
                    IOobject
                    (
                        arrayName,
                        "",
                        obj
                    ),
                    size
                )
            );
            readBlock(inFile, fieldVals().size(), fieldVals());
            regIOobject::store(fieldVals);
        }
        break;

        case VTK_STRING:
        {
            if (debug)
            {
                Info<< "Reading strings:" << size << endl;
            }
            autoPtr<stringIOList> fieldVals
            (
                new stringIOList
                (
                    IOobject
                    (
                        arrayName,
                        "",
                        obj
                    ),
                    size
                )
            );
            // Consume current line.
            inFile.getLine(fieldVals()[0]);
            // Read without parsing
            forAll(fieldVals(), i)
            {
                inFile.getLine(fieldVals()[i]);
            }
            regIOobject::store(fieldVals);
        }
        break;

        default:
        {
            IOWarningIn("vtkUnstructuredReader::extractCells(..)", inFile)
                << "Unhandled type " << vtkDataTypeNames[dataType] << endl
                << "Skipping " << size
                << " words." << endl;
            scalarField fieldVals;
            readBlock(inFile, size, fieldVals);
        }
        break;
    }
}


// Reads fields, stores them on the objectRegistry. Returns a list of
// read fields
Foam::wordList Foam::vtkUnstructuredReader::readFieldArray
(
    ISstream& inFile,
    objectRegistry& obj,
    const label wantedSize
) const
{
    DynamicList<word> fields;

    word dataName(inFile);
    if (debug)
    {
        Info<< "dataName:" << dataName << endl;
    }
    label numArrays(readLabel(inFile));
    if (debug)
    {
        Pout<< "numArrays:" << numArrays << endl;
    }
    for (label i = 0; i < numArrays; i++)
    {
        word arrayName(inFile);
        label numComp(readLabel(inFile));
        label numTuples(readLabel(inFile));
        word dataType(inFile);

        if (debug)
        {
            Info<< "Reading field " << arrayName
                << " of " << numTuples << " tuples of rank " << numComp << endl;
        }

        if (wantedSize != -1 && numTuples != wantedSize)
        {
            FatalIOErrorIn("vtkUnstructuredReader::readFieldArray(..)", inFile)
                << "Expected " << wantedSize << " tuples but only have "
                << numTuples << exit(FatalIOError);
        }

        readField
        (
            inFile,
            obj,
            arrayName,
            dataType,
            numTuples*numComp
        );
        fields.append(arrayName);
    }
    return fields.shrink();
}


Foam::objectRegistry& Foam::vtkUnstructuredReader::selectRegistry
(
    const parseMode readMode
)
{
    if (readMode == CELL_DATA)
    {
        return cellData_;
    }
    else if (readMode == POINT_DATA)
    {
        return pointData_;
    }
    else
    {
        return otherData_;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::vtkUnstructuredReader::vtkUnstructuredReader
(
    const objectRegistry& obr,
    ISstream& inFile
)
:
    cellData_(IOobject("cellData", obr)),
    pointData_(IOobject("pointData", obr)),
    otherData_(IOobject("otherData", obr))
{
    read(inFile);
}


void Foam::vtkUnstructuredReader::read(ISstream& inFile)
{
    inFile.getLine(header_);
    if (debug)
    {
        Info<< "Header   : " << header_ << endl;
    }
    inFile.getLine(title_);
    if (debug)
    {
        Info<< "Title    : " << title_ << endl;
    }
    inFile.getLine(dataType_);
    if (debug)
    {
        Info<< "dataType : " << dataType_ << endl;
    }

    parseMode readMode = NOMODE;
    label wantedSize = -1;


    // Temporary storage for vertices of cells.
    labelList cellVerts;

    while (inFile.good())
    {
        word tag(inFile);

        if (!inFile.good())
        {
            break;
        }

        if (debug)
        {
            Info<< "line:" << inFile.lineNumber()
                << " tag:" << tag << endl;
        }

        if (tag == "DATASET")
        {
            word geomType(inFile);
            if (debug)
            {
                Info<< "geomType : " << geomType << endl;
            }
            readMode = parseModeNames[geomType];
            wantedSize = -1;
        }
        else if (tag == "POINTS")
        {
            label nPoints(readLabel(inFile));
            points_.setSize(nPoints);    ///3);
            if (debug)
            {
                Info<< "Reading " << nPoints << " numbers representing "
                    << points_.size() << " coordinates." << endl;
            }

            word primitiveTag(inFile);
            if (primitiveTag != "float")
            {
                FatalIOErrorIn("vtkUnstructuredReader::read(..)", inFile)
                    << "Expected 'float' entry but found "
                    << primitiveTag
                    << exit(FatalIOError);
            }
            forAll(points_, i)
            {
                inFile >> points_[i].x() >> points_[i].y() >> points_[i].z();
            }
        }
        else if (tag == "CELLS")
        {
            label nCells(readLabel(inFile));
            label nNumbers(readLabel(inFile));
            if (debug)
            {
                Info<< "Reading " << nCells << " cells or faces." << endl;
            }
            readBlock(inFile, nNumbers, cellVerts);
        }
        else if (tag == "CELL_TYPES")
        {
            label nCellTypes(readLabel(inFile));

            labelList cellTypes;
            readBlock(inFile, nCellTypes, cellTypes);

            if (cellTypes.size() > 0 && cellVerts.size() == 0)
            {
                FatalIOErrorIn("vtkUnstructuredReader::read(..)", inFile)
                    << "Found " << cellTypes.size()
                    << " cellTypes but no cells."
                    << exit(FatalIOError);
            }

            extractCells(inFile, cellTypes, cellVerts);
            cellVerts.clear();
        }
        else if (tag == "LINES")
        {
            label nLines(readLabel(inFile));
            label nNumbers(readLabel(inFile));
            if (debug)
            {
                Info<< "Reading " << nLines << " lines." << endl;
            }
            labelList lineVerts;
            readBlock(inFile, nNumbers, lineVerts);
            lines_.setSize(nLines);
            label elemI = 0;
            forAll(lines_, lineI)
            {
                labelList& f = lines_[lineI];
                f.setSize(lineVerts[elemI++]);
                forAll(f, fp)
                {
                    f[fp] = lineVerts[elemI++];
                }
            }
        }
        else if (tag == "POLYGONS")
        {
            // If in polydata mode

            label nFaces(readLabel(inFile));
            label nNumbers(readLabel(inFile));
            if (debug)
            {
                Info<< "Reading " << nFaces << " faces." << endl;
            }
            labelList faceVerts;
            readBlock(inFile, nNumbers, faceVerts);
            faces_.setSize(nFaces);
            label elemI = 0;
            forAll(faces_, faceI)
            {
                face& f = faces_[faceI];
                f.setSize(faceVerts[elemI++]);
                forAll(f, fp)
                {
                    f[fp] = faceVerts[elemI++];
                }
            }
        }
        else if (tag == "POINT_DATA")
        {
            readMode = POINT_DATA;
            wantedSize = points_.size();

            label nPoints(readLabel(inFile));
            if (nPoints != wantedSize)
            {
                FatalIOErrorIn("vtkUnstructuredReader::read(..)", inFile)
                    << "Reading POINT_DATA : expected " << wantedSize
                    << " but read " << nPoints << exit(FatalIOError);
            }
        }
        else if (tag == "CELL_DATA")
        {
            readMode = CELL_DATA;
            wantedSize = cells_.size()+faces_.size();

            label nCells(readLabel(inFile));
            if (nCells != wantedSize)
            {
                FatalIOErrorIn("vtkUnstructuredReader::read(..)", inFile)
                    << "Reading CELL_DATA : expected "
                    << wantedSize
                    << " but read " << nCells << exit(FatalIOError);
            }
        }
        else if (tag == "FIELD")
        {
            // wantedSize already set according to type we expected to read.
            readFieldArray(inFile, selectRegistry(readMode), wantedSize);
        }
        else if (tag == "SCALARS")
        {
            string line;
            inFile.getLine(line);
            IStringStream is(line);
            word dataName(is);
            word dataType(is);
            //label numComp(readLabel(inFile));

            if (debug)
            {
                Info<< "Reading scalar " << dataName
                    << " of type " << dataType
                    << " from lookup table" << endl;
            }

            word lookupTableTag(inFile);
            if (lookupTableTag != "LOOKUP_TABLE")
            {
                FatalIOErrorIn("vtkUnstructuredReader::read(..)", inFile)
                    << "Expected tag LOOKUP_TABLE but read "
                    << lookupTableTag
                    << exit(FatalIOError);
            }

            word lookupTableName(inFile);

            readField
            (
                inFile,
                selectRegistry(readMode),
                dataName,
                dataType,
                wantedSize//*numComp
            );
        }
        else if (tag == "VECTORS" || tag == "NORMALS")
        {
            string line;
            inFile.getLine(line);
            IStringStream is(line);
            word dataName(is);
            word dataType(is);
            if (debug)
            {
                Info<< "Reading vector " << dataName
                    << " of type " << dataType << endl;
            }

            objectRegistry& reg = selectRegistry(readMode);

            readField
            (
                inFile,
                reg,
                dataName,
                dataType,
                3*wantedSize
            );

            if (vtkDataTypeNames[dataType] == VTK_FLOAT)
            {
                objectRegistry::iterator iter = reg.find(dataName);
                scalarField s(*dynamic_cast<const scalarField*>(iter()));
                reg.erase(iter);
                autoPtr<vectorIOField> fieldVals
                (
                    new vectorIOField
                    (
                        IOobject
                        (
                            dataName,
                            "",
                            reg
                        ),
                        s.size()/3
                    )
                );

                label elemI = 0;
                forAll(fieldVals(), i)
                {
                    fieldVals()[i].x() = s[elemI++];
                    fieldVals()[i].y() = s[elemI++];
                    fieldVals()[i].z() = s[elemI++];
                }
                regIOobject::store(fieldVals);
            }
        }
        else
        {
            FatalIOErrorIn("vtkUnstructuredReader::read(..)", inFile)
                << "Unsupported tag "
                << tag << exit(FatalIOError);
        }
    }

    if (debug)
    {
        Info<< "Read points:" << points_.size()
            << " cellShapes:" << cells_.size()
            << " faces:" << faces_.size()
            << " lines:" << lines_.size()
            << nl << endl;

        Info<< "Cell fields:" << endl;
        printFieldStats<vectorIOField>(cellData_);
        printFieldStats<scalarIOField>(cellData_);
        printFieldStats<labelIOField>(cellData_);
        printFieldStats<stringIOList>(cellData_);
        Info<< nl << endl;

        Info<< "Point fields:" << endl;
        printFieldStats<vectorIOField>(pointData_);
        printFieldStats<scalarIOField>(pointData_);
        printFieldStats<labelIOField>(pointData_);
        printFieldStats<stringIOList>(pointData_);
        Info<< nl << endl;

        Info<< "Other fields:" << endl;
        printFieldStats<vectorIOField>(otherData_);
        printFieldStats<scalarIOField>(otherData_);
        printFieldStats<labelIOField>(otherData_);
        printFieldStats<stringIOList>(otherData_);
    }
}


template<class Type>
void Foam::vtkUnstructuredReader::printFieldStats
(
    const objectRegistry& obj
) const
{
    wordList fieldNames(obj.names(Type::typeName));

    if (fieldNames.size() > 0)
    {
        Info<< "Read " << fieldNames.size() << " " << Type::typeName
            << " fields:" << endl;
        Info<< "Size\tName" << nl
            << "----\t----" << endl;
        forAll(fieldNames, i)
        {
            Info<< obj.lookupObject<Type>(fieldNames[i]).size()
                << "\t" << fieldNames[i]
                << endl;
        }
        Info<< endl;
    }
}


// ************************************************************************* //

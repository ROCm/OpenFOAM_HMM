/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015 OpenCFD Ltd
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

#include "IOmanip.H"
#include "IFstream.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "ensightPartFaces.H"
#include "ensightPTraits.H"
#include "OStringStream.H"
#include "regExp.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::fileName Foam::ensightSurfaceWriter::writeUncollated
(
    const fileName& outputDir,
    const fileName& surfaceName,
    const pointField& points,
    const faceList& faces,
    const word& fieldName,
    const Field<Type>& values,
    const bool isNodeValues,
    const bool verbose
) const
{
    if (!isDir(outputDir/fieldName))
    {
        mkDir(outputDir/fieldName);
    }

    // const scalar timeValue = Foam::name(this->mesh().time().timeValue());
    const scalar timeValue = 0.0;

    OFstream osCase(outputDir/fieldName/surfaceName + ".case");
    ensightGeoFile osGeom
    (
        outputDir/fieldName/surfaceName + ".0000.mesh",
        writeFormat_
    );
    ensightFile osField
    (
        outputDir/fieldName/surfaceName + ".0000." + fieldName,
        writeFormat_
    );

    if (verbose)
    {
        Info<< "Writing case file to " << osCase.name() << endl;
    }

    osCase
        << "FORMAT" << nl
        << "type: ensight gold" << nl
        << nl
        << "GEOMETRY" << nl
        << "model:        1     " << osGeom.name().name() << nl
        << nl
        << "VARIABLE" << nl
        << ensightPTraits<Type>::typeName << " per "
        << word(isNodeValues ? "node:" : "element:") << setw(10) << 1
        << "       " << fieldName
        << "       " << surfaceName.c_str() << ".****." << fieldName << nl
        << nl
        << "TIME" << nl
        << "time set:                      1" << nl
        << "number of steps:               1" << nl
        << "filename start number:         0" << nl
        << "filename increment:            1" << nl
        << "time values:" << nl
        << timeValue << nl
        << nl;

    ensightPartFaces ensPart(0, osGeom.name().name(), points, faces, true);
    osGeom << ensPart;

    // Write field
    osField.writeKeyword(ensightPTraits<Type>::typeName);
    ensPart.writeField(osField, values, isNodeValues);

    return osCase.name();
}


template<class Type>
Foam::fileName Foam::ensightSurfaceWriter::writeCollated
(
    const fileName& outputDir,
    const fileName& surfaceName,
    const pointField& points,
    const faceList& faces,
    const word& fieldName,
    const Field<Type>& values,
    const bool isNodeValues,
    const bool verbose
) const
{

    const fileName baseDir = outputDir.path()/surfaceName;
    const fileName timeDir = outputDir.name();

    if (!isDir(baseDir))
    {
        mkDir(baseDir);
    }

    const fileName meshFile(baseDir/surfaceName + ".0000.mesh");
    const scalar timeValue = readScalar(IStringStream(timeDir)());
    label timeIndex = 0;


    // Do case file
    {
        dictionary dict;
        scalarList times;
        bool stateChanged = false;

        if (isFile(baseDir/"fieldsDict"))
        {
            IFstream is(baseDir/"fieldsDict");
            if (is.good() && dict.read(is))
            {
                dict.lookup("times") >> times;
                const scalar timeValue = readScalar(IStringStream(timeDir)());
                label index = findLower(times, timeValue);
                timeIndex = index+1;
            }
        }


        // Update stored times list
        times.setSize(timeIndex+1, -1);

        if (times[timeIndex] != timeValue)
        {
            stateChanged = true;
        }
        times[timeIndex] = timeValue;


        // Add my information to dictionary
        {
            dict.set("times", times);
            if (dict.found("fields"))
            {
                dictionary& fieldsDict = dict.subDict("fields");
                if (!fieldsDict.found(fieldName))
                {
                    dictionary fieldDict;
                    fieldDict.set("type", ensightPTraits<Type>::typeName);
                    fieldsDict.set(fieldName, fieldDict);

                    stateChanged = true;
                }
            }
            else
            {
                dictionary fieldDict;
                fieldDict.set("type", ensightPTraits<Type>::typeName);

                dictionary fieldsDict;
                fieldsDict.set(fieldName, fieldDict);

                dict.set("fields", fieldsDict);

                stateChanged = true;
            }
        }


        if (stateChanged)
        {
            if (verbose)
            {
                Info<< "Writing state file to fieldsDict" << endl;
            }
            OFstream os(baseDir/"fieldsDict");
            os << dict;


            OFstream osCase(baseDir/surfaceName + ".case");

            if (verbose)
            {
                Info<< "Writing case file to " << osCase.name() << endl;
            }

            osCase
                << "FORMAT" << nl
                << "type: ensight gold" << nl
                << nl
                << "GEOMETRY" << nl
                << "model:        1     " << meshFile.name() << nl
                << nl
                << "VARIABLE" << nl;
            const dictionary& fieldsDict = dict.subDict("fields");
            forAllConstIter(dictionary, fieldsDict, iter)
            {
                const word& fieldName = iter().keyword();
                const word fieldType(iter().dict().lookup("type"));

                osCase
                    << fieldType << " per "
                    << word(isNodeValues ? "node:" : "element:")
                    << setw(10) << 1
                    << setw(15) << fieldName
                    << "       " << surfaceName.c_str() << ".****." << fieldName
                    << nl;
            }
            osCase << nl;

            osCase
                << "TIME" << nl
                << "time set:                      1" << nl
                << "number of steps:               " << timeIndex+1 << nl
                << "filename start number:         0" << nl
                << "filename increment:            1" << nl
                << "time values:" << nl;
            forAll(times, timeI)
            {
                osCase << setw(12) << times[timeI] << " ";

                if (timeI != 0 && (timeI % 6) == 0)
                {
                    osCase << nl;
                }
            }
            osCase << nl;
        }
    }


    // Write geometry
    ensightPartFaces ensPart(0, meshFile.name(), points, faces, true);
    if (!exists(meshFile))
    {
        if (verbose)
        {
            Info<< "Writing mesh file to " << meshFile.name() << endl;
        }
        ensightGeoFile osGeom(meshFile, writeFormat_);
        osGeom << ensPart;
    }


    // Get string representation
    string timeString;
    {
        OStringStream os;
        os.stdStream().fill('0');
        os << setw(4) << timeIndex;
        timeString = os.str();
    }

    // Write field
    ensightFile osField
    (
        baseDir/surfaceName + "." + timeString + "." + fieldName,
        writeFormat_
    );
    if (verbose)
    {
        Info<< "Writing field file to " << osField.name() << endl;
    }
    osField.writeKeyword(ensightPTraits<Type>::typeName);
    ensPart.writeField(osField, values, isNodeValues);

    return baseDir/surfaceName + ".case";
}


template<class Type>
Foam::fileName Foam::ensightSurfaceWriter::writeTemplate
(
    const fileName& outputDir,
    const fileName& surfaceName,
    const pointField& points,
    const faceList& faces,
    const word& fieldName,
    const Field<Type>& values,
    const bool isNodeValues,
    const bool verbose
) const
{
    if (collateTimes_)
    {
        return writeCollated
        (
            outputDir,
            surfaceName,
            points,
            faces,
            fieldName,
            values,
            isNodeValues,
            verbose
        );
    }
    else
    {
        return writeUncollated
        (
            outputDir,
            surfaceName,
            points,
            faces,
            fieldName,
            values,
            isNodeValues,
            verbose
        );
    }
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015-2016 OpenCFD Ltd.
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
    const ensight::FileName surfName(surfaceName);
    const ensight::VarName  varName(fieldName);

    // use variable name as sub-directory for results
    // eg, something like this:
    // - VAR1/SURF1.case
    // - VAR1/SURF1.0000.mesh
    // - VAR1/SURF1.0001.VAR1
    // - VAR1/SURF2.case
    // - VAR1/SURF2.0000.mesh
    // - VAR1/SURF2.0001.VAR1
    // and
    // - VAR2/SURF1.case
    // - VAR2/SURF1.0000.mesh
    // - VAR2/SURF1.0001.VAR2

    const fileName baseDir = outputDir/varName;
    const fileName timeDir = outputDir.name();

    if (!isDir(baseDir))
    {
        mkDir(baseDir);
    }

    // const scalar timeValue = Foam::name(this->mesh().time().timeValue());
    const scalar timeValue = readScalar(IStringStream(timeDir)());

    OFstream osCase(baseDir/surfName + ".case");
    ensightGeoFile osGeom
    (
        baseDir,
        surfName + ".00000000.mesh",
        writeFormat_
    );
    ensightFile osField
    (
        baseDir,
        surfName + ".00000000." + varName,
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
        << "model:  1   " << osGeom.name().name() << nl
        << nl
        << "VARIABLE" << nl
        << ensightPTraits<Type>::typeName << " per "
        << word(isNodeValues ? "node:" : "element:")
        << setw(3) << 1
        << setw(15) << varName
        << "   " << surfName.c_str() << ".********." << varName << nl
        << nl
        << "TIME" << nl
        << "time set:               1" << nl
        << "number of steps:        1" << nl
        << "filename start number:  0" << nl
        << "filename increment:     1" << nl
        << "time values:" << nl
        << "    " << timeValue
        << nl << nl << "# end" << nl;

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
    const ensight::FileName surfName(surfaceName);
    const ensight::VarName  varName(fieldName);

    // use surface name as sub-directory for results
    // eg, something like this:
    // - SURF1/SURF1.case
    // - SURF1/SURF1.0000.mesh
    // - SURF1/SURF1/data/0000/VAR1
    // - SURF1/SURF1/data/0000/VAR2
    // and
    // - SURF2/SURF2.case
    // - SURF2/SURF2.0000.mesh
    // - SURF2/SURF2/data/0000/VAR1
    // - SURF2/SURF2/data/0000/VAR2

    const fileName baseDir = outputDir.path()/surfName;
    const fileName timeDir = outputDir.name();

    if (!isDir(baseDir))
    {
        mkDir(baseDir);
    }

    // surfName already validated
    const fileName meshFile(baseDir/surfName + ".000000.mesh");
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
                    fieldDict.set("name", varName); // ensight variable name

                    fieldsDict.set(fieldName, fieldDict);

                    stateChanged = true;
                }
            }
            else
            {
                dictionary fieldDict;
                fieldDict.set("type", ensightPTraits<Type>::typeName);
                fieldDict.set("name", varName); // ensight variable name

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
            {
                OFstream os(baseDir/"fieldsDict");
                os << "// summary of ensight times/fields" << nl << nl;
                dict.write(os, false);
            }

            OFstream osCase(baseDir/surfName + ".case");

            if (verbose)
            {
                Info<< "Writing case file to " << osCase.name() << endl;
            }

            osCase
                << "FORMAT" << nl
                << "type: ensight gold" << nl
                << nl
                << "GEOMETRY" << nl
                << "model:  1   " << meshFile.name() << nl
                << nl
                << "VARIABLE" << nl;

            const dictionary& fieldsDict = dict.subDict("fields");
            forAllConstIter(dictionary, fieldsDict, iter)
            {
                const dictionary& subDict = iter().dict();
                const word fieldType(subDict.lookup("type"));
                const word varName = subDict.lookupOrDefault
                (
                    "name",
                    iter().keyword() // fieldName as fallback
                );

                osCase
                    << fieldType << " per "
                    << word(isNodeValues ? "node:" : "element:")
                    << setw(3)  << 1
                    << setw(15) << varName
                    << "   data/******/" << varName
                    << nl;
            }
            osCase << nl;

            osCase
                << "TIME" << nl
                << "time set:               1" << nl
                << "number of steps:        " << timeIndex+1 << nl
                << "filename start number:  0" << nl
                << "filename increment:     1" << nl
                << "time values:" << nl;

            label count = 0;
            forAll(times, timeI)
            {
                osCase << ' ' << setw(12) << times[timeI];

                if (++count % 6 == 0)
                {
                    osCase << nl;
                }
            }
            osCase << nl << nl << "# end" << nl;
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
        // use two-argument form for path-name to avoid validating the base-dir
        ensightGeoFile osGeom(meshFile.path(), meshFile.name(), writeFormat_);
        osGeom << ensPart;
    }


    // Get string representation
    string timeString;
    {
        OStringStream os;
        os.stdStream().fill('0');
        os << setw(6) << timeIndex;
        timeString = os.str();
    }

    fileName dataDir = baseDir/"data"/timeString;

    // as per mkdir -p "data/000000"
    mkDir(dataDir);

    // Write field
    ensightFile osField
    (
        dataDir,
        varName,
        writeFormat_
    );
    if (verbose)
    {
        Info<< "Writing field file to " << osField.name() << endl;
    }
    osField.writeKeyword(ensightPTraits<Type>::typeName);
    ensPart.writeField(osField, values, isNodeValues);

    // place a timestamp in the directory for future reference
    {
        OFstream timeStamp(dataDir/"time");
        timeStamp
            << "#   timestep time" << nl
            << dataDir.name() << " " << timeValue << nl;
    }

    return baseDir/surfName + ".case";
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

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2010-2010 OpenCFD Ltd.
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

#include "ensightSurfaceWriter.H"

#include "OFstream.H"
#include "OSspecific.H"
#include "IOmanip.H"
#include "ensightGeoFile.H"
#include "ensightPartNonMeshFaces.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

namespace Foam
{

    // Write scalarField in ensight format
    template<>
    void Foam::ensightSurfaceWriter<Foam::scalar>::writeData
    (
        Ostream& os,
        const Field<Foam::scalar>& values
    )
    {
        forAll(values, i)
        {
            os << values[i] << nl;
        }
    }


    // Write booField in ensight format
    template<>
    void Foam::ensightSurfaceWriter<bool>::writeData
    (
        Ostream& os,
        const Field<bool>& values
    )
    {
        forAll(values, i)
        {
            os << values[i] << nl;
        }
    }
}


// Write generic field in ensight format
template<class Type>
void Foam::ensightSurfaceWriter<Type>::writeData
(
    Ostream& os,
    const Field<Type>& values
)
{
    for (direction cmpt = 0; cmpt < vector::nComponents; cmpt++)
    {
        scalarField v = values.component(cmpt);
        forAll(v, i)
        {
            os << v[i] << nl;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
template<class Type>
Foam::ensightSurfaceWriter<Type>::ensightSurfaceWriter()
:
    surfaceWriter<Type>()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::ensightSurfaceWriter<Type>::~ensightSurfaceWriter()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::ensightSurfaceWriter<Type>::write
(
    const fileName& outputDir,
    const fileName& surfaceName,
    const pointField& points,
    const faceList& faces,
    const bool verbose
) const
{
    if (!isDir(outputDir))
    {
        mkDir(outputDir);
    }

    //const scalar timeValue = Foam::name(this->mesh().time().timeValue());
    const scalar timeValue = 0.0;


    OFstream caseStr(outputDir/surfaceName + ".case");
    ensightGeoFile geomStr
    (
        outputDir/surfaceName + ".000.mesh",
        IOstream::ASCII
    );

    if (verbose)
    {
        Info<< "Writing case file to " << caseStr.name() << endl;
    }

    caseStr
        << "FORMAT" << nl
        << "type: ensight gold" << nl
        << nl
        << "GEOMETRY" << nl
        << "model:        1     " << geomStr.name().name() << nl
        << nl
        << "TIME" << nl
        << "time set:                      1" << nl
        << "number of steps:               1" << nl
        << "filename start number:         0" << nl
        << "filename increment:            1" << nl
        << "time values:" << nl
        << timeValue << nl
        << nl;

    ensightPartNonMeshFaces faceWriter(0, geomStr.name().name(), faces, points);
    faceWriter.writeGeometry(geomStr);
}


template<class Type>
void Foam::ensightSurfaceWriter<Type>::write
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

    //const scalar timeValue = Foam::name(this->mesh().time().timeValue());
    const scalar timeValue = 0.0;


    OFstream caseStr(outputDir/fieldName/surfaceName + ".case");
    ensightGeoFile geomStr
    (
        outputDir/fieldName/surfaceName + ".000.mesh",
        IOstream::ASCII
    );
    ensightFile fieldStr
    (
        outputDir/fieldName/surfaceName + ".000." + fieldName,
        IOstream::ASCII
    );

    if (verbose)
    {
        Info<< "Writing case file to " << caseStr.name() << endl;
    }

    caseStr
        << "FORMAT" << nl
        << "type: ensight gold" << nl
        << nl
        << "GEOMETRY" << nl
        << "model:        1     " << geomStr.name().name() << nl
        << nl
        << "VARIABLE" << nl;
    if (isNodeValues)
    {
        caseStr
            << pTraits<Type>::typeName << " per node:" << setw(10) << 1
            << "       " << fieldName
            << "       " << surfaceName.c_str() << ".***." << fieldName << nl;
    }
    else
    {
        caseStr
            << pTraits<Type>::typeName << " per element:" << setw(10) << 1
            << "       " << fieldName
            << "       " << surfaceName.c_str() << ".***." << fieldName << nl;
    }

    caseStr
        << nl
        << "TIME" << nl
        << "time set:                      1" << nl
        << "number of steps:               1" << nl
        << "filename start number:         0" << nl
        << "filename increment:            1" << nl
        << "time values:" << nl
        << timeValue << nl
        << nl;

    ensightPartNonMeshFaces faceWriter(0, geomStr.name().name(), faces, points);
    faceWriter.writeGeometry(geomStr);

    // Write field
    fieldStr
        << pTraits<Type>::typeName << nl
        << "part" << nl
        << setw(10) << 1 << nl;

    if (isNodeValues)
    {
        fieldStr << "coordinates" << nl;
        writeData(fieldStr, values);
    }
    else
    {
        //faceWriter.writeField(fieldStr, values);
        forAll(faceWriter.elementTypes(), elemI)
        {
            if (faceWriter.elemLists()[elemI].size())
            {
                fieldStr.writeKeyword(faceWriter.elementTypes()[elemI]);
                writeData
                (
                    fieldStr,
                    Field<Type>
                    (
                        values,
                        faceWriter.elemLists()[elemI]
                    )
                );
            }
        }
    }
}


// ************************************************************************* //

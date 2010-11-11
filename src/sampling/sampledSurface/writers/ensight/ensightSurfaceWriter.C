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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::ensightSurfaceWriter<Type>::binShapes
(
    const pointField& points,
    const faceList& faces,
    DynamicList<label>& tris,
    DynamicList<label>& quads,
    DynamicList<label>& polys
)
{
    tris.setCapacity(faces.size());
    quads.setCapacity(faces.size());
    polys.setCapacity(faces.size());

    forAll(faces, i)
    {
        const face& f = faces[i];
        if (f.size() == 3)
        {
            tris.append(i);
        }
        else if (f.size() == 4)
        {
            quads.append(i);
        }
        else
        {
            polys.append(i);
        }
    }
}


template<class Type>
void Foam::ensightSurfaceWriter<Type>::writeGeometry
(
    const fileName& surfaceName,
    Ostream& os,
    const pointField& points,
    const faceList& faces,
    const DynamicList<label>& tris,
    const DynamicList<label>& quads,
    const DynamicList<label>& polys
)
{
    os  << "EnSight Geometry File" << nl
        << (string("written by OpenFOAM-") + Foam::FOAMversion).c_str() << nl
        << "node id assign" << nl
        << "element id assign" << nl
        << "part" << nl
        << setw(10) << 1 << nl
        << surfaceName.c_str() << nl;

    os  << "coordinates" << nl
        << setw(10) << points.size() << nl;
    for (direction cmpt = 0; cmpt < vector::nComponents; cmpt++)
    {
        scalarField v = points.component(cmpt);
        forAll(v, i)
        {
            os  << v[i] << nl;
        }
    }


    if (tris.size())
    {
        os  << "tria3" << nl
            << setw(10) << tris.size() << nl;
        forAll(tris, i)
        {
            const face& f = faces[tris[i]];
            forAll(f, fp)
            {
                os  << setw(10) << f[fp]+1;
            }
            os  << nl;
        }
    }
    if (quads.size())
    {
        os  << "quad4" << nl
            << setw(10) << quads.size() << nl;
        forAll(quads, i)
        {
            const face& f = faces[quads[i]];
            forAll(f, fp)
            {
                os  << setw(10) << f[fp]+1;
            }
            os  << nl;
        }
    }
    if (polys.size())
    {
        os  << "nsided" << nl
            << setw(10) << polys.size() << nl;
        forAll(polys, i)
        {
            const face& f = faces[polys[i]];
            forAll(f, fp)
            {
                os  << setw(10) << f[fp]+1;
            }
            os  << nl;
        }
    }
}


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
    OFstream geomStr(outputDir/surfaceName + ".***.mesh");

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

    DynamicList<label> tris;
    DynamicList<label> quads;
    DynamicList<label> polys;
    binShapes(points, faces, tris, quads, polys);

    writeGeometry(surfaceName, geomStr, points, faces, tris, quads, polys);
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
    OFstream geomStr(outputDir/fieldName/surfaceName + ".000.mesh");
    OFstream fieldStr(outputDir/fieldName/surfaceName + ".000." + fieldName);

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

    DynamicList<label> tris;
    DynamicList<label> quads;
    DynamicList<label> polys;
    binShapes(points, faces, tris, quads, polys);

    writeGeometry(surfaceName, geomStr, points, faces, tris, quads, polys);

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
        if (tris.size())
        {
            fieldStr << "tria3" << nl;
            writeData(fieldStr, Field<Type>(values, tris));
        }
        if (quads.size())
        {
            fieldStr << "quad4" << nl;
            writeData(fieldStr, Field<Type>(values, quads));
        }
        if (polys.size())
        {
            fieldStr << "nsided" << nl;
            writeData(fieldStr, Field<Type>(values, polys));
        }
    }
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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

#include "ensightSetWriter.H"
#include "coordSet.H"
#include "IOmanip.H"
#include "ensightGeoFile.H"
#include "ensightPTraits.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::ensightSetWriter<Type>::ensightSetWriter()
:
    writer<Type>()
{}


template<class Type>
Foam::ensightSetWriter<Type>::ensightSetWriter(const dictionary& dict)
:
    writer<Type>(dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::fileName Foam::ensightSetWriter<Type>::getFileName
(
    const coordSet& points,
    const wordList& valueSetNames
) const
{
    return
        this->getBaseName(points, valueSetNames)
      //+ '_'
      //+ ensightPTraits<Type>::typeName
      + ".case";
}


template<class Type>
void Foam::ensightSetWriter<Type>::write
(
    const coordSet& points,
    const wordList& valueSetNames,
    const List<const Field<Type>*>& valueSets,
    Ostream& os
) const
{
    const fileName base(os.name().lessExt());
    const fileName meshFile(base + ".mesh");

    // Write .case file
    os  << "FORMAT" << nl
        << "type: ensight gold" << nl
        << nl
        << "GEOMETRY" << nl
        << "model:        1     " << meshFile.name().c_str() << nl
        << nl
        << "VARIABLE"
        << nl;

    for (const word& valueName : valueSetNames)
    {
        fileName dataFile(base + ".***." + valueName);

        os.setf(ios_base::left);
        os  << ensightPTraits<Type>::typeName
            << " per node:            1       "
            << setw(15) << valueName
            << " " << dataFile.name().c_str()
            << nl;
    }
    os  << nl
        << "TIME" << nl
        << "time set:                      1" << nl
        << "number of steps:               1" << nl
        << "filename start number:         0" << nl
        << "filename increment:            1" << nl
        << "time values:" << nl
        << "0.00000e+00" << nl;

    // Write .mesh file
    {
        string desc("Written by OpenFOAM");
        OFstream os(meshFile);
        os.setf(ios_base::scientific, ios_base::floatfield);
        os.precision(5);

        os  << "Ensight Geometry File" << nl
            << desc.c_str() << nl
            << "node id assign" << nl
            << "element id assign" << nl
            << "part" << nl
            << setw(10) << 1 << nl
            << "internalMesh" << nl
            << "coordinates" << nl
            << setw(10) << points.size() << nl;

        for (direction cmpt = 0; cmpt < vector::nComponents; ++cmpt)
        {
            for (const point& p : points)
            {
                const float comp = narrowFloat(p[cmpt]);
                os  << setw(12) << comp << nl;
            }
        }
        os  << "point" << nl
            << setw(10) << points.size() << nl;
        forAll(points, pointi)
        {
            os  << setw(10) << pointi+1 << nl;
        }
    }

    // Write data files
    forAll(valueSetNames, seti)
    {
        const word& valueName = valueSetNames[seti];
        const Field<Type>& fld = *(valueSets[seti]);

        fileName dataFile(base + ".000." + valueName);
        OFstream os(dataFile);
        os.setf(ios_base::scientific, ios_base::floatfield);
        os.precision(5);

        os  << ensightPTraits<Type>::typeName << nl
            << "part" << nl
            << setw(10) << 1 << nl
            << "coordinates" << nl;

        for (direction d=0; d < pTraits<Type>::nComponents; ++d)
        {
            const direction cmpt = ensightPTraits<Type>::componentOrder[d];

            for (const Type& val : fld)
            {
                const float comp = narrowFloat(component(val, cmpt));
                os  << setw(12) << comp << nl;
            }
        }
    }
}


template<class Type>
void Foam::ensightSetWriter<Type>::write
(
    const bool writeTracks,
    const List<scalarField>& times,
    const PtrList<coordSet>& tracks,
    const wordList& valueSetNames,
    const List<List<Field<Type>>>& valueSets,
    Ostream& os
) const
{
    const fileName base(os.name().lessExt());
    const fileName meshFile(base + ".mesh");

    // Write .case file
    os  << "FORMAT" << nl
        << "type: ensight gold" << nl
        << nl
        << "GEOMETRY" << nl
        << "model:        1     " << meshFile.name().c_str() << nl
        << nl
        << "VARIABLE"
        << nl;

    for (const word& valueName : valueSetNames)
    {
        fileName dataFile(base + ".***." + valueName);

        os.setf(ios_base::left);
        os  << ensightPTraits<Type>::typeName
            << " per node:            1       "
            << setw(15) << valueName
            << " " << dataFile.name().c_str()
            << nl;
    }
    os  << nl
        << "TIME" << nl
        << "time set:                      1" << nl
        << "number of steps:               1" << nl
        << "filename start number:         0" << nl
        << "filename increment:            1" << nl
        << "time values:" << nl
        << "0.00000e+00" << nl;

    // Write .mesh file
    {
        string desc("Written by OpenFOAM");
        OFstream os(meshFile);
        os.setf(ios_base::scientific, ios_base::floatfield);
        os.precision(5);
        os  << "Ensight Geometry File" << nl
            << desc.c_str() << nl
            << "node id assign" << nl
            << "element id assign" << nl;

        forAll(tracks, tracki)
        {
            const coordSet& points = tracks[tracki];

            os  << "part" << nl
                << setw(10) << tracki+1 << nl
                << "internalMesh" << nl
                << "coordinates" << nl
                << setw(10) << points.size() << nl;

            for (direction cmpt = 0; cmpt < vector::nComponents; ++cmpt)
            {
                for (const point& p : points)
                {
                    const float comp = narrowFloat(p[cmpt]);
                    os  << setw(12) << comp << nl;
                }
            }

            if (writeTracks)
            {
                os  << "bar2" << nl
                    << setw(10) << points.size()-1 << nl;
                for (label i = 0; i < points.size()-1; i++)
                {
                    os  << setw(10) << i+1
                        << setw(10) << i+2
                        << nl;
                }
            }
        }
    }


    // Write data files
    forAll(valueSetNames, seti)
    {
        const word& valueName = valueSetNames[seti];
        const List<Field<Type>>& fieldVals = valueSets[seti];

        fileName dataFile(base + ".000." + valueName);
        OFstream os(dataFile);
        os.setf(ios_base::scientific, ios_base::floatfield);
        os.precision(5);
        {
            os  << ensightPTraits<Type>::typeName << nl;

            forAll(fieldVals, tracki)
            {
                const Field<Type>& fld = fieldVals[tracki];

                os  << "part" << nl
                    << setw(10) << tracki+1 << nl
                    << "coordinates" << nl;

                for (direction d=0; d < pTraits<Type>::nComponents; ++d)
                {
                    const direction cmpt =
                        ensightPTraits<Type>::componentOrder[d];

                    for (const Type& val : fld)
                    {
                        const float comp = narrowFloat(component(val, cmpt));
                        os  << setw(12) << comp << nl;
                    }
                }
            }
        }
    }
}


// ************************************************************************* //

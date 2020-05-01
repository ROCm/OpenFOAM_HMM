/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2020 OpenCFD Ltd.
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

#include "volFields.H"
#include "meshStructure.H"
#include "globalIndex.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
bool Foam::functionObjects::columnAverage::columnAverageField
(
    const word& fieldName
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    const fieldType* fldPtr = cfindObject<fieldType>(fieldName);

    if (fldPtr)
    {
        const fieldType& fld = *fldPtr;

        const word resultName(averageName(fieldName));

        fieldType* resPtr = obr_.getObjectPtr<fieldType>(resultName);

        if (!resPtr)
        {
            resPtr = new fieldType
            (
                IOobject
                (
                    resultName,
                    fld.mesh().time().timeName(),
                    fld.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                fld
            );
            obr_.objectRegistry::store(resPtr);
        }
        fieldType& res = *resPtr;

        const meshStructure& ms = meshAddressing(fld.mesh());
        if (globalFaces_().empty())
        {
            return false;
        }

        const labelList& cellToPatchFace = ms.cellToPatchFaceAddressing();

        // Brute force: collect per-global-patchface on all processors
        Field<Type> regionField(globalFaces_().size(), Zero);
        labelList regionCount(globalFaces_().size(), 0);

        forAll(cellToPatchFace, celli)
        {
            const label regioni = cellToPatchFace[celli];
            regionField[regioni] += fld[celli];
            regionCount[regioni]++;
        }

        // Global sum
        Pstream::listCombineGather(regionField, plusEqOp<Type>());
        Pstream::listCombineScatter(regionField);
        Pstream::listCombineGather(regionCount, plusEqOp<label>());
        Pstream::listCombineScatter(regionCount);

        forAll(regionField, regioni)
        {
            regionField[regioni] /= regionCount[regioni];
        }

        // And send result back
        forAll(cellToPatchFace, celli)
        {
            const label regioni = cellToPatchFace[celli];
            res[celli] = regionField[regioni];
        }
        res.correctBoundaryConditions();

        return true;
    }

    return false;
}


// ************************************************************************* //

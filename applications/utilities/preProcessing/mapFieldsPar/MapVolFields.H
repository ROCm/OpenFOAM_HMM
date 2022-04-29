/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2018-2021 OpenCFD Ltd.
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

#ifndef MapConsistentVolFields_H
#define MapConsistentVolFields_H

#include "GeometricField.H"
#include "meshToMesh.H"
#include "IOobjectList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
void evaluateConstraintTypes(GeometricField<Type, fvPatchField, volMesh>& fld)
{
    auto& fldBf = fld.boundaryFieldRef();

    const UPstream::commsTypes commsType(UPstream::defaultCommsType);

    if
    (
        commsType == UPstream::commsTypes::blocking
     || commsType == UPstream::commsTypes::nonBlocking
    )
    {
        const label startOfRequests = UPstream::nRequests();

        forAll(fldBf, patchi)
        {
            fvPatchField<Type>& tgtField = fldBf[patchi];

            if
            (
                tgtField.type() == tgtField.patch().patch().type()
             && polyPatch::constraintType(tgtField.patch().patch().type())
            )
            {
                tgtField.initEvaluate(commsType);
            }
        }

        // Wait for outstanding requests
        if
        (
            UPstream::parRun()
         && commsType == UPstream::commsTypes::nonBlocking
        )
        {
            UPstream::waitRequests(startOfRequests);
        }

        forAll(fldBf, patchi)
        {
            fvPatchField<Type>& tgtField = fldBf[patchi];

            if
            (
                tgtField.type() == tgtField.patch().patch().type()
             && polyPatch::constraintType(tgtField.patch().patch().type())
            )
            {
                tgtField.evaluate(commsType);
            }
        }
    }
    else if (commsType == UPstream::commsTypes::scheduled)
    {
        const lduSchedule& patchSchedule =
            fld.mesh().globalData().patchSchedule();

        for (const auto& schedEval : patchSchedule)
        {
            const label patchi = schedEval.patch;

            fvPatchField<Type>& tgtField = fldBf[patchi];

            if
            (
                tgtField.type() == tgtField.patch().patch().type()
             && polyPatch::constraintType(tgtField.patch().patch().type())
            )
            {
                if (schedEval.init)
                {
                    tgtField.initEvaluate(commsType);
                }
                else
                {
                    tgtField.evaluate(commsType);
                }
            }
        }
    }
}


template<class Type, class CombineOp>
void MapVolFields
(
    const IOobjectList& objects,
    const wordRes& selectedFields,
    const meshToMesh& interp,
    const CombineOp& cop
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    const fvMesh& meshSource = static_cast<const fvMesh&>(interp.srcRegion());
    const fvMesh& meshTarget = static_cast<const fvMesh&>(interp.tgtRegion());

    // Available fields, sorted order
    const wordList fieldNames =
    (
        selectedFields.empty()
      ? objects.sortedNames<fieldType>()
      : objects.sortedNames<fieldType>(selectedFields)
    );

    for (const word& fieldName : fieldNames)
    {
        const fieldType fieldSource(*(objects[fieldName]), meshSource, false);

        IOobject targetIO
        (
            fieldName,
            meshTarget.time().timeName(),
            meshTarget,
            IOobject::MUST_READ
        );

        if (targetIO.typeHeaderOk<fieldType>(true))
        {
            Info<< "    interpolating onto existing field "
                << fieldName << endl;

            fieldType fieldTarget(targetIO, meshTarget, false);

            interp.mapSrcToTgt(fieldSource, cop, fieldTarget);

            evaluateConstraintTypes(fieldTarget);

            fieldTarget.write();
        }
        else
        {
            Info<< "    creating new field "
                << fieldName << endl;

            targetIO.readOpt(IOobject::NO_READ);

            tmp<fieldType> tfieldTarget
            (
                interp.mapSrcToTgt(fieldSource, cop)
            );

            fieldType fieldTarget(targetIO, tfieldTarget);

            evaluateConstraintTypes(fieldTarget);

            fieldTarget.write();
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

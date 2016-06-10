/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2016 OpenFOAM Foundation
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

#include "fvcCellReduce.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "extrapolatedCalculatedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class CombineOp>
tmp<GeometricField<Type, fvPatchField, volMesh>> cellReduce
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& ssf,
    const CombineOp& cop,
    const Type& nullValue
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> volFieldType;

    const fvMesh& mesh = ssf.mesh();

    tmp<volFieldType> tresult
    (
        new volFieldType
        (
            IOobject
            (
                "cellReduce(" + ssf.name() + ')',
                ssf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<Type>("initialValue", ssf.dimensions(), nullValue),
            extrapolatedCalculatedFvPatchField<Type>::typeName
        )
    );

    volFieldType& result = tresult.ref();

    const labelUList& own = mesh.owner();
    const labelUList& nbr = mesh.neighbour();

    // Internal field
    const Field<Type>& iField = ssf.internalField();
    forAll(iField, faceI)
    {
        label cellOwn = own[faceI];
        cop(result[cellOwn], iField[faceI]);

        label cellNbr = nbr[faceI];
        cop(result[cellNbr], iField[faceI]);
    }

    // Boundary field
    forAll(ssf.boundaryField(), patchI)
    {
        const fvsPatchField<Type>& pf = ssf.boundaryField()[patchI];
        const label start = pf.patch().start();

        forAll(pf, i)
        {
            label faceI = start + i;
            label cellI = own[faceI];
            cop(result[cellI], pf[i]);
        }
    }

    result.correctBoundaryConditions();

    return tresult;
}


template<class Type, class CombineOp>
tmp<GeometricField<Type, fvPatchField, volMesh>> cellReduce
(
    const tmp<GeometricField<Type, fvsPatchField, surfaceMesh>&> tssf,
    const CombineOp& cop,
    const Type& nullValue
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh>>
        tvf(cellReduce(cop, tssf, nullValue));

    tssf.clear();

    return tvf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

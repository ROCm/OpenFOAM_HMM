/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015 OpenCFD Ltd.
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

#include "gaussConvectionScheme.H"
#include "blendedSchemeBase.H"
#include "fvcCellReduce.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::blendingFactor::calc()
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    if (!obr_.foundObject<fieldType>(fieldName_))
    {
        return;
    }

    const fvMesh& mesh = refCast<const fvMesh>(obr_);

    const fieldType& field = mesh.lookupObject<fieldType>(fieldName_);

    const word divScheme("div(" + phiName_ + ',' + fieldName_ + ')');
    ITstream& its = mesh.divScheme(divScheme);

    const surfaceScalarField& phi =
        mesh.lookupObject<surfaceScalarField>(phiName_);

    tmp<fv::convectionScheme<Type> > cs =
        fv::convectionScheme<Type>::New(mesh, phi, its);

    const fv::gaussConvectionScheme<Type>& gcs =
        refCast<const fv::gaussConvectionScheme<Type> >(cs());

    const surfaceInterpolationScheme<Type>& interpScheme =
        gcs.interpScheme();

    if (!isA<blendedSchemeBase<Type> >(interpScheme))
    {
        FatalErrorInFunction
            << interpScheme.typeName << " is not a blended scheme"
            << exit(FatalError);
    }

    // Retrieve the face-based blending factor
    const blendedSchemeBase<Type>& blendedScheme =
        refCast<const blendedSchemeBase<Type> >(interpScheme);
    const surfaceScalarField factorf(blendedScheme.blendingFactor(field));

    // Convert into vol field whose values represent the local face minima
    // Note: factor applied to 1st scheme, and (1-factor) to 2nd scheme
    volScalarField& indicator =
        const_cast<volScalarField&>
        (
            obr_.lookupObject<volScalarField>(resultName_)
        );
    indicator = 1 - fvc::cellReduce(factorf, minEqOp<scalar>(), GREAT);
    indicator.correctBoundaryConditions();

    // Generate scheme statistics
    label nCellsScheme1 = 0;
    label nCellsScheme2 = 0;
    label nCellsBlended = 0;
    forAll(indicator, cellI)
    {
        scalar i = indicator[cellI];

        if (i < tolerance_)
        {
            nCellsScheme1++;
        }
        else if (i > (1 - tolerance_))
        {
            nCellsScheme2++;
        }
        else
        {
            nCellsBlended++;
        }
    }

    reduce(nCellsScheme1, sumOp<label>());
    reduce(nCellsScheme2, sumOp<label>());
    reduce(nCellsBlended, sumOp<label>());

    if (log_) Info
        << type() << " " << name_ << " output:" << nl
        << "    scheme 1 cells :  " << nCellsScheme1 << nl
        << "    scheme 2 cells :  " << nCellsScheme2 << nl
        << "    blended cells  :  " << nCellsBlended << nl
        << endl;

    file()
        << obr_.time().time().value()
        << token::TAB << nCellsScheme1
        << token::TAB << nCellsScheme2
        << token::TAB << nCellsBlended
        << endl;
}


// ************************************************************************* //

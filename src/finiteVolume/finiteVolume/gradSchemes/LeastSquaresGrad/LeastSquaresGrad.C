/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

#include "LeastSquaresGrad.H"
#include "LeastSquaresVectors.H"
#include "gaussGrad.H"
#include "fvMesh.H"
#include "volMesh.H"
#include "zeroGradientFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class Stencil>
Foam::tmp
<
    Foam::GeometricField
    <
        typename Foam::outerProduct<Foam::vector, Type>::type,
        Foam::fvPatchField,
    Foam::volMesh
    >
>
Foam::fv::LeastSquaresGrad<Type, Stencil>::calcGrad
(
    const GeometricField<Type, fvPatchField, volMesh>& vtf,
    const word& name
) const
{
    typedef typename outerProduct<vector, Type>::type GradType;

    const fvMesh& mesh = vtf.mesh();

    // Get reference to least square vectors
    const LeastSquaresVectors<Stencil>& lsv = LeastSquaresVectors<Stencil>::New
    (
        mesh
    );

    tmp<GeometricField<GradType, fvPatchField, volMesh> > tlsGrad
    (
        new GeometricField<GradType, fvPatchField, volMesh>
        (
            IOobject
            (
                name,
                vtf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<GradType>
            (
                "zero",
                vtf.dimensions()/dimLength,
                pTraits<GradType>::zero
            ),
            zeroGradientFvPatchField<GradType>::typeName
        )
    );
    GeometricField<GradType, fvPatchField, volMesh>& lsGrad = tlsGrad();
    Field<GradType>& lsGradIf = lsGrad;

    const extendedCentredCellToCellStencil& stencil = lsv.stencil();
    const List<List<vector> >& lsvs = lsv.vectors();

    // lsGrad = stencil.weightedSum(vtf, lsv.vectors());

    List<List<Type> > stencilVtf;
    extendedCellToFaceStencil::collectData
    (
        stencil.map(),
        stencil.stencil(),
        vtf,
        stencilVtf
    );

    forAll(lsGradIf, celli)
    {
        const List<Type>& stfc = stencilVtf[celli];
        const List<vector>& lsvc = lsvs[celli];

        forAll(stfc, i)
        {
            lsGradIf[celli] += lsvc[i]*stfc[i];
        }
    }

    lsGrad.correctBoundaryConditions();

    gaussGrad<Type>::correctBoundaryConditions(vtf, lsGrad);

    return tlsGrad;
}


// ************************************************************************* //

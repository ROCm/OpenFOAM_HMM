/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2007-2019 PCOpt/NTUA
                            | Copyright (C) 2013-2019 FOSS GP
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

#include "emptyFvPatch.H"
#include "adjointBoundaryCondition.H"
#include "adjointSolverManager.H"
#include "HashTable.H"
#include "surfaceInterpolationScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
tmp
<
    Field<typename Foam::outerProduct<Foam::vector, Type>::type>
>
adjointBoundaryCondition::computePatchGrad(word name)
{
    // Return field
    typedef typename outerProduct<vector, Type>::type GradType;
    auto tresGrad = tmp<Field<GradType>>::New(patch_.size(), Zero);
    auto& resGrad = tresGrad.ref();

    const labelList& faceCells = patch_.faceCells();
    const fvMesh& mesh = patch_.boundaryMesh().mesh();
    const cellList& cells = mesh.cells();

    // Go through the surfaceInterpolation scheme defined in gradSchemes for
    // consistency
    const GeometricField<Type, fvPatchField, volMesh>& field =
        mesh.lookupObject<volVectorField>(name);

    // Gives problems when grad(AdjointVar) is computed using a limited scheme,
    // since it is not possible to know a priori how many words to expect in the
    // stream.
    // Interpolation scheme is now read through interpolation schemes.
    /*
    word  gradSchemeName       ("grad(" + name + ')');
    Istream& is = mesh.gradScheme(gradSchemeName);
    word schemeData(is);
    */

    tmp<surfaceInterpolationScheme<Type>> tinterpScheme
    (
        surfaceInterpolationScheme<Type>::New
        (
            mesh,
            mesh.interpolationScheme("interpolate(" + name + ")")
        )
    );

    GeometricField<Type, fvsPatchField, surfaceMesh> surfField
    (
        tinterpScheme().interpolate(field)
    );

    // Auxiliary fields
    const surfaceVectorField& Sf = mesh.Sf();
    tmp<vectorField> tnf = patch_.nf();
    const vectorField& nf = tnf();
    const scalarField& V = mesh.V();
    const labelUList& owner = mesh.owner();

    // Compute grad value of cell adjacent to the boundary
    forAll(faceCells, fI)
    {
        const label cI = faceCells[fI];
        const cell& cellI = cells[cI];
        for (const label faceI : cellI) // global face numbering
        {
            label patchID = mesh.boundaryMesh().whichPatch(faceI);
            if (patchID == -1) //face is internal
            {
                const label own = owner[faceI];
                tensor flux = Sf[faceI]*surfField[faceI];
                if (cI == own)
                {
                    resGrad[fI] += flux;
                }
                else
                {
                    resGrad[fI] -= flux;
                }
            }
            else  // Face is boundary. Covers coupled patches as well
            {
                if (!isA<emptyFvPatch>(mesh.boundary()[patchID]))
                {
                    const fvPatch& patchForFlux = mesh.boundary()[patchID];
                    const label boundaryFaceI = faceI - patchForFlux.start();
                    const vectorField& Sfb = Sf.boundaryField()[patchID];
                    resGrad[fI] +=
                        Sfb[boundaryFaceI]
                       *surfField.boundaryField()[patchID][boundaryFaceI];
                }
            }
        }
        resGrad[fI] /= V[cI];
    }

    // This has concluded the computation of the grad at the cell next to the
    // boundary. We now need to compute the grad at the boundary face
    const fvPatchField<Type>& bField = field.boundaryField()[patch_.index()];
    resGrad = nf*bField.snGrad() + (resGrad - nf*(nf & resGrad));

    return tresGrad;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
adjointBoundaryCondition::adjointBoundaryCondition
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const word& solverName
)
:
    patch_(p),
    managerName_("objectiveManager" + solverName),
    adjointSolverName_(solverName),
    simulationType_("incompressible"),
    boundaryContrPtr_(nullptr),
    addATCUaGradUTerm_(nullptr)
{
    // Set the boundaryContribution pointer
    setBoundaryContributionPtr();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

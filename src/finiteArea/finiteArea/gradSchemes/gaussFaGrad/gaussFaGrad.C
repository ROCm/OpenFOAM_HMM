/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2016-2017 Wikki Ltd
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

#include "gaussFaGrad.H"
#include "facGrad.H"
#include "areaFields.H"
#include "edgeFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fa
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp
<
    GeometricField
    <
        typename outerProduct<vector, Type>::type, faPatchField, areaMesh
    >
>
gaussGrad<Type>::grad
(
    const GeometricField<Type, faPatchField, areaMesh>& vsf
) const
{
    typedef typename outerProduct<vector, Type>::type GradType;

    tmp<GeometricField<GradType, faPatchField, areaMesh>> tgGrad
    (
        fac::edgeIntegrate
        (
            vsf.mesh().Le()
           *tinterpScheme_().interpolate(vsf)
        )
    );

    GeometricField<GradType, faPatchField, areaMesh>& gGrad = tgGrad.ref();

    // Removed for consistencty.  Matthias Rauter, 6/Dec/2016
    gGrad.correctBoundaryConditions();

    gGrad.rename("grad(" + vsf.name() + ')');
    correctBoundaryConditions(vsf, gGrad);

    return tgGrad;
}


template<class Type>
void gaussGrad<Type>::correctBoundaryConditions
(
    const GeometricField<Type, faPatchField, areaMesh>& vsf,
    GeometricField
    <
        typename outerProduct<vector, Type>::type, faPatchField, areaMesh
    >& gGrad
)
{
    forAll(vsf.boundaryField(), patchI)
    {
        if (!vsf.boundaryField()[patchI].coupled())
        {
            const vectorField m
            (
                vsf.mesh().Le().boundaryField()[patchI]/
                vsf.mesh().magLe().boundaryField()[patchI]
            );

            gGrad.boundaryFieldRef()[patchI] += m*
            (
                vsf.boundaryField()[patchI].snGrad()
              - (m & gGrad.boundaryField()[patchI])
            );
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fa

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

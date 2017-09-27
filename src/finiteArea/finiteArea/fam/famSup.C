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

#include "areaFields.H"
#include "edgeFields.H"
#include "faMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<faMatrix<Type>>
Su
(
    const GeometricField<Type, faPatchField, areaMesh>& su,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    const faMesh& mesh = vf.mesh();

    tmp<faMatrix<Type>> tfam
    (
        new faMatrix<Type>
        (
            vf,
            dimArea*su.dimensions()
        )
    );
    faMatrix<Type>& fam = tfam.ref();

    fam.source() -= mesh.S()*su.internalField();

    return tfam;
}

template<class Type>
tmp<faMatrix<Type>>
Su
(
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tsu,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    tmp<faMatrix<Type>> tfam = fam::Su(tsu(), vf);
    tsu.clear();
    return tfam;
}


template<class Type>
tmp<faMatrix<Type>>
Sp
(
    const areaScalarField& sp,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    const faMesh& mesh = vf.mesh();

    tmp<faMatrix<Type>> tfam
    (
        new faMatrix<Type>
        (
            vf,
            dimArea*sp.dimensions()*vf.dimensions()
        )
    );
    faMatrix<Type>& fam = tfam.ref();

    fam.diag() += mesh.S()*sp.internalField();

    return tfam;
}

template<class Type>
tmp<faMatrix<Type>>
Sp
(
    const tmp<areaScalarField>& tsp,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    tmp<faMatrix<Type>> tfam = fam::Sp(tsp(), vf);
    tsp.clear();
    return tfam;
}


template<class Type>
tmp<faMatrix<Type>>
Sp
(
    const dimensionedScalar& sp,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    const faMesh& mesh = vf.mesh();

    tmp<faMatrix<Type>> tfam
    (
        new faMatrix<Type>
        (
            vf,
            dimArea*sp.dimensions()*vf.dimensions()
        )
    );
    faMatrix<Type>& fam = tfam.ref();

    fam.diag() += mesh.S()*sp.value();

    return tfam;
}


template<class Type>
tmp<faMatrix<Type>>
SuSp
(
    const areaScalarField& sp,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    const faMesh& mesh = vf.mesh();

    tmp<faMatrix<Type>> tfam
    (
        new faMatrix<Type>
        (
            vf,
            dimArea*sp.dimensions()*vf.dimensions()
        )
    );
    faMatrix<Type>& fam = tfam.ref();

    fam.diag() += mesh.S()*max(sp.internalField(), scalar(0));

    fam.source() -= mesh.S()*min(sp.internalField(), scalar(0))
        *vf.internalField();

    return tfam;
}

template<class Type>
tmp<faMatrix<Type>>
SuSp
(
    const tmp<areaScalarField>& tsp,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    tmp<faMatrix<Type>> tfam = fam::SuSp(tsp(), vf);
    tsp.clear();
    return tfam;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

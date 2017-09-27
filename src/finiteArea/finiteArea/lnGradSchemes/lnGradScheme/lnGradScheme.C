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

#include "fa.H"
#include "lnGradScheme.H"
#include "areaFields.H"
#include "edgeFields.H"
#include "HashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fa
{

// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

template<class Type>
tmp<lnGradScheme<Type>> lnGradScheme<Type>::New
(
    const faMesh& mesh,
    Istream& schemeData
)
{
    if (fa::debug)
    {
        InfoInFunction
            << "constructing lnGradScheme<Type>"
            << endl;
    }

    if (schemeData.eof())
    {
        FatalIOErrorInFunction(schemeData)
            << "Discretisation scheme not specified" << nl << nl
            << "Valid schemes are :" << endl
            << MeshConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    const word schemeName(schemeData);

    typename MeshConstructorTable::iterator constructorIter =
        MeshConstructorTablePtr_->find(schemeName);

    if (!constructorIter.found())
    {
        FatalIOErrorInFunction(schemeData)
            << "Unknown discretisation scheme "
            << schemeName << nl << nl
            << "Valid schemes are :" << nl
            << MeshConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    return constructorIter()(mesh, schemeData);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
lnGradScheme<Type>::~lnGradScheme()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, faePatchField, edgeMesh>>
lnGradScheme<Type>::lnGrad
(
    const GeometricField<Type, faPatchField, areaMesh>& vf,
    const tmp<edgeScalarField>& tdeltaCoeffs,
    const word& lnGradName
)
{
    const faMesh& mesh = vf.mesh();

    // construct GeometricField<Type, faePatchField, edgeMesh>
    tmp<GeometricField<Type, faePatchField, edgeMesh>> tssf
    (
        new GeometricField<Type, faePatchField, edgeMesh>
        (
            IOobject
            (
                lnGradName + "("+vf.name()+')',
                vf.instance(),
                vf.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            vf.dimensions()*tdeltaCoeffs().dimensions()
        )
    );
    GeometricField<Type, faePatchField, edgeMesh>& ssf = tssf.ref();

    // set reference to difference factors array
    const scalarField& deltaCoeffs = tdeltaCoeffs().internalField();

    // owner/neighbour addressing
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    forAll(owner, faceI)
    {
        ssf[faceI] =
            deltaCoeffs[faceI]*(vf[neighbour[faceI]] - vf[owner[faceI]]);
    }

    typename GeometricField<Type, faePatchField, edgeMesh>::Boundary& ssfb =
        ssf.boundaryFieldRef();

    forAll(vf.boundaryField(), patchI)
    {
        ssfb[patchI] = vf.boundaryField()[patchI].snGrad();
    }

    return tssf;
}


template<class Type>
tmp<GeometricField<Type, faePatchField, edgeMesh>>
lnGradScheme<Type>::lnGrad
(
    const GeometricField<Type, faPatchField, areaMesh>& vf
) const
{
    tmp<GeometricField<Type, faePatchField, edgeMesh>> tsf =
        lnGrad(vf, deltaCoeffs(vf));

    if (corrected())
    {
        tsf.ref() += correction(vf);
    }

    return tsf;
}


template<class Type>
tmp<GeometricField<Type, faePatchField, edgeMesh>>
lnGradScheme<Type>::lnGrad
(
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tvf
) const
{
    tmp<GeometricField<Type, faePatchField, edgeMesh>> tinterpVf =
        lnGrad(tvf());
    tvf.clear();
    return tinterpVf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fa

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

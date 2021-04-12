/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2021 OpenCFD Ltd.
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
#include "surfaceFields.H"
#include "fixedValueFvPatchField.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::DMDModels::STDMD::filterIndexed
(
    List<Type>& lst,
    const UList<label>& indices
)
{
    // Elements within [a, b]
    List<Type> lstWithin(indices.size());

    // Copy if frequency of element is within [a, b]
    label j = 0;
    for (const auto& i : indices)
    {
        lstWithin[j] = lst[i];
        ++j;
    }
    lst.transfer(lstWithin);
}


template<class MatrixType>
void Foam::DMDModels::STDMD::filterIndexed
(
    MatrixType& mat,
    const UList<label>& indices
)
{
    // Elements within [a, b]
    MatrixType matWithin(labelPair(mat.m(), indices.size()));

    // Copy if frequency of element is within [a, b]
    label j = 0;
    for (const auto& i : indices)
    {
        matWithin.subColumn(j) = mat.subColumn(i);
        ++j;
    }
    mat.transfer(matWithin);
}


template<class Type>
bool Foam::DMDModels::STDMD::modes()
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> SurfaceFieldType;

    if (mesh_.foundObject<VolFieldType>(fieldName_))
    {
        return calcModes<VolFieldType>();
    }
    else if (mesh_.foundObject<SurfaceFieldType>(fieldName_))
    {
        return calcModes<SurfaceFieldType>();
    }

    return false;
}


template<class GeoFieldType>
bool Foam::DMDModels::STDMD::calcModes()
{
    typedef typename GeoFieldType::value_type Type;

    // Resize the number of output variables to "nModes" if requested
    if (magsi_.size() > nModes_)
    {
        magsi_.resize(nModes_);
    }

    // Compute and write modes one by one
    const RMatrix primitiveMode(Qupper_*RxInv_);
    Qupper_.clear();
    RxInv_.clear();

    label modei = 0;
    for (const label i : magsi_)
    {
        GeoFieldType modeRe
        (
            IOobject
            (
                "modeRe_" + std::to_string(modei)
                + "_" + fieldName_ + "_" + name_,
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensioned<Type>(dimless, Zero),
            fixedValueFvPatchField<Type>::typeName
        );

        GeoFieldType modeIm
        (
            IOobject
            (
                "modeIm_" + std::to_string(modei)
                + "_" + fieldName_ + "_" + name_,
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensioned<Type>(dimless, Zero),
            fixedValueFvPatchField<Type>::typeName
        );

        if (modeRe.size() != 0 && !empty_)
        {
            if (patch_.empty())
            {
                auto& inModeRe = modeRe.primitiveFieldRef();
                auto& inModeIm = modeIm.primitiveFieldRef();

                calcMode(inModeRe, inModeIm, primitiveMode, i);
            }
            else
            {
                const label patchi = mesh_.boundaryMesh().findPatchID(patch_);

                auto& bfModeRe = modeRe.boundaryFieldRef()[patchi];
                auto& bfModeIm = modeIm.boundaryFieldRef()[patchi];

                calcMode(bfModeRe, bfModeIm, primitiveMode, i);
            }
        }

        modeRe.write();
        modeIm.write();
        ++modei;
    }

    return true;
}


template<class GeoFieldType>
typename std::enable_if
<
    std::is_same<Foam::scalar, typename GeoFieldType::value_type>::value,
    void
>::type Foam::DMDModels::STDMD::calcMode
(
    GeoFieldType& modeRe,
    GeoFieldType& modeIm,
    const RMatrix& primitiveMode,
    const label i
)
{
    const label fieldSize = modeRe.size();

    for (label p = 0; p < primitiveMode.m(); ++p)
    {
        complex mode(Zero);
        for (label q = 0; q < evecs_.m(); ++q)
        {
            mode += primitiveMode(p, q)*evecs_(q, i);
        }
        label p1 = p%fieldSize;
        modeRe[p1] = mode.real();
        modeIm[p1] = mode.imag();
    }
}


template<class GeoFieldType>
typename std::enable_if
<
    !std::is_same<Foam::scalar, typename GeoFieldType::value_type>::value,
    void
>::type Foam::DMDModels::STDMD::calcMode
(
    GeoFieldType& modeRe,
    GeoFieldType& modeIm,
    const RMatrix& primitiveMode,
    const label i
)
{
    const label fieldSize = modeRe.size();

    for (label p = 0; p < primitiveMode.m(); ++p)
    {
        complex mode(Zero);
        for (label q = 0; q < evecs_.m(); ++q)
        {
            mode += primitiveMode(p, q)*evecs_(q, i);
        }
        label p1 = p%fieldSize;
        label p2 = p/fieldSize;
        modeRe[p1][p2] = mode.real();
        modeIm[p1][p2] = mode.imag();
    }
}


// ************************************************************************* //

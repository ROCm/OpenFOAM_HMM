/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "exprOps.H"
#include "FieldOps.H"
#include "surfaceInterpolate.H"
#include "volPointInterpolation.H"
#include "interpolatePointToCell.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::expressions::volumeExpr::parseDriver::setInternalFieldResult
(
    const Field<Type>& fld
)
{
    if (isLogical_)
    {
        // Eg, volScalarField -> volLogicalField
        resultType_.replace("Scalar", "Logical");

        Field<bool> bools(fld.size());
        FieldOps::assign(bools, fld, expressions::boolOp<Type>());

        this->result().setResult(std::move(bools), this->isPointData());
    }
    else
    {
        // Deep copy
        this->result().setResult(fld, this->isPointData());
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::expressions::volumeExpr::parseDriver::setResult
(
    GeometricField<Type, fvPatchField, volMesh>* ptr,
    bool logical
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    resultField_.reset(nullptr);

    // Characteristics
    resultType_ = pTraits<fieldType>::typeName;
    isLogical_ = logical;
    fieldGeoType_ = VOLUME_DATA;

    // Assign dimensions
    if (hasDimensions_ && !logical)
    {
        ptr->dimensions().reset(resultDimensions_);
    }

    setInternalFieldResult(ptr->primitiveField());

    // Take ownership
    resultField_.reset(ptr);
}


template<class Type>
void Foam::expressions::volumeExpr::parseDriver::setResult
(
    GeometricField<Type, fvsPatchField, surfaceMesh>* ptr,
    bool logical
)
{
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> fieldType;

    resultField_.reset(nullptr);

    // Characteristics
    resultType_ = pTraits<fieldType>::typeName;
    isLogical_ = logical;
    fieldGeoType_ = FACE_DATA;

    // Assign dimensions
    if (hasDimensions_ && !logical)
    {
        ptr->dimensions().reset(resultDimensions_);
    }

    setInternalFieldResult(ptr->primitiveField());

    // Take ownership
    resultField_.reset(ptr);
}


template<class Type>
void Foam::expressions::volumeExpr::parseDriver::setResult
(
    GeometricField<Type, pointPatchField, pointMesh>* ptr,
    bool logical
)
{
    typedef GeometricField<Type, pointPatchField, pointMesh> fieldType;

    resultField_.reset(nullptr);

    // Characteristics
    resultType_ = pTraits<fieldType>::typeName;
    isLogical_ = logical;
    fieldGeoType_ = POINT_DATA;

    // Assign dimensions
    if (hasDimensions_ && !logical)
    {
        ptr->dimensions().reset(resultDimensions_);
    }

    setInternalFieldResult(ptr->primitiveField());

    // Take ownership
    resultField_.reset(ptr);
}


template<class GeomField>
const GeomField*
Foam::expressions::volumeExpr::parseDriver::isResultType() const
{
    return dynamic_cast<const GeomField*>(resultField_.get());
}


template<class GeomField>
const GeomField*
Foam::expressions::volumeExpr::parseDriver::isResultType
(
    bool logical,
    bool dieOnNull
) const
{
    const regIOobject* ptr = resultField_.get();

    if (dieOnNull && !ptr)
    {
        FatalErrorInFunction
            << "No result available. Requested "
            << pTraits<GeomField>::typeName << nl
            << exit(FatalError);
    }

    if (isLogical_ == logical)
    {
        return dynamic_cast<const GeomField*>(ptr);
    }

    return nullptr;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>>
Foam::expressions::volumeExpr::parseDriver::getVolField
(
    const word& fldName,
    bool getOldTime
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    return this->getOrReadField<fieldType>
    (
        fldName,
        true,  // mandatory
        getOldTime
    );
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::expressions::volumeExpr::parseDriver::getSurfaceField
(
    const word& fldName,
    bool getOldTime
)
{
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> fieldType;

    return this->getOrReadField<fieldType>
    (
        fldName,
        true,  // mandatory
        getOldTime
    );
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::pointPatchField, Foam::pointMesh>>
Foam::expressions::volumeExpr::parseDriver::getPointField
(
    const word& fldName,
    bool getOldTime
)
{
    typedef GeometricField<Type, pointPatchField, pointMesh> fieldType;

    return this->getOrReadPointField<fieldType>
    (
        fldName,
        true,  // mandatory
        getOldTime
    );
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>>
Foam::expressions::volumeExpr::parseDriver::newVolField
(
    const Type& val
) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    return fieldType::New
    (
        word("constant.") + word(pTraits<Type>::typeName),
        mesh(),
        dimensioned<Type>("", dimless, val)
    );
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::expressions::volumeExpr::parseDriver::newSurfaceField
(
    const Type& val
) const
{
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> fieldType;

    return fieldType::New
    (
        word("constant.") + word(pTraits<Type>::typeName),
        mesh(),
        dimensioned<Type>("", dimless, val)
    );
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::pointPatchField, Foam::pointMesh>>
Foam::expressions::volumeExpr::parseDriver::newPointField
(
    const Type& val
) const
{
    typedef GeometricField<Type, pointPatchField, pointMesh> fieldType;

    return fieldType::New
    (
        word("constant.") + word(pTraits<Type>::typeName),
        pointMesh::New(mesh()),
        dimensioned<Type>("", dimless, val)
    );
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::expressions::volumeExpr::parseDriver::cellToFace
(
    const GeometricField<Type, fvPatchField, volMesh>& field
) const
{
    return fvc::interpolate(field);
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::pointPatchField, Foam::pointMesh>>
Foam::expressions::volumeExpr::parseDriver::cellToPoint
(
    const GeometricField<Type, fvPatchField, volMesh>& field
) const
{
    volPointInterpolation interp(this->mesh());
    return interp.interpolate(field);
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type,Foam::fvPatchField,Foam::volMesh>>
Foam::expressions::volumeExpr::parseDriver::pointToCell
(
    const GeometricField<Type, pointPatchField, pointMesh>& field
) const
{
    auto tresult = newVolField<Type>();
    auto& result = tresult.ref();

    forAll(result,celli)
    {
        result[celli] = interpolatePointToCell(field, celli);
    }

    return tresult;
}


// ************************************************************************* //

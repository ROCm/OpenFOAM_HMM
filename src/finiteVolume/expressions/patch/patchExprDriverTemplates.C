/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
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

#include "primitivePatchInterpolation.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::expressions::patchExpr::parseDriver::getVariableIfAvailable
(
    const word& name
) const
{
    bool isPointVal = false;
    bool isUniformVal = false;

    tmp<Field<Type>> tfield;

    if (hasVariable(name) && variable(name).isType<Type>())
    {
        const expressions::exprResult& var = variable(name);

        isPointVal = var.isPointValue();
        isUniformVal = var.isUniform();

        tfield = var.cref<Type>().clone();
    }
    else if (isGlobalVariable<Type>(name, false))
    {
        const expressions::exprResult& var = lookupGlobal(name);

        isUniformVal = var.isUniform();

        tfield = var.cref<Type>().clone();
    }

    if (tfield.valid())
    {
        const label fldLen = tfield().size();
        const label len = (isPointVal ? this->pointSize() : this->size());

        if (returnReduce((fldLen == len), andOp<bool>()))
        {
            return tfield;
        }

        if (!isUniformVal)
        {
            WarningInFunction
                << "Variable " << name
                << " does not fit the size and is not a uniform value." << nl
                << "Using average value" << endl;
        }

        return tmp<Field<Type>>::New(this->size(), gAverage(tfield));
    }

    return tfield;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::expressions::patchExpr::parseDriver::getVolField(const word& name)
{
    return getField<Type>(name);
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::expressions::patchExpr::parseDriver::getSurfaceField(const word& name)
{
    return getField<Type>(name);
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::expressions::patchExpr::parseDriver::getPointField(const word& name)
{
    return getField<Type>(name);
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::expressions::patchExpr::parseDriver::getField(const word& name)
{
    tmp<Field<Type>> tfield = getVariableIfAvailable<Type>(name);

    if (tfield.valid())
    {
        return tfield;
    }


    typedef GeometricField<Type, fvPatchField, volMesh> vfieldType;
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> sfieldType;
    typedef GeometricField<Type, pointPatchField, pointMesh> pfieldType;

    const objectRegistry& obr = this->mesh().thisDb();

    const vfieldType* vfield = obr.findObject<vfieldType>(name);
    const sfieldType* sfield = obr.findObject<sfieldType>(name);
    const pfieldType* pfield = obr.findObject<pfieldType>(name);

    // Local, temporary storage
    tmp<vfieldType> t_vfield;
    tmp<sfieldType> t_sfield;
    tmp<pfieldType> t_pfield;

    if (searchFiles() && !vfield && !sfield && !pfield)
    {
        const word fldType = this->getTypeOfField(name);

        if (fldType == vfieldType::typeName)
        {
            t_vfield = this->readAndRegister<vfieldType>(name, mesh());
            vfield = t_vfield.get();
        }
        else if (fldType == sfieldType::typeName)
        {
            t_sfield = this->readAndRegister<sfieldType>(name, mesh());
            sfield = t_sfield.get();
        }
        else if (fldType == pfieldType::typeName)
        {
            t_pfield = this->readAndRegister<pfieldType>
            (
                name,
                pointMesh::New(mesh())
            );
            pfield = t_pfield.get();
        }
     }

    const label patchIndex = patch_.index();

    if (vfield)
    {
        return tmp<Field<Type>>::New
        (
            vfield->boundaryField()[patchIndex]
        );
    }

    if (sfield)
    {
        return tmp<Field<Type>>::New
        (
            sfield->boundaryField()[patchIndex]
        );
    }

    if (pfield)
    {
        return pfield->boundaryField()[patchIndex].patchInternalField();
    }


    FatalErrorInFunction
        << "No field '" << name << "' of type "
        << pTraits<Type>::typeName << nl << nl
        << vfieldType::typeName << " Fields: "
        << flatOutput(obr.sortedNames<vfieldType>()) << nl
        << sfieldType::typeName << " Fields: "
        << flatOutput(obr.sortedNames<sfieldType>()) << nl
        << pfieldType::typeName << " Fields: "
        << flatOutput(obr.sortedNames<pfieldType>()) << nl
        << exit(FatalError);

    return tmp<Field<Type>>::New();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::expressions::patchExpr::parseDriver::faceToPoint
(
    const Field<Type>& field
) const
{
    primitivePatchInterpolation interp(patch_.patch());

    return interp.pointToFaceInterpolate(field);
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::expressions::patchExpr::parseDriver::pointToFace
(
    const Field<Type>& field
) const
{
    primitivePatchInterpolation interp(patch_.patch());

    return interp.faceToPointInterpolate(field);
}


// ************************************************************************* //

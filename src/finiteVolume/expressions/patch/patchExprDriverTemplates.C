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

#include "primitivePatchInterpolation.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::expressions::patchExpr::parseDriver::getVariableIfAvailable
(
    const word& name
) const
{
    bool hasPointData = false;

    refPtr<expressions::exprResult> tvar;

    if (hasVariable(name) && variable(name).isType<Type>())
    {
        tvar.cref(variable(name));
        hasPointData = tvar().isPointData();
    }
    else if (isGlobalVariable<Type>(name))
    {
        tvar.cref(lookupGlobal(name));
    }


    if (tvar.valid())
    {
        const auto& var = tvar.cref();
        const Field<Type>& vals = var.cref<Type>();

        const label len = (hasPointData ? this->pointSize() : this->size());

        if (returnReduce((vals.size() == len), andOp<bool>()))
        {
            // Return a copy of the field
            return tmp<Field<Type>>::New(vals);
        }

        if (!var.isUniform())
        {
            WarningInFunction
                << "Variable " << name
                << " is nonuniform and does not fit the size"
                << ". Using average" << endl;
        }

        return tmp<Field<Type>>::New(this->size(), gAverage(vals));
    }

    return nullptr;
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

    const objectRegistry& obr = this->mesh().thisDb();
    const label patchIndex = patch_.index();


    // Field types

    typedef GeometricField<Type, fvPatchField, volMesh> vfieldType;
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> sfieldType;
    typedef GeometricField<Type, pointPatchField, pointMesh> pfieldType;

    // Local, temporary storage and/or lookup values
    bool found = false;
    tmp<vfieldType> vfield;
    tmp<sfieldType> sfield;
    tmp<pfieldType> pfield;

    for (int checki = 0; !found && checki < 2; ++checki)
    {
        // Check 0: object context (first)
        // Check 1: regular objectRegistry
        const regIOobject* ioptr =
        (
            (checki == 0)
          ? exprDriver::cfindContextIOobject(name)
          : obr.cfindIOobject(name)
        );
        if (!ioptr) continue;

        if (!found)
        {
            vfield.cref(dynamic_cast<const vfieldType*>(ioptr));
            found = vfield.valid();
        }
        if (!found)
        {
            sfield.cref(dynamic_cast<const sfieldType*>(ioptr));
            found = sfield.valid();
        }
        if (!found)
        {
            pfield.cref(dynamic_cast<const pfieldType*>(ioptr));
            found = pfield.valid();
        }
    }


    // Finally, search files if necessary (and permitted)
    if (!found && searchFiles())
    {
        const word fldType = this->getTypeOfField(name);

        if (fldType == vfieldType::typeName)
        {
            vfield = this->readAndRegister<vfieldType>(name, mesh());
        }
        else if (fldType == sfieldType::typeName)
        {
            sfield = this->readAndRegister<sfieldType>(name, mesh());
        }
        else if (fldType == pfieldType::typeName)
        {
            pfield = this->readAndRegister<pfieldType>
            (
                name,
                pointMesh::New(mesh())
            );
        }
    }


    if (vfield.valid())
    {
        return tmp<Field<Type>>::New
        (
            vfield().boundaryField()[patchIndex]
        );
    }
    if (sfield.valid())
    {
        return tmp<Field<Type>>::New
        (
            sfield().boundaryField()[patchIndex]
        );
    }
    if (pfield.valid())
    {
        return pfield().boundaryField()[patchIndex].patchInternalField();
    }


    FatalErrorInFunction
        << "No field '" << name << "' of type "
        << pTraits<Type>::typeName << nl << nl;

    FatalError
        << vfieldType::typeName << " Fields: "
        << flatOutput(obr.sortedNames<vfieldType>()) << nl;

    FatalError
        << sfieldType::typeName << " Fields: "
        << flatOutput(obr.sortedNames<sfieldType>()) << nl;

    FatalError
        << pfieldType::typeName << " Fields: "
        << flatOutput(obr.sortedNames<pfieldType>()) << nl;

    FatalError
        << exit(FatalError);

    return tmp<Field<Type>>::New();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::expressions::patchExpr::parseDriver::patchInternalField
(
    const word& name
)
{
    tmp<Field<Type>> tfield = getVariableIfAvailable<Type>(name);

    if (tfield.valid())
    {
        return tfield;
    }

    const objectRegistry& obr = this->mesh().thisDb();
    const label patchIndex = patch_.index();


    // Field types

    typedef GeometricField<Type, fvPatchField, volMesh> vfieldType;
    typedef GeometricField<Type, pointPatchField, pointMesh> pfieldType;

    // Local, temporary storage and/or lookup values
    bool found = false;
    tmp<vfieldType> vfield;
    tmp<pfieldType> pfield;

    for (int checki = 0; !found && checki < 2; ++checki)
    {
        // Check 0: object context (first)
        // Check 1: regular objectRegistry
        const regIOobject* ioptr =
        (
            (checki == 0)
          ? exprDriver::cfindContextIOobject(name)
          : obr.cfindIOobject(name)
        );
        if (!ioptr) continue;

        if (!found)
        {
            vfield.cref(dynamic_cast<const vfieldType*>(ioptr));
            found = vfield.valid();
        }
        if (!found)
        {
            pfield.cref(dynamic_cast<const pfieldType*>(ioptr));
            found = pfield.valid();
        }
    }


    // Finally, search files if necessary (and permitted)
    if (!found && searchFiles())
    {
        const word fldType = this->getTypeOfField(name);

        if (fldType == vfieldType::typeName)
        {
            vfield = this->readAndRegister<vfieldType>(name, mesh());
        }
        else if (fldType == pfieldType::typeName)
        {
            pfield = this->readAndRegister<pfieldType>
            (
                name,
                pointMesh::New(mesh())
            );
        }
    }


    if (vfield.valid())
    {
        return vfield().boundaryField()[patchIndex].patchInternalField();
    }
    if (pfield.valid())
    {
        return pfield().boundaryField()[patchIndex].patchInternalField();
    }


    FatalErrorInFunction
        << "No field '" << name << "' of type "
        << pTraits<Type>::typeName << nl << nl;

    FatalError
        << vfieldType::typeName << " Fields: "
        << flatOutput(obr.sortedNames<vfieldType>()) << nl;

    FatalError
        << pfieldType::typeName << " Fields: "
        << flatOutput(obr.sortedNames<pfieldType>()) << nl;

    FatalError
        << exit(FatalError);

    return tmp<Field<Type>>::New();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::expressions::patchExpr::parseDriver::patchNeighbourField
(
    const word& name
)
{
    tmp<Field<Type>> tfield = getVariableIfAvailable<Type>(name);

    if (tfield.valid())
    {
        return tfield;
    }

    const objectRegistry& obr = this->mesh().thisDb();
    const label patchIndex = patch_.index();


    // Field types

    typedef GeometricField<Type, fvPatchField, volMesh> vfieldType;

    // Local, temporary storage and/or lookup values
    bool found = false;
    tmp<vfieldType> vfield;

    for (int checki = 0; !found && checki < 2; ++checki)
    {
        // Check 0: object context (first)
        // Check 1: regular objectRegistry
        const regIOobject* ioptr =
        (
            (checki == 0)
          ? exprDriver::cfindContextIOobject(name)
          : obr.cfindIOobject(name)
        );
        if (!ioptr) continue;

        if (!found)
        {
            vfield.cref(dynamic_cast<const vfieldType*>(ioptr));
            found = vfield.valid();
        }
    }


    // Finally, search files if necessary (and permitted)
    if (!found && searchFiles())
    {
        const word fldType = this->getTypeOfField(name);

        if (fldType == vfieldType::typeName)
        {
            vfield = this->readAndRegister<vfieldType>(name, mesh());
        }
    }


    if (vfield.valid())
    {
        return vfield().boundaryField()[patchIndex].patchNeighbourField();
    }


    FatalErrorInFunction
        << "No field '" << name << "' of type "
        << pTraits<Type>::typeName << nl << nl;

    FatalError
        << vfieldType::typeName << " Fields: "
        << flatOutput(obr.sortedNames<vfieldType>()) << nl;

    FatalError
        << exit(FatalError);

    return tmp<Field<Type>>::New();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::expressions::patchExpr::parseDriver::patchNormalField
(
    const word& name
)
{
    tmp<Field<Type>> tfield = getVariableIfAvailable<Type>(name);

    if (tfield.valid())
    {
        return tfield;
    }

    const objectRegistry& obr = this->mesh().thisDb();
    const label patchIndex = patch_.index();


    // Field types

    typedef GeometricField<Type, fvPatchField, volMesh> vfieldType;

    // Local, temporary storage and/or lookup values
    bool found = false;
    tmp<vfieldType> vfield;

    for (int checki = 0; !found && checki < 2; ++checki)
    {
        // Check 0: object context (first)
        // Check 1: regular objectRegistry
        const regIOobject* ioptr =
        (
            (checki == 0)
          ? exprDriver::cfindContextIOobject(name)
          : obr.cfindIOobject(name)
        );
        if (!ioptr) continue;

        if (!found)
        {
            vfield.cref(dynamic_cast<const vfieldType*>(ioptr));
            found = vfield.valid();
        }
    }


    // Finally, search files if necessary (and permitted)
    if (!found && searchFiles())
    {
        const word fldType = this->getTypeOfField(name);

        if (fldType == vfieldType::typeName)
        {
            vfield = this->readAndRegister<vfieldType>(name, mesh());
        }
    }


    if (vfield.valid())
    {
        return vfield().boundaryField()[patchIndex].snGrad();
    }


    FatalErrorInFunction
        << "No field '" << name << "' of type "
        << pTraits<Type>::typeName << nl << nl;

    FatalError
        << vfieldType::typeName << " Fields: "
        << flatOutput(obr.sortedNames<vfieldType>()) << nl;

    FatalError
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

    return interp.faceToPointInterpolate(field);
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::expressions::patchExpr::parseDriver::pointToFace
(
    const Field<Type>& field
) const
{
    primitivePatchInterpolation interp(patch_.patch());

    return interp.pointToFaceInterpolate(field);
}


// ************************************************************************* //

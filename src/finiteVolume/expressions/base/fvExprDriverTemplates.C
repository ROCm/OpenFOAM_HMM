/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2010-2018 Bernhard Gschaider <bgschaid@hfd-research.com>
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

#include "surfaceMesh.H"
#include "fvsPatchField.H"
#include "pointPatchField.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
bool Foam::expressions::fvExprDriver::isGlobalVariable
(
    const word& name,
    bool isPointVal,
    label expectedSize
) const
{
    DebugInfo
        << "Looking for global" << (isPointVal ? " point" : "")
        << " field name:" << name;

    const exprResult& result = lookupGlobal(name);

    DebugInfo
        << " - found (" << result.valueType() << ' ' << result.isPointValue() << ')';


    bool good = (result.isType<Type>() && result.isPointValue(isPointVal));

    // Do size checking if requested
    if (good && expectedSize >= 0)
    {
        good = (result.size() == expectedSize);
        reduce(good, andOp<bool>());

        if (debug && !good)
        {
            Info<< " size is";
        }
    }

    DebugInfo << (good ? " good" : " bad") << endl;

    return good;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::expressions::fvExprDriver::getVariable
(
    const word& name,
    label expectedSize,
    const bool mandatory
) const
{
    tmp<Field<Type>> tresult;

    bool isSingleValue = false;

    if (hasVariable(name) && variable(name).isType<Type>())
    {
        isSingleValue = variable(name).isUniform();
        tresult = variable(name).cref<Type>().clone();
    }
    else if (isGlobalVariable<Type>(name, false))
    {
        const exprResult& var = lookupGlobal(name);

        isSingleValue = var.isUniform();

        tresult = var.cref<Type>().clone();
    }

    if (tresult.valid())
    {
        if
        (
            expectedSize < 0
         || returnReduce((tresult->size() == expectedSize), andOp<bool>())
        )
        {
            return tresult;
        }

        if (!isSingleValue)
        {
            WarningInFunction
                << "Variable " << name
                << " is not a single value and does not fit the size "
                << expectedSize << ". Using average" << endl;
        }

        return tmp<Field<Type>>::New(expectedSize, gAverage(tresult()));
    }

    if (mandatory)
    {
        FatalErrorInFunction
            << "Variable (" << name << ") not found." << nl
            << exit(FatalError);
    }

    return nullptr;
}


template<class Type>
bool Foam::expressions::fvExprDriver::foundField
(
    const word& name
) const
{
    if (debug)
    {
        Info<< "fvExprDriver::foundField. Name: " << name
            << " Type: " << Type::typeName
            << " registry:" << searchInMemory()
            << " disk:" << searchFiles() << endl;
    }

    // if (std::is_void<Type>::value) ...

    if (searchInMemory())
    {
        const regIOobject* ioptr =
            this->mesh().findObject<regIOobject>(name);

        if (this->mesh().foundObject<Type>(name))
        {
            if (debug)
            {
                Info<< "Found registered: " << name << endl;
            }
            return true;
        }

        if (debug)
        {
            Info<< "Registered " << name;

            if (ioptr)
            {
                Info<< " type:" << ioptr->headerClassName();
            }
            Info<< ", not type:" << Type::typeName << nl;
        }
    }


    if (searchFiles() && getTypeOfField(name) == Type::typeName)
    {
        if (debug)
        {
            Info<< "Found file: " << name << nl;
        }
        return true;
    }
    else
    {
        if (debug)
        {
            Info<< name << " not found" << endl;
        }
    }

    return false;
}


template<class Type>
bool Foam::expressions::fvExprDriver::isField
(
    const word& name,
    bool isPointVal,
    label
) const
{
    if (debug)
    {
        Info<< "fvExprDriver::isField <" << name << '>' << endl;
    }

    typedef GeometricField<Type, fvPatchField, volMesh> vfieldType;
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> sfieldType;
    typedef GeometricField<Type, pointPatchField, pointMesh> pfieldType;

    return
    (
        isPointVal
      ? this->foundField<pfieldType>(name)
      :
        (
            this->foundField<vfieldType>(name)
         || this->foundField<sfieldType>(name)
        )
    );
}


template<class GeomField, class Mesh>
Foam::tmp<GeomField> Foam::expressions::fvExprDriver::getOrReadFieldImpl
(
    const word& name,
    const Mesh& meshRef,
    bool mandatory,
    bool getOldTime
)
{
    typedef typename GeomField::value_type Type;

    if (debug)
    {
        Info<< "fvExprDriver::getOrReadField <" << name
            << "> Type: " << GeomField::typeName << endl;
    }

    tmp<GeomField> tfield;

    if
    (
        (hasVariable(name) && variable(name).isType<Type>())
     || isGlobalVariable<Type>(name, false)
    )
    {
        if (debug)
        {
            Info<< "Getting " << name << " from variables" << endl;
        }

        if (debug)
        {
            Info<< "Creating field " << name << " of type "
                << GeomField::typeName << nl;
        }

        tfield.reset
        (
            GeomField::New(name, meshRef, dimensioned<Type>(Zero))
        );

        GeomField& fld = tfield.ref();

        if (debug)
        {
            Info<< "New field: " << name << " ownedByRegistry"
                << fld.ownedByRegistry() << endl;
        }

        Field<Type> vals;

        if (hasVariable(name) && variable(name).isType<Type>())
        {
            vals = variable(name).cref<Type>();
        }
        else
        {
            vals = lookupGlobal(name).cref<Type>();
        }

        if (debug)
        {
            Pout<< "sizes: " << vals.size() << ' ' << fld.size() << endl;
        }

        if (returnReduce((vals.size() == fld.size()), andOp<bool>()))
        {
            fld.primitiveFieldRef() = vals;
        }
        else
        {
            Type avg = gAverage(vals);

            bool noWarn = false;

            if (!noWarn)
            {
                MinMax<Type> range = gMinMax(vals);

                if (range.mag() > SMALL)
                {
                    WarningInFunction
                        << "The min/max ranges differ " << range
                        << " - using average " << avg << nl;
                }
            }

            fld.primitiveFieldRef() = avg;
        }

        correctField(fld);

        return tfield;
    }


    const objectRegistry& obr = meshRef.thisDb();

    if (searchInMemory() && obr.foundObject<GeomField>(name))
    {
        if (debug)
        {
            Info<< "Retrieve registered: " << name << nl;
        }

        const GeomField& origFld = obr.lookupObject<GeomField>(name);

        // Avoid shadowing the original object

        tfield.reset
        (
            GeomField::New(name + "_exprDriverCopy", origFld)
        );

        if (getOldTime)
        {
            if (debug)
            {
                Info<< "Getting oldTime of " << name << " has "
                    << origFld.nOldTimes() << endl;
            }

            if (!origFld.nOldTimes() && this->prevIterIsOldTime())
            {
                if (debug)
                {
                    Info<< "No oldTime, using previous iteration" << endl;
                }
                tfield.ref().oldTime() = origFld.prevIter();
            }
        }
    }
    else if (searchFiles() && getTypeOfField(name) == GeomField::typeName)
    {
        if (debug)
        {
            Info<< "Reading " << name << " from disc" << endl;
        }

        tfield.reset
        (
            this->readAndRegister<GeomField>(name, meshRef)
        );
        // oldTime automatically read
    }

    if (debug)
    {
        Info<< "field: valid()=" << tfield.valid() << endl;
    }

    if (tfield.valid())
    {
        GeomField& fld = tfield.ref();

        if (debug)
        {
            Info<< "Valid " << name << " found. Removing dimensions" << nl;
        }

        fld.dimensions().clear();

        if (fld.nOldTimes())
        {
            if (debug)
            {
                Info<< "Removing dimensions of oldTime of " << name
                    << " has " << fld.nOldTimes() << nl;
            }

            // Switch dimension checking off
            const int oldDebug = dimensionSet::debug;
            dimensionSet::debug = 0;

            // go through ALL old times
            GeomField* fp = &(fld);

            while (fp->nOldTimes())
            {
                fp = &(fp->oldTime());
                fp->dimensions().clear();
            }

            // Restore old value of dimension checking
            dimensionSet::debug = oldDebug;
        }
    }
    else if (mandatory)
    {
        FatalErrorInFunction
            << "Could not find field " << name
            << " in registry or on file-system" << nl
            << exit(FatalError);
    }

    return tfield;
}


template<class T>
Foam::autoPtr<T> Foam::expressions::fvExprDriver::getTopoSet
(
    const fvMesh& mesh,
    const word& name,
    SetOrigin& origin
) const
{
    // Avoid possible name clashes
    const word regName = name + "RegisteredNameFor" + T::typeName;

    if (debug)
    {
        Info<< "Looking for " << T::typeName << " named " << name;

        Info<< " or registered as " << regName << " with mesh "
            << "Caching:" << cacheSets()
            << " Found:" << (mesh.foundObject<T>(name))
            << " Found registered:" << mesh.foundObject<T>(regName)
            << endl;
    }


    origin = SetOrigin::INVALID;
    autoPtr<T> setPtr;

    if
    (
        !cacheSets()
     ||
        (
            !mesh.thisDb().foundObject<T>(regName)
         && !mesh.thisDb().foundObject<T>(name)
        )
    )
    {
        if (debug)
        {
            Info<< "Constructing new " << T::typeName << ' ' << name << nl;

            if (debug > 1)
            {
                Pout<< mesh.thisDb().names();
            }
        }

        origin = SetOrigin::FILE;
        setPtr.reset(new T(mesh, name, IOobject::MUST_READ));

        if (cacheSets())
        {
            if (debug)
            {
                Info<< "Registering a copy of " << name << " with mesh" << nl;
            }

            autoPtr<T> toCache(new T(mesh, regName, *setPtr));
            toCache->store(toCache);
        }
    }
    else
    {
        const T* ptr = mesh.thisDb().findObject<T>(name);

        if (ptr)
        {
            if (debug)
            {
                Info<< "Getting existing " << name << endl;
            }

            origin = SetOrigin::MEMORY;
            setPtr.reset(new T(mesh, name, *ptr));
        }
        else
        {
            if (debug)
            {
                Info<< "Getting existing " << regName << endl;
            }

            origin = SetOrigin::CACHE;
            setPtr.reset(new T(mesh, name, mesh.lookupObject<T>(regName)));
        }
    }


    return setPtr;
}


template<class T>
bool Foam::expressions::fvExprDriver::updateSet
(
    autoPtr<T>& setPtr,
    const word& name,
    SetOrigin origin
) const
{
    const label oldSize = setPtr->size();

    bool updated = false;
    const polyMesh& mesh = dynamic_cast<const polyMesh&>(setPtr->db());

    if (debug)
    {
        Info<< "UpdateSet: " << setPtr->name() << " Id: " << name
            << " Origin: " << int(origin) << endl;
    }

    switch (origin)
    {
        case SetOrigin::FILE:
        {
            IOobject header
            (
                name,
                mesh.time().timeName(),
                polyMesh::meshSubDir/"sets",
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            );

            if (header.typeHeaderOk<T>())
            {
                if (debug)
                {
                    Pout<< "Rereading from "
                        << header.localFilePath(T::typeName) << endl;
                }
                setPtr.reset(new T(header));
                updated = true;
            }
            break;
        }

        case SetOrigin::NEW:
        case SetOrigin::MEMORY:
        case SetOrigin::CACHE:
        {
            if (origin == SetOrigin::NEW)
            {
                WarningInFunction
                    << "State NEW shouldn't exist"
                    << endl;
            }

            word sName = name;

            const T* ptr = mesh.thisDb().findObject<T>(name);

            if (ptr)
            {
                if (debug)
                {
                    Info<< "Found " << name
                        << " and rereading it" << endl;
                }

                setPtr.reset(new T(mesh, name, *ptr));
            }
            else
            {
                FatalErrorInFunction
                    << name << " Not found" << endl
                    << "In registry: " << mesh.thisDb().names() << endl
                    << exit(FatalError);
            }
            updated = true;
            break;
        }

        case INVALID:
        {
            FatalErrorInFunction
                << T::typeName << ' ' << name << " is invalid" << endl
                << exit(FatalError);
            break;
        }

        default:
        {
            if (debug)
            {
                Info<< "Origin " << int(origin) << " not implemented" << endl;
            }
            break;
        }
    }

    if (debug)
    {
        Pout<< name << " old size " << oldSize << " new: "
            << setPtr->size() << endl;
    }

    return updated;
}


// ************************************************************************* //

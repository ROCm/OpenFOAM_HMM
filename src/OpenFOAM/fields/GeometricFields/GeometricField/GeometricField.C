/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2015-2021 OpenCFD Ltd.
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

#include "GeometricField.H"
#include "Time.H"
#include "demandDrivenData.H"
#include "dictionary.H"
#include "localIOdictionary.H"
#include "data.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define checkField(gf1, gf2, op)                                               \
if ((gf1).mesh() != (gf2).mesh())                                              \
{                                                                              \
    FatalErrorInFunction                                                       \
        << "different mesh for fields "                                        \
        << (gf1).name() << " and " << (gf2).name()                             \
        << " during operation " <<  op                                         \
        << abort(FatalError);                                                  \
}


// * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::readFields
(
    const dictionary& dict
)
{
    Internal::readField(dict, "internalField");

    boundaryField_.readField(*this, dict.subDict("boundaryField"));

    Type refLevel;

    if (dict.readIfPresent("referenceLevel", refLevel))
    {
        Field<Type>::operator+=(refLevel);

        forAll(boundaryField_, patchi)
        {
            boundaryField_[patchi] == boundaryField_[patchi] + refLevel;
        }
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::readFields()
{
    const localIOdictionary dict
    (
        IOobject
        (
            this->name(),
            this->instance(),
            this->local(),
            this->db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        ),
        typeName
    );

    this->close();

    readFields(dict);
}


template<class Type, template<class> class PatchField, class GeoMesh>
bool Foam::GeometricField<Type, PatchField, GeoMesh>::readIfPresent()
{
    if
    (
        this->readOpt() == IOobject::MUST_READ
     || this->readOpt() == IOobject::MUST_READ_IF_MODIFIED
    )
    {
        WarningInFunction
            << "read option IOobject::MUST_READ or MUST_READ_IF_MODIFIED"
            << " suggests that a read constructor for field " << this->name()
            << " would be more appropriate." << endl;
    }
    else if
    (
        this->readOpt() == IOobject::READ_IF_PRESENT
     && this->template typeHeaderOk<GeometricField<Type, PatchField, GeoMesh>>
        (
            true
        )
    )
    {
        readFields();

        // Check compatibility between field and mesh
        if (this->size() != GeoMesh::size(this->mesh()))
        {
            FatalIOErrorInFunction(this->readStream(typeName))
                << "   number of field elements = " << this->size()
                << " number of mesh elements = "
                << GeoMesh::size(this->mesh())
                << exit(FatalIOError);
        }

        readOldTimeIfPresent();

        return true;
    }

    return false;
}


template<class Type, template<class> class PatchField, class GeoMesh>
bool Foam::GeometricField<Type, PatchField, GeoMesh>::readOldTimeIfPresent()
{
    // Read the old time field if present
    IOobject field0
    (
        this->name()  + "_0",
        this->time().timeName(),
        this->db(),
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE,
        this->registerObject()
    );

    if
    (
        field0.template typeHeaderOk<GeometricField<Type, PatchField, GeoMesh>>
        (
            true
        )
    )
    {
        DebugInFunction
            << "Reading old time level for field" << nl << this->info() << endl;

        field0Ptr_ = new GeometricField<Type, PatchField, GeoMesh>
        (
            field0,
            this->mesh()
        );

        // Ensure the old time field oriented flag is set to the parent's state
        // Note: required for backwards compatibility in case of restarting from
        // an old run where the oriented state may not have been set
        field0Ptr_->oriented() = this->oriented();

        field0Ptr_->timeIndex_ = timeIndex_ - 1;

        if (!field0Ptr_->readOldTimeIfPresent())
        {
            field0Ptr_->oldTime();
        }

        return true;
    }

    return false;
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const IOobject& io,
    const Mesh& mesh,
    const dimensionSet& ds,
    const word& patchFieldType
)
:
    Internal(io, mesh, ds, false),
    timeIndex_(this->time().timeIndex()),
    field0Ptr_(nullptr),
    fieldPrevIterPtr_(nullptr),
    boundaryField_(mesh.boundary(), *this, patchFieldType)
{
    DebugInFunction
        << "Creating temporary" << nl << this->info() << endl;

    readIfPresent();
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const IOobject& io,
    const Mesh& mesh,
    const dimensionSet& ds,
    const wordList& patchFieldTypes,
    const wordList& actualPatchTypes
)
:
    Internal(io, mesh, ds, false),
    timeIndex_(this->time().timeIndex()),
    field0Ptr_(nullptr),
    fieldPrevIterPtr_(nullptr),
    boundaryField_(mesh.boundary(), *this, patchFieldTypes, actualPatchTypes)
{
    DebugInFunction
        << "Creating temporary" << nl << this->info() << endl;

    readIfPresent();
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const IOobject& io,
    const Mesh& mesh,
    const dimensioned<Type>& dt,
    const word& patchFieldType
)
:
    Internal(io, mesh, dt, false),
    timeIndex_(this->time().timeIndex()),
    field0Ptr_(nullptr),
    fieldPrevIterPtr_(nullptr),
    boundaryField_(mesh.boundary(), *this, patchFieldType)
{
    DebugInFunction
        << "Creating temporary" << nl << this->info() << endl;

    boundaryField_ == dt.value();

    readIfPresent();
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const IOobject& io,
    const Mesh& mesh,
    const dimensioned<Type>& dt,
    const wordList& patchFieldTypes,
    const wordList& actualPatchTypes
)
:
    Internal(io, mesh, dt, false),
    timeIndex_(this->time().timeIndex()),
    field0Ptr_(nullptr),
    fieldPrevIterPtr_(nullptr),
    boundaryField_(mesh.boundary(), *this, patchFieldTypes, actualPatchTypes)
{
    DebugInFunction
        << "Creating temporary" << nl << this->info() << endl;

    boundaryField_ == dt.value();

    readIfPresent();
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const IOobject& io,
    const Internal& diField,
    const PtrList<PatchField<Type>>& ptfl
)
:
    Internal(io, diField),
    timeIndex_(this->time().timeIndex()),
    field0Ptr_(nullptr),
    fieldPrevIterPtr_(nullptr),
    boundaryField_(this->mesh().boundary(), *this, ptfl)
{
    DebugInFunction
        << "Copy construct from components" << nl << this->info() << endl;

    readIfPresent();
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const IOobject& io,
    const Mesh& mesh,
    const dimensionSet& ds,
    const Field<Type>& iField,
    const word& patchFieldType
)
:
    Internal(io, mesh, ds, iField),
    timeIndex_(this->time().timeIndex()),
    field0Ptr_(nullptr),
    fieldPrevIterPtr_(nullptr),
    boundaryField_(mesh.boundary(), *this, patchFieldType)
{
    DebugInFunction
        << "Copy construct from internal field" << nl << this->info() << endl;

    readIfPresent();
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const IOobject& io,
    const Mesh& mesh,
    const dimensionSet& ds,
    Field<Type>&& iField,
    const word& patchFieldType
)
:
    Internal(io, mesh, ds, std::move(iField)),
    timeIndex_(this->time().timeIndex()),
    field0Ptr_(nullptr),
    fieldPrevIterPtr_(nullptr),
    boundaryField_(mesh.boundary(), *this, patchFieldType)
{
    DebugInFunction
        << "Move construct from internal field" << nl << this->info() << endl;

    readIfPresent();
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const IOobject& io,
    const Mesh& mesh,
    const dimensionSet& ds,
    const Field<Type>& iField,
    const PtrList<PatchField<Type>>& ptfl
)
:
    Internal(io, mesh, ds, iField),
    timeIndex_(this->time().timeIndex()),
    field0Ptr_(nullptr),
    fieldPrevIterPtr_(nullptr),
    boundaryField_(mesh.boundary(), *this, ptfl)
{
    DebugInFunction
        << "Copy construct from components" << nl << this->info() << endl;

    readIfPresent();
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const IOobject& io,
    const Mesh& mesh,
    const bool readOldTime
)
:
    Internal(io, mesh, dimless, false),
    timeIndex_(this->time().timeIndex()),
    field0Ptr_(nullptr),
    fieldPrevIterPtr_(nullptr),
    boundaryField_(mesh.boundary())
{
    readFields();

    // Check compatibility between field and mesh

    if (this->size() != GeoMesh::size(this->mesh()))
    {
        FatalIOErrorInFunction(this->readStream(typeName))
            << "   number of field elements = " << this->size()
            << " number of mesh elements = " << GeoMesh::size(this->mesh())
            << exit(FatalIOError);
    }

    if (readOldTime)
    {
        readOldTimeIfPresent();
    }

    DebugInFunction
        << "Finishing read-construction" << nl << this->info() << endl;
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const IOobject& io,
    const Mesh& mesh,
    const dictionary& dict
)
:
    Internal(io, mesh, dimless, false),
    timeIndex_(this->time().timeIndex()),
    field0Ptr_(nullptr),
    fieldPrevIterPtr_(nullptr),
    boundaryField_(mesh.boundary())
{
    readFields(dict);

    // Check compatibility between field and mesh

    if (this->size() != GeoMesh::size(this->mesh()))
    {
        FatalIOErrorInFunction(dict)
            << "   number of field elements = " << this->size()
            << " number of mesh elements = " << GeoMesh::size(this->mesh())
            << exit(FatalIOError);
    }

    DebugInFunction
        << "Finishing dictionary-construct" << nl << this->info() << endl;
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const GeometricField<Type, PatchField, GeoMesh>& gf
)
:
    Internal(gf),
    timeIndex_(gf.timeIndex()),
    field0Ptr_(nullptr),
    fieldPrevIterPtr_(nullptr),
    boundaryField_(*this, gf.boundaryField_)
{
    DebugInFunction
        << "Copy construct" << nl << this->info() << endl;

    if (gf.field0Ptr_)
    {
        field0Ptr_ = new GeometricField<Type, PatchField, GeoMesh>
        (
            *gf.field0Ptr_
        );
    }

    this->writeOpt(IOobject::NO_WRITE);
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const tmp<GeometricField<Type, PatchField, GeoMesh>>& tgf
)
:
    Internal(tgf.constCast(), tgf.movable()),
    timeIndex_(tgf().timeIndex()),
    field0Ptr_(nullptr),
    fieldPrevIterPtr_(nullptr),
    boundaryField_(*this, tgf().boundaryField_)
{
    DebugInFunction
        << "Constructing from tmp" << nl << this->info() << endl;

    this->writeOpt(IOobject::NO_WRITE);

    tgf.clear();
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const IOobject& io,
    const GeometricField<Type, PatchField, GeoMesh>& gf
)
:
    Internal(io, gf),
    timeIndex_(gf.timeIndex()),
    field0Ptr_(nullptr),
    fieldPrevIterPtr_(nullptr),
    boundaryField_(*this, gf.boundaryField_)
{
    DebugInFunction
        << "Copy construct, resetting IO params" << nl
        << this->info() << endl;

    if (!readIfPresent() && gf.field0Ptr_)
    {
        field0Ptr_ = new GeometricField<Type, PatchField, GeoMesh>
        (
            io.name() + "_0",
            *gf.field0Ptr_
        );
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const IOobject& io,
    const tmp<GeometricField<Type, PatchField, GeoMesh>>& tgf
)
:
    Internal(io, tgf.constCast(), tgf.movable()),
    timeIndex_(tgf().timeIndex()),
    field0Ptr_(nullptr),
    fieldPrevIterPtr_(nullptr),
    boundaryField_(*this, tgf().boundaryField_)
{
    DebugInFunction
        << "Constructing from tmp resetting IO params" << nl
        << this->info() << endl;

    tgf.clear();

    readIfPresent();
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const word& newName,
    const GeometricField<Type, PatchField, GeoMesh>& gf
)
:
    Internal(newName, gf),
    timeIndex_(gf.timeIndex()),
    field0Ptr_(nullptr),
    fieldPrevIterPtr_(nullptr),
    boundaryField_(*this, gf.boundaryField_)
{
    DebugInFunction
        << "Copy construct, resetting name" << nl
        << this->info() << endl;

    if (!readIfPresent() && gf.field0Ptr_)
    {
        field0Ptr_ = new GeometricField<Type, PatchField, GeoMesh>
        (
            newName + "_0",
            *gf.field0Ptr_
        );
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const word& newName,
    const tmp<GeometricField<Type, PatchField, GeoMesh>>& tgf
)
:
    Internal(newName, tgf.constCast(), tgf.movable()),
    timeIndex_(tgf().timeIndex()),
    field0Ptr_(nullptr),
    fieldPrevIterPtr_(nullptr),
    boundaryField_(*this, tgf().boundaryField_)
{
    DebugInFunction
        << "Constructing from tmp resetting name" << nl
        << this->info() << endl;

    tgf.clear();
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const IOobject& io,
    const GeometricField<Type, PatchField, GeoMesh>& gf,
    const word& patchFieldType
)
:
    Internal(io, gf),
    timeIndex_(gf.timeIndex()),
    field0Ptr_(nullptr),
    fieldPrevIterPtr_(nullptr),
    boundaryField_(this->mesh().boundary(), *this, patchFieldType)
{
    DebugInFunction
        << "Copy construct, resetting IO params" << nl
        << this->info() << endl;

    boundaryField_ == gf.boundaryField_;

    if (!readIfPresent() && gf.field0Ptr_)
    {
        field0Ptr_ = new GeometricField<Type, PatchField, GeoMesh>
        (
            io.name() + "_0",
            *gf.field0Ptr_
        );
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const IOobject& io,
    const GeometricField<Type, PatchField, GeoMesh>& gf,
    const wordList& patchFieldTypes,
    const wordList& actualPatchTypes
)
:
    Internal(io, gf),
    timeIndex_(gf.timeIndex()),
    field0Ptr_(nullptr),
    fieldPrevIterPtr_(nullptr),
    boundaryField_
    (
        this->mesh().boundary(),
        *this,
        patchFieldTypes,
        actualPatchTypes
    )
{
    DebugInFunction
        << "Copy construct, resetting IO params and patch types" << nl
        << this->info() << endl;

    boundaryField_ == gf.boundaryField_;

    if (!readIfPresent() && gf.field0Ptr_)
    {
        field0Ptr_ = new GeometricField<Type, PatchField, GeoMesh>
        (
            io.name() + "_0",
            *gf.field0Ptr_
        );
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const IOobject& io,
    const GeometricField<Type, PatchField, GeoMesh>& gf,
    const labelList& patchIDs,
    const word& patchFieldType
)
:
    Internal(io, gf),
    timeIndex_(gf.timeIndex()),
    field0Ptr_(nullptr),
    fieldPrevIterPtr_(nullptr),
    boundaryField_(*this, gf.boundaryField_, patchIDs, patchFieldType)
{
    DebugInFunction
        << "Copy construct, resetting IO params and setting patchFieldType "
        << "for patchIDs" << nl
        << this->info() << endl;

    if (!readIfPresent() && gf.field0Ptr_)
    {
        field0Ptr_ = new GeometricField<Type, PatchField, GeoMesh>
        (
            io.name() + "_0",
            *gf.field0Ptr_
        );
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::GeometricField
(
    const IOobject& io,
    const tmp<GeometricField<Type, PatchField, GeoMesh>>& tgf,
    const wordList& patchFieldTypes,
    const wordList& actualPatchTypes
)
:
    Internal(io, tgf.constCast(), tgf.movable()),
    timeIndex_(tgf().timeIndex()),
    field0Ptr_(nullptr),
    fieldPrevIterPtr_(nullptr),
    boundaryField_
    (
        this->mesh().boundary(),
        *this,
        patchFieldTypes,
        actualPatchTypes
    )
{
    DebugInFunction
        << "Constructing from tmp resetting IO params and patch types" << nl
        << this->info() << endl;

    boundaryField_ == tgf().boundaryField_;

    tgf.clear();
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::tmp<Foam::GeometricField<Type, PatchField, GeoMesh>>
Foam::GeometricField<Type, PatchField, GeoMesh>::clone() const
{
    return tmp<GeometricField<Type, PatchField, GeoMesh>>::New(*this);
}


// * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::~GeometricField()
{
    deleteDemandDrivenData(field0Ptr_);
    deleteDemandDrivenData(fieldPrevIterPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
typename
Foam::GeometricField<Type, PatchField, GeoMesh>::Internal&
Foam::GeometricField<Type, PatchField, GeoMesh>::ref
(
    const bool updateAccessTime
)
{
    if (updateAccessTime)
    {
        this->setUpToDate();
        storeOldTimes();
    }
    return *this;
}


template<class Type, template<class> class PatchField, class GeoMesh>
typename
Foam::GeometricField<Type, PatchField, GeoMesh>::Internal::FieldType&
Foam::GeometricField<Type, PatchField, GeoMesh>::primitiveFieldRef
(
    const bool updateAccessTime
)
{
    if (updateAccessTime)
    {
        this->setUpToDate();
        storeOldTimes();
    }
    return *this;
}


template<class Type, template<class> class PatchField, class GeoMesh>
typename
Foam::GeometricField<Type, PatchField, GeoMesh>::Boundary&
Foam::GeometricField<Type, PatchField, GeoMesh>::boundaryFieldRef
(
    const bool updateAccessTime
)
{
    if (updateAccessTime)
    {
        this->setUpToDate();
        storeOldTimes();
    }
    return boundaryField_;
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::storeOldTimes() const
{
    if
    (
        field0Ptr_
     && timeIndex_ != this->time().timeIndex()
     && !this->name().ends_with("_0")
    )
    {
        storeOldTime();
        timeIndex_ = this->time().timeIndex();
    }

    // Correct time index
    //timeIndex_ = this->time().timeIndex();
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::storeOldTime() const
{
    if (field0Ptr_)
    {
        field0Ptr_->storeOldTime();

        DebugInFunction
            << "Storing old time field for field" << nl << this->info() << endl;

        *field0Ptr_ == *this;
        field0Ptr_->timeIndex_ = timeIndex_;

        if (field0Ptr_->field0Ptr_)
        {
            field0Ptr_->writeOpt(this->writeOpt());
        }
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::label Foam::GeometricField<Type, PatchField, GeoMesh>::nOldTimes() const
{
    if (field0Ptr_)
    {
        return field0Ptr_->nOldTimes() + 1;
    }

    return 0;
}


template<class Type, template<class> class PatchField, class GeoMesh>
const Foam::GeometricField<Type, PatchField, GeoMesh>&
Foam::GeometricField<Type, PatchField, GeoMesh>::oldTime() const
{
    if (!field0Ptr_)
    {
        field0Ptr_ = new GeometricField<Type, PatchField, GeoMesh>
        (
            IOobject
            (
                this->name() + "_0",
                this->time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                this->registerObject()
            ),
            *this
        );

        if (debug)
        {
            InfoInFunction
                << "created old time field " << field0Ptr_->info() << endl;

            if (debug&2)
            {
                error::printStack(Info);
            }
        }
    }
    else
    {
        storeOldTimes();
    }

    return *field0Ptr_;
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>&
Foam::GeometricField<Type, PatchField, GeoMesh>::oldTime()
{
    static_cast<const GeometricField<Type, PatchField, GeoMesh>&>(*this)
        .oldTime();

    return *field0Ptr_;
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::storePrevIter() const
{
    if (!fieldPrevIterPtr_)
    {
        DebugInFunction
            << "Allocating previous iteration field" << nl
            << this->info() << endl;

        fieldPrevIterPtr_ = new GeometricField<Type, PatchField, GeoMesh>
        (
            this->name() + "PrevIter",
            *this
        );
    }
    else
    {
        *fieldPrevIterPtr_ == *this;
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
const Foam::GeometricField<Type, PatchField, GeoMesh>&
Foam::GeometricField<Type, PatchField, GeoMesh>::prevIter() const
{
    if (!fieldPrevIterPtr_)
    {
        FatalErrorInFunction
            << "previous iteration field" << endl << this->info() << endl
            << "  not stored."
            << "  Use field.storePrevIter() at start of iteration."
            << abort(FatalError);
    }

    return *fieldPrevIterPtr_;
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::
correctBoundaryConditions()
{
    this->setUpToDate();
    storeOldTimes();
    boundaryField_.evaluate();
}


template<class Type, template<class> class PatchField, class GeoMesh>
bool Foam::GeometricField<Type, PatchField, GeoMesh>::needReference() const
{
    // Search all boundary conditions, if any are
    // fixed-value or mixed (Robin) do not set reference level for solution.

    bool needRef = true;

    forAll(boundaryField_, patchi)
    {
        if (boundaryField_[patchi].fixesValue())
        {
            needRef = false;
            break;
        }
    }

    reduce(needRef, andOp<bool>());

    return needRef;
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::relax(const scalar alpha)
{
    DebugInFunction
        << "Relaxing" << nl << this->info() << " by " << alpha << endl;

    operator==(prevIter() + alpha*(*this - prevIter()));
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::relax()
{
    word name = this->name();

    if
    (
        this->mesh().data::template getOrDefault<bool>
        (
            "finalIteration",
            false
        )
    )
    {
        name += "Final";
    }

    if (this->mesh().relaxField(name))
    {
        relax(this->mesh().fieldRelaxationFactor(name));
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::word Foam::GeometricField<Type, PatchField, GeoMesh>::select
(
    bool final
) const
{
    if (final)
    {
        return this->name() + "Final";
    }

    return this->name();
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::writeMinMax
(
    Ostream& os
) const
{
    MinMax<Type> range = Foam::minMax(*this).value();

    os  << "min/max(" << this->name() << ") = "
        << range.min() << ", " << range.max() << endl;
}


template<class Type, template<class> class PatchField, class GeoMesh>
bool Foam::GeometricField<Type, PatchField, GeoMesh>::
writeData(Ostream& os) const
{
    os << *this;
    return os.good();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
Foam::tmp<Foam::GeometricField<Type, PatchField, GeoMesh>>
Foam::GeometricField<Type, PatchField, GeoMesh>::T() const
{
    auto tresult = tmp<GeometricField<Type, PatchField, GeoMesh>>::New
    (
        IOobject
        (
            this->name() + ".T()",
            this->instance(),
            this->db()
        ),
        this->mesh(),
        this->dimensions()
    );

    Foam::T(tresult.ref().primitiveFieldRef(), primitiveField());
    Foam::T(tresult.ref().boundaryFieldRef(), boundaryField());

    return tresult;
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::tmp
<
    Foam::GeometricField
    <
        typename Foam::GeometricField<Type, PatchField, GeoMesh>::cmptType,
        PatchField,
        GeoMesh
    >
>
Foam::GeometricField<Type, PatchField, GeoMesh>::component
(
    const direction d
) const
{
    auto tresult = tmp<GeometricField<cmptType, PatchField, GeoMesh>>::New
    (
        IOobject
        (
            this->name() + ".component(" + Foam::name(d) + ')',
            this->instance(),
            this->db()
        ),
        this->mesh(),
        this->dimensions()
    );

    Foam::component(tresult.ref().primitiveFieldRef(), primitiveField(), d);
    Foam::component(tresult.ref().boundaryFieldRef(), boundaryField(), d);

    return tresult;
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::replace
(
    const direction d,
    const GeometricField
    <
        typename GeometricField<Type, PatchField, GeoMesh>::cmptType,
        PatchField,
        GeoMesh
     >& gcf
)
{
    primitiveFieldRef().replace(d, gcf.primitiveField());
    boundaryFieldRef().replace(d, gcf.boundaryField());
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::replace
(
    const direction d,
    const dimensioned<cmptType>& ds
)
{
    primitiveFieldRef().replace(d, ds.value());
    boundaryFieldRef().replace(d, ds.value());
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::min
(
    const dimensioned<Type>& dt
)
{
    Foam::min(primitiveFieldRef(), primitiveField(), dt.value());
    Foam::min(boundaryFieldRef(), boundaryField(), dt.value());
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::max
(
    const dimensioned<Type>& dt
)
{
    Foam::max(primitiveFieldRef(), primitiveField(), dt.value());
    Foam::max(boundaryFieldRef(), boundaryField(), dt.value());
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::clip
(
    const dimensioned<MinMax<Type>>& range
)
{
    Foam::clip(primitiveFieldRef(), primitiveField(), range.value());
    Foam::clip(boundaryFieldRef(), boundaryField(), range.value());
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::clip
(
    const dimensioned<Type>& minVal,
    const dimensioned<Type>& maxVal
)
{
    MinMax<Type> range(minVal.value(), maxVal.value());

    Foam::clip(primitiveFieldRef(), primitiveField(), range);
    Foam::clip(boundaryFieldRef(), boundaryField(), range);
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::maxMin
(
    const dimensioned<Type>& minVal,
    const dimensioned<Type>& maxVal
)
{
    this->clip(minVal, maxVal);
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::negate()
{
    primitiveFieldRef().negate();
    boundaryFieldRef().negate();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::operator=
(
    const GeometricField<Type, PatchField, GeoMesh>& gf
)
{
    if (this == &gf)
    {
        return;  // Self-assignment is a no-op
    }

    checkField(*this, gf, "=");

    // Only assign field contents not ID

    ref() = gf();
    boundaryFieldRef() = gf.boundaryField();
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::operator=
(
    const tmp<GeometricField<Type, PatchField, GeoMesh>>& tgf
)
{
    const auto& gf = tgf();

    if (this == &gf)
    {
        return;  // Self-assignment is a no-op
    }

    checkField(*this, gf, "=");

    // Only assign field contents not ID

    this->dimensions() = gf.dimensions();
    this->oriented() = gf.oriented();

    if (tgf.movable())
    {
        // Transfer storage from the tmp
        primitiveFieldRef().transfer(tgf.constCast().primitiveFieldRef());
    }
    else
    {
        primitiveFieldRef() = gf.primitiveField();
    }

    boundaryFieldRef() = gf.boundaryField();

    tgf.clear();
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::operator=
(
    const dimensioned<Type>& dt
)
{
    ref() = dt;
    boundaryFieldRef() = dt.value();
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::operator==
(
    const tmp<GeometricField<Type, PatchField, GeoMesh>>& tgf
)
{
    const auto& gf = tgf();

    checkField(*this, gf, "==");

    // Only assign field contents not ID

    ref() = gf();
    boundaryFieldRef() == gf.boundaryField();

    tgf.clear();
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::operator==
(
    const dimensioned<Type>& dt
)
{
    ref() = dt;
    boundaryFieldRef() == dt.value();
}


#define COMPUTED_ASSIGNMENT(TYPE, op)                                          \
                                                                               \
template<class Type, template<class> class PatchField, class GeoMesh>          \
void Foam::GeometricField<Type, PatchField, GeoMesh>::operator op              \
(                                                                              \
    const GeometricField<TYPE, PatchField, GeoMesh>& gf                        \
)                                                                              \
{                                                                              \
    checkField(*this, gf, #op);                                                \
                                                                               \
    ref() op gf();                                                             \
    boundaryFieldRef() op gf.boundaryField();                                  \
}                                                                              \
                                                                               \
template<class Type, template<class> class PatchField, class GeoMesh>          \
void Foam::GeometricField<Type, PatchField, GeoMesh>::operator op              \
(                                                                              \
    const tmp<GeometricField<TYPE, PatchField, GeoMesh>>& tgf                  \
)                                                                              \
{                                                                              \
    operator op(tgf());                                                        \
    tgf.clear();                                                               \
}                                                                              \
                                                                               \
template<class Type, template<class> class PatchField, class GeoMesh>          \
void Foam::GeometricField<Type, PatchField, GeoMesh>::operator op              \
(                                                                              \
    const dimensioned<TYPE>& dt                                                \
)                                                                              \
{                                                                              \
    ref() op dt;                                                               \
    boundaryFieldRef() op dt.value();                                          \
}

COMPUTED_ASSIGNMENT(Type, +=)
COMPUTED_ASSIGNMENT(Type, -=)
COMPUTED_ASSIGNMENT(scalar, *=)
COMPUTED_ASSIGNMENT(scalar, /=)

#undef COMPUTED_ASSIGNMENT


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const GeometricField<Type, PatchField, GeoMesh>& gf
)
{
    gf().writeData(os, "internalField");
    os  << nl;
    gf.boundaryField().writeEntry("boundaryField", os);

    os.check(FUNCTION_NAME);
    return os;
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const tmp<GeometricField<Type, PatchField, GeoMesh>>& tgf
)
{
    os << tgf();
    tgf.clear();
    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#undef checkField

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "GeometricFieldNew.C"
#include "GeometricBoundaryField.C"
#include "GeometricFieldFunctions.C"

// ************************************************************************* //

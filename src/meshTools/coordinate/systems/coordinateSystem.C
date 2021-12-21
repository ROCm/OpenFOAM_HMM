/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2021 OpenCFD Ltd.
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

#include "coordinateSystem.H"
#include "cartesianCS.H"
#include "IOstream.H"
#include "axesRotation.H"
#include "identityRotation.H"
#include "transform.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(coordinateSystem, 0);
    defineRunTimeSelectionTable(coordinateSystem, dictionary);
    defineRunTimeSelectionTable(coordinateSystem, registry);
}

Foam::coordinateSystem Foam::coordinateSystem::dummy_(nullptr);


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{
    //- Is it cartesian?
    //  For output, can treat the base class as Cartesian too,
    //  since it defaults to cartesian on input.
    static inline bool isCartesian(const word& modelType)
    {
        return
        (
            modelType == coordinateSystem::typeName_()
         || modelType == coordSystem::cartesian::typeName_()
        );
    }

} // End namespace Foam


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::coordinateSystem::assign(const dictionary& dict)
{
    dict.readEntry("origin", origin_);

    note_.clear();
    dict.readIfPresent("note", note_);

    // Non-recursive, no pattern search for "rotation"
    // or "coordinateRotation" (older) sub-dictionary.
    // Don't warn about older naming for now (OCT-2018)

    const auto finder = dict.csearchCompat
    (
        "rotation", {{"coordinateRotation", -1806}},
        keyType::LITERAL
    );

    if (finder.isDict())
    {
        spec_ = coordinateRotation::New(finder.dict());
    }
    else if (finder.good() && (finder->stream().peek().isWord("none")))
    {
        spec_.reset(new coordinateRotations::identity());
    }
    else
    {
        // Fall through to expecting e1/e2/e3 specification in the dictionary
        spec_.reset(new coordinateRotations::axes(dict));
    }

    rot_ = spec_->R();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coordinateSystem::coordinateSystem(std::nullptr_t)
:
    spec_(),
    origin_(Zero),
    rot_(sphericalTensor::I),
    name_(),
    note_()
{}


Foam::coordinateSystem::coordinateSystem()
:
    spec_(new coordinateRotations::identity()),
    origin_(Zero),
    rot_(sphericalTensor::I),
    name_(),
    note_()
{}


Foam::coordinateSystem::coordinateSystem(const coordinateRotation& crot)
:
    coordinateSystem(word::null, point::zero, crot)
{}


Foam::coordinateSystem::coordinateSystem(coordinateRotation&& crot)
:
    coordinateSystem(word::null, point::zero, std::move(crot))
{}


Foam::coordinateSystem::coordinateSystem(const coordinateSystem& csys)
:
    spec_(csys.spec_.clone()),
    origin_(csys.origin_),
    rot_(csys.rot_),
    name_(csys.name_),
    note_(csys.note_)
{}


Foam::coordinateSystem::coordinateSystem(coordinateSystem&& csys)
:
    spec_(std::move(csys.spec_)),
    origin_(std::move(csys.origin_)),
    rot_(std::move(csys.rot_)),
    name_(std::move(csys.name_)),
    note_(std::move(csys.note_))
{}


Foam::coordinateSystem::coordinateSystem(autoPtr<coordinateSystem>&& csys)
:
    coordinateSystem(nullptr)
{
    if (csys)
    {
        // Has valid autoPtr - move.
        coordinateSystem::operator=(std::move(*csys));
        csys.clear();
    }
    else
    {
        // No valid autoPtr - treat like identity
        spec_.reset(new coordinateRotations::identity());
    }
}


Foam::coordinateSystem::coordinateSystem
(
    const word& name,
    const coordinateSystem& csys
)
:
    spec_(csys.spec_.clone()),
    origin_(csys.origin_),
    rot_(csys.rot_),
    name_(name),
    note_(csys.note_)
{}


Foam::coordinateSystem::coordinateSystem
(
    const point& origin,
    const coordinateRotation& crot
)
:
    coordinateSystem(word::null, origin, crot)
{}


Foam::coordinateSystem::coordinateSystem
(
    const word& name,
    const point& origin,
    const coordinateRotation& crot
)
:
    spec_(crot.clone()),
    origin_(origin),
    rot_(spec_->R()),
    name_(name),
    note_()
{}


Foam::coordinateSystem::coordinateSystem
(
    const point& origin,
    const vector& axis,
    const vector& dirn
)
:
    coordinateSystem(word::null, origin, axis, dirn)
{}


Foam::coordinateSystem::coordinateSystem
(
    const word& name,
    const point& origin,
    const vector& axis,
    const vector& dirn
)
:
    spec_(new coordinateRotations::axes(axis, dirn)),
    origin_(origin),
    rot_(spec_->R()),
    name_(name),
    note_()
{}


Foam::coordinateSystem::coordinateSystem
(
    const word& name,
    const dictionary& dict
)
:
    spec_(nullptr),
    origin_(Zero),
    rot_(sphericalTensor::I),
    name_(name),
    note_()
{
    assign(dict);
}


Foam::coordinateSystem::coordinateSystem(const dictionary& dict)
:
    coordinateSystem(word::null, dict)
{}


Foam::coordinateSystem::coordinateSystem
(
    const dictionary& dict,
    const word& dictName
)
:
    coordinateSystem(nullptr)
{
    if (dictName.size())
    {
        assign(dict.subDict(dictName));
    }
    else
    {
        assign(dict);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::coordinateSystem::clear()
{
    spec_->clear();
    origin_ = Zero;
    rot_ = sphericalTensor::I;
    note_.clear();
}


Foam::tensor Foam::coordinateSystem::R(const point& global) const
{
    return rot_;
}


Foam::tmp<Foam::tensorField> Foam::coordinateSystem::R
(
    const UList<point>& global
) const
{
    return rotationsImpl(global);
}


Foam::tmp<Foam::tensorField> Foam::coordinateSystem::R
(
    const pointUIndList& global
) const
{
    return rotationsImpl(global);
}


Foam::point Foam::coordinateSystem::transformPoint
(
    const point& localCart
) const
{
    return Foam::transform(rot_, localCart) + origin_;
}


Foam::point Foam::coordinateSystem::invTransformPoint
(
    const point& global
) const
{
    return Foam::invTransform(rot_, global - origin_);
}


Foam::vector Foam::coordinateSystem::localToGlobal
(
    const vector& local,
    bool translate
) const
{
    if (translate)
    {
        return this->transform(local) + origin_;
    }

    return this->transform(local);
}


Foam::tmp<Foam::vectorField> Foam::coordinateSystem::localToGlobal
(
    const vectorField& local,
    bool translate
) const
{
    if (translate)
    {
        return this->transform(local) + origin_;
    }

    return this->transform(local);
}


Foam::vector Foam::coordinateSystem::globalToLocal
(
    const vector& global,
    bool translate
) const
{
    if (translate)
    {
        return this->invTransform(global - origin_);
    }

    return this->invTransform(global);
}


Foam::tmp<Foam::vectorField> Foam::coordinateSystem::globalToLocal
(
    const vectorField& global,
    bool translate
) const
{
    if (translate)
    {
        return this->invTransform(global - origin_);
    }

    return this->invTransform(global);
}


void Foam::coordinateSystem::rotation(autoPtr<coordinateRotation>&& crot)
{
    spec_.reset(std::move(crot));
    if (spec_)
    {
        rot_ = spec_->R();
    }
    else
    {
        rot_ = sphericalTensor::I;
    }
}


void Foam::coordinateSystem::write(Ostream& os) const
{
    if (!valid())
    {
        return;
    }

    // Suppress output of type for Cartesian
    if (!isCartesian(type()))
    {
        os << type() << ' ';
    }

    os << "origin: " << origin_ << ' ';
    spec_->write(os);
}


void Foam::coordinateSystem::writeEntry(const word& keyword, Ostream& os) const
{
    if (!valid())
    {
        return;
    }

    const bool subDict = !keyword.empty();

    if (subDict)
    {
        os.beginBlock(keyword);

        // Suppress output of type for Cartesian
        if (!isCartesian(type()))
        {
            os.writeEntry<word>("type", type());
        }

        if (note_.size())
        {
            // The 'note' is optional
            os.writeEntry("note", note_);
        }
    }

    os.writeEntry("origin", origin_);

    spec_->writeEntry("rotation", os);

    if (subDict)
    {
        os.endBlock();
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::coordinateSystem::operator=(const coordinateSystem& csys)
{
    name_ = csys.name_;
    note_ = csys.note_;
    origin_ = csys.origin_;

    // Some extra safety
    if (csys.spec_)
    {
        rotation(csys.spec_.clone());
    }
    else
    {
        spec_.reset(new coordinateRotations::identity());
        rot_ = sphericalTensor::I;
    }
}


void Foam::coordinateSystem::operator=(coordinateSystem&& csys)
{
    name_ = std::move(csys.name_);
    note_ = std::move(csys.note_);
    spec_ = std::move(csys.spec_);
    origin_ = csys.origin_;
    rot_ = csys.rot_;
}


void Foam::coordinateSystem::operator=(const autoPtr<coordinateSystem>& csys)
{
    coordinateSystem::operator=(*csys);
}


void Foam::coordinateSystem::operator=(autoPtr<coordinateSystem>&& csys)
{
    coordinateSystem::operator=(std::move(*csys));
    csys.clear();
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

bool Foam::operator!=(const coordinateSystem& a, const coordinateSystem& b)
{
    return
    (
        a.type() != b.type()
     || a.origin() != b.origin()
     || a.R() != b.R()
    );
}


Foam::Ostream& Foam::operator<<(Ostream& os, const coordinateSystem& csys)
{
    csys.write(os);
    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //

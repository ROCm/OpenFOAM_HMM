/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "fieldValue.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fieldValue, 0);

    defineTemplateTypeNameAndDebug(IOList<vector>, 0);
    defineTemplateTypeNameAndDebug(IOList<sphericalTensor>, 0);
    defineTemplateTypeNameAndDebug(IOList<symmTensor>, 0);
    defineTemplateTypeNameAndDebug(IOList<tensor>, 0);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fieldValue::updateMesh(const mapPolyMesh&)
{
    // Do nothing
}


void Foam::fieldValue::movePoints(const Field<point>&)
{
    // Do nothing
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fieldValue::fieldValue
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    name_(name),
    obr_(obr),
    active_(true),
    log_(false),
    sourceName_(dict.lookup("sourceName")),
    fields_(dict.lookup("fields")),
    valueOutput_(dict.lookup("valueOutput"))
{
    // Only active if obr is an fvMesh
    if (isA<fvMesh>(obr_))
    {
        read(dict);
    }
    else
    {
        WarningIn
        (
            "fieldValue::fieldValue"
            "("
                "const word&, "
                "const objectRegistry&, "
                "const dictionary&, "
                "const bool"
            ")"
        )   << "No fvMesh available, deactivating."
            << nl << endl;
        active_ = false;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fieldValue::~fieldValue()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::word& Foam::fieldValue::name() const
{
    return name_;
}


const Foam::objectRegistry& Foam::fieldValue::obr() const
{
    return obr_;
}


bool Foam::fieldValue::active() const
{
    return active_;
}


const Foam::Switch& Foam::fieldValue::log() const
{
    return log_;
}


const Foam::word& Foam::fieldValue::sourceName() const
{
    return sourceName_;
}


const Foam::wordList& Foam::fieldValue::fields() const
{
    return fields_;
}


const Foam::Switch& Foam::fieldValue::valueOutput() const
{
    return valueOutput_;
}


const Foam::fvMesh& Foam::fieldValue::mesh() const
{
    return refCast<const fvMesh>(obr_);
}


void Foam::fieldValue::read(const dictionary& dict)
{
    if (active_)
    {
        log_ = dict.lookupOrDefault<Switch>("log", false);
        dict.lookup("fields") >> fields_;
        dict.lookup("valueOutput") >> valueOutput_;
    }
}


void Foam::fieldValue::execute()
{
    // Do nothing
}


void Foam::fieldValue::end()
{
    // Do nothing
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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

Description

\*---------------------------------------------------------------------------*/

#include "sampledSurface.H"
#include "polyMesh.H"
#include "demandDrivenData.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sampledSurface, 0);
    defineRunTimeSelectionTable(sampledSurface, word);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void Foam::sampledSurface::clearGeom() const
{
    deleteDemandDrivenData(SfPtr_);
    deleteDemandDrivenData(magSfPtr_);
    deleteDemandDrivenData(CfPtr_);
}



void Foam::sampledSurface::makeSf() const
{
    // It is an error to recalculate if the pointer is already set
    if (SfPtr_)
    {
        FatalErrorIn("Foam::sampledSurface::makeSf()")
            << "face area vectors already exist"
            << abort(FatalError);
    }

    const faceList& sampleFaces = faces();
    SfPtr_ = new vectorField(sampleFaces.size());

    vectorField& values = *SfPtr_;

    forAll(sampleFaces, faceI)
    {
        values[faceI] = sampleFaces[faceI].normal(points());
    }
}


void Foam::sampledSurface::makeMagSf() const
{
    // It is an error to recalculate if the pointer is already set
    if (magSfPtr_)
    {
        FatalErrorIn("Foam::sampledSurface::makeMagSf()")
            << "mag face areas already exist"
            << abort(FatalError);
    }

    const faceList& sampleFaces = faces();
    magSfPtr_ = new scalarField(sampleFaces.size());

    scalarField& values = *magSfPtr_;

    forAll(sampleFaces, faceI)
    {
        values[faceI] = sampleFaces[faceI].mag(points());
    }
}


void Foam::sampledSurface::makeCf() const
{
    // It is an error to recalculate if the pointer is already set
    if (CfPtr_)
    {
        FatalErrorIn("Foam::sampledSurface::makeCf()")
            << "face centres already exist"
            << abort(FatalError);
    }

    const faceList& sampleFaces = faces();
    CfPtr_ = new vectorField(sampleFaces.size());

    vectorField& values = *CfPtr_;

    forAll(sampleFaces, faceI)
    {
        values[faceI] = sampleFaces[faceI].centre(points());
    }
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //


Foam::autoPtr<Foam::sampledSurface>
Foam::sampledSurface::New
(
    const word& sampleType,
    const polyMesh& mesh,
    const dictionary& dict
)
{
    wordConstructorTable::iterator cstrIter =
        wordConstructorTablePtr_
            ->find(sampleType);

    if (cstrIter == wordConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "sampledSurface::New(const word&, "
            "const polyMesh&, const dictionary&)"
        )   << "Unknown sample type " << sampleType
            << endl << endl
            << "Valid sample types : " << endl
            << wordConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<sampledSurface>
    (
        cstrIter()
        (
            mesh,
            dict
        )
    );
}


bool Foam::sampledSurface::getBool
(
    const dictionary& dict,
    const word& key,
    const bool defaultVal
)
{
    if (dict.found(key))
    {
        return readBool(dict.lookup(key));
    }
    else
    {
        return defaultVal;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from mesh, name
Foam::sampledSurface::sampledSurface
(
    const polyMesh& mesh,
    const word& name,
    const bool triangulate
)
:
    mesh_(mesh),
    name_(name),
    triangulate_(triangulate),
    interpolate_(false),
    SfPtr_(NULL),
    magSfPtr_(NULL),
    CfPtr_(NULL)
{}


// Construct from dictionary
Foam::sampledSurface::sampledSurface
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    name_(type()),
    triangulate_(getBool(dict, "triangulate", true)),
    interpolate_(getBool(dict, "interpolate", false)),
    SfPtr_(NULL),
    magSfPtr_(NULL),
    CfPtr_(NULL)
{
    if (dict.found("name"))
    {
        dict.lookup("name") >> name_;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sampledSurface::~sampledSurface()
{
    clearGeom();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::vectorField& Foam::sampledSurface::Sf() const
{
    if (!SfPtr_)
    {
        makeSf();
    }

    return *SfPtr_;
}

const Foam::scalarField& Foam::sampledSurface::magSf() const
{
    if (!magSfPtr_)
    {
        makeMagSf();
    }

    return *magSfPtr_;
}

const Foam::vectorField& Foam::sampledSurface::Cf() const
{
    if (!CfPtr_)
    {
        makeCf();
    }

    return *CfPtr_;
}


// do not project scalar - just copy values
Foam::tmp<Foam::Field<Foam::scalar> >
Foam::sampledSurface::project(const Field<scalar>& field) const
{
    tmp<Field<scalar> > tRes(new Field<scalar>(faces().size()));
    Field<scalar>& res = tRes();

    forAll(faces(), faceI)
    {
        res[faceI] = field[faceI];
    }

    return tRes;
}

Foam::tmp<Foam::Field<Foam::scalar> >
Foam::sampledSurface::project(const Field<vector>& field) const
{
    tmp<Field<scalar> > tRes(new Field<scalar>(faces().size()));
    project(tRes(), field);
    return tRes;
}

Foam::tmp<Foam::Field<Foam::vector> >
Foam::sampledSurface::project(const Field<sphericalTensor>& field) const
{
    tmp<Field<vector> > tRes(new Field<vector>(faces().size()));
    project(tRes(), field);
    return tRes;
}


Foam::tmp<Foam::Field<Foam::vector> >
Foam::sampledSurface::project(const Field<symmTensor>& field) const
{
    tmp<Field<vector> > tRes(new Field<vector>(faces().size()));
    project(tRes(), field);
    return tRes;
}


Foam::tmp<Foam::Field<Foam::vector> >
Foam::sampledSurface::project(const Field<tensor>& field) const
{
    tmp<Field<vector> > tRes(new Field<vector>(faces().size()));
    project(tRes(), field);
    return tRes;
}



void Foam::sampledSurface::print(Ostream& os) const
{
    os  << type();
}

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream &os, const sampledSurface& s)
{
    s.print(os);
    os.check("Ostream& operator<<(Ostream&, const sampledSurface&");
    return os;
}

// ************************************************************************* //

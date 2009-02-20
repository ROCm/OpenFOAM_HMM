/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

#include "timeActivatedExplicitSource.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<>
const char* Foam::NamedEnum
<
    Foam::timeActivatedExplicitSource::volumeType,
    2
>::names[] =
{
    "specific",
    "absolute"
};

const Foam::NamedEnum<Foam::timeActivatedExplicitSource::volumeType, 2>
Foam::timeActivatedExplicitSource::volumeTypeNames_;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::timeActivatedExplicitSource::updateCellSet()
{
    cellSelector_->applyToSet(topoSetSource::NEW, selectedCellSet_);

    Info<< "    " << sourceName_ << ": selected "
        << returnReduce(selectedCellSet_.size(), sumOp<label>())
        << " cells" << nl << endl;

    V_ = scalarField(selectedCellSet_.size(), 1.0);
    if (volumeType_ == vtAbsolute)
    {
        label i = 0;
        forAllConstIter(cellSet, selectedCellSet_, iter)
        {
            V_[i++] = mesh_.V()[iter.key()];
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::timeActivatedExplicitSource::timeActivatedExplicitSource
(
    const word& sourceName,
    const fvMesh& mesh,
    const dimensionSet& dims
)
:
    IOdictionary
    (
        IOobject
        (
            sourceName + "Properties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    sourceName_(sourceName),
    mesh_(mesh),
    runTime_(mesh.time()),
    dimensions_(dims),
    volumeType_(volumeTypeNames_.read(lookup("volumeType"))),
    timeStart_(readScalar(lookup("timeStart"))),
    duration_(readScalar(lookup("duration"))),
    onValue_(readScalar(lookup("onValue"))),
    offValue_(readScalar(lookup("offValue"))),
    currentValue_(0.0),
    V_(0),
    cellSource_(lookup("cellSource")),
    cellSelector_
    (
        topoSetSource::New
        (
            cellSource_,
            mesh,
            subDict(cellSource_ + "Coeffs")
        )
    ),
    selectedCellSet_
    (
        mesh,
        sourceName + "SourceCellSet",
        mesh.nCells()/10 + 1  // Reasonable size estimate.
    )
{
    // Create the cell set
    updateCellSet();

    // Initialise the value
    update();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::timeActivatedExplicitSource::timeStart() const
{
    return timeStart_;
}


Foam::scalar Foam::timeActivatedExplicitSource::duration() const
{
    return duration_;
}


const Foam::dimensionedScalar
Foam::timeActivatedExplicitSource::currentValue() const
{
    return dimensionedScalar
    (
        sourceName_,
        dimensions_,
        currentValue_
    );
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh> >
Foam::timeActivatedExplicitSource::Su() const
{
    tmp<DimensionedField<scalar, volMesh> > tSource
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                sourceName_ + "Su",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimensions_, 0.0)
        )
    );

    DimensionedField<scalar, volMesh>& sourceField = tSource();

    label i = 0;
    forAllConstIter(cellSet, selectedCellSet_, iter)
    {
        sourceField[iter.key()] = currentValue_/V_[i++];
    }

    return tSource;
}


bool Foam::timeActivatedExplicitSource::read()
{
    if (regIOobject::read())
    {
        volumeType_ = volumeTypeNames_.read(lookup("volumeType"));
        lookup("timeStart") >> duration_;
        lookup("duration") >> duration_;
        lookup("onValue") >> onValue_;
        lookup("offValue") >> offValue_;
        lookup("cellSource") >> cellSource_;
        cellSelector_ =
            topoSetSource::New
            (
                cellSource_,
                mesh_,
                subDict(cellSource_ + "Coeffs")
            );
        updateCellSet();

        return true;
    }
    else
    {
        return false;
    }
}


void Foam::timeActivatedExplicitSource::update()
{
    // Set the source value
    if
    (
        (runTime_.time().value() >= timeStart_)
     && (runTime_.time().value() <= timeStart_ + duration_)
    )
    {
        currentValue_ = onValue_;
    }
    else
    {
        currentValue_ = offValue_;
    }

    // Update the cell set if the mesh is changing
    if (mesh_.changing())
    {
        updateCellSet();
    }
}


// ************************************************************************* //

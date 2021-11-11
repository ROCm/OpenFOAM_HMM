/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2015-2020 OpenCFD Ltd.
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

#include "sampledSets.H"
#include "dictionary.H"
#include "Time.H"
#include "volFields.H"
#include "volPointInterpolation.H"
#include "mapPolyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sampledSets, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        sampledSets,
        dictionary
    );
}

bool Foam::sampledSets::verbose_ = false;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::sampledSets::combineSampledSets
(
    PtrList<coordSet>& masterSampledSets,
    labelListList& indexSets
)
{
    // Combine sampleSets from processors. Sort by curveDist. Return
    // ordering in indexSets.
    // Note: only master results are valid

    masterSampledSets_.clear();
    masterSampledSets_.setSize(size());
    indexSets_.setSize(size());

    const PtrList<sampledSet>& sampledSets = *this;

    forAll(sampledSets, setI)
    {
        labelList segments;
        masterSampledSets.set
        (
            setI,
            sampledSets[setI].gather(indexSets[setI], segments)
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sampledSets::sampledSets
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    functionObjects::regionFunctionObject(name, runTime, dict),
    PtrList<sampledSet>(),
    mesh_(refCast<const fvMesh>(obr_)),
    loadFromFiles_(false),
    outputPath_(fileName::null),
    searchEngine_(mesh_),
    interpolationScheme_(word::null),
    writeFormat_(word::null),
    writeFormatOptions_(dict.subOrEmptyDict("formatOptions"))
{
    outputPath_ =
    (
        mesh_.time().globalPath()/functionObject::outputPrefix/name
    );

    if (mesh_.name() != polyMesh::defaultRegion)
    {
        outputPath_ /= mesh_.name();
    }

    outputPath_.clean();  // Remove unneeded ".."

    read(dict);
}


Foam::sampledSets::sampledSets
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    functionObjects::regionFunctionObject(name, obr, dict),
    PtrList<sampledSet>(),
    mesh_(refCast<const fvMesh>(obr)),
    loadFromFiles_(loadFromFiles),
    outputPath_(fileName::null),
    searchEngine_(mesh_),
    interpolationScheme_(word::null),
    writeFormat_(word::null),
    writeFormatOptions_(dict.subOrEmptyDict("formatOptions"))
{
    outputPath_ =
    (
        mesh_.time().globalPath()/functionObject::outputPrefix/name
    );

    if (mesh_.name() != polyMesh::defaultRegion)
    {
        outputPath_ /= mesh_.name();
    }

    outputPath_.clean();  // Remove unneeded ".."

    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::sampledSets::verbose(const bool on)
{
    bool old(verbose_);
    verbose_ = on;
    return old;
}


bool Foam::sampledSets::execute()
{
    return true;
}


bool Foam::sampledSets::write()
{
    if (size())
    {
        const label nFields = classifyFields();

        if (Pstream::master())
        {
            if (debug)
            {
                Pout<< "timeName = " << mesh_.time().timeName() << nl
                    << "scalarFields    " << scalarFields_ << nl
                    << "vectorFields    " << vectorFields_ << nl
                    << "sphTensorFields " << sphericalTensorFields_ << nl
                    << "symTensorFields " << symmTensorFields_ <<nl
                    << "tensorFields    " << tensorFields_ <<nl;
            }

            if (nFields)
            {
                if (debug)
                {
                    Pout<< "Creating directory "
                        << outputPath_/mesh_.time().timeName()
                        << nl << endl;
                }

                mkDir(outputPath_/mesh_.time().timeName());
            }
            else
            {
                Info<< "No fields to sample" << endl;
            }
        }

        if (nFields)
        {
            sampleAndWrite(scalarFields_);
            sampleAndWrite(vectorFields_);
            sampleAndWrite(sphericalTensorFields_);
            sampleAndWrite(symmTensorFields_);
            sampleAndWrite(tensorFields_);
        }
    }

    return true;
}


bool Foam::sampledSets::read(const dictionary& dict)
{
    dict_ = dict;

    if (dict_.found("sets"))
    {
        dict_.readEntry("fields", fieldSelection_);
        clearFieldGroups();

        dict.readEntry("interpolationScheme", interpolationScheme_);
        dict.readEntry("setFormat", writeFormat_);

        PtrList<sampledSet> newList
        (
            dict_.lookup("sets"),
            sampledSet::iNew(mesh_, searchEngine_)
        );
        transfer(newList);
        combineSampledSets(masterSampledSets_, indexSets_);

        if (this->size())
        {
            Info<< "Reading set description:" << nl;
            forAll(*this, setI)
            {
                Info<< "    " << operator[](setI).name() << nl;
            }
            Info<< endl;
        }
    }

    if (Pstream::master() && debug)
    {
        Pout<< "sample fields:" << fieldSelection_ << nl
            << "sample sets:" << nl << "(" << nl;

        forAll(*this, setI)
        {
            Pout<< "  " << operator[](setI) << endl;
        }
        Pout<< ")" << endl;
    }

    return true;
}


void Foam::sampledSets::correct()
{
    if (dict_.found("sets"))
    {
        searchEngine_.correct();

        PtrList<sampledSet> newList
        (
            dict_.lookup("sets"),
            sampledSet::iNew(mesh_, searchEngine_)
        );
        transfer(newList);
        combineSampledSets(masterSampledSets_, indexSets_);
    }
}


void Foam::sampledSets::updateMesh(const mapPolyMesh& mpm)
{
    if (&mpm.mesh() == &mesh_)
    {
        correct();
    }
}


void Foam::sampledSets::movePoints(const polyMesh& mesh)
{
    if (&mesh == &mesh_)
    {
        correct();
    }
}


void Foam::sampledSets::readUpdate(const polyMesh::readUpdateState state)
{
    if (state != polyMesh::UNCHANGED)
    {
        correct();
    }
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
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

#include "sampledIsoSurface.H"
#include "dictionary.H"
#include "volFields.H"
#include "volPointInterpolation.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sampledIsoSurface, 0);
    addNamedToRunTimeSelectionTable(sampledSurface, sampledIsoSurface, word, isoSurface);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::sampledIsoSurface::getIsoFields() const
{
    const fvMesh& fvm = static_cast<const fvMesh&>(mesh());

    word pointFldName = "volPointInterpolate(" + isoField_ + ')';

    // Get volField
    // ~~~~~~~~~~~~

    if (fvm.foundObject<volScalarField>(isoField_))
    {
        if (debug)
        {
            Info<< "sampledIsoSurface::getIsoField() : lookup "
                << isoField_ << endl;
        }
        storedVolFieldPtr_.clear();
        volFieldPtr_ = &fvm.lookupObject<volScalarField>(isoField_);
    }
    else
    {
        // Bit of a hack. Read field and store.

        if (debug)
        {
            Info<< "sampledIsoSurface::getIsoField() : checking "
                << isoField_ << " for same time " << fvm.time().timeName()
                << endl;
        }

        if
        (
           !storedVolFieldPtr_.valid()
         || (fvm.time().timeName() != storedVolFieldPtr_().instance())
        )
        {
            if (debug)
            {
                Info<< "sampledIsoSurface::getIsoField() : reading "
                    << isoField_ << " from time " << fvm.time().timeName()
                    << endl;
            }

            storedVolFieldPtr_.reset
            (
                new volScalarField
                (
                    IOobject
                    (
                        isoField_,
                        fvm.time().timeName(),
                        fvm,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE,
                        false
                    ),
                    fvm
                )
            );
            volFieldPtr_ = storedVolFieldPtr_.operator->();

            // Interpolate to get pointField

            if (debug)
            {
                Info<< "sampledIsoSurface::getIsoField() : interpolating "
                    << pointFldName << endl;
            }

            storedPointFieldPtr_.reset
            (
                volPointInterpolation::New(fvm).interpolate(*volFieldPtr_).ptr()
            );
            pointFieldPtr_ = storedPointFieldPtr_.operator->();
        }
    }



    // Get pointField
    // ~~~~~~~~~~~~~~

    if (fvm.foundObject<pointScalarField>(pointFldName))
    {
        if (debug)
        {
            Info<< "sampledIsoSurface::getIsoField() : lookup "
                << pointFldName << endl;
        }
        pointFieldPtr_ = &fvm.lookupObject<pointScalarField>(pointFldName);
    }
    else
    {
        // Not in registry. Interpolate.

        if (debug)
        {
            Info<< "sampledIsoSurface::getIsoField() : checking interpolate "
                << isoField_ << " for same time " << fvm.time().timeName()
                << endl;
        }

        if
        (
           !storedPointFieldPtr_.valid()
         || (fvm.time().timeName() != storedPointFieldPtr_().instance())
        )
        {
            if (debug)
            {
                Info<< "sampledIsoSurface::getIsoField() : interpolating "
                    << pointFldName << endl;
            }

            storedPointFieldPtr_.reset
            (
                volPointInterpolation::New(fvm).interpolate(*volFieldPtr_).ptr()
            );
            pointFieldPtr_ = storedPointFieldPtr_.operator->();
        }
    }

    if (debug)
    {
        Info<< "sampledIsoSurface::getIsoField() : volField "
            << volFieldPtr_->name() << " min:" << min(*volFieldPtr_).value()
            << " max:" << max(*volFieldPtr_).value() << endl;
        Info<< "sampledIsoSurface::getIsoField() : pointField "
            << pointFieldPtr_->name()
            << " min:" << gMin(pointFieldPtr_->internalField())
            << " max:" << gMax(pointFieldPtr_->internalField()) << endl;
    }
}


void Foam::sampledIsoSurface::createGeometry() const
{
    const fvMesh& fvm = static_cast<const fvMesh&>(mesh());

    if (fvm.time().timeIndex() != storedTimeIndex_)
    {
        storedTimeIndex_ = fvm.time().timeIndex();

        getIsoFields();

        // Clear any stored topo
        surfPtr_.clear();
        facesPtr_.clear();

        if (average_)
        {
            //- From point field and interpolated cell.
            volScalarField cellAvg
            (
                IOobject
                (
                    "cellAvg",
                    fvm.time().timeName(),
                    fvm.time(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                fvm,
                dimensionedScalar("zero", dimless, scalar(0.0))
            );
            labelField nPointCells(fvm.nCells(), 0);
            {
                for (label pointI = 0; pointI < fvm.nPoints(); pointI++)
                {
                    const labelList& pCells = fvm.pointCells(pointI);

                    forAll(pCells, i)
                    {
                        label cellI = pCells[i];

                        cellAvg[cellI] += (*pointFieldPtr_)[pointI];
                        nPointCells[cellI]++;
                    }
                }
            }
            forAll(cellAvg, cellI)
            {
                cellAvg[cellI] /= nPointCells[cellI];
            }
            // Give value to calculatedFvPatchFields
            cellAvg.correctBoundaryConditions();

            surfPtr_.reset
            (
                new isoSurface
                (
                    cellAvg,
                    *pointFieldPtr_,
                    isoVal_,
                    regularise_
                )
            );
        }
        else
        {
            //- Direct from cell field and point field.
            surfPtr_.reset
            (
                new isoSurface
                (
                    *volFieldPtr_,
                    *pointFieldPtr_,
                    isoVal_,
                    regularise_
                )
            );
        }


        if (debug)
        {
            Pout<< "sampledIsoSurface::createGeometry() : constructed iso:"
                << nl
                << "    regularise     : " << regularise_ << nl
                << "    average        : " << average_ << nl
                << "    isoField       : " << isoField_ << nl
                << "    isoValue       : " << isoVal_ << nl
                << "    points         : " << points().size() << nl
                << "    tris           : " << surface().size() << nl
                << "    cut cells      : " << surface().meshCells().size()
                << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sampledIsoSurface::sampledIsoSurface
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    sampledSurface(name, mesh, dict),
    isoField_(dict.lookup("isoField")),
    isoVal_(readScalar(dict.lookup("isoValue"))),
    regularise_(dict.lookupOrDefault("regularise", true)),
    average_(dict.lookupOrDefault("average", false)),
    zoneName_(word::null),
    surfPtr_(NULL),
    facesPtr_(NULL),
    storedTimeIndex_(-1),
    storedVolFieldPtr_(NULL),
    volFieldPtr_(NULL),
    storedPointFieldPtr_(NULL),
    pointFieldPtr_(NULL)
{
    if (!sampledSurface::interpolate())
    {
        FatalErrorIn
        (
            "sampledIsoSurface::sampledIsoSurface"
            "(const word&, const polyMesh&, const dictionary&)"
        )   << "Non-interpolated iso surface not supported since triangles"
            << " span across cells." << exit(FatalError);
    }

//    label zoneId = -1;
//    if (dict.readIfPresent("zone", zoneName_))
//    {
//        zoneId = mesh.cellZones().findZoneID(zoneName_);
//        if (debug && zoneId < 0)
//        {
//            Info<< "cellZone \"" << zoneName_
//                << "\" not found - using entire mesh"
//                << endl;
//        }
//    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sampledIsoSurface::~sampledIsoSurface()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::sampledIsoSurface::correct(const bool meshChanged)
{
    // Only change of mesh changes plane - zone restriction gets lost
    if (meshChanged)
    {
        surfPtr_.clear();
        facesPtr_.clear();
    }
}


Foam::tmp<Foam::scalarField>
Foam::sampledIsoSurface::sample
(
    const volScalarField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::vectorField>
Foam::sampledIsoSurface::sample
(
    const volVectorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::sphericalTensorField>
Foam::sampledIsoSurface::sample
(
    const volSphericalTensorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::symmTensorField>
Foam::sampledIsoSurface::sample
(
    const volSymmTensorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::tensorField>
Foam::sampledIsoSurface::sample
(
    const volTensorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::scalarField>
Foam::sampledIsoSurface::interpolate
(
    const interpolation<scalar>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::vectorField>
Foam::sampledIsoSurface::interpolate
(
    const interpolation<vector>& interpolator
) const
{
    return interpolateField(interpolator);
}

Foam::tmp<Foam::sphericalTensorField>
Foam::sampledIsoSurface::interpolate
(
    const interpolation<sphericalTensor>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::symmTensorField>
Foam::sampledIsoSurface::interpolate
(
    const interpolation<symmTensor>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::tensorField>
Foam::sampledIsoSurface::interpolate
(
    const interpolation<tensor>& interpolator
) const
{
    return interpolateField(interpolator);
}


void Foam::sampledIsoSurface::print(Ostream& os) const
{
    os  << "sampledIsoSurface: " << name() << " :"
        << "  field:" << isoField_
        << "  value:" << isoVal_;
        //<< "  faces:" << faces().size()       // note: possibly no geom yet
        //<< "  points:" << points().size();
}


// ************************************************************************* //

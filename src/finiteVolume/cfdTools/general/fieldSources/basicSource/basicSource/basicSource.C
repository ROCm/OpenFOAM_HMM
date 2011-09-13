/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "basicSource.H"
#include "fvMesh.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(basicSource, 0);
    defineRunTimeSelectionTable(basicSource, dictionary);


    // * * * * * * * * * * Static Member Functions * * * * * * * * * * * * * //

    template<> const char* NamedEnum
    <
        basicSource::selectionModeType,
        4
        >::names[] =
    {
        "points",
        "cellSet",
        "cellZone",
        "all"
    };

    const NamedEnum<basicSource::selectionModeType, 4>
        basicSource::selectionModeTypeNames_;
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::basicSource::setSelection(const dictionary& dict)
{
    switch (selectionMode_)
    {
        case smPoints:
        {
            dict.lookup("points") >> points_;
            break;
        }
        case smCellSet:
        {
            dict.lookup("cellSet") >> cellSetName_;
            break;
        }
        case smCellZone:
        {
            dict.lookup("cellZone") >> cellSetName_;
            break;
        }
        case smAll:
        {
            break;
        }
        default:
        {
            FatalErrorIn
            (
                "basicSource::setSelection(const dictionary&)"
            )   << "Unknown selectionMode "
                << selectionModeTypeNames_[selectionMode_]
                << ". Valid selectionMode types are" << selectionModeTypeNames_
                << exit(FatalError);
        }
    }
}


void Foam::basicSource::setCellSet()
{
    Info<< incrIndent << indent << "Source: " << name_ << endl;
    switch (selectionMode_)
    {
        case smPoints:
        {
            Info<< indent << "- selecting cells using points" << endl;

            labelHashSet selectedCells;

            forAll(points_, i)
            {
                label cellI = mesh_.findCell(points_[i]);
                if (cellI >= 0)
                {
                    selectedCells.insert(cellI);
                }

                label globalCellI = returnReduce(cellI, maxOp<label>());
                if (globalCellI < 0)
                {
                    WarningIn("TimeActivatedExplicitSource<Type>::setCellIds()")
                        << "Unable to find owner cell for point " << points_[i]
                        << endl;
                }
            }

            cells_ = selectedCells.toc();

            break;
        }
        case smCellSet:
        {
            Info<< indent << "- selecting cells using cellSet "
                << cellSetName_ << endl;

            cellSet selectedCells(mesh_, cellSetName_);
            cells_ = selectedCells.toc();

            break;
        }
        case smCellZone:
        {
            Info<< indent << "- selecting cells using cellZone "
                << cellSetName_ << endl;
            label zoneID = mesh_.cellZones().findZoneID(cellSetName_);
            if (zoneID == -1)
            {
                FatalErrorIn("basicSource<Type>::setCellIds()")
                    << "Cannot find cellZone " << cellSetName_ << endl
                    << "Valid cellZones are " << mesh_.cellZones().names()
                    << exit(FatalError);
            }
            cells_ = mesh_.cellZones()[zoneID];

            break;
        }
        case smAll:
        {
            Info<< indent << "- selecting all cells" << endl;
            cells_ = identity(mesh_.nCells());

            break;
        }
        default:
        {
            FatalErrorIn("basicSource<Type>::setCellIds()")
                << "Unknown selectionMode "
                << selectionModeTypeNames_[selectionMode_]
                << ". Valid selectionMode types are" << selectionModeTypeNames_
                << exit(FatalError);
        }
    }

    // Set volume information
    V_ = 0.0;
    forAll(cells_, i)
    {
        V_ += mesh_.V()[cells_[i]];
    }
    reduce(V_, sumOp<scalar>());

    Info<< indent << "- selected "
        << returnReduce(cells_.size(), sumOp<label>())
        << " cell(s) with volume " << V_ << nl << decrIndent << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::basicSource::basicSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    name_(name),
    mesh_(mesh),
    dict_(dict),
    active_(readBool(dict_.lookup("active"))),
    timeStart_(readScalar(dict_.lookup("timeStart"))),
    duration_(readScalar(dict_.lookup("duration"))),
    selectionMode_
    (
        selectionModeTypeNames_.read(dict_.lookup("selectionMode"))
    ),
    cellSetName_("none"),
    V_(0.0)
{
    setSelection(dict_);

    setCellSet();
}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::basicSource> Foam::basicSource::New
(
    const word& name,
    const dictionary& dict,
    const fvMesh& mesh
)
{
    word modelType(dict.lookup("type"));

    Info<< "Selecting model type " << modelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "basicSource::New(const volVectorField&, "
            "const surfaceScalarField&, transportModel&)"
        )   << "Unknown Model type " << modelType
            << nl << nl
            << "Valid model types are :" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<basicSource>(cstrIter()(name, modelType, dict, mesh));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::basicSource::isActive()
{
    if
    (
        active_
     && (mesh_.time().value() >= timeStart_)
     && (mesh_.time().value() <= timeEnd())
    )
    {
        // Update the cell set if the mesh is changing
        if (mesh_.changing())
        {
            setCellSet();
        }
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //

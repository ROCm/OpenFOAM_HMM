/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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

#include "patchCellsSource.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "boundarySourcePatch.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(patchCellsSource, 0);
    addToRunTimeSelectionTable(option, patchCellsSource, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::patchCellsSource::patchCellsSource
(
    const word& sourceName,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fv::option(sourceName, modelType, dict, mesh),
    curTimeIndex_(-1),
    UName_(coeffs_.getOrDefault<word>("U", "none")),
    hName_(coeffs_.getOrDefault<word>("he", "none")),
    speciesName_(coeffs_.getOrDefault<word>("species", "none"))
{
    label nFields = 0;
    if (UName_ != "none")
    {
        nFields++;
    }
    if (hName_ != "none")
    {
        nFields++;
    }
    if (speciesName_ != "none")
    {
        nFields++;
    }

    if (nFields > 1)
    {
        FatalErrorInFunction
            << "patchCellsSource : "
            << " cannot be used for more than one field."
            << exit(FatalError);
    }

    fieldNames_.resize(1);
    if (speciesName_ != "none")
    {
        fieldNames_[0] = speciesName_;
    }

    if (UName_ != "none")
    {
        fieldNames_[0] = UName_;
    }

    if (hName_ != "none")
    {
        fieldNames_[0] = hName_;
    }

    fv::option::resetApplied();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::patchCellsSource::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    if (debug)
    {
        Info<< type() << ": applying source to " << eqn.psi().name() << endl;
    }

    if (curTimeIndex_ == mesh_.time().timeIndex())
    {
        return;
    }

    volScalarField* psiPtr;

    // If source applied to he, we need to loop over T for BC's
    if (hName_ != "none")
    {
        psiPtr = mesh_.getObjectPtr<volScalarField>("T");
    }
    else
    {
        auto* psi =
            mesh_.getObjectPtr<volScalarField>(eqn.psi().name());

        psiPtr = psi;
    }

    const volScalarField::Boundary& psib = psiPtr->boundaryField();

    volScalarField mDot
    (
        IOobject
        (
            "mDot",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(eqn.dimensions()/dimVolume, Zero)
    );

    forAll(psib, patchi)
    {
        if (isA<boundarySourcePatch>(psib[patchi]))
        {
            const boundarySourcePatch& pp =
                refCast<const boundarySourcePatch>(psib[patchi]);

            const labelUList& fc = mesh_.boundary()[patchi].faceCells();

            tmp<scalarField> tsb = pp.patchSource();
            const scalarField& sb = tsb.cref();

            forAll(fc, facei)
            {
                const label celli = fc[facei];
                mDot[celli] += sb[facei];
            }
        }
    }
    eqn += mDot;

    curTimeIndex_ = mesh_.time().timeIndex();

    if (debug)
    {
         Info<< " Field source rate min/max : "
             << gMin(mDot) << " / " << gMax(mDot) << endl;
    }
}


bool Foam::fv::patchCellsSource::read(const dictionary& dict)
{
    if (!fv::option::read(dict))
    {
        return false;
    }

    return true;
}


// ************************************************************************* //

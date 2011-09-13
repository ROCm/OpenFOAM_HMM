/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*----------------------------------------------------------------------------*/

#include "radialActuationDiskSource.H"
#include "geometricOneField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(radialActuationDiskSource, 0);
    addToRunTimeSelectionTable
    (
        basicSource,
        radialActuationDiskSource,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radialActuationDiskSource::radialActuationDiskSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    actuationDiskSource(name, modelType, dict, mesh),
    dict_(dict.subDict(modelType + "Coeffs")),
    coeffs_()
{
    dict_.lookup("coeffs") >> coeffs_;
    Info<< "    - creating radial actuation disk zone: "
        << this->name() << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::radialActuationDiskSource::addSu(fvMatrix<vector>& UEqn)
{
    bool compressible = false;
    if (UEqn.dimensions() == dimensionSet(1, 1, -2, 0, 0))
    {
        compressible = true;
    }

    const scalarField& cellsV = this->mesh().V();
    vectorField& Usource = UEqn.source();
    const vectorField& U = UEqn.psi();

    if (V() > VSMALL)
    {
        if (compressible)
        {
            addRadialActuationDiskAxialInertialResistance
            (
                Usource,
                cells_,
                cellsV,
                this->mesh().lookupObject<volScalarField>("rho"),
                U
            );
        }
        else
        {
            addRadialActuationDiskAxialInertialResistance
            (
                Usource,
                cells_,
                cellsV,
                geometricOneField(),
                U
            );
        }
    }
}


void Foam::radialActuationDiskSource::writeData(Ostream& os) const
{
    actuationDiskSource::writeData(os);
}


bool Foam::radialActuationDiskSource::read(const dictionary& dict)
{
    if (basicSource::read(dict))
    {
        const dictionary& sourceDict = dict.subDict(name());
        const dictionary& subDictCoeffs =
            sourceDict.subDict(typeName + "Coeffs");
        subDictCoeffs.readIfPresent("diskDir", diskDir_);
        subDictCoeffs.readIfPresent("Cp", Cp_);
        subDictCoeffs.readIfPresent("Ct", Ct_);
        subDictCoeffs.readIfPresent("diskArea", diskArea_);
        subDictCoeffs.lookup("coeffs") >> coeffs_;
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //

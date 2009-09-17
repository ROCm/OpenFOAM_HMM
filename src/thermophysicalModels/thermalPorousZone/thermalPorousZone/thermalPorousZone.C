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

\*----------------------------------------------------------------------------*/

#include "thermalPorousZone.H"
#include "fvMesh.H"
#include "fvMatrices.H"
#include "basicThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermalPorousZone::thermalPorousZone
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    porousZone(name, mesh, dict),
    T_("T", dimTemperature, -GREAT)
{
    if (const dictionary* dictPtr = dict.subDictPtr("thermalModel"))
    {
        word thermalModel(dictPtr->lookup("type"));

        if (thermalModel == "fixedTemperature")
        {
            dictPtr->lookup("T") >> T_;
        }
        else
        {
            FatalIOErrorIn
            (
                "thermalPorousZone::thermalPorousZone"
                "("
                "const word& name, "
                "const fvMesh& mesh, "
                "const dictionary& dict"
                ")",
                *dictPtr
            )   << "thermalModel " << thermalModel << " is not supported" << nl
                << "    Supported thermalModels are: fixedTemperature"
                << exit(FatalIOError);
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::thermalPorousZone::addEnthalpySource
(
    const basicThermo& thermo,
    const volScalarField& rho,
    fvScalarMatrix& hEqn
) const
{
    if (zoneId() == -1 || T_.value() < 0.0)
    {
        return;
    }

    const labelList& cells = mesh().cellZones()[zoneId()];
    const scalarField& V = mesh().V();
    scalarField& hDiag = hEqn.diag();
    scalarField& hSource = hEqn.source();

    scalarField hZone = thermo.h(scalarField(cells.size(), T_.value()), cells);
    scalar rate = 1e6;

    forAll (cells, i)
    {
        hDiag[cells[i]] += rate*V[cells[i]]*rho[cells[i]];
        hSource[cells[i]] += rate*V[cells[i]]*rho[cells[i]]*hZone[i];
    }
}


// ************************************************************************* //

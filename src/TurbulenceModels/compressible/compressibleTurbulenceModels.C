/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

#include "CompressibleTurbulenceModel.H"
#include "laminar.H"
#include "RASModel.H"
#include "kEpsilon.H"
#include "fluidThermo.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
    typedef TurbulenceModel
    <
        geometricOneField,
        volScalarField,
        compressibleTurbulenceModel,
        fluidThermo
    > baseCompressibleFluidThermoTurbulenceModel;

    defineTemplateRunTimeSelectionTable
    (
        baseCompressibleFluidThermoTurbulenceModel,
        dictionary
    );


    typedef CompressibleTurbulenceModel<fluidThermo>
        compressibleFluidThermoTurbulenceModel;

    typedef laminar<compressibleFluidThermoTurbulenceModel> compressibleLaminar;

    defineNamedTemplateTypeNameAndDebug(compressibleLaminar, 0);

    addToRunTimeSelectionTable
    (
        baseCompressibleFluidThermoTurbulenceModel,
        compressibleLaminar,
        dictionary
    );


    typedef RASModel<compressibleFluidThermoTurbulenceModel>
        compressibleRASModel;

    defineNamedTemplateTypeNameAndDebug(compressibleRASModel, 0);

    defineTemplateRunTimeSelectionTable(compressibleRASModel, dictionary);

    addToRunTimeSelectionTable
    (
        baseCompressibleFluidThermoTurbulenceModel,
        compressibleRASModel,
        dictionary
    );


    namespace RASModels
    {
        typedef kEpsilon<compressibleFluidThermoTurbulenceModel>
            compressibleKEpsilon;

        defineNamedTemplateTypeNameAndDebug(compressibleKEpsilon, 0);

        addToRunTimeSelectionTable
        (
            compressibleRASModel,
            compressibleKEpsilon,
            dictionary
        );
    }
}


// ************************************************************************* //

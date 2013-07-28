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

#include "PhaseIncompressibleTurbulenceModel.H"
#include "laminar.H"
#include "RASModel.H"
#include "kEpsilon.H"
#include "LaheyKEpsilon.H"
#include "continuousGasKEpsilon.H"
#include "kineticTheoryModel.H"
#include "phasePressureModel.H"
#include "phaseModel.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
    typedef TurbulenceModel
    <
        volScalarField,
        geometricOneField,
        incompressibleTurbulenceModel,
        phaseModel
    > basePhaseIncompressibleTransportTurbulenceModel;

    defineTemplateRunTimeSelectionTable
    (
        basePhaseIncompressibleTransportTurbulenceModel,
        dictionary
    );

    typedef PhaseIncompressibleTurbulenceModel<phaseModel>
        incompressibleTransportTurbulenceModel;

    typedef laminar<incompressibleTransportTurbulenceModel>
        incompressibleLaminar;

    defineNamedTemplateTypeNameAndDebug(incompressibleLaminar, 0);

    addToRunTimeSelectionTable
    (
        basePhaseIncompressibleTransportTurbulenceModel,
        incompressibleLaminar,
        dictionary
    );


    typedef RASModel<incompressibleTransportTurbulenceModel>
        incompressibleRASModel;

    defineNamedTemplateTypeNameAndDebug(incompressibleRASModel, 0);

    defineTemplateRunTimeSelectionTable(incompressibleRASModel, dictionary);

    addToRunTimeSelectionTable
    (
        basePhaseIncompressibleTransportTurbulenceModel,
        incompressibleRASModel,
        dictionary
    );

    namespace RASModels
    {
        typedef kEpsilon<incompressibleTransportTurbulenceModel>
            incompressiblekEpsilon;

        defineNamedTemplateTypeNameAndDebug(incompressiblekEpsilon, 0);

        addToRunTimeSelectionTable
        (
            incompressibleRASModel,
            incompressiblekEpsilon,
            dictionary
        );
    }

    namespace RASModels
    {
        typedef LaheyKEpsilon<incompressibleTransportTurbulenceModel>
            incompressibleLaheyKEpsilon;

        defineNamedTemplateTypeNameAndDebug(incompressibleLaheyKEpsilon, 0);

        addToRunTimeSelectionTable
        (
            incompressibleRASModel,
            incompressibleLaheyKEpsilon,
            dictionary
        );
    }

    namespace RASModels
    {
        typedef continuousGasKEpsilon<incompressibleTransportTurbulenceModel>
            incompressiblecontinuousGasKEpsilon;

        defineNamedTemplateTypeNameAndDebug
        (
            incompressiblecontinuousGasKEpsilon,
            0
        );

        addToRunTimeSelectionTable
        (
            incompressibleRASModel,
            incompressiblecontinuousGasKEpsilon,
            dictionary
        );
    }
}


namespace Foam
{
    typedef PhaseIncompressibleTurbulenceModel<phaseModel>
        incompressibleTransportTurbulenceModel;

    typedef RASModel<incompressibleTransportTurbulenceModel>
        incompressibleRASModel;

    defineTypeNameAndDebug(kineticTheoryModel, 0);

    addToRunTimeSelectionTable
    (
        incompressibleRASModel,
        kineticTheoryModel,
        dictionary
    );
}


namespace Foam
{
    typedef PhaseIncompressibleTurbulenceModel<phaseModel>
        incompressibleTransportTurbulenceModel;

    typedef RASModel<incompressibleTransportTurbulenceModel>
        incompressibleRASModel;

    defineTypeNameAndDebug(phasePressureModel, 0);

    addToRunTimeSelectionTable
    (
        incompressibleRASModel,
        phasePressureModel,
        dictionary
    );
}


// ************************************************************************* //

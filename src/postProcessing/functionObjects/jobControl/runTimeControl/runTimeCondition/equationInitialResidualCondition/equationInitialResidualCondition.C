/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015 OpenCFD Ltd.
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

#include "equationInitialResidualCondition.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(equationInitialResidualCondition, 0);
    addToRunTimeSelectionTable
    (
        runTimeCondition,
        equationInitialResidualCondition,
        dictionary
    );

    template<>
    const char* Foam::NamedEnum
    <
        equationInitialResidualCondition::operatingMode,
        2
    >::names[] =
    {
        "minimum",
        "maximum"
    };

    const NamedEnum<Foam::equationInitialResidualCondition::operatingMode, 2>
        Foam::equationInitialResidualCondition::operatingModeNames;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::equationInitialResidualCondition::equationInitialResidualCondition
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    functionObjectState& state
)
:
    runTimeCondition(name, obr, dict, state),
    fieldNames_(dict.lookup("fields")),
    value_(readScalar(dict.lookup("value"))),
    timeStart_(dict.lookupOrDefault("timeStart", -GREAT)),
    mode_(operatingModeNames.read(dict.lookup("mode")))
{
    if (!fieldNames_.size())
    {
        WarningInFunction
            << "No fields supplied: deactivating" << endl;

        active_ = false;
    }
}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

Foam::equationInitialResidualCondition::
~equationInitialResidualCondition()
{}


// * * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * //

bool Foam::equationInitialResidualCondition::apply()
{
    bool satisfied = false;

    if (!active_)
    {
        return true;
    }

    if ((obr_.time().timeIndex() < 3) || (obr_.time().value() < timeStart_))
    {
        // Do not start checking until reached start time
        return false;
    }

    const fvMesh& mesh = refCast<const fvMesh>(obr_);
    const dictionary& solverDict = mesh.solverPerformanceDict();

    List<scalar> result(fieldNames_.size(), -VGREAT);

    forAll(fieldNames_, fieldI)
    {
        const word& fieldName = fieldNames_[fieldI];

        if (solverDict.found(fieldName))
        {
            const List<solverPerformance> sp(solverDict.lookup(fieldName));
            const scalar residual = sp.first().initialResidual();
            result[fieldI] = residual;

            switch (mode_)
            {
                case omMin:
                {
                    if (residual < value_)
                    {
                        satisfied = true;
                    }
                    break;
                }
                case omMax:
                {
                    if (residual > value_)
                    {
                        satisfied = true;
                    }
                    break;
                }
                default:
                {
                    FatalErrorInFunction
                        << "Unhandled enumeration "
                        << operatingModeNames[mode_]
                        << abort(FatalError);
                }
            }
        }
    }

    bool valid = false;
    forAll(result, i)
    {
        if (result[i] < 0)
        {
            WarningInFunction
                << "Initial residual data not found for field "
                << fieldNames_[i] << endl;
        }
        else
        {
            valid = true;
        }
    }

    if (!valid)
    {
        WarningInFunction
            << "Initial residual data not found for any fields: "
            << "deactivating" << endl;

        active_ = false;
    }

    if (satisfied && valid)
    {
        if (log_)
        {
            Info<< type() << ": " << name_
                << ": satisfied using threshold value: " << value_ << nl;
        }

        forAll(result, resultI)
        {
            if (result[resultI] > 0)
            {
                if (log_)
                {
                    Info<< "    field: " << fieldNames_[resultI]
                        << ", residual: " << result[resultI] << nl;
                }
            }
        }
        if (log_) Info<< endl;
    }

    return satisfied;
}


void Foam::equationInitialResidualCondition::write()
{
    // do nothing
}


// ************************************************************************* //

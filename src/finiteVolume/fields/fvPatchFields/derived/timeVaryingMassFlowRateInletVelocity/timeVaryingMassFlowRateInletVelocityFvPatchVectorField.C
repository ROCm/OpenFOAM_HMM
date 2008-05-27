/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2006-07 OpenCFD Ltd.
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

#include "timeVaryingMassFlowRateInletVelocityFvPatchVectorField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"
#include "Time.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::
timeVaryingMassFlowRateInletVelocityFvPatchVectorField::
timeVaryingMassFlowRateInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    massFlowRateInletVelocityFvPatchVectorField(p, iF)
{}


Foam::
timeVaryingMassFlowRateInletVelocityFvPatchVectorField::
timeVaryingMassFlowRateInletVelocityFvPatchVectorField
(
    const timeVaryingMassFlowRateInletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    massFlowRateInletVelocityFvPatchVectorField(ptf, p, iF, mapper),
    timeDataFile_(ptf.timeDataFile_),
    timeSeries_(ptf.timeBounding())
{}


Foam::
timeVaryingMassFlowRateInletVelocityFvPatchVectorField::
timeVaryingMassFlowRateInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    massFlowRateInletVelocityFvPatchVectorField(p, iF, dict),
    timeDataFile_(dict.lookup("timeDataFile")),
    timeSeries_(word(dict.lookup("timeBounding")))
{}


Foam::
timeVaryingMassFlowRateInletVelocityFvPatchVectorField::
timeVaryingMassFlowRateInletVelocityFvPatchVectorField
(
    const timeVaryingMassFlowRateInletVelocityFvPatchVectorField& ptf
)
:
    massFlowRateInletVelocityFvPatchVectorField(ptf),
    timeDataFile_(ptf.timeDataFile_),
    timeSeries_(ptf.timeBounding())
{}


Foam::
timeVaryingMassFlowRateInletVelocityFvPatchVectorField::
timeVaryingMassFlowRateInletVelocityFvPatchVectorField
(
    const timeVaryingMassFlowRateInletVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    massFlowRateInletVelocityFvPatchVectorField(ptf, iF),
    timeDataFile_(ptf.timeDataFile_),
    timeSeries_(ptf.timeBounding())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar
Foam::timeVaryingMassFlowRateInletVelocityFvPatchVectorField::
currentValue()
{
    if (timeSeries_.size() == 0)
    {
        fileName fName(timeDataFile_);
        fName.expand();

        if (fName.size() == 0)
        {
            FatalErrorIn
            (
                "timeVaryingMassFlowRateInletVelocity"
                "::currentValue()"
            )   << "timeDataFile not specified for Patch "
                << this->patch().name()
                << exit(FatalError);
        }
        else
        {
            // relative path
            if (fName[0] != '/')
            {
                fName = this->db().path()/fName;
            }

            // just in case we change the interface to timeSeries
            word boundType = timeBounding();

            IFstream(fName)() >> timeSeries_;
            timeSeries_.bounding(boundType);

            // be a bit paranoid and check that the list is okay
            timeSeries_.check();
        }

        if (timeSeries_.size() == 0)
        {
            FatalErrorIn
            (
                "timeVaryingMassFlowRateInletVelocity"
                "::currentValue()"
            )   << "empty time series for Patch "
                << this->patch().name()
                << exit(FatalError);
        }
    }

    return timeSeries_(this->db().time().timeOutputValue());
}


void Foam::
timeVaryingMassFlowRateInletVelocityFvPatchVectorField::
updateCoeffs()
{
    if (updated())
    {
        return;
    }

    massFlowRate() = currentValue();
    massFlowRateInletVelocityFvPatchVectorField::updateCoeffs();
}


void Foam::
timeVaryingMassFlowRateInletVelocityFvPatchVectorField::
write(Ostream& os) const
{
    massFlowRateInletVelocityFvPatchVectorField::write(os);
    os.writeKeyword("timeDataFile")
        << timeDataFile_ << token::END_STATEMENT << nl;
    os.writeKeyword("timeBounding")
        << timeBounding() << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchVectorField,
       timeVaryingMassFlowRateInletVelocityFvPatchVectorField
   );
}


// ************************************************************************* //

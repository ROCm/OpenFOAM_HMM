/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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

#include "timeVaryingUniformFixedValueFvPatchField.H"
#include "Time.H"
#include "Tuple2.H"
#include "IFstream.H"
#include "interpolateXY.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
timeVaryingUniformFixedValueFvPatchField<Type>::
timeVaryingUniformFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF)
{}


template<class Type>
timeVaryingUniformFixedValueFvPatchField<Type>::
timeVaryingUniformFixedValueFvPatchField
(
    const timeVaryingUniformFixedValueFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper),
    timeDataFileName_(ptf.timeDataFileName_)
{}


template<class Type>
timeVaryingUniformFixedValueFvPatchField<Type>::
timeVaryingUniformFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF),
    timeDataFileName_(fileName(dict.lookup("timeDataFileName")).expand())
{
    if (dict.found("value"))
    {
        fvPatchField<Type>::operator==(Field<Type>("value", dict, p.size()));
    }
    else
    {
        updateCoeffs();
    }
}


template<class Type>
timeVaryingUniformFixedValueFvPatchField<Type>::
timeVaryingUniformFixedValueFvPatchField
(
    const timeVaryingUniformFixedValueFvPatchField<Type>& ptf
)
:
    fixedValueFvPatchField<Type>(ptf),
    timeDataFileName_(ptf.timeDataFileName_)
{}


template<class Type>
timeVaryingUniformFixedValueFvPatchField<Type>::
timeVaryingUniformFixedValueFvPatchField
(
    const timeVaryingUniformFixedValueFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
    timeDataFileName_(ptf.timeDataFileName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void timeVaryingUniformFixedValueFvPatchField<Type>::checkTable()
{
    const Time& tm = this->db().time();

    if (times_.size() == 0)
    {
        if (timeDataFileName_.size() == 0)
        {
            FatalErrorIn
            (
                "timeVaryingUniformFixedValueFvPatchField<Type>"
                "::checkTable()"
            )   << "timeDataFileName not specified for Patch "
                << this->patch().name()
                << exit(FatalError);
        }
        else
        {
            IFstream str(timeDataFileName_);

            List<Tuple2<scalar, Type> > timeValues(str);

            times_.setSize(timeValues.size());
            values_.setSize(timeValues.size());

            forAll(timeValues, i)
            {
                times_[i] = timeValues[i].first();
                values_[i] = timeValues[i].second();
            }
        }
    }

    if (tm.value() < times_[0])
    {
        WarningIn
        (
            "timeVaryingUniformFixedValueFvPatchField<Type>::checkTable()"
        )   << "current time (" << tm.value()
            << ") is less than the minimum in the data table ("
            << times_[0] << ')' << endl
            << "    Continuing with the value for the smallest time"
            << endl;
    }

    if (tm.value() > times_[times_.size()-1])
    {
        WarningIn
        (
            "timeVaryingUniformFixedValueFvPatchField<Type>::checkTable()"
        )   << "current time (" << tm.value()
            << ") is greater than the maximum in the data table ("
            << times_[times_.size()-1] << ')' << endl
            << "    Continuing with the value for the largest time"
            << endl;
    }
}


template<class Type>
void timeVaryingUniformFixedValueFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    checkTable();

    this->operator==
    (
        interpolateXY
        (
            this->db().time().value(),
            times_,
            values_
        )
    );

    fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void timeVaryingUniformFixedValueFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    os.writeKeyword("timeDataFileName")
        << timeDataFileName_ << token::END_STATEMENT << nl;
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

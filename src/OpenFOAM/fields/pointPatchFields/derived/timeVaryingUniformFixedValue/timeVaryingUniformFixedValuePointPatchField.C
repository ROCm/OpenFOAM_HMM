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

#include "timeVaryingUniformFixedValuePointPatchField.H"
#include "Time.H"
#include "pointMesh.H"
#include "Tuple2.H"
#include "IFstream.H"
#include "interpolateXY.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * //

template<class Type>
timeVaryingUniformFixedValuePointPatchField<Type>::
timeVaryingUniformFixedValuePointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF
)
:
    fixedValuePointPatchField<Type>(p, iF)
{}


template<class Type>
timeVaryingUniformFixedValuePointPatchField<Type>::
timeVaryingUniformFixedValuePointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const dictionary& dict
)
:
    fixedValuePointPatchField<Type>(p, iF, dict),
    timeDataFileName_(fileName(dict.lookup("timeDataFileName")).expand())
{}


template<class Type>
timeVaryingUniformFixedValuePointPatchField<Type>::
timeVaryingUniformFixedValuePointPatchField
(
    const timeVaryingUniformFixedValuePointPatchField<Type>& ptf,
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    fixedValuePointPatchField<Type>(ptf, p, iF, mapper),
    timeDataFileName_(ptf.timeDataFileName_)
{}


template<class Type>
timeVaryingUniformFixedValuePointPatchField<Type>::
timeVaryingUniformFixedValuePointPatchField
(
    const timeVaryingUniformFixedValuePointPatchField<Type>& ptf
)
:
    fixedValuePointPatchField<Type>(ptf),
    timeDataFileName_(ptf.timeDataFileName_)
{}


template<class Type>
timeVaryingUniformFixedValuePointPatchField<Type>::
timeVaryingUniformFixedValuePointPatchField
(
    const timeVaryingUniformFixedValuePointPatchField<Type>& ptf,
    const DimensionedField<Type, pointMesh>& iF
)
:
    fixedValuePointPatchField<Type>(ptf, iF),
    timeDataFileName_(ptf.timeDataFileName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void timeVaryingUniformFixedValuePointPatchField<Type>::checkTable()
{
    const Time& tm = this->db().time();

    if (times_.size() == 0)
    {
        if (timeDataFileName_.size() == 0)
        {
            FatalErrorIn
            (
                "timeVaryingUniformFixedValuePointPatchField<Type>"
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

    if (tm.value() < min(times_))
    {
        WarningIn
        (
            "timeVaryingUniformFixedValuePointPatchField<Type>::checkTable()"
        )   << "current time (" << tm.value()
            << ") is less than the minimum in the data table ("
            << min(times_) << ')' << endl
            << "    Continuing with the value for the smallest time"
            << endl;
    }

    if (tm.value() > max(times_))
    {
        WarningIn
        (
            "timeVaryingUniformFixedValuePointPatchField<Type>::checkTable()"
        )   << "current time (" << tm.value()
            << ") is greater than the maximum in the data table ("
            << max(times_) << ')' << endl
            << "    Continuing with the value for the largest time"
            << endl;
    }
}


template<class Type>
void timeVaryingUniformFixedValuePointPatchField<Type>::updateCoeffs()
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

    fixedValuePointPatchField<Type>::updateCoeffs();
}


template<class Type>
void timeVaryingUniformFixedValuePointPatchField<Type>::
write(Ostream& os) const
{
    fixedValuePointPatchField<Type>::write(os);
    os.writeKeyword("timeDataFileName")
        << timeDataFileName_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

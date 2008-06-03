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
#include "IFstream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::timeVaryingUniformFixedValueFvPatchField<Type>::
timeVaryingUniformFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF)
{}


template<class Type>
Foam::timeVaryingUniformFixedValueFvPatchField<Type>::
timeVaryingUniformFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF),
    timeDataFile_(dict.lookup("timeDataFile")),
    timeSeries_(word(dict.lookup("timeBounding")))
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
Foam::timeVaryingUniformFixedValueFvPatchField<Type>::
timeVaryingUniformFixedValueFvPatchField
(
    const timeVaryingUniformFixedValueFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper),
    timeDataFile_(ptf.timeDataFile_),
    timeSeries_(ptf.timeBounding())
{}


template<class Type>
Foam::timeVaryingUniformFixedValueFvPatchField<Type>::
timeVaryingUniformFixedValueFvPatchField
(
    const timeVaryingUniformFixedValueFvPatchField<Type>& ptf
)
:
    fixedValueFvPatchField<Type>(ptf),
    timeDataFile_(ptf.timeDataFile_),
    timeSeries_(ptf.timeBounding())
{}


template<class Type>
Foam::timeVaryingUniformFixedValueFvPatchField<Type>::
timeVaryingUniformFixedValueFvPatchField
(
    const timeVaryingUniformFixedValueFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
    timeDataFile_(ptf.timeDataFile_),
    timeSeries_(ptf.timeBounding())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Type Foam::timeVaryingUniformFixedValueFvPatchField<Type>::
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
                "timeVaryingUniformFixedValueFvPatchField"
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
                "timeVaryingUniformFixedValueFvPatchField"
                "::currentValue()"
            )   << "empty time series for Patch "
                << this->patch().name()
                << exit(FatalError);
        }
    }

    return timeSeries_(this->db().time().timeOutputValue());
}


template<class Type>
void Foam::timeVaryingUniformFixedValueFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    fvPatchField<Type>::operator==(currentValue());
    fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::timeVaryingUniformFixedValueFvPatchField<Type>::write
(
    Ostream& os
) const
{
    fvPatchField<Type>::write(os);
    os.writeKeyword("timeDataFile")
        << timeDataFile_ << token::END_STATEMENT << nl;
    os.writeKeyword("timeBounding")
        << timeBounding() << token::END_STATEMENT << nl;
    this->writeEntry("value", os);
}


// ************************************************************************* //

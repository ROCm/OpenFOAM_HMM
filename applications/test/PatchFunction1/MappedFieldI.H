/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenCFD Ltd.
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

#include "MappedField.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
inline Foam::tmp<Foam::Field<Type>>
Foam::PatchFunction1Types::MappedField<Type>::value
(
    const scalar x
) const
{
    const polyMesh& mesh = this->patch_.boundaryMesh().mesh();
    checkTable();

    tmp<Field<Type>> tfld(new Field<Type>(this->patch_.size()));
    Field<Type>& fld = tfld.ref();
    Type wantedAverage;

    if (endSampleTime_ == -1)
    {
        // Only start value
        if (debug)
        {
            Pout<< "updateCoeffs : Sampled, non-interpolated values"
                << " from start time:"
                << sampleTimes_[startSampleTime_].name() << nl;
        }

        fld = startSampledValues_;
        wantedAverage = startAverage_;
    }
    else
    {
        scalar start = sampleTimes_[startSampleTime_].value();
        scalar end = sampleTimes_[endSampleTime_].value();

        scalar s = (mesh.time().value() - start)/(end - start);

        if (debug)
        {
            Pout<< "updateCoeffs : Sampled, interpolated values"
                << " between start time:"
                << sampleTimes_[startSampleTime_].name()
                << " and end time:" << sampleTimes_[endSampleTime_].name()
                << " with weight:" << s << endl;
        }

        fld = ((1 - s)*startSampledValues_ + s*endSampledValues_);
        wantedAverage = (1 - s)*startAverage_ + s*endAverage_;
    }

    // Enforce average. Either by scaling (if scaling factor > 0.5) or by
    // offsetting.
    if (setAverage_)
    {
        const scalarField magSf(mag(this->patch_.faceAreas()));

        Type averagePsi = gSum(magSf*fld)/gSum(magSf);

        if (debug)
        {
            Pout<< "updateCoeffs :"
                << " actual average:" << averagePsi
                << " wanted average:" << wantedAverage
                << endl;
        }

        if (mag(averagePsi) < VSMALL)
        {
            // Field too small to scale. Offset instead.
            const Type offset = wantedAverage - averagePsi;
            if (debug)
            {
                Pout<< "updateCoeffs :"
                    << " offsetting with:" << offset << endl;
            }
            fld += offset;
        }
        else
        {
            const scalar scale = mag(wantedAverage)/mag(averagePsi);

            if (debug)
            {
                Pout<< "updateCoeffs :"
                    << " scaling with:" << scale << endl;
            }
            fld *= scale;
        }
    }

    // Apply offset to mapped values
    if (offset_.valid())
    {
        const scalar t = mesh.time().timeOutputValue();
        fld += offset_->value(t);
    }

    if (debug)
    {
        Pout<< "updateCoeffs : set fixedValue to min:" << gMin(fld)
            << " max:" << gMax(fld)
            << " avg:" << gAverage(fld) << endl;
    }

    return this->transform(tfld);
}


template<class Type>
inline Foam::tmp<Foam::Field<Type>>
Foam::PatchFunction1Types::MappedField<Type>::integrate
(
    const scalar x1,
    const scalar x2
) const
{
    NotImplemented;
    return tmp<Field<Type>>(nullptr);
}


// ************************************************************************* //

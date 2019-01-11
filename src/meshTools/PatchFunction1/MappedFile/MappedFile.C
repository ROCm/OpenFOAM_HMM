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

#include "polyMesh.H"
#include "IFstream.H"
#include "AverageField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::PatchFunction1Types::MappedFile<Type>::MappedFile
(
    const polyPatch& pp,
    const word& type,
    const word& entryName,
    const dictionary& dict,
    const bool faceValues
)
:
    PatchFunction1<Type>(pp, entryName, dict, faceValues),
    dictConstructed_(true),
    fieldTableName_(entryName),
    setAverage_(dict.lookupOrDefault("setAverage", false)),
    perturb_(dict.lookupOrDefault("perturb", 1e-5)),
    pointsName_(dict.lookupOrDefault<word>("points", "points")),
    mapMethod_
    (
        dict.lookupOrDefault<word>
        (
            "mapMethod",
            "planarInterpolation"
        )
    ),
    mapperPtr_(nullptr),
    sampleTimes_(0),
    startSampleTime_(-1),
    startSampledValues_(0),
    startAverage_(Zero),
    endSampleTime_(-1),
    endSampledValues_(0),
    endAverage_(Zero),
    offset_()
{
    if (dict.found("offset"))
    {
        offset_ = Function1<Type>::New("offset", dict);
    }

    if
    (
        mapMethod_ != "planarInterpolation"
     && mapMethod_ != "nearest"
    )
    {
        FatalIOErrorInFunction(dict)
            << "mapMethod should be one of 'planarInterpolation'"
            << ", 'nearest'" << exit(FatalIOError);
    }

    dict.readIfPresent("fieldTable", fieldTableName_);
}


template<class Type>
Foam::PatchFunction1Types::MappedFile<Type>::MappedFile
(
    const polyPatch& pp,
    const word& entryName,
    const dictionary& dict,
    const word& fieldTableName,
    const bool faceValues
)
:
    PatchFunction1<Type>(pp, entryName, dict, faceValues),
    dictConstructed_(false),
    fieldTableName_(fieldTableName),
    setAverage_(dict.lookupOrDefault("setAverage", false)),
    perturb_(dict.lookupOrDefault("perturb", 1e-5)),
    pointsName_(dict.lookupOrDefault<word>("points", "points")),
    mapMethod_
    (
        dict.lookupOrDefault<word>
        (
            "mapMethod",
            "planarInterpolation"
        )
    ),
    mapperPtr_(nullptr),
    sampleTimes_(0),
    startSampleTime_(-1),
    startSampledValues_(0),
    startAverage_(Zero),
    endSampleTime_(-1),
    endSampledValues_(0),
    endAverage_(Zero),
    offset_()
{
    if (dict.found("offset"))
    {
        offset_ = Function1<Type>::New("offset", dict);
    }

    if
    (
        mapMethod_ != "planarInterpolation"
     && mapMethod_ != "nearest"
    )
    {
        FatalIOErrorInFunction(dict)
            << "mapMethod should be one of 'planarInterpolation'"
            << ", 'nearest'" << exit(FatalIOError);
    }
}


template<class Type>
Foam::PatchFunction1Types::MappedFile<Type>::MappedFile
(
    const MappedFile<Type>& ut
)
:
    PatchFunction1<Type>(ut),
    dictConstructed_(ut.dictConstructed_),
    fieldTableName_(ut.fieldTableName_),
    setAverage_(ut.setAverage_),
    perturb_(ut.perturb_),
    pointsName_(ut.pointsName_),
    mapMethod_(ut.mapMethod_),
    mapperPtr_(ut.mapperPtr_.clone()),
    sampleTimes_(ut.sampleTimes_),
    startSampleTime_(ut.startSampleTime_),
    startSampledValues_(ut.startSampledValues_),
    startAverage_(ut.startAverage_),
    endSampleTime_(ut.endSampleTime_),
    endSampledValues_(ut.endSampledValues_),
    endAverage_(ut.endAverage_),
    offset_(ut.offset_.clone())
{}


template<class Type>
Foam::PatchFunction1Types::MappedFile<Type>::MappedFile
(
    const MappedFile<Type>& ut,
    const polyPatch& pp
)
:
    PatchFunction1<Type>(ut, pp),
    dictConstructed_(ut.dictConstructed_),
    fieldTableName_(ut.fieldTableName_),
    setAverage_(ut.setAverage_),
    perturb_(ut.perturb_),
    pointsName_(ut.pointsName_),
    mapMethod_(ut.mapMethod_),
    mapperPtr_(ut.mapperPtr_.clone()),
    sampleTimes_(ut.sampleTimes_),
    startSampleTime_(ut.startSampleTime_),
    startSampledValues_(ut.startSampledValues_),
    startAverage_(ut.startAverage_),
    endSampleTime_(ut.endSampleTime_),
    endSampledValues_(ut.endSampledValues_),
    endAverage_(ut.endAverage_),
    offset_(ut.offset_.clone())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::PatchFunction1Types::MappedFile<Type>::autoMap
(
    const FieldMapper& mapper
)
{
    PatchFunction1<Type>::autoMap(mapper);

    if (startSampledValues_.size())
    {
        startSampledValues_.autoMap(mapper);
        endSampledValues_.autoMap(mapper);
    }
    // Clear interpolator
    mapperPtr_.clear();
    startSampleTime_ = -1;
    endSampleTime_ = -1;
}


template<class Type>
void Foam::PatchFunction1Types::MappedFile<Type>::rmap
(
    const PatchFunction1<Type>& pf1,
    const labelList& addr
)
{
    PatchFunction1<Type>::rmap(pf1, addr);

    const PatchFunction1Types::MappedFile<Type>& tiptf =
        refCast<const PatchFunction1Types::MappedFile<Type>>(pf1);

    startSampledValues_.rmap(tiptf.startSampledValues_, addr);
    endSampledValues_.rmap(tiptf.endSampledValues_, addr);

    // Clear interpolator
    mapperPtr_.clear();
    startSampleTime_ = -1;
    endSampleTime_ = -1;
}


template<class Type>
void Foam::PatchFunction1Types::MappedFile<Type>::checkTable
(
    const scalar t
) const
{
    const polyMesh& mesh = this->patch_.boundaryMesh().mesh();

    // Initialise
    if (mapperPtr_.empty())
    {
        // Reread values and interpolate
        fileName samplePointsFile
        (
            mesh.time().path()
           /mesh.time().caseConstant()
           /"boundaryData"
           /this->patch_.name()
           /pointsName_
        );

        pointField samplePoints((IFstream(samplePointsFile)()));

        DebugInfo
            << "Read " << samplePoints.size() << " sample points from "
            << samplePointsFile << endl;


        // tbd: run-time selection
        bool nearestOnly =
        (
           !mapMethod_.empty()
         && mapMethod_ != "planarInterpolation"
        );

        // Allocate the interpolator
        if (this->faceValues_)
        {
            mapperPtr_.reset
            (
                new pointToPointPlanarInterpolation
                (
                    samplePoints,
                    this->localPosition(this->patch_.faceCentres()),
                    perturb_,
                    nearestOnly
                )
            );
        }
        else
        {
            mapperPtr_.reset
            (
                new pointToPointPlanarInterpolation
                (
                    samplePoints,
                    this->localPosition(this->patch_.localPoints()),
                    perturb_,
                    nearestOnly
                )
            );
        }



        // Read the times for which data is available
        const fileName samplePointsDir = samplePointsFile.path();
        sampleTimes_ = Time::findTimes(samplePointsDir);

        DebugInfo
            << "In directory "
            << samplePointsDir << " found times "
            << pointToPointPlanarInterpolation::timeNames(sampleTimes_)
            << endl;
    }


    // Find current time in sampleTimes
    label lo = -1;
    label hi = -1;

    bool foundTime = mapperPtr_().findTime
    (
        sampleTimes_,
        startSampleTime_,
        t,  //mesh.time().value(),
        lo,
        hi
    );

    if (!foundTime)
    {
        FatalErrorInFunction
            << "Cannot find starting sampling values for index "
            << t << nl
            << "Have sampling values for "
            << pointToPointPlanarInterpolation::timeNames(sampleTimes_) << nl
            << "In directory "
            <<  mesh.time().constant()/"boundaryData"/this->patch_.name()
            << "\n    on patch " << this->patch_.name()
            << " of field " << fieldTableName_
            << exit(FatalError);
    }


    // Update sampled data fields.

    if (lo != startSampleTime_)
    {
        startSampleTime_ = lo;

        if (startSampleTime_ == endSampleTime_)
        {
            // No need to reread since are end values
            if (debug)
            {
                Pout<< "checkTable : Setting startValues to (already read) "
                    << "boundaryData"
                      /this->patch_.name()
                      /sampleTimes_[startSampleTime_].name()
                    << endl;
            }
            startSampledValues_ = endSampledValues_;
            startAverage_ = endAverage_;
        }
        else
        {
            if (debug)
            {
                Pout<< "checkTable : Reading startValues from "
                    << "boundaryData"
                      /this->patch_.name()
                      /sampleTimes_[lo].name()
                    << endl;
            }


            // Reread values and interpolate
            fileName valsFile
            (
                mesh.time().path()
               /mesh.time().caseConstant()
               /"boundaryData"
               /this->patch_.name()
               /sampleTimes_[startSampleTime_].name()
               /fieldTableName_
            );

            Field<Type> vals;

            if (setAverage_)
            {
                AverageField<Type> avals((IFstream(valsFile)()));
                vals = avals;
                startAverage_ = avals.average();
            }
            else
            {
                IFstream(valsFile)() >> vals;
            }

            if (vals.size() != mapperPtr_().sourceSize())
            {
                FatalErrorInFunction
                    << "Number of values (" << vals.size()
                    << ") differs from the number of points ("
                    <<  mapperPtr_().sourceSize()
                    << ") in file " << valsFile << exit(FatalError);
            }

            startSampledValues_ = mapperPtr_().interpolate(vals);
        }
    }

    if (hi != endSampleTime_)
    {
        endSampleTime_ = hi;

        if (endSampleTime_ == -1)
        {
            // endTime no longer valid. Might as well clear endValues.
            if (debug)
            {
                Pout<< "checkTable : Clearing endValues" << endl;
            }
            endSampledValues_.clear();
        }
        else
        {
            if (debug)
            {
                Pout<< "checkTable : Reading endValues from "
                    << "boundaryData"
                      /this->patch_.name()
                      /sampleTimes_[endSampleTime_].name()
                    << endl;
            }

            // Reread values and interpolate
            fileName valsFile
            (
                mesh.time().path()
               /mesh.time().caseConstant()
               /"boundaryData"
               /this->patch_.name()
               /sampleTimes_[endSampleTime_].name()
               /fieldTableName_
            );

            Field<Type> vals;

            if (setAverage_)
            {
                AverageField<Type> avals((IFstream(valsFile)()));
                vals = avals;
                endAverage_ = avals.average();
            }
            else
            {
                IFstream(valsFile)() >> vals;
            }

            if (vals.size() != mapperPtr_().sourceSize())
            {
                FatalErrorInFunction
                    << "Number of values (" << vals.size()
                    << ") differs from the number of points ("
                    <<  mapperPtr_().sourceSize()
                    << ") in file " << valsFile << exit(FatalError);
            }

            endSampledValues_ = mapperPtr_().interpolate(vals);
        }
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::PatchFunction1Types::MappedFile<Type>::value
(
    const scalar x
) const
{
    checkTable(x);

    auto tfld = tmp<Field<Type>>::New(startSampledValues_.size());
    Field<Type>& fld = tfld.ref();
    Type wantedAverage;

    if (endSampleTime_ == -1)
    {
        // Only start value
        if (debug)
        {
            Pout<< "MappedFile<Type>::value : Sampled, non-interpolated values"
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

        scalar s = (x - start)/(end - start);

        if (debug)
        {
            Pout<< "MappedFile<Type>::value : Sampled, interpolated values"
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
        Type averagePsi;
        if (this->faceValues_)
        {
            const scalarField magSf(mag(this->patch_.faceAreas()));
            averagePsi = gSum(magSf*fld)/gSum(magSf);
        }
        else
        {
            averagePsi = gAverage(fld);
        }

        if (debug)
        {
            Pout<< "MappedFile<Type>::value :"
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
                Pout<< "MappedFile<Type>::value :"
                    << " offsetting with:" << offset << endl;
            }
            fld += offset;
        }
        else
        {
            const scalar scale = mag(wantedAverage)/mag(averagePsi);

            if (debug)
            {
                Pout<< "MappedFile<Type>::value :"
                    << " scaling with:" << scale << endl;
            }
            fld *= scale;
        }
    }

    // Apply offset to mapped values
    if (offset_.valid())
    {
        fld += offset_->value(x);
    }

    if (debug)
    {
        Pout<< "MappedFile<Type>::value : set fixedValue to min:" << gMin(fld)
            << " max:" << gMax(fld)
            << " avg:" << gAverage(fld) << endl;
    }

    return this->transform(tfld);
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::PatchFunction1Types::MappedFile<Type>::integrate
(
    const scalar x1,
    const scalar x2
) const
{
    NotImplemented;
    return nullptr;
}


template<class Type>
void Foam::PatchFunction1Types::MappedFile<Type>::writeData
(
    Ostream& os
) const
{
    PatchFunction1<Type>::writeData(os);

    // Check if field name explicitly provided (e.g. through timeVaryingMapped
    // bc)
    if (dictConstructed_)
    {
        os.writeEntry(this->name(), type());

        os.writeEntryIfDifferent
        (
            "fieldTable",
            this->name(),
            fieldTableName_
        );
    }

    if (setAverage_)
    {
        os.writeEntry("setAverage", setAverage_);
    }

    os.writeEntryIfDifferent<scalar>("perturb", 1e-5, perturb_);

    os.writeEntryIfDifferent<word>("points", "points", pointsName_);

    os.writeEntryIfDifferent<word>
    (
        "mapMethod",
        "planarInterpolation",
        mapMethod_
    );

    if (offset_.valid())
    {
        offset_->writeData(os);
    }
}


// ************************************************************************* //

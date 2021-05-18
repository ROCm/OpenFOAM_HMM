/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2021 OpenCFD Ltd.
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
#include "rawIOField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::PatchFunction1Types::MappedFile<Type>::MappedFile
(
    const polyPatch& pp,
    const word& redirectType,
    const word& entryName,
    const dictionary& dict,
    const bool faceValues
)
:
    PatchFunction1<Type>(pp, entryName, dict, faceValues),
    dictConstructed_(true),
    setAverage_(dict.getOrDefault("setAverage", false)),
    fieldTableName_(entryName),
    perturb_(dict.getOrDefault<scalar>("perturb", 1e-5)),
    pointsName_(dict.getOrDefault<word>("points", "points")),
    mapMethod_
    (
        dict.getOrDefault<word>
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
    offset_(Function1<Type>::NewIfPresent("offset", dict))
{
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
    setAverage_(dict.getOrDefault("setAverage", false)),
    fieldTableName_(fieldTableName),
    perturb_(dict.getOrDefault<scalar>("perturb", 1e-5)),
    pointsName_(dict.getOrDefault<word>("points", "points")),
    mapMethod_
    (
        dict.getOrDefault<word>
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
    offset_(Function1<Type>::NewIfPresent("offset", dict))
{
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
    const MappedFile<Type>& rhs
)
:
    MappedFile<Type>(rhs, rhs.patch())
{}


template<class Type>
Foam::PatchFunction1Types::MappedFile<Type>::MappedFile
(
    const MappedFile<Type>& rhs,
    const polyPatch& pp
)
:
    PatchFunction1<Type>(rhs, pp),
    dictConstructed_(rhs.dictConstructed_),
    setAverage_(rhs.setAverage_),
    fieldTableName_(rhs.fieldTableName_),
    perturb_(rhs.perturb_),
    pointsName_(rhs.pointsName_),
    mapMethod_(rhs.mapMethod_),
    mapperPtr_(rhs.mapperPtr_.clone()),
    sampleTimes_(rhs.sampleTimes_),
    startSampleTime_(rhs.startSampleTime_),
    startSampledValues_(rhs.startSampledValues_),
    startAverage_(rhs.startAverage_),
    endSampleTime_(rhs.endSampleTime_),
    endSampledValues_(rhs.endSampledValues_),
    endAverage_(rhs.endAverage_),
    offset_(rhs.offset_.clone())
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
    }

    if (endSampledValues_.size())
    {
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

    if (tiptf.startSampledValues_.size())
    {
        startSampledValues_.setSize(this->size());
        startSampledValues_.rmap(tiptf.startSampledValues_, addr);
    }

    if (tiptf.endSampledValues_.size())
    {
        endSampledValues_.setSize(this->size());
        endSampledValues_.rmap(tiptf.endSampledValues_, addr);
    }

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
    const Time& time = mesh.time();

    // Initialise
    if (!mapperPtr_)
    {
        // Reread values and interpolate
        const fileName samplePointsFile
        (
            time.globalPath()
           /time.constant()
           /"boundaryData"
           /this->patch_.name()
           /pointsName_
        );

        IOobject io
        (
            samplePointsFile,   // absolute path
            time,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false,              // no need to register
            true                // is global object (currently not used)
        );

        // Read data
        const rawIOField<point> samplePoints(io, false);

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
        if (this->faceValues())
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
            <<  time.constant()/"boundaryData"/this->patch_.name()
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
            const fileName valsFile
            (
                time.globalPath()
               /time.constant()
               /"boundaryData"
               /this->patch_.name()
               /sampleTimes_[startSampleTime_].name()
               /fieldTableName_
            );

            IOobject io
            (
                valsFile,   // absolute path
                time,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false,              // no need to register
                true                // is global object (currently not used)
            );

            const rawIOField<Type> vals(io, setAverage_);
            if (setAverage_)
            {
                startAverage_ = vals.average();
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
                time.globalPath()
               /time.constant()
               /"boundaryData"
               /this->patch_.name()
               /sampleTimes_[endSampleTime_].name()
               /fieldTableName_
            );

            IOobject io
            (
                valsFile,   // absolute path
                time,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false,              // no need to register
                true                // is global object (currently not used)
            );

            const rawIOField<Type> vals(io, setAverage_);
            if (setAverage_)
            {
                endAverage_ = vals.average();
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
    auto& fld = tfld.ref();
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
        if (this->faceValues())
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
    if (offset_)
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
void Foam::PatchFunction1Types::MappedFile<Type>::writeEntries
(
    Ostream& os
) const
{
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

    if (offset_)
    {
        offset_->writeData(os);
    }
}


template<class Type>
void Foam::PatchFunction1Types::MappedFile<Type>::writeData
(
    Ostream& os
) const
{
    PatchFunction1<Type>::writeData(os);

    // Check if field name explicitly provided
    // (e.g. through timeVaryingMapped bc)
    if (dictConstructed_)
    {
        os.writeEntry(this->name(), type());

        os.writeEntryIfDifferent
        (
            "fieldTable",
            this->name(),
            fieldTableName_
        );

        os.beginBlock(word(this->name() + "Coeffs"));
        writeEntries(os);
        os.endBlock();
    }
    else
    {
        writeEntries(os);
    }
}


// ************************************************************************* //

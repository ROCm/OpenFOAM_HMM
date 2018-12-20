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
Foam::PatchFunction1Types::MappedField<Type>::MappedField
(
    const polyPatch& pp,
    const word& entryName,
    const dictionary& dict
)
:
    PatchFunction1<Type>(pp, entryName, dict),
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
Foam::PatchFunction1Types::MappedField<Type>::MappedField
(
    const MappedField<Type>& ut
)
:
    PatchFunction1<Type>(ut),
    fieldTableName_(ut.fieldTableName_),
    setAverage_(ut.setAverage_),
    perturb_(ut.perturb_),
    pointsName_(ut.pointsName_),
    mapMethod_(ut.mapMethod_),
    mapperPtr_(nullptr),
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
Foam::PatchFunction1Types::MappedField<Type>::MappedField
(
    const MappedField<Type>& ut,
    const polyPatch& pp
)
:
    PatchFunction1<Type>(ut, pp),
    fieldTableName_(ut.fieldTableName_),
    setAverage_(ut.setAverage_),
    perturb_(ut.perturb_),
    pointsName_(ut.pointsName_),
    mapMethod_(ut.mapMethod_),
    mapperPtr_(nullptr),
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
void Foam::PatchFunction1Types::MappedField<Type>::autoMap
(
    const FieldMapper& mapper
)
{
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
void Foam::PatchFunction1Types::MappedField<Type>::rmap
(
    const PatchFunction1<Type>& pf1,
    const labelList& addr
)
{
    const PatchFunction1Types::MappedField<Type>& tiptf =
        refCast<const PatchFunction1Types::MappedField<Type>>(pf1);

    startSampledValues_.rmap(tiptf.startSampledValues_, addr);
    endSampledValues_.rmap(tiptf.endSampledValues_, addr);

    // Clear interpolator
    mapperPtr_.clear();
    startSampleTime_ = -1;
    endSampleTime_ = -1;
}


template<class Type>
void Foam::PatchFunction1Types::MappedField<Type>::checkTable() const
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
            << " Read " << samplePoints.size() << " sample points from "
            << samplePointsFile << endl;


        // tbd: run-time selection
        bool nearestOnly =
        (
           !mapMethod_.empty()
         && mapMethod_ != "planarInterpolation"
        );

        // Allocate the interpolator
        mapperPtr_.reset
        (
            new pointToPointPlanarInterpolation
            (
                samplePoints,
                this->patch_.faceCentres(),
                perturb_,
                nearestOnly
            )
        );

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
        mesh.time().value(),
        lo,
        hi
    );

    if (!foundTime)
    {
        FatalErrorInFunction
            << "Cannot find starting sampling values for current time "
            << mesh.time().value() << nl
            << "Have sampling values for times "
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
void Foam::PatchFunction1Types::MappedField<Type>::writeData
(
    Ostream& os
) const
{
    PatchFunction1<Type>::writeData(os);
    //os  << token::END_STATEMENT << nl;
//    uniformValuePtr_->writeData(os);
    //os  << endl;
}


// ************************************************************************* //

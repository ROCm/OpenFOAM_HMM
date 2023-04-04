/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2023 OpenCFD Ltd.
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

#include "parProfiling.H"
#include "addToRunTimeSelectionTable.H"
#include "Pstream.H"
#include "PstreamReduceOps.H"
#include "profilingPstream.H"
#include "Tuple2.H"
#include "FixedList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(parProfiling, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        parProfiling,
        dictionary
    );

} // End namespace functionObject
} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::parProfiling::parProfiling
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    functionObject(name),
    detailLevel_(0)
{
    dict.readIfPresent("detail", detailLevel_);
    profilingPstream::enable();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::parProfiling::~parProfiling()
{
    profilingPstream::disable();
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Loop over all values (with striding) and extract the value at given index
template<class Type>
inline static void extractValues
(
    UList<Type>& result,
    const int index,
    const UList<Type>& allValues
)
{
    if (result.empty())
    {
        return;
    }

    const label numProc = result.size();
    const Type* values = allValues.cbegin();
    const label stride = allValues.size() / numProc;

    if (!values || !stride)
    {
        result = Type(0);
        return;
    }

    for (label proci = 0; proci < numProc; ++proci, values += stride)
    {
        result[proci] = values[index];
    }
}


// Loop over all values (with striding) and extract combined value
// using the given unary function
template<class Type, class Extract>
inline static void extractValues
(
    UList<Type>& result,
    const UList<Type>& allValues,
    const Extract& extract
)
{
    if (result.empty())
    {
        return;
    }

    const label numProc = result.size();
    const Type* values = allValues.cbegin();
    const label stride = allValues.size() / numProc;

    if (!values || !stride)
    {
        result = Type(0);
        return;
    }

    for (label proci = 0; proci < numProc; ++proci, values += stride)
    {
        result[proci] = extract(values);
    }
}


inline static void printTimingDetail(const UList<double>& values)
{
    const label numProc = values.size();

    if (numProc)
    {
        Info<< indent << "    times   " << numProc << '(';

        for (label proci = 0; proci < numProc; ++proci)
        {
            if (proci) Info<< ' ';
            Info<< values[proci];
        }

        Info<< ')' << nl;
    }
}


inline static void printTimingDetail(const UList<uint64_t>& values)
{
    const label numProc = values.size();

    if (numProc)
    {
        // Output via std::ostream to avoid conversion to Foam::label
        // that Ostream performs

        auto& os = Info.stdStream();

        Info<< indent << "    counts  " << numProc << '(';

        for (label proci = 0; proci < numProc; ++proci)
        {
            if (proci) os << ' ';
            os << values[proci];
        }

        Info<< ')' << nl;
    }
}

} // End namespace Foam


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::functionObjects::parProfiling::report()
{
    const label numProc = (UPstream::parRun() ? UPstream::nProcs() : 1);

    if (!profilingPstream::active() || numProc < 2)
    {
        return;
    }

    // Use mpiGather on all values and perform the combinations
    // and statistics locally. This reduces the overall number of MPI
    // calls. For detailed output we need this information anyhow.

    // NB: profilingPstream uses a FixedList for timings(), counts()
    // so the sizes are guaranteed to be consistent and identical
    // everywhere.

    List<double> allTimes;
    List<uint64_t> allCounts;

    // Avoid disturbing the counts
    profilingPstream::suspend();

    {
        // The timings
        const auto& procTimes = profilingPstream::times();

        if (Pstream::master())
        {
            allTimes.resize(numProc * procTimes.size());
        }

        UPstream::mpiGather
        (
            procTimes.cdata_bytes(), // Send
            procTimes.size_bytes(),  // Num send per proc
            allTimes.data_bytes(),   // Recv
            procTimes.size_bytes(),  // Num recv per proc
            UPstream::commWorld()
        );
    }

    if (detailLevel_ > 1)
    {
        // The counts
        const auto& procCounts = profilingPstream::counts();

        if (Pstream::master())
        {
            allCounts.resize(numProc * procCounts.size());
        }

        UPstream::mpiGather
        (
            procCounts.cdata_bytes(), // Send
            procCounts.size_bytes(),  // Num send per proc
            allCounts.data_bytes(),   // Recv
            procCounts.size_bytes(),  // Num recv per proc
            UPstream::commWorld()
        );
    }

    profilingPstream::resume();


    // (Time, Processor) for each of: min/max/sum(avg)
    typedef FixedList<Tuple2<double, int>, 3> statData;

    // Extract min/max/average
    auto calcStats = [](const UList<double>& data) -> statData
    {
        statData stats;
        stats = Tuple2<double, int>((data.empty() ? 0 : data[0]), 0);

        const label np = data.size();
        for (label proci = 1; proci < np; ++proci)
        {
            Tuple2<double, int> tup(data[proci], proci);

            // 0: min, 1: max, 2: total(avg)
            if (stats[0].first() > tup.first()) stats[0] = tup;
            if (stats[1].first() < tup.first()) stats[1] = tup;
            stats[2].first() += tup.first();
        }

        // From total -> average value
        if (np) { stats[2].first() /= np; }

        return stats;
    };


    const auto printTimingStats =
        [&](Ostream& os, const char* tag, const statData& stats)
        {
            os  << indent << tag << ": avg = " << stats[2].first()
                << ", min = " << stats[0].first()
                << " (proc " << stats[0].second() << ')'
                << ", max = " << stats[1].first()
                << " (proc " << stats[1].second() << ')'
                << nl;
        };


    if (Pstream::master())
    {
        statData stats;
        List<double> extractedTimes(numProc);
        List<uint64_t> extractedCounts;

        if (detailLevel_ > 1)
        {
            extractedCounts.resize(numProc);
        }

        Info<< type() << ':' << nl
            << incrIndent;

        // Total times
        {
            extractValues
            (
                extractedTimes,
                allTimes,
                [=](const double values[])
                {
                    double total = 0;
                    for (unsigned i = 0; i < profilingPstream::nCategories; ++i)
                    {
                        total += values[i];
                    }
                    return total;
                }
            );
            stats = calcStats(extractedTimes);

            printTimingStats(Info(), "total     ", stats);
            if (detailLevel_ > 0) printTimingDetail(extractedTimes);
        }

        // all-all
        {
            const int index = int(profilingPstream::ALL_TO_ALL);

            extractValues(extractedTimes, index, allTimes);
            extractValues(extractedCounts, index, allCounts);
            stats = calcStats(extractedTimes);

            printTimingStats(Info(), "all-all   ", stats);
            if (detailLevel_ > 0) printTimingDetail(extractedTimes);
            if (detailLevel_ > 1) printTimingDetail(extractedCounts);
        }

        // broadcast
        {
            const int index = int(profilingPstream::BROADCAST);

            extractValues(extractedTimes, index, allTimes);
            extractValues(extractedCounts, index, allCounts);
            stats = calcStats(extractedTimes);

            printTimingStats(Info(), "broadcast ", stats);
            if (detailLevel_ > 0) printTimingDetail(extractedTimes);
            if (detailLevel_ > 1) printTimingDetail(extractedCounts);
        }

        // probe
        {
            const int index = int(profilingPstream::PROBE);

            extractValues(extractedTimes, index, allTimes);
            extractValues(extractedCounts, index, allCounts);
            stats = calcStats(extractedTimes);

            printTimingStats(Info(), "probe     ", stats);
            if (detailLevel_ > 0) printTimingDetail(extractedTimes);
            if (detailLevel_ > 1) printTimingDetail(extractedCounts);
        }

        // Reduce/scatter times
        {
            // const int index = int(profilingPstream::REDUCE);

            extractValues
            (
                extractedTimes,
                allTimes,
                [=](const double values[])
                {
                    return
                    (
                        values[profilingPstream::REDUCE]
                      + values[profilingPstream::GATHER]
                      + values[profilingPstream::SCATTER]
                    );
                }
            );
            extractValues
            (
                extractedCounts,
                allCounts,
                [=](const uint64_t values[])
                {
                    return
                    (
                        values[profilingPstream::REDUCE]
                      + values[profilingPstream::GATHER]
                      + values[profilingPstream::SCATTER]
                    );
                }
            );
            stats = calcStats(extractedTimes);

            printTimingStats(Info(), "reduce    ", stats);
            if (detailLevel_ > 0) printTimingDetail(extractedTimes);
            if (detailLevel_ > 1) printTimingDetail(extractedCounts);
        }

        // request
        {
            const int index = int(profilingPstream::REQUEST);

            extractValues(extractedTimes, index, allTimes);
            extractValues(extractedCounts, index, allCounts);
            stats = calcStats(extractedTimes);

            printTimingStats(Info(), "request   ", stats);

            if (detailLevel_ > 0) printTimingDetail(extractedTimes);
            if (detailLevel_ > 1) printTimingDetail(extractedCounts);
        }

        // wait
        {
            const int index = int(profilingPstream::WAIT);

            extractValues(extractedTimes, index, allTimes);
            extractValues(extractedCounts, index, allCounts);
            stats = calcStats(extractedTimes);

            printTimingStats(Info(), "wait      ", stats);

            if (detailLevel_ > 0) printTimingDetail(extractedTimes);
            if (detailLevel_ > 1) printTimingDetail(extractedCounts);
        }

        Info<< decrIndent;
    }
}


bool Foam::functionObjects::parProfiling::execute()
{
    report();
    return true;
}


bool Foam::functionObjects::parProfiling::write()
{
    return true;
}


bool Foam::functionObjects::parProfiling::end()
{
    profilingPstream::disable();
    return true;
}


// ************************************************************************* //

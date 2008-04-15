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


#include "Pstream.H"
#include "ListListOps.H"
#include "combineSampleSets.H"
#include "SortableList.H"
#include "IndirectList.H"

using namespace Foam;

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

// Combine points from all processors. Sort by curveDist and produce
// index list. Valid result only on master processor.
void combineSampleSets
(
    const PtrList<sampleSet>& sampleSets,
    PtrList<coordSet>& masterSampleSets,
    labelListList& indexSets
)
{
    forAll(sampleSets, setI)
    {
        const sampleSet& samplePts = sampleSets[setI];

        // Collect data from all processors
        List<List<point> > gatheredPts(Pstream::nProcs());
        gatheredPts[Pstream::myProcNo()] = samplePts;
        Pstream::gatherList(gatheredPts);

        List<labelList> gatheredSegments(Pstream::nProcs());
        gatheredSegments[Pstream::myProcNo()] = samplePts.segments();
        Pstream::gatherList(gatheredSegments);

        List<scalarList> gatheredDist(Pstream::nProcs());
        gatheredDist[Pstream::myProcNo()] = samplePts.curveDist();
        Pstream::gatherList(gatheredDist);


        // Combine processor lists into one big list.
        List<point> allPts
        (
            ListListOps::combine<List<point> >
            (
                gatheredPts, accessOp<List<point> >()
            )
        );
        labelList allSegments
        (
            ListListOps::combine<labelList>
            (
                gatheredSegments, accessOp<labelList>()
            )
        );
        scalarList allCurveDist
        (
            ListListOps::combine<scalarList>
            (
                gatheredDist, accessOp<scalarList>()
            )
        );

        // Sort curveDist and use to fill masterSamplePts
        SortableList<scalar> sortedDist(allCurveDist);
        indexSets[setI] = sortedDist.indices();

        // Get reference point (note: only master has all points)
        point refPt;

        if (allPts.size() > 0)
        {
            refPt = samplePts.getRefPoint(allPts);
        }
        else
        {
            refPt = vector::zero;
        }


        masterSampleSets.set
        (
            setI,
            new coordSet
            (
                samplePts.name(),
                samplePts.axis(),
                IndirectList<point>(allPts, indexSets[setI]),
                refPt
            )
        );
    }
}


// ************************************************************************* //

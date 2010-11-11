/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2010 OpenCFD Ltd.
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

#include "hierarchicalDensityWeightedStochastic.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(hierarchicalDensityWeightedStochastic, 0);
addToRunTimeSelectionTable
(
    initialPointsMethod,
    hierarchicalDensityWeightedStochastic,
    dictionary
);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::hierarchicalDensityWeightedStochastic::writeOBJ
(
    const treeBoundBox& bb,
    fileName name
) const
{
    OFstream str(name + ".obj");

    Info<< "Writing " << str.name() << endl;

    label vertI = 0;

    pointField bbPoints(bb.points());

    label pointVertI = vertI;

    forAll(bbPoints, i)
    {
        meshTools::writeOBJ(str, bbPoints[i]);
        vertI++;
    }

    forAll(treeBoundBox::edges, i)
    {
        const edge& e = treeBoundBox::edges[i];

        str << "l " << e[0] + pointVertI + 1
            << ' ' << e[1] + pointVertI + 1
            << nl;
    }
}


void Foam::hierarchicalDensityWeightedStochastic::recurseAndFill
(
    std::vector<Vb::Point>& initialPoints,
    const treeBoundBox& bb,
    label levelLimit,
    word recursionName
) const
{
    const conformationSurfaces& geometry = cvMesh_.geometryToConformTo();

    for (direction i = 0; i < 8; i++)
    {
        treeBoundBox subBB = bb.subBbox(i);

        word newName = recursionName + "_" + Foam::name(i);

        if (geometry.overlaps(subBB))
        {
            if (levelLimit > 0)
            {
                recurseAndFill
                (
                    initialPoints,
                    subBB,
                    levelLimit - 1,
                    newName
                );
            }
            else
            {
                // writeOBJ
                // (
                //     subBB,
                //     word(newName + "_overlap")
                // );

                if (!fillBox(initialPoints, subBB, true))
                {
                    recurseAndFill
                    (
                        initialPoints,
                        subBB,
                        levelLimit - 1,
                        newName
                    );
                }
            }
        }
        else if (geometry.inside(subBB.midpoint()))
        {
            // writeOBJ
            // (
            //     subBB,
            //     newName + "_inside"
            // );

            if (!fillBox(initialPoints, subBB, false))
            {
                recurseAndFill
                (
                    initialPoints,
                    subBB,
                    levelLimit - 1,
                    newName
                );
            }
        }
        else
        {
            // writeOBJ
            // (
            //     subBB,
            //     newName + "_outside"
            // );
        }
    }
}


bool Foam::hierarchicalDensityWeightedStochastic::fillBox
(
    std::vector<Vb::Point>& initialPoints,
    const treeBoundBox& bb,
    bool overlapping
) const
{
    const conformationSurfaces& geometry = cvMesh_.geometryToConformTo();

    Random& rnd = cvMesh_.rndGen();

    unsigned int initialSize = initialPoints.size();

    scalar maxCellSize = -GREAT;

    scalar minCellSize = GREAT;

    scalar maxDensity = 1/pow3(minCellSize);

    scalar volumeAdded = 0.0;

    const point& min = bb.min();

    vector span = bb.span();

    scalar totalVolume = bb.volume();

    label trialPoints = 0;

    bool wellInside = false;

    if (!overlapping)
    {
        // Check the nearest point on the surface to the box, if it is far
        // enough away, then the surface sampling of the box can be skipped.
        // Checking if the nearest piece of surface is at least 1.5*bb.span away
        // from the bb.midpoint.

        pointIndexHit surfHit;
        label hitSurface;

        geometry.findSurfaceNearest
        (
            bb.midpoint(),
            2.25*magSqr(span),
            surfHit,
            hitSurface
        );

        if (surfHit.hit())
        {
            // Info<< "box wellInside, no need to sample surface." << endl;

            wellInside = true;
        }
    }

    if (!overlapping && !wellInside)
    {
        // If this is an inside box then then it is possible to fill points very
        // close to the boundary, to prevent this, check the corners and sides
        // of the box so ensure that they are "wellInside".  If not, set as an
        // overlapping box.

        pointField corners(bb.points());

        scalarField cornerSizes(8, 0.0);

        Field<bool> insideCorners(8, true);

        cornerSizes = cvMesh_.cellSizeControl().cellSize
        (
            corners,
            List<bool>(8, false)
        );

        insideCorners = geometry.wellInside
        (
            corners,
            minimumSurfaceDistanceCoeffSqr_*sqr(cornerSizes)
        );

        forAll(insideCorners, i)
        {
            // Use the sizes to improve the min/max cell size estimate
            scalar& s = cornerSizes[i];

            if (s > maxCellSize)
            {
                maxCellSize = s;
            }

            if (s < minCellSize)
            {
                minCellSize = max(s, minCellSizeLimit_);
            }

            if (maxCellSize/minCellSize > maxSizeRatio_)
            {
                // Info<< "Abort fill at corner sample stage, maxSizeRatio = "
                //     << maxCellSize/minCellSize << endl;

                return false;
            }

            if (!insideCorners[i])
            {
                // If one or more corners is not "wellInside", then treat this
                // as an overlapping box.

                // Info<< "Inside box found to have some non-wellInside "
                //     << "corners, using overlapping fill."
                //     << endl;

                overlapping = true;

                break;
            }
        }

        if (!overlapping)
        {
            vector delta = span/(surfRes_ - 1);

            label nLine = 6*(surfRes_ - 2);

            pointField linePoints(nLine, vector::zero);

            scalarField lineSizes(nLine, 0.0);

            Field<bool> insideLines(nLine, true);

            for (label i = 0; i < surfRes_; i++)
            {
                label lPI = 0;

                for (label j = 1; j < surfRes_ - 1 ; j++)
                {
                    linePoints[lPI++] =
                        min
                      + vector(0, delta.y()*i, delta.z()*j);

                    linePoints[lPI++] =
                        min
                      + vector
                        (
                            delta.x()*(surfRes_ - 1),
                            delta.y()*i,
                            delta.z()*j
                        );

                    linePoints[lPI++] =
                        min
                      + vector(delta.x()*j, 0, delta.z()*i);

                    linePoints[lPI++] =
                        min
                      + vector
                        (
                            delta.x()*j,
                            delta.y()*(surfRes_ - 1),
                            delta.z()*i
                        );

                    linePoints[lPI++] =
                        min
                      + vector(delta.x()*i, delta.y()*j, 0);

                    linePoints[lPI++] =
                        min
                      + vector
                        (
                            delta.x()*i,
                            delta.y()*j,
                            delta.z()*(surfRes_ - 1)
                        );
                }

                lineSizes = cvMesh_.cellSizeControl().cellSize
                (
                    linePoints,
                    List<bool>(nLine, false)
                );

                insideLines = geometry.wellInside
                (
                    linePoints,
                    minimumSurfaceDistanceCoeffSqr_*sqr(lineSizes)
                );

                forAll(insideLines, i)
                {
                    // Use the sizes to improve the min/max cell size estimate
                    scalar& s = lineSizes[i];

                    if (s > maxCellSize)
                    {
                        maxCellSize = s;
                    }

                    if (s < minCellSize)
                    {
                        minCellSize = max(s, minCellSizeLimit_);
                    }

                    if (maxCellSize/minCellSize > maxSizeRatio_)
                    {
                        // Info<< "Abort fill at surface sample stage, "
                        //     << "maxSizeRatio = "
                        //     << maxCellSize/minCellSize << endl;

                        return false;
                    }

                    if (!insideLines[i])
                    {
                        // If one or more surface points is not "wellInside",
                        // then treat this as an overlapping box.
                        overlapping = true;

                        // Info<< "Inside box found to have some non-"
                        //     << "wellInside surface points, using "
                        //     << "overlapping fill."
                        //     << endl;

                        break;
                    }
                }
            }
        }
    }

    if (overlapping)
    {
        // Sample the box to find an estimate of the min size, and a volume
        // estimate when overlapping == true.

        pointField samplePoints
        (
            volRes_*volRes_*volRes_,
            vector::zero
        );

        vector delta = span/volRes_;

        label pI = 0;

        for (label i = 0; i < volRes_; i++)
        {
            for (label j = 0; j < volRes_; j++)
            {
                for (label k = 0; k < volRes_; k++)
                {
                    // Perturb the points to avoid creating degenerate positions
                    // in the Delaunay tessellation.

                    samplePoints[pI++] =
                        min
                      + vector
                        (
                            delta.x()*(i + 0.5 + 0.1*(rnd.scalar01() - 0.5)),
                            delta.y()*(j + 0.5 + 0.1*(rnd.scalar01() - 0.5)),
                            delta.z()*(k + 0.5 + 0.1*(rnd.scalar01() - 0.5))
                        );
                }
            }
        }

        // Randomise the order of the points to (potentially) improve the speed
        // of assessing the density ratio, and prevent a box being filled from a
        // corner when only some these points are required.
        shuffle(samplePoints);

        scalarField sampleSizes = cvMesh_.cellSizeControl().cellSize
        (
            samplePoints,
            List<bool>(samplePoints.size(), false)
        );

        Field<bool> insidePoints = geometry.wellInside
        (
            samplePoints,
            minimumSurfaceDistanceCoeffSqr_*sqr(sampleSizes)
        );

        label nInside = 0;

        forAll(insidePoints, i)
        {
            if (insidePoints[i])
            {
                nInside++;

                scalar& s = sampleSizes[i];

                if (s > maxCellSize)
                {
                    maxCellSize = s;
                }

                if (s < minCellSize)
                {
                    minCellSize = max(s, minCellSizeLimit_);
                }

                if (maxCellSize/minCellSize > maxSizeRatio_)
                {
                    // Info<< "Abort fill at sample stage, maxSizeRatio = "
                    //     << maxCellSize/minCellSize << endl;

                    return false;
                }
            }
        }

        if (nInside == 0)
        {
            // Info<< "No sample points found inside box" << endl;

            return true;
        }

        // Info<< scalar(nInside)/scalar(samplePoints.size())
        //     << " full overlapping box" << endl;

        totalVolume *= scalar(nInside)/scalar(samplePoints.size());

        // Info<< "Total volume to fill = " << totalVolume << endl;

        // Using the sampledPoints as the first test locations as they are
        // randomly shuffled, but unfiormly sampling space and have wellInside
        // and size data already

        maxDensity = 1/pow3(max(minCellSize, SMALL));

        forAll(insidePoints, i)
        {
            if (insidePoints[i])
            {
                trialPoints++;

                point p = samplePoints[i];

                scalar localSize = sampleSizes[i];

                scalar localDensity = 1/pow3(localSize);

                // No need to look at max/min cell size here, already handled
                // by sampling

                // Accept possible placements proportional to the relative
                // local density

                // TODO - is there a lot of cost in the 1/density calc?  Could
                // assess on
                //    (1/maxDensity)/(1/localDensity) = minVolume/localVolume
                if (localDensity/maxDensity > rnd.scalar01())
                {
                    scalar localVolume = 1/localDensity;

                    if (volumeAdded + localVolume > totalVolume)
                    {
                        // Add the final box with a probability of to the ratio
                        // of the remaining volume to the volume to be added,
                        // i.e. insert a box of volume 0.5 into a remaining
                        // volume of 0.1 20% of the time.
                        scalar addProbability =
                           (totalVolume - volumeAdded)/localVolume;

                        scalar r = rnd.scalar01();

                        // Info<< "totalVolume " << totalVolume << nl
                        //     << "volumeAdded " << volumeAdded << nl
                        //     << "localVolume " << localVolume << nl
                        //     << "addProbability " << addProbability << nl
                        //     << "random " << r
                        //     << endl;

                        if (addProbability > r)
                        {
                            // Do not place this volume, finished filling this
                            // box

                            // Info<< "Final volume probability break accept"
                            //     << endl;

                            initialPoints.push_back
                            (
                                Vb::Point(p.x(), p.y(), p.z())
                            );

                            volumeAdded += localVolume;
                        }

                        break;
                    }

                    initialPoints.push_back(Vb::Point(p.x(), p.y(), p.z()));

                    volumeAdded += localVolume;
                }
            }
        }
    }

    if (volumeAdded < totalVolume)
    {
        // Info<< "Adding random points, remaining volume "
        //     << totalVolume - volumeAdded
        //     << endl;

        maxDensity = 1/pow3(max(minCellSize, SMALL));

        while (true)
        {
            trialPoints++;

            point p =
            min
            + vector
            (
                span.x()*rnd.scalar01(),
                span.y()*rnd.scalar01(),
                span.z()*rnd.scalar01()
            );

            scalar localSize = cvMesh_.cellSizeControl().cellSize(p);

            if (localSize > maxCellSize)
            {
                maxCellSize = localSize;
            }

            if (localSize < minCellSize)
            {
                minCellSize = max(localSize, minCellSizeLimit_);

                localSize = minCellSize;

                // 1/(minimum cell size)^3, gives the maximum permissible point
                // density
                maxDensity = 1/pow3(max(minCellSize, SMALL));
            }

            if (maxCellSize/minCellSize > maxSizeRatio_)
            {
                // Info<< "Abort fill at random fill stage, maxSizeRatio = "
                //     << maxCellSize/minCellSize << endl;

                // Discard any points already filled into this box by setting
                // size of initialPoints back to its starting value
                initialPoints.resize(initialSize);

                return false;
            }

            scalar localDensity = 1/pow3(max(localSize, SMALL));

            // Accept possible placements proportional to the relative local
            // density
            if (localDensity/maxDensity > rnd.scalar01())
            {
                bool insidePoint = false;

                if (!overlapping)
                {
                    insidePoint = true;
                }
                else
                {
                    // Determine if the point is "wellInside" the domain
                    insidePoint = geometry.wellInside
                    (
                        p,
                        minimumSurfaceDistanceCoeffSqr_*sqr(localSize)
                    );
                }

                if (insidePoint)
                {
                    scalar localVolume = 1/localDensity;

                    if (volumeAdded + localVolume > totalVolume)
                    {
                        // Add the final box with a probability of to the ratio
                        // of the remaining volume to the volume to be added,
                        // i.e. insert a box of volume 0.5 into a remaining
                        // volume of 0.1 20% of the time.
                        scalar addProbability =
                            (totalVolume - volumeAdded)/localVolume;

                        scalar r = rnd.scalar01();

                        // Info<< "totalVolume " << totalVolume << nl
                        //     << "volumeAdded " << volumeAdded << nl
                        //     << "localVolume " << localVolume << nl
                        //     << "addProbability " << addProbability << nl
                        //     << "random " << r
                        //     << endl;

                        if (addProbability > r)
                        {
                            // Do not place this volume, finished filling this
                            // box

                            // Info<< "Final volume probability break accept"
                            //     << endl;

                            initialPoints.push_back
                            (
                                Vb::Point(p.x(), p.y(), p.z())
                            );

                            volumeAdded += localVolume;
                        }

                        break;
                    }

                    initialPoints.push_back(Vb::Point(p.x(), p.y(), p.z()));

                    volumeAdded += localVolume;
                }
            }
        }
    }

    globalTrialPoints_ += trialPoints;

    // Info<< trialPoints
    //     << " locations queried, " << initialPoints.size() - initialSize
    //     << " points placed, ("
    //     << scalar(initialPoints.size() - initialSize)/scalar(trialPoints)
    //     << " success rate)." << nl
    //     << "minCellSize " << minCellSize
    //     << ", maxCellSize " << maxCellSize
    //     << ", ratio " << maxCellSize/minCellSize
    //     << nl << endl;

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

hierarchicalDensityWeightedStochastic::hierarchicalDensityWeightedStochastic
(
    const dictionary& initialPointsDict,
    const conformalVoronoiMesh& cvMesh
)
:
    initialPointsMethod(typeName, initialPointsDict, cvMesh),
    globalTrialPoints_(0),
    minCellSizeLimit_
    (
        detailsDict().lookupOrDefault<scalar>("minCellSizeLimit", 0.0)
    ),
    maxLevels_(readLabel(detailsDict().lookup("maxLevels"))),
    maxSizeRatio_(readScalar(detailsDict().lookup("maxSizeRatio"))),
    volRes_(readLabel(detailsDict().lookup("sampleResolution"))),
    surfRes_
    (
        detailsDict().lookupOrDefault<label>("surfaceSampleResolution", volRes_)
    )
{
    if (maxSizeRatio_ <= 1.0)
    {
        maxSizeRatio_ = 2.0;

        WarningIn
        (
            "hierarchicalDensityWeightedStochastic::"
            "hierarchicalDensityWeightedStochastic"
            "("
                "const dictionary& initialPointsDict,"
                "const conformalVoronoiMesh& cvMesh"
            ")"
        )
            << "The maxSizeRatio must be greater than one to be sensible, "
            << "setting to " << maxSizeRatio_
            << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

std::vector<Vb::Point>
hierarchicalDensityWeightedStochastic::initialPoints() const
{
    const conformationSurfaces& geometry = cvMesh_.geometryToConformTo();

    treeBoundBox hierBB = treeBoundBox(geometry.bounds()).extend
    (
        cvMesh_.rndGen(),
        1E-4
    );

    std::vector<Vb::Point> initialPoints;

    recurseAndFill
    (
        initialPoints,
        hierBB,
        maxLevels_,
        "recursionBox"
    );

    Info<< nl << "    " << typeName << nl
        << "        " << initialPoints.size() << " points placed" << nl
        << "        " << globalTrialPoints_ << " locations queried" << nl
        << "        " << scalar(initialPoints.size())/scalar(globalTrialPoints_)
        << " success rate"
        << endl;

    return initialPoints;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

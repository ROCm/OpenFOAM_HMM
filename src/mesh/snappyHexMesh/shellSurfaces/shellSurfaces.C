/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015-2017 OpenCFD Ltd.
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

#include "shellSurfaces.H"
#include "searchableSurface.H"
#include "boundBox.H"
#include "triSurfaceMesh.H"
#include "refinementSurfaces.H"
#include "searchableSurfaces.H"
#include "orientedSurface.H"
#include "pointIndexHit.H"
#include "volumeType.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::shellSurfaces::refineMode
>
Foam::shellSurfaces::refineModeNames_
{
    { refineMode::INSIDE, "inside" },
    { refineMode::OUTSIDE, "outside" },
    { refineMode::DISTANCE, "distance" },
};


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::shellSurfaces::setAndCheckLevels
(
    const label shellI,
    const List<Tuple2<scalar, label>>& distLevels
)
{
    if (modes_[shellI] != DISTANCE && distLevels.size() != 1)
    {
        FatalErrorInFunction
            << "For refinement mode "
            << refineModeNames_[modes_[shellI]]
            << " specify only one distance+level."
            << " (its distance gets discarded)"
            << exit(FatalError);
    }
    // Extract information into separate distance and level
    distances_[shellI].setSize(distLevels.size());
    levels_[shellI].setSize(distLevels.size());

    forAll(distLevels, j)
    {
        distances_[shellI][j] = distLevels[j].first();
        levels_[shellI][j] = distLevels[j].second();

        // Check in incremental order
        if (j > 0)
        {
            if
            (
                (distances_[shellI][j] <= distances_[shellI][j-1])
             || (levels_[shellI][j] > levels_[shellI][j-1])
            )
            {
                FatalErrorInFunction
                    << "For refinement mode "
                    << refineModeNames_[modes_[shellI]]
                    << " : Refinement should be specified in order"
                    << " of increasing distance"
                    << " (and decreasing refinement level)." << endl
                    << "Distance:" << distances_[shellI][j]
                    << " refinementLevel:" << levels_[shellI][j]
                    << exit(FatalError);
            }
        }
    }

    const searchableSurface& shell = allGeometry_[shells_[shellI]];

    if (modes_[shellI] == DISTANCE)
    {
        Info<< "Refinement level according to distance to "
            << shell.name() << endl;
        forAll(levels_[shellI], j)
        {
            Info<< "    level " << levels_[shellI][j]
                << " for all cells within " << distances_[shellI][j]
                << " metre." << endl;
        }
    }
    else
    {
        if (!shell.hasVolumeType())
        {
            FatalErrorInFunction
                << "Shell " << shell.name()
                << " does not support testing for "
                << refineModeNames_[modes_[shellI]] << endl
                << "Probably it is not closed."
                << exit(FatalError);
        }

        if (modes_[shellI] == INSIDE)
        {
            Info<< "Refinement level " << levels_[shellI][0]
                << " for all cells inside " << shell.name() << endl;
        }
        else
        {
            Info<< "Refinement level " << levels_[shellI][0]
                << " for all cells outside " << shell.name() << endl;
        }
    }
}


void Foam::shellSurfaces::checkGapLevels
(
    const dictionary& shellDict,
    const label shellI,
    const List<FixedList<label, 3>>& levels
)
{
    const searchableSurface& shell = allGeometry_[shells_[shellI]];

    forAll(levels, regionI)
    {
        const FixedList<label, 3>& info = levels[regionI];

        if (info[2] > 0)
        {
            if (modes_[shellI] == DISTANCE)
            {
                FatalIOErrorInFunction(shellDict)
                    << "'gapLevel' specification cannot be used with mode "
                    << refineModeNames_[DISTANCE]
                    << " for shell " << shell.name()
                    << exit(FatalIOError);
            }
        }
    }

    // Hardcode for region 0
    if (levels[0][0] > 0)
    {
        Info<< "Refinement level up to " << levels[0][2]
            << " for all cells in gaps for shell "
            << shell.name() << endl;
    }
}



// Specifically orient triSurfaces using a calculated point outside.
// Done since quite often triSurfaces not of consistent orientation which
// is (currently) necessary for sideness calculation
void Foam::shellSurfaces::orient()
{
    // Determine outside point.
    boundBox overallBb = boundBox::invertedBox;

    bool hasSurface = false;

    forAll(shells_, shellI)
    {
        const searchableSurface& s = allGeometry_[shells_[shellI]];

        if (modes_[shellI] != DISTANCE && isA<triSurfaceMesh>(s))
        {
            const triSurfaceMesh& shell = refCast<const triSurfaceMesh>(s);

            if (shell.triSurface::size())
            {
                hasSurface = true;
                // Assume surface is compact!
                overallBb.add(shell.points());
            }
        }
    }

    if (hasSurface)
    {
        const point outsidePt = overallBb.max() + overallBb.span();

        //Info<< "Using point " << outsidePt << " to orient shells" << endl;

        forAll(shells_, shellI)
        {
            const searchableSurface& s = allGeometry_[shells_[shellI]];

            if (modes_[shellI] != DISTANCE && isA<triSurfaceMesh>(s))
            {
                triSurfaceMesh& shell = const_cast<triSurfaceMesh&>
                (
                    refCast<const triSurfaceMesh>(s)
                );

                // Flip surface so outsidePt is outside.
                bool anyFlipped = orientedSurface::orient
                (
                    shell,
                    outsidePt,
                    true
                );

                if (anyFlipped)
                {
                    // orientedSurface will have done a clearOut of the surface.
                    // we could do a clearout of the triSurfaceMeshes::trees()
                    // but these aren't affected by orientation
                    // (except for cached
                    // sideness which should not be set at this point.
                    // !!Should check!)

                    Info<< "shellSurfaces : Flipped orientation of surface "
                        << s.name()
                        << " so point " << outsidePt << " is outside." << endl;
                }
            }
        }
    }
}


// Find maximum level of a shell.
void Foam::shellSurfaces::findHigherLevel
(
    const pointField& pt,
    const label shellI,
    labelList& maxLevel
) const
{
    const labelList& levels = levels_[shellI];

    if (modes_[shellI] == DISTANCE)
    {
        // Distance mode.

        const scalarField& distances = distances_[shellI];

        // Collect all those points that have a current maxLevel less than
        // (any of) the shell. Also collect the furthest distance allowable
        // to any shell with a higher level.

        pointField candidates(pt.size());
        labelList candidateMap(pt.size());
        scalarField candidateDistSqr(pt.size());
        label candidateI = 0;

        forAll(maxLevel, pointi)
        {
            forAllReverse(levels, levelI)
            {
                if (levels[levelI] > maxLevel[pointi])
                {
                    candidates[candidateI] = pt[pointi];
                    candidateMap[candidateI] = pointi;
                    candidateDistSqr[candidateI] = sqr(distances[levelI]);
                    candidateI++;
                    break;
                }
            }
        }
        candidates.setSize(candidateI);
        candidateMap.setSize(candidateI);
        candidateDistSqr.setSize(candidateI);

        // Do the expensive nearest test only for the candidate points.
        List<pointIndexHit> nearInfo;
        allGeometry_[shells_[shellI]].findNearest
        (
            candidates,
            candidateDistSqr,
            nearInfo
        );

        // Update maxLevel
        forAll(nearInfo, i)
        {
            if (nearInfo[i].hit())
            {
                // Check which level it actually is in.
                label minDistI = findLower
                (
                    distances,
                    mag(nearInfo[i].hitPoint()-candidates[i])
                );

                label pointi = candidateMap[i];

                // pt is inbetween shell[minDistI] and shell[minDistI+1]
                maxLevel[pointi] = levels[minDistI+1];
            }
        }
    }
    else
    {
        // Inside/outside mode

        // Collect all those points that have a current maxLevel less than the
        // shell.

        pointField candidates(pt.size());
        labelList candidateMap(pt.size());
        label candidateI = 0;

        forAll(maxLevel, pointi)
        {
            if (levels[0] > maxLevel[pointi])
            {
                candidates[candidateI] = pt[pointi];
                candidateMap[candidateI] = pointi;
                candidateI++;
            }
        }
        candidates.setSize(candidateI);
        candidateMap.setSize(candidateI);

        // Do the expensive nearest test only for the candidate points.
        List<volumeType> volType;
        allGeometry_[shells_[shellI]].getVolumeType(candidates, volType);

        forAll(volType, i)
        {
            label pointi = candidateMap[i];

            if
            (
                (
                    modes_[shellI] == INSIDE
                 && volType[i] == volumeType::INSIDE
                )
             || (
                    modes_[shellI] == OUTSIDE
                 && volType[i] == volumeType::OUTSIDE
                )
            )
            {
                maxLevel[pointi] = levels[0];
            }
        }
    }
}


void Foam::shellSurfaces::findHigherGapLevel
(
    const pointField& pt,
    const labelList& ptLevel,
    const label shellI,
    labelList& gapShell,
    List<FixedList<label, 3>>& gapInfo,
    List<volumeType>& gapMode
) const
{
    //TBD: hardcoded for region 0 information
    const FixedList<label, 3>& info = extendedGapLevel_[shellI][0];
    volumeType mode = extendedGapMode_[shellI][0];

    if (info[2] == 0)
    {
        return;
    }


    // Collect all those points that have a current maxLevel less than the
    // shell.

    labelList candidateMap(pt.size());
    label candidateI = 0;

    forAll(ptLevel, pointI)
    {
        if (ptLevel[pointI] >= info[1] && ptLevel[pointI] < info[2])
        {
            candidateMap[candidateI++] = pointI;
        }
    }
    candidateMap.setSize(candidateI);

    // Do the expensive nearest test only for the candidate points.
    List<volumeType> volType;
    allGeometry_[shells_[shellI]].getVolumeType
    (
        pointField(pt, candidateMap),
        volType
    );

    forAll(volType, i)
    {
        label pointI = candidateMap[i];

        bool isInside = (volType[i] == volumeType::INSIDE);

        if
        (
            (
                (modes_[shellI] == INSIDE && isInside)
             || (modes_[shellI] == OUTSIDE && !isInside)
            )
         && info[2] > gapInfo[pointI][2]
        )
        {
            gapShell[pointI] = shellI;
            gapInfo[pointI] = info;
            gapMode[pointI] = mode;
        }
    }
}


void Foam::shellSurfaces::findLevel
(
    const pointField& pt,
    const label shellI,
    labelList& minLevel,
    labelList& shell
) const
{
    const labelList& levels = levels_[shellI];

    if (modes_[shellI] == DISTANCE)
    {
        // Distance mode.

        const scalarField& distances = distances_[shellI];

        // Collect all those points that have a current level equal/greater
        // (any of) the shell. Also collect the furthest distance allowable
        // to any shell with a higher level.

        pointField candidates(pt.size());
        labelList candidateMap(pt.size());
        scalarField candidateDistSqr(pt.size());
        label candidateI = 0;

        forAll(shell, pointI)
        {
            if (shell[pointI] == -1)
            {
                forAllReverse(levels, levelI)
                {
                    if (levels[levelI] <= minLevel[pointI])
                    {
                        candidates[candidateI] = pt[pointI];
                        candidateMap[candidateI] = pointI;
                        candidateDistSqr[candidateI] = sqr(distances[levelI]);
                        candidateI++;
                        break;
                    }
                }
            }
        }
        candidates.setSize(candidateI);
        candidateMap.setSize(candidateI);
        candidateDistSqr.setSize(candidateI);

        // Do the expensive nearest test only for the candidate points.
        List<pointIndexHit> nearInfo;
        allGeometry_[shells_[shellI]].findNearest
        (
            candidates,
            candidateDistSqr,
            nearInfo
        );

        // Update maxLevel
        forAll(nearInfo, i)
        {
            if (nearInfo[i].hit())
            {
                // Check which level it actually is in.
                label minDistI = findLower
                (
                    distances,
                    mag(nearInfo[i].hitPoint()-candidates[i])
                );

                label pointI = candidateMap[i];

                // pt is inbetween shell[minDistI] and shell[minDistI+1]
                shell[pointI] = shellI;
                minLevel[pointI] = levels[minDistI+1];
            }
        }
    }
    else
    {
        // Inside/outside mode

        // Collect all those points that have a current maxLevel less than the
        // shell.

        pointField candidates(pt.size());
        labelList candidateMap(pt.size());
        label candidateI = 0;

        forAll(shell, pointI)
        {
            if (shell[pointI] == -1 && levels[0] <= minLevel[pointI])
            {
                candidates[candidateI] = pt[pointI];
                candidateMap[candidateI] = pointI;
                candidateI++;
            }
        }
        candidates.setSize(candidateI);
        candidateMap.setSize(candidateI);

        // Do the expensive nearest test only for the candidate points.
        List<volumeType> volType;
        allGeometry_[shells_[shellI]].getVolumeType(candidates, volType);

        forAll(volType, i)
        {
            if
            (
                (
                    modes_[shellI] == INSIDE
                 && volType[i] == volumeType::INSIDE
                )
             || (
                    modes_[shellI] == OUTSIDE
                 && volType[i] == volumeType::OUTSIDE
                )
            )
            {
                label pointI = candidateMap[i];
                shell[pointI] = shellI;
                minLevel[pointI] = levels[0];
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::shellSurfaces::shellSurfaces
(
    const searchableSurfaces& allGeometry,
    const dictionary& shellsDict
)
:
    allGeometry_(allGeometry)
{
    // Wilcard specification : loop over all surfaces and try to find a match.

    // Count number of shells.
    label shellI = 0;
    forAll(allGeometry.names(), geomI)
    {
        const word& geomName = allGeometry_.names()[geomI];

        if (shellsDict.found(geomName))
        {
            shellI++;
        }
    }


    // Size lists
    shells_.setSize(shellI);
    modes_.setSize(shellI);
    distances_.setSize(shellI);
    levels_.setSize(shellI);

    extendedGapLevel_.setSize(shellI);
    extendedGapMode_.setSize(shellI);

    FixedList<label, 3> nullGapLevel;
    nullGapLevel[0] = 0;
    nullGapLevel[1] = 0;
    nullGapLevel[2] = 0;


    HashSet<word> unmatchedKeys(shellsDict.toc());
    shellI = 0;

    forAll(allGeometry_.names(), geomI)
    {
        const word& geomName = allGeometry_.names()[geomI];

        const entry* ePtr = shellsDict.lookupEntryPtr(geomName, false, true);

        if (ePtr)
        {
            const dictionary& dict = ePtr->dict();
            unmatchedKeys.erase(ePtr->keyword());

            shells_[shellI] = geomI;
            modes_[shellI] = refineModeNames_.lookup("mode", dict);

            // Read pairs of distance+level
            setAndCheckLevels(shellI, dict.lookup("levels"));



            // Gap specification
            // ~~~~~~~~~~~~~~~~~


            // Shell-wide gap level specification
            const searchableSurface& surface = allGeometry_[geomI];
            const wordList& regionNames = surface.regions();

            FixedList<label, 3> gapSpec
            (
                dict.lookupOrDefault
                (
                    "gapLevel",
                    nullGapLevel
                )
            );
            extendedGapLevel_[shellI].setSize(regionNames.size());
            extendedGapLevel_[shellI] = gapSpec;

            volumeType gapModeSpec
            (
                volumeType::names
                [
                    dict.lookupOrDefault<word>
                    (
                        "gapMode",
                        volumeType::names[volumeType::MIXED]
                    )
                ]
            );
            extendedGapMode_[shellI].setSize(regionNames.size());
            extendedGapMode_[shellI] = gapModeSpec;


            // Override on a per-region basis?

            if (dict.found("regions"))
            {
                const dictionary& regionsDict = dict.subDict("regions");
                forAll(regionNames, regionI)
                {
                    if (regionsDict.found(regionNames[regionI]))
                    {
                        // Get the dictionary for region
                        const dictionary& regionDict = regionsDict.subDict
                        (
                            regionNames[regionI]
                        );
                        FixedList<label, 3> gapSpec
                        (
                            regionDict.lookupOrDefault
                            (
                                "gapLevel",
                                nullGapLevel
                            )
                        );
                        extendedGapLevel_[shellI][regionI] = gapSpec;

                        volumeType gapModeSpec
                        (
                            volumeType::names
                            [
                                regionDict.lookupOrDefault<word>
                                (
                                    "gapMode",
                                    volumeType::names[volumeType::MIXED]
                                )
                            ]
                        );
                        extendedGapMode_[shellI][regionI] = gapModeSpec;
                    }
                }
            }


            checkGapLevels(dict, shellI, extendedGapLevel_[shellI]);

            shellI++;
        }
    }

    if (unmatchedKeys.size() > 0)
    {
        IOWarningInFunction
        (
            shellsDict
        )   << "Not all entries in refinementRegions dictionary were used."
            << " The following entries were not used : "
            << unmatchedKeys.sortedToc()
            << endl;
    }

    // Orient shell surfaces before any searching is done. Note that this
    // only needs to be done for inside or outside. Orienting surfaces
    // constructs lots of addressing which we want to avoid.
    orient();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Highest shell level
Foam::label Foam::shellSurfaces::maxLevel() const
{
    label overallMax = 0;
    forAll(levels_, shellI)
    {
        overallMax = max(overallMax, max(levels_[shellI]));
    }
    return overallMax;
}


Foam::labelList Foam::shellSurfaces::maxGapLevel() const
{
    labelList surfaceMax(extendedGapLevel_.size(), 0);

    forAll(extendedGapLevel_, shelli)
    {
        const List<FixedList<label, 3>>& levels = extendedGapLevel_[shelli];
        forAll(levels, i)
        {
            surfaceMax[shelli] = max(surfaceMax[shelli], levels[i][2]);
        }
    }
    return surfaceMax;
}


void Foam::shellSurfaces::findHigherLevel
(
    const pointField& pt,
    const labelList& ptLevel,
    labelList& maxLevel
) const
{
    // Maximum level of any shell. Start off with level of point.
    maxLevel = ptLevel;

    forAll(shells_, shelli)
    {
        findHigherLevel(pt, shelli, maxLevel);
    }
}


void Foam::shellSurfaces::findHigherGapLevel
(
    const pointField& pt,
    const labelList& ptLevel,
    labelList& gapShell,
    List<FixedList<label, 3>>& gapInfo,
    List<volumeType>& gapMode
) const
{
    gapShell.setSize(pt.size());
    gapShell = -1;

    FixedList<label, 3> nullGapLevel;
    nullGapLevel[0] = 0;
    nullGapLevel[1] = 0;
    nullGapLevel[2] = 0;

    gapInfo.setSize(pt.size());
    gapInfo = nullGapLevel;

    gapMode.setSize(pt.size());
    gapMode = volumeType::MIXED;

    forAll(shells_, shelli)
    {
        findHigherGapLevel(pt, ptLevel, shelli, gapShell, gapInfo, gapMode);
    }
}


void Foam::shellSurfaces::findHigherGapLevel
(
    const pointField& pt,
    const labelList& ptLevel,
    List<FixedList<label, 3>>& gapInfo,
    List<volumeType>& gapMode
) const
{
    labelList gapShell;
    findHigherGapLevel(pt, ptLevel, gapShell, gapInfo, gapMode);
}


void Foam::shellSurfaces::findLevel
(
    const pointField& pt,
    const labelList& ptLevel,
    labelList& shell
) const
{
    shell.setSize(pt.size());
    shell = -1;

    labelList minLevel(ptLevel);

    forAll(shells_, shelli)
    {
        findLevel(pt, shelli, minLevel, shell);
    }
}


// ************************************************************************* //

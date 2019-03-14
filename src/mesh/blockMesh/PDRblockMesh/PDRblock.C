/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
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

#include "PDRblock.H"
#include "ListOps.H"
#include "emptyPolyPatch.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{
    //- Calculate geometric expansion factor from expansion ratio
    inline scalar calcGexp(const scalar expRatio, const label nDiv)
    {
        return nDiv > 1 ? pow(expRatio, 1.0/(nDiv - 1)) : 0.0;
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::PDRblock::readGridControl
(
    const direction cmpt,
    const dictionary& dict,
    const scalar scaleFactor
)
{
    //- The begin/end nodes for each segment
    scalarList knots;

    // The number of division per segment
    labelList divisions;

    // The expansion ratio per segment
    scalarList expansion;

    if (verbose_)
    {
        Info<< "Reading grid control for "
            << vector::componentNames[cmpt] << " direction" << nl;
    }

    // Points
    // ~~~~~~

    dict.readEntry("points", knots);

    const label nSegments = (knots.size()-1);

    if (nSegments < 1)
    {
        FatalIOErrorInFunction(dict)
            << "Must define at least two control points:"
            << " in " << vector::componentNames[cmpt]
            << " direction" << nl
            << exit(FatalIOError);
    }

    // Apply scaling
    if (scaleFactor > 0)
    {
        for (scalar& pt : knots)
        {
            pt *= scaleFactor;
        }
    }

    // Do points increase monotonically?
    {
        const scalar& minVal = knots.first();

        for (label pointi = 1; pointi < knots.size(); ++pointi)
        {
            if (knots[pointi] <= minVal)
            {
                FatalIOErrorInFunction(dict)
                    << "Points are not monotonically increasing:"
                    << " in " << vector::componentNames[cmpt]
                    << " direction" << nl
                    << exit(FatalIOError);
            }
        }
    }


    if (verbose_)
    {
        Info<< "    points : " << flatOutput(knots) << nl;
    }


    // Divisions
    // ~~~~~~~~~

    dict.readEntry("nCells", divisions);

    label nTotalDivs = 0;
    for (const label ndiv : divisions)
    {
        nTotalDivs += ndiv;

        if (ndiv < 1)
        {
            FatalIOErrorInFunction(dict)
                << "Negative or zero divisions:"
                << " in " << vector::componentNames[cmpt]
                << " direction" << nl
                << exit(FatalIOError);
        }
    }

    if (divisions.size() != nSegments)
    {
        FatalIOErrorInFunction(dict)
            << "Mismatch in number of divisions and number of points:"
            << " in " << vector::componentNames[cmpt]
            << " direction" << nl
            << exit(FatalIOError);
    }
    else if (!nTotalDivs)
    {
        FatalIOErrorInFunction(dict)
            << "No grid divisions at all:"
            << " in " << vector::componentNames[cmpt]
            << " direction" << nl
            << exit(FatalIOError);
    }


    if (verbose_)
    {
        Info<< "    nCells : " << flatOutput(divisions) << nl;
    }


    // Expansion ratios
    // ~~~~~~~~~~~~~~~~

    dict.readIfPresent("ratios", expansion);

    // Correct expansion ratios - negative is the same as inverse.
    for (scalar& expRatio : expansion)
    {
        if (expRatio < 0)
        {
            expRatio = 1.0/(-expRatio);
        }
    }

    if (expansion.size() && expansion.size() != nSegments)
    {
        FatalIOErrorInFunction(dict)
            << "Mismatch in number of expansion ratios and divisions:"
            << " in " << vector::componentNames[cmpt]
            << " direction" << nl
            << exit(FatalIOError);
    }

    if (expansion.empty())
    {
        expansion.resize(nSegments, scalar(1));

        if (verbose_)
        {
            Info<< "Warning: 'ratios' not specified, using constant spacing"
                << nl;
        }
    }
    else
    {
        if (verbose_)
        {
            Info<< "    ratios : " << flatOutput(expansion) << nl;
        }
    }



    // Define the grid points

    scalarList& gridPoint = grid_[cmpt];
    gridPoint.resize(nTotalDivs+1);

    label start = 0;

    for (label segmenti=0; segmenti < nSegments; ++segmenti)
    {
        const label nDiv = divisions[segmenti];
        const scalar expRatio = expansion[segmenti];

        SubList<scalar> subPoint(gridPoint, nDiv+1, start);
        start += nDiv;

        subPoint[0] = knots[segmenti];
        subPoint[nDiv] = knots[segmenti+1];

        const scalar dist = (subPoint.last() - subPoint.first());

        if (equal(expRatio, 1.0))
        {
            for (label i=1; i < nDiv; ++i)
            {
                subPoint[i] = (subPoint[0] + (dist * i)/nDiv);
            }
        }
        else
        {
            // Geometric expansion factor from the expansion ratio
            const scalar expFact = calcGexp(expRatio, nDiv);

            for (label i=1; i < nDiv; ++i)
            {
                subPoint[i] =
                (
                    subPoint[0]
                  + dist * (1.0 - pow(expFact, i))/(1.0 - pow(expFact, nDiv))
                );
            }
        }
    }
}


void Foam::PDRblock::readBoundary(const dictionary& dict)
{
    if (verbose_)
    {
        Info<< "Reading boundary entries" << nl;
    }

    PtrList<entry> patchEntries;

    if (dict.found("boundary"))
    {
        PtrList<entry> newEntries(dict.lookup("boundary"));
        patchEntries.transfer(newEntries);
    }

    patches_.clear();
    patches_.resize(patchEntries.size());

    // Hex cell has 6 sides
    const labelRange validFaceId(0, 6);

    // Track which sides have been handled
    labelHashSet handled;

    forAll(patchEntries, patchi)
    {
        if (!patchEntries.set(patchi))
        {
            continue;
        }

        const entry& e = patchEntries[patchi];

        if (!e.isDict())
        {
            FatalIOErrorInFunction(dict)
                << "Entry " << e
                << " in boundary section is not a valid dictionary."
                << exit(FatalIOError);
        }

        const dictionary& patchDict = e.dict();

        const word& patchName = e.keyword();

        const word patchType(patchDict.get<word>("type"));

        labelList faceLabels(patchDict.get<labelList>("faces"));

        labelHashSet patchFaces(faceLabels);

        if (faceLabels.size() != patchFaces.size())
        {
            WarningInFunction
                << "Patch: " << patchName
                << ": Ignoring duplicate face ids" << nl;
        }

        label nErased = patchFaces.filterKeys(validFaceId);
        if (nErased)
        {
            WarningInFunction
                << "Patch: " << patchName << ": Ignoring " << nErased
                << " faces with invalid identifiers" << nl;
        }

        nErased = patchFaces.erase(handled);
        if (nErased)
        {
            WarningInFunction
                << "Patch: " << patchName << ": Ignoring " << nErased
                << " faces, which were already handled" << nl;
        }

        // Mark these faces as having been handled
        handled += patchFaces;


        // Save information for later access during mesh creation.

        patches_.set(patchi, new boundaryEntry());

        boundaryEntry& bentry = patches_[patchi];

        bentry.name_ = patchName;
        bentry.type_ = patchType;
        bentry.size_ = 0;
        bentry.faces_ = patchFaces.sortedToc();

        // Count patch faces
        for (const label shapeFacei : bentry.faces_)
        {
            bentry.size_ += nBoundaryFaces(shapeFacei);
        }
    }


    // Deal with unhandled faces

    labelHashSet missed(identity(6));
    missed.erase(handled);

    if (missed.size())
    {
        patches_.append(new boundaryEntry());

        boundaryEntry& bentry = patches_.last();

        bentry.name_ = "defaultFaces";
        bentry.type_ = emptyPolyPatch::typeName;
        bentry.size_ = 0;
        bentry.faces_ = missed.sortedToc();

        // Count patch faces
        for (const label shapeFacei : bentry.faces_)
        {
            bentry.size_ += nBoundaryFaces(shapeFacei);
        }


        // Use name/type from "defaultPatch" entry if available
        // - be generous with error handling

        const dictionary* defaultEntry = dict.findDict("defaultPatch");
        if (defaultEntry)
        {
            defaultEntry->readIfPresent("name", bentry.name_);
            defaultEntry->readIfPresent("type", bentry.type_);
        }
    }

    if (verbose_)
    {
        label patchi = 0;
        for (const boundaryEntry& bentry : patches_)
        {
            Info<< "    patch: " << patchi
                << " (size: " << bentry.size_
                << " type: " << bentry.type_
                << ") name: " << bentry.name_
                << " faces: " << flatOutput(bentry.faces_) << nl;

            ++patchi;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PDRblock::PDRblock(const dictionary& dict, bool verbose)
:
    PDRblock()
{
    verbose_ = verbose;

    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::PDRblock::read(const dictionary& dict)
{
    // Mark no scaling with invalid value
    const scalar scaleFactor = dict.lookupOrDefault<scalar>("scale", -1);

    // Grid controls
    readGridControl(0, dict.subDict("x"), scaleFactor);
    readGridControl(1, dict.subDict("y"), scaleFactor);
    readGridControl(2, dict.subDict("z"), scaleFactor);

    // Adjust i-j-k addressing
    sizes().x() = grid_.x().nCells();
    sizes().y() = grid_.y().nCells();
    sizes().z() = grid_.z().nCells();

    // Adjust boundBox
    bounds_.min().x() = grid_.x().first();
    bounds_.min().y() = grid_.y().first();
    bounds_.min().z() = grid_.z().first();

    bounds_.max().x() = grid_.x().last();
    bounds_.max().y() = grid_.y().last();
    bounds_.max().z() = grid_.z().last();


    // Boundaries
    readBoundary(dict);

    return true;
}


Foam::labelVector Foam::PDRblock::findCell(const point& pt) const
{
    // Out-of-bounds is handled explicitly, for efficiency and consistency,
    // but principally to ensure that findLower() returns a valid
    // result when the point is to the right of the bounds.

    labelVector pos(-1,-1,-1);
    if (bounds_.contains(pt))
    {
        for (direction cmpt=0; cmpt < labelVector::nComponents; ++cmpt)
        {
            // Binary search
            pos[cmpt] = findLower(grid_[cmpt], pt[cmpt]);
        }
    }

    return pos;
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2018 OpenCFD Ltd.
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

#include "searchableSphere.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(searchableSphere, 0);
    addToRunTimeSelectionTable
    (
        searchableSurface,
        searchableSphere,
        dict
    );
    addNamedToRunTimeSelectionTable
    (
        searchableSurface,
        searchableSphere,
        dict,
        sphere
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::pointIndexHit Foam::searchableSphere::findNearest
(
    const point& sample,
    const scalar nearestDistSqr
) const
{
    pointIndexHit info(false, sample, -1);

    const vector n(sample - origin_);
    scalar magN = mag(n);

    if (nearestDistSqr >= sqr(magN - radius_))
    {
        if (magN < ROOTVSMALL)
        {
            info.rawPoint() = origin_ + vector(1,0,0)*radius_;
        }
        else
        {
            info.rawPoint() = origin_ + n/magN*radius_;
        }
        info.setHit();
        info.setIndex(0);
    }

    return info;
}


// From Graphics Gems - intersection of sphere with ray
void Foam::searchableSphere::findLineAll
(
    const point& start,
    const point& end,
    pointIndexHit& near,
    pointIndexHit& far
) const
{
    near.setMiss();
    far.setMiss();

    vector dir(end-start);
    scalar magSqrDir = magSqr(dir);

    if (magSqrDir > ROOTVSMALL)
    {
        const vector toCentre(origin_ - start);
        scalar magSqrToCentre = magSqr(toCentre);

        dir /= Foam::sqrt(magSqrDir);

        scalar v = (toCentre & dir);

        scalar disc = sqr(radius_) - (magSqrToCentre - sqr(v));

        if (disc >= 0)
        {
            scalar d = Foam::sqrt(disc);

            scalar nearParam = v-d;

            if (nearParam >= 0 && sqr(nearParam) <= magSqrDir)
            {
                near.setHit();
                near.setPoint(start + nearParam*dir);
                near.setIndex(0);
            }

            scalar farParam = v+d;

            if (farParam >= 0 && sqr(farParam) <= magSqrDir)
            {
                far.setHit();
                far.setPoint(start + farParam*dir);
                far.setIndex(0);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::searchableSphere::searchableSphere
(
    const IOobject& io,
    const point& origin,
    const scalar radius
)
:
    searchableSurface(io),
    origin_(origin),
    radius_(radius)
{
    bounds() = boundBox
    (
        origin_ - radius_*vector::one,
        origin_ + radius_*vector::one
    );
}


Foam::searchableSphere::searchableSphere
(
    const IOobject& io,
    const dictionary& dict
)
:
    searchableSphere
    (
        io,
        dict.getCompat<vector>("origin", {{"centre", -1806}}),
        dict.get<scalar>("radius")
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::searchableSphere::overlaps(const boundBox& bb) const
{
    return bb.overlaps(origin_, sqr(radius_));
}


const Foam::wordList& Foam::searchableSphere::regions() const
{
    if (regions_.empty())
    {
        regions_.resize(1);
        regions_.first() = "region0";
    }
    return regions_;
}



void Foam::searchableSphere::boundingSpheres
(
    pointField& centres,
    scalarField& radiusSqr
) const
{
    centres.resize(1);
    centres[0] = origin_;

    radiusSqr.resize(1);
    radiusSqr[0] = Foam::sqr(radius_);

    // Add a bit to make sure all points are tested inside
    radiusSqr += Foam::sqr(SMALL);
}


void Foam::searchableSphere::findNearest
(
    const pointField& samples,
    const scalarField& nearestDistSqr,
    List<pointIndexHit>& info
) const
{
    info.setSize(samples.size());

    forAll(samples, i)
    {
        info[i] = findNearest(samples[i], nearestDistSqr[i]);
    }
}


void Foam::searchableSphere::findLine
(
    const pointField& start,
    const pointField& end,
    List<pointIndexHit>& info
) const
{
    info.setSize(start.size());

    pointIndexHit b;

    forAll(start, i)
    {
        // Pick nearest intersection. If none intersected take second one.
        findLineAll(start[i], end[i], info[i], b);
        if (!info[i].hit() && b.hit())
        {
            info[i] = b;
        }
    }
}


void Foam::searchableSphere::findLineAny
(
    const pointField& start,
    const pointField& end,
    List<pointIndexHit>& info
) const
{
    info.setSize(start.size());

    pointIndexHit b;

    forAll(start, i)
    {
        // Discard far intersection
        findLineAll(start[i], end[i], info[i], b);
        if (!info[i].hit() && b.hit())
        {
            info[i] = b;
        }
    }
}


void Foam::searchableSphere::findLineAll
(
    const pointField& start,
    const pointField& end,
    List<List<pointIndexHit>>& info
) const
{
    info.setSize(start.size());

    forAll(start, i)
    {
        pointIndexHit near, far;
        findLineAll(start[i], end[i], near, far);

        if (near.hit())
        {
            if (far.hit())
            {
                info[i].setSize(2);
                info[i][0] = near;
                info[i][1] = far;
            }
            else
            {
                info[i].setSize(1);
                info[i][0] = near;
            }
        }
        else
        {
            if (far.hit())
            {
                info[i].setSize(1);
                info[i][0] = far;
            }
            else
            {
                info[i].clear();
            }
        }
    }
}


void Foam::searchableSphere::getRegion
(
    const List<pointIndexHit>& info,
    labelList& region
) const
{
    region.setSize(info.size());
    region = 0;
}


void Foam::searchableSphere::getNormal
(
    const List<pointIndexHit>& info,
    vectorField& normal
) const
{
    normal.setSize(info.size());
    normal = Zero;

    forAll(info, i)
    {
        if (info[i].hit())
        {
            normal[i] = normalised(info[i].hitPoint() - origin_);
        }
        else
        {
            // Set to what?
        }
    }
}


void Foam::searchableSphere::getVolumeType
(
    const pointField& points,
    List<volumeType>& volType
) const
{
    volType.setSize(points.size());

    const scalar rad2 = sqr(radius_);

    forAll(points, pointi)
    {
        const point& pt = points[pointi];

        volType[pointi] =
        (
            (magSqr(pt - origin_) <= rad2)
          ? volumeType::INSIDE : volumeType::OUTSIDE
        );
    }
}


// ************************************************************************* //

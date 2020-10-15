/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2018-2020 OpenCFD Ltd.
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


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

// General handling
namespace Foam
{

// Dictionary entry with single scalar or vector quantity
inline static vector getRadius(const word& name, const dictionary& dict)
{
    if (token(dict.lookup(name)).isNumber())
    {
        return vector::uniform(dict.get<scalar>(name));
    }

    return dict.get<vector>(name);
}


// Test point for negative components, return the sign-changes
inline static unsigned getOctant(const point& p)
{
    unsigned octant = 0;

    if (p.x() < 0) { octant |= 1; }
    if (p.y() < 0) { octant |= 2; }
    if (p.z() < 0) { octant |= 4; }

    return octant;
}


// Apply sign-changes to point
inline static void applyOctant(point& p, unsigned octant)
{
    if (octant & 1) { p.x() = -p.x(); }
    if (octant & 2) { p.y() = -p.y(); }
    if (octant & 4) { p.z() = -p.z(); }
}

} // End namespace Foam


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

inline Foam::searchableSphere::componentOrder
Foam::searchableSphere::getOrdering(const vector& radii)
{
    #ifdef FULLDEBUG
    for (direction cmpt = 0; cmpt < vector::nComponents; ++cmpt)
    {
        if (radii[cmpt] <= 0)
        {
            FatalErrorInFunction
                << "Radii must be positive, non-zero: " << radii << endl
                << exit(FatalError);
        }
    }
    #endif

    std::array<direction, 3> idx{0, 1, 2};

    // Reverse sort by magnitude (largest first...)
    // Radii are positive (checked above, or just always true)
    std::stable_sort
    (
        idx.begin(),
        idx.end(),
        [&](direction a, direction b){ return radii[a] > radii[b]; }
    );

    componentOrder order{idx[0], idx[1], idx[2], shapeType::GENERAL};

    if (equal(radii[order.major], radii[order.minor]))
    {
        order.shape = shapeType::SPHERE;
    }
    else if (equal(radii[order.major], radii[order.mezzo]))
    {
        order.shape = shapeType::OBLATE;
    }
    else if (equal(radii[order.mezzo], radii[order.minor]))
    {
        order.shape = shapeType::PROLATE;
    }

    return order;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::pointIndexHit Foam::searchableSphere::findNearest
(
    const point& sample,
    const scalar nearestDistSqr
) const
{
    pointIndexHit info(false, sample, -1);

    // Handle special cases first

    // if (order_.shape == shapeType::SPHERE)
    if (true)
    {
        // Point relative to origin, simultaneously the normal on the sphere
        const vector n(sample - origin_);
        const scalar magN = mag(n);

        // It is a sphere, all radii are identical

        if (nearestDistSqr >= sqr(magN - radii_[0]))
        {
            info.setPoint
            (
                origin_
              + (magN < ROOTVSMALL ? vector(radii_[0],0,0) : (radii_[0]*n/magN))
            );
            info.setHit();
            info.setIndex(0);
        }

        return info;
    }

    //[code]

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

    // if (order_.shape == shapeType::SPHERE)
    if (true)
    {
        vector dir(end-start);
        const scalar magSqrDir = magSqr(dir);

        if (magSqrDir > ROOTVSMALL)
        {
            dir /= Foam::sqrt(magSqrDir);

            const vector relStart(start - origin_);

            const scalar v = -(relStart & dir);

            const scalar disc = sqr(radius()) - (magSqr(relStart) - sqr(v));

            if (disc >= 0)
            {
                const scalar d = Foam::sqrt(disc);

                const scalar nearParam = v - d;
                const scalar farParam = v + d;

                if (nearParam >= 0 && sqr(nearParam) <= magSqrDir)
                {
                    near.setHit();
                    near.setPoint(start + nearParam*dir);
                    near.setIndex(0);
                }

                if (farParam >= 0 && sqr(farParam) <= magSqrDir)
                {
                    far.setHit();
                    far.setPoint(start + farParam*dir);
                    far.setIndex(0);
                }
            }
        }
    }

    //[code]
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::searchableSphere::searchableSphere
(
    const IOobject& io,
    const point& origin,
    const scalar radius
)
:
    searchableSphere(io, origin, vector::uniform(radius))
{}


Foam::searchableSphere::searchableSphere
(
    const IOobject& io,
    const point& origin,
    const vector& radii
)
:
    searchableSurface(io),
    origin_(origin),
    // radii_(radii),
    radii_(vector::uniform(cmptMax(radii))) /* Transition */,
    order_{getOrdering(radii_)}
{
    bounds().min() = (centre() - radii_);
    bounds().max() = (centre() + radii_);
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
        getRadius("radius", dict)
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::point Foam::searchableSphere::surfacePoint
(
    const scalar theta,
    const scalar phi
) const
{
    return point
    (
        origin_.x() + radii_.x() * cos(theta)*sin(phi),
        origin_.y() + radii_.y() * sin(theta)*sin(phi),
        origin_.z() + radii_.z() * cos(phi)
    );
}


Foam::vector Foam::searchableSphere::surfaceNormal
(
    const scalar theta,
    const scalar phi
) const
{
    // Normal is (x0/r0^2, x1/r1^2, x2/r2^2)

    return vector
    (
        cos(theta)*sin(phi) / radii_.x(),
        sin(theta)*sin(phi) / radii_.y(),
        cos(phi) / radii_.z()
    ).normalise();
}


bool Foam::searchableSphere::overlaps(const boundBox& bb) const
{
    // if (order_.shape == shapeType::SPHERE)
    if (true)
    {
        return bb.overlaps(origin_, sqr(radius()));
    }

    if (!bb.valid())
    {
        return false;
    }

    //[code]

    return true;
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
    radiusSqr.resize(1);

    centres[0] = origin_;
    radiusSqr[0] = Foam::sqr(radius());

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
    info.resize(samples.size());

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
    info.resize(start.size());

    pointIndexHit b;

    forAll(start, i)
    {
        // Pick nearest intersection.
        // If none intersected take second one.

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
    info.resize(start.size());

    pointIndexHit b;

    forAll(start, i)
    {
        // Pick nearest intersection.
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
    info.resize(start.size());

    forAll(start, i)
    {
        pointIndexHit near, far;

        findLineAll(start[i], end[i], near, far);

        if (near.hit())
        {
            if (far.hit())
            {
                info[i].resize(2);
                info[i][0] = near;
                info[i][1] = far;
            }
            else
            {
                info[i].resize(1);
                info[i][0] = near;
            }
        }
        else
        {
            if (far.hit())
            {
                info[i].resize(1);
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
    region.resize(info.size());
    region = 0;
}


void Foam::searchableSphere::getNormal
(
    const List<pointIndexHit>& info,
    vectorField& normal
) const
{
    normal.resize(info.size());

    forAll(info, i)
    {
        if (info[i].hit())
        {
            normal[i] = normalised(info[i].point() - origin_);
        }
        else
        {
            // Set to what?
            normal[i] = Zero;
        }
    }
}


void Foam::searchableSphere::getVolumeType
(
    const pointField& points,
    List<volumeType>& volType
) const
{
    volType.resize(points.size());

    // if (order_.shape == shapeType::SPHERE)
    if (true)
    {
        const scalar rad2 = sqr(radius());

        forAll(points, pointi)
        {
            const point& p = points[pointi];

            volType[pointi] =
            (
                (magSqr(p - origin_) <= rad2)
              ? volumeType::INSIDE : volumeType::OUTSIDE
            );
        }
    }

    //[code]
}


// ************************************************************************* //

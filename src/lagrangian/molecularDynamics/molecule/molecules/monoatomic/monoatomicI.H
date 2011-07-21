/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2011 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::monoatomic::constantProperties::constantProperties()
:
    sites_(),
    electrostaticSites_(-1),
    pairPotSites_(-1),
    mass_(0)
{}


inline Foam::monoatomic::constantProperties::constantProperties
(
    const dictionary& dict,
    const List<label>& siteIds
)
:
    sites_(),
    electrostaticSites_(-1),
    pairPotSites_(-1),
    mass_(0)
{
    if (siteIds.size() != 1)
    {
        FatalErrorIn
        (
            "Foam::monoatomic::constantProperties::constantProperties"
            "("
                "const dictionary& dict, "
                "const List<label>& siteIds"
            ")"
        )
            << "monoatomic, single site only, given: " << dict
            << nl << abort(FatalError);
    }

    FixedList<word, 1> siteIdNames = dict.lookup("siteIds");

    FixedList<vector, 1> siteReferencePositions
    (
        dict.lookup("siteReferencePositions")
    );

    FixedList<scalar, 1> siteMasses(dict.lookup("siteMasses"));

    FixedList<scalar, 1> siteCharges(dict.lookup("siteCharges"));

    FixedList<word, 1> pairPotentialIds(dict.lookup("pairPotentialSiteIds"));

    constPropSite site = sites_[0];

    site = constPropSite
    (
        siteReferencePositions[0],
        siteMasses[0],
        siteCharges[0],
        siteIds[0],
        siteIdNames[0],
        (findIndex(pairPotentialIds, siteIdNames[0]) != -1),    // pair
        (mag(siteCharges[0]) > VSMALL)                          // charge
    );

    mass_ = site.siteMass();

    if (site.pairPotentialSite())
    {
        pairPotSites_[0] = 0;
    }

    if (site.electrostaticSite())
    {
        electrostaticSites_[0] = 0;
    }

    if (!site.pairPotentialSite() && !site.electrostaticSite())
    {
        WarningIn
        (
            "Foam::monoatomic::constantProperties::constantProperties"
            "("
            "const dictionary& dict"
            ")"
        )
        << siteIdNames[0] << " is a non-interacting site." << endl;
    }

    // Single site monoatomic - no rotational motion.
    site.siteReferencePosition() = vector::zero;
}


inline Foam::monoatomic::monoatomic
(
    const polyMesh& mesh,
    const vector& position,
    const label cellI,
    const label tetFaceI,
    const label tetPtI,
    const scalar temperature,
    const vector& bulkVelocity,
    const vector& specialPosition,
    const constantProperties& cP,
    monoatomic::trackingData& td,
    const label special,
    const label id
)
:
    particle(mesh, position, cellI, tetFaceI, tetPtI),
    v_(bulkVelocity),
    a_(vector::zero),
    specialPosition_(specialPosition),
    potentialEnergy_(0.0),
    rf_(tensor::zero),
    special_(special),
    id_(id),
    siteForces_(vector::zero),
    sitePositions_()
{
    setSitePositions(cP);

    v_ += td.cloud().equipartitionLinearVelocity(temperature, cP.mass());
}


inline Foam::monoatomic::monoatomic
(
    const polyMesh& mesh,
    const vector& position,
    const label cellI,
    const label tetFaceI,
    const label tetPtI,
    const vector& v,
    const vector& a,
    const vector& specialPosition,
    const constantProperties& cP,
    const label special,
    const label id

)
:
    particle(mesh, position, cellI, tetFaceI, tetPtI),
    v_(v),
    a_(a),
    specialPosition_(specialPosition),
    potentialEnergy_(0.0),
    rf_(tensor::zero),
    special_(special),
    id_(id),
    siteForces_(vector::zero),
    sitePositions_()
{
    setSitePositions(cP);
}

// * * * * * * * constantProperties Member Functions * * * * * * * * * * * * //

inline const Foam::FixedList<Foam::constPropSite, 1>&
Foam::monoatomic::constantProperties::sites() const
{
    return sites_;
}


inline const Foam::FixedList<Foam::label, 1>&
Foam::monoatomic::constantProperties::pairPotSites() const
{
    return pairPotSites_;
}


inline const Foam::FixedList<Foam::label, 1>&
Foam::monoatomic::constantProperties::electrostaticSites() const
{
    return electrostaticSites_;
}


inline Foam::label
Foam::monoatomic::constantProperties::degreesOfFreedom() const
{
    return 3;
}


inline Foam::scalar Foam::monoatomic::constantProperties::mass() const
{
    return mass_;
}


inline Foam::label Foam::monoatomic::constantProperties::nSites() const
{
    return 1;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

inline const Foam::vector& Foam::monoatomic::v() const
{
    return v_;
}


inline Foam::vector& Foam::monoatomic::v()
{
    return v_;
}


inline const Foam::vector& Foam::monoatomic::a() const
{
    return a_;
}


inline Foam::vector& Foam::monoatomic::a()
{
    return a_;
}


inline const Foam::FixedList<Foam::vector, 1>&
Foam::monoatomic::siteForces() const
{
    return siteForces_;
}


inline Foam::FixedList<Foam::vector, 1>& Foam::monoatomic::siteForces()
{
    return siteForces_;
}


inline const Foam::FixedList<Foam::vector, 1>&
Foam::monoatomic::sitePositions() const
{
    return sitePositions_;
}


inline Foam::FixedList<Foam::vector, 1>& Foam::monoatomic::sitePositions()
{
    return sitePositions_;
}


inline const Foam::vector& Foam::monoatomic::specialPosition() const
{
    return specialPosition_;
}


inline Foam::vector& Foam::monoatomic::specialPosition()
{
    return specialPosition_;
}


inline Foam::scalar Foam::monoatomic::potentialEnergy() const
{
    return potentialEnergy_;
}


inline Foam::scalar& Foam::monoatomic::potentialEnergy()
{
    return potentialEnergy_;
}


inline const Foam::tensor& Foam::monoatomic::rf() const
{
    return rf_;
}


inline Foam::tensor& Foam::monoatomic::rf()
{
    return rf_;
}


inline Foam::label Foam::monoatomic::special() const
{
    return special_;
}


inline bool Foam::monoatomic::tethered() const
{
    return special_ == SPECIAL_TETHERED;
}


inline Foam::label Foam::monoatomic::id() const
{
    return id_;
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::molecule::constantProperties::constantProperties()
:
    siteReferencePositions_(List<vector>(0)),
    siteCharges_(List<scalar>(0)),
    siteIds_(List<label>(0)),
    momentOfInertia_(diagTensor(0, 0, 0)),
    mass_(0)
{}


inline Foam::molecule::constantProperties::constantProperties
(
    const List<vector>& siteReferencePositions,
    const List<scalar>& siteCharges,
    const List<label>& siteIds,
    const diagTensor& momentOfInertia,
    scalar mass
)
:
    siteReferencePositions_(siteReferencePositions),
    siteCharges_(siteCharges),
    siteIds_(siteIds),
    momentOfInertia_(momentOfInertia),
    mass_(mass)
{
    checkSiteListSizes();
}


inline Foam::molecule::constantProperties::constantProperties
(
    const dictionary& dict
)
:
    siteReferencePositions_(dict.lookup("siteReferencePositions")),
    siteCharges_(dict.lookup("siteCharges")),
    siteIds_(List<word>(dict.lookup("siteIds")).size(), -1),
    momentOfInertia_(dict.lookup("momentOfInertia")),
    mass_(readScalar(dict.lookup("mass")))
{
    checkSiteListSizes();
}


inline Foam::molecule::molecule
(
    const Cloud<molecule>& c,
    const vector& position,
    const label celli,
    const tensor& R,
    const vector& v,
    const vector& a,
    const vector& omega,
    const vector& alpha,
    const vector& specialPosition,
    const constantProperties& constProps,
    const label special,
    const label id

)
:
    Particle<molecule>(c, position, celli),
    R_(R),
    v_(v),
    a_(a),
    omega_(omega),
    alpha_(alpha),
    siteForces_(constProps.nSites(), vector::zero),
    sitePositions_(constProps.nSites()),
    specialPosition_(specialPosition),
    potentialEnergy_(0.0),
    rf_(tensor::zero),
    special_(special),
    id_(id)
{
    setSitePositions(constProps);
}


// * * * constantProperties Private Member Functions * * * * * * * * * * * * //


inline void Foam::molecule::constantProperties::checkSiteListSizes() const
{
    if
    (
        siteIds_.size() != siteReferencePositions_.size()
     || siteIds_.size() != siteCharges_.size()
    )
    {
        FatalErrorIn("molecule::constantProperties::checkSiteListSizes")
            << "Sizes of site id, charge and "
            << "referencePositions are not the same. "
            << nl << abort(FatalError);
    }
}


// * * * * * * * constantProperties Member Functions * * * * * * * * * * * * //

inline const Foam::List<Foam::vector>&
    Foam::molecule::constantProperties::siteReferencePositions() const
{
    return siteReferencePositions_;
}


inline const Foam::List<Foam::scalar>&
    Foam::molecule::constantProperties::siteCharges() const
{
    return siteCharges_;
}


inline const Foam::List<Foam::label>&
    Foam::molecule::constantProperties::siteIds() const
{
    return siteIds_;
}


inline Foam::List<Foam::label>&
    Foam::molecule::constantProperties::siteIds()
{
    return siteIds_;
}


inline const Foam::diagTensor&
    Foam::molecule::constantProperties::momentOfInertia() const
{
    return momentOfInertia_;
}


inline Foam::scalar Foam::molecule::constantProperties::mass() const
{
    return mass_;
}


inline Foam::label Foam::molecule::constantProperties::nSites() const
{
    return siteIds_.size();

}

// * * * * * * * * * * *  trackData Member Functions * * * * * * * * * * * * //

inline Foam::moleculeCloud& Foam::molecule::trackData::molCloud()
{
    return molCloud_;
}


inline Foam::label Foam::molecule::trackData::part() const
{
    return part_;
}


// * * * * * * * * * * * * molecule Member Functions * * * * * * * * * * * * //

inline const Foam::tensor& Foam::molecule::R() const
{
    return R_;
}


inline Foam::tensor& Foam::molecule::R()
{
    return R_;
}


inline const Foam::vector& Foam::molecule::v() const
{
    return v_;
}


inline Foam::vector& Foam::molecule::v()
{
    return v_;
}


inline const Foam::vector& Foam::molecule::a() const
{
    return a_;
}


inline Foam::vector& Foam::molecule::a()
{
    return a_;
}


inline const Foam::vector& Foam::molecule::omega() const
{
    return omega_;
}


inline Foam::vector& Foam::molecule::omega()
{
    return omega_;
}


inline const Foam::vector& Foam::molecule::alpha() const
{
    return alpha_;
}


inline Foam::vector& Foam::molecule::alpha()
{
    return alpha_;
}


inline const Foam::List<Foam::vector>& Foam::molecule::siteForces() const
{
    return siteForces_;
}


inline Foam::List<Foam::vector>& Foam::molecule::siteForces()
{
    return siteForces_;
}


inline const Foam::List<Foam::vector>& Foam::molecule::sitePositions() const
{
    return sitePositions_;
}


inline Foam::List<Foam::vector>& Foam::molecule::sitePositions()
{
    return sitePositions_;
}


inline const Foam::vector& Foam::molecule::specialPosition() const
{
    return specialPosition_;
}


inline Foam::vector& Foam::molecule::specialPosition()
{
    return specialPosition_;
}


inline Foam::scalar Foam::molecule::potentialEnergy() const
{
    return potentialEnergy_;
}


inline Foam::scalar& Foam::molecule::potentialEnergy()
{
    return potentialEnergy_;
}


inline const Foam::tensor& Foam::molecule::rf() const
{
    return rf_;
}


inline Foam::tensor& Foam::molecule::rf()
{
    return rf_;
}


inline Foam::label Foam::molecule::special() const
{
    return special_;
}


inline bool Foam::molecule::tethered() const
{
    return special_ == SPECIAL_TETHERED;
}


inline Foam::label Foam::molecule::id() const
{
    return id_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //

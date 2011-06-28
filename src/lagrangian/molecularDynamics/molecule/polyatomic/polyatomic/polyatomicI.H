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

inline Foam::polyatomic::constantProperties::constantProperties()
:
    siteReferencePositions_(Field<vector>(0)),
    siteMasses_(List<scalar>(0)),
    siteCharges_(List<scalar>(0)),
    siteIds_(List<label>(0)),
    pairPotentialSites_(List<bool>(false)),
    electrostaticSites_(List<bool>(false)),
    momentOfInertia_(diagTensor(0, 0, 0)),
    mass_(0)
{}


inline Foam::polyatomic::constantProperties::constantProperties
(
    const dictionary& dict
)
:
    siteReferencePositions_(dict.lookup("siteReferencePositions")),
    siteMasses_(dict.lookup("siteMasses")),
    siteCharges_(dict.lookup("siteCharges")),
    siteIds_(List<word>(dict.lookup("siteIds")).size(), -1),
    pairPotentialSites_(),
    electrostaticSites_(),
    momentOfInertia_(),
    mass_()
{
    checkSiteListSizes();

    setInteracionSiteBools
    (
        List<word>(dict.lookup("siteIds")),
        List<word>(dict.lookup("pairPotentialSiteIds"))
    );

    mass_ = sum(siteMasses_);

    vector centreOfMass(vector::zero);

    // Calculate the centre of mass of the body and subtract it from each
    // position

    forAll(siteReferencePositions_, i)
    {
        centreOfMass += siteReferencePositions_[i]*siteMasses_[i];
    }

    centreOfMass /= mass_;

    siteReferencePositions_ -= centreOfMass;

    if (siteIds_.size() == 1)
    {
        // Single site polyatomic - no rotational motion.

        siteReferencePositions_[0] = vector::zero;

        momentOfInertia_ = diagTensor(-1, -1, -1);
    }
    else if (linearMoleculeTest())
    {
        // Linear polyatomic.

        Info<< nl << "Linear polyatomic." << endl;

        vector dir = siteReferencePositions_[1] - siteReferencePositions_[0];

        dir /= mag(dir);

        tensor Q = rotationTensor(dir, vector(1,0,0));

        siteReferencePositions_ = (Q & siteReferencePositions_);

        // The rotation was around the centre of mass but remove any
        // components that have crept in due to floating point errors

        centreOfMass = vector::zero;

        forAll(siteReferencePositions_, i)
        {
            centreOfMass += siteReferencePositions_[i]*siteMasses_[i];
        }

        centreOfMass /= mass_;

        siteReferencePositions_ -= centreOfMass;

        diagTensor momOfInertia = diagTensor::zero;

        forAll(siteReferencePositions_, i)
        {
            const vector& p(siteReferencePositions_[i]);

            momOfInertia +=
                siteMasses_[i]*diagTensor(0, p.x()*p.x(), p.x()*p.x());
        }

        momentOfInertia_ = diagTensor
         (
             -1,
             momOfInertia.yy(),
             momOfInertia.zz()
         );
    }
    else
    {
        // Fully 6DOF polyatomic

        // Calculate the inertia tensor in the current orientation

        tensor momOfInertia(tensor::zero);

        forAll(siteReferencePositions_, i)
        {
            const vector& p(siteReferencePositions_[i]);

            momOfInertia += siteMasses_[i]*tensor
            (
                p.y()*p.y() + p.z()*p.z(), -p.x()*p.y(), -p.x()*p.z(),
                -p.y()*p.x(), p.x()*p.x() + p.z()*p.z(), -p.y()*p.z(),
                -p.z()*p.x(), -p.z()*p.y(), p.x()*p.x() + p.y()*p.y()
            );
        }

        if (eigenValues(momOfInertia).x() < VSMALL)
        {
            FatalErrorIn("polyatomic::constantProperties::constantProperties")
                << "An eigenvalue of the inertia tensor is zero, the polyatomic "
                << dict.name() << " is not a valid 6DOF shape."
                << nl << abort(FatalError);
        }

        // Normalise the inertia tensor magnitude to avoid SMALL numbers in the
        // components causing problems

        momOfInertia /= eigenValues(momOfInertia).x();

        tensor e = eigenVectors(momOfInertia);

        // Calculate the transformation between the principle axes to the
        // global axes

        tensor Q =
            vector(1,0,0)*e.x() + vector(0,1,0)*e.y() + vector(0,0,1)*e.z();

        // Transform the site positions

        siteReferencePositions_ = (Q & siteReferencePositions_);

        // Recalculating the moment of inertia with the new site positions

        // The rotation was around the centre of mass but remove any
        // components that have crept in due to floating point errors

        centreOfMass = vector::zero;

        forAll(siteReferencePositions_, i)
        {
            centreOfMass += siteReferencePositions_[i]*siteMasses_[i];
        }

        centreOfMass /= mass_;

        siteReferencePositions_ -= centreOfMass;

        // Calculate the moment of inertia in the principle component
        // reference frame

        momOfInertia = tensor::zero;

        forAll(siteReferencePositions_, i)
        {
            const vector& p(siteReferencePositions_[i]);

            momOfInertia += siteMasses_[i]*tensor
            (
                p.y()*p.y() + p.z()*p.z(), -p.x()*p.y(), -p.x()*p.z(),
                -p.y()*p.x(), p.x()*p.x() + p.z()*p.z(), -p.y()*p.z(),
                -p.z()*p.x(), -p.z()*p.y(), p.x()*p.x() + p.y()*p.y()
            );
        }

        momentOfInertia_ = diagTensor
        (
            momOfInertia.xx(),
            momOfInertia.yy(),
            momOfInertia.zz()
        );
    }
}


inline Foam::polyatomic::polyatomic
(
    const polyMesh& mesh,
    const vector& position,
    const label cellI,
    const label tetFaceI,
    const label tetPtI,
    const tensor& Q,
    const vector& v,
    const vector& a,
    const vector& pi,
    const vector& tau,
    const vector& specialPosition,
    const constantProperties& constProps,
    const label special,
    const label id

)
:
    particle(mesh, position, cellI, tetFaceI, tetPtI),
    Q_(Q),
    v_(v),
    a_(a),
    pi_(pi),
    tau_(tau),
    specialPosition_(specialPosition),
    potentialEnergy_(0.0),
    rf_(tensor::zero),
    special_(special),
    id_(id),
    siteForces_(constProps.nSites(), vector::zero),
    sitePositions_(constProps.nSites())
{
    setSitePositions(constProps);
}


// * * * * * * * constantProperties Private Member Functions * * * * * * * * //

inline void Foam::polyatomic::constantProperties::checkSiteListSizes() const
{
    if
    (
        siteIds_.size() != siteReferencePositions_.size()
     || siteIds_.size() != siteCharges_.size()
    )
    {
        FatalErrorIn("polyatomic::constantProperties::checkSiteListSizes")
            << "Sizes of site id, charge and "
            << "referencePositions are not the same. "
            << nl << abort(FatalError);
    }
}


inline void Foam::polyatomic::constantProperties::setInteracionSiteBools
(
    const List<word>& siteIds,
    const List<word>& pairPotSiteIds
)
{
    pairPotentialSites_.setSize(siteIds_.size());

    electrostaticSites_.setSize(siteIds_.size());

    forAll(siteIds_, i)
    {
        const word& id(siteIds[i]);

        pairPotentialSites_[i] = (findIndex(pairPotSiteIds, id) > -1);

        electrostaticSites_[i] = (mag(siteCharges_[i]) > VSMALL);
    }
}


inline bool Foam::polyatomic::constantProperties::linearMoleculeTest() const
{
    if (siteIds_.size() == 2)
    {
        return true;
    }

    vector refDir = siteReferencePositions_[1] - siteReferencePositions_[0];

    refDir /= mag(refDir);

    for
    (
        label i = 2;
        i < siteReferencePositions_.size();
        i++
    )
    {
        vector dir = siteReferencePositions_[i] - siteReferencePositions_[i-1];

        dir /= mag(dir);

        if (mag(refDir & dir) < 1 - SMALL)
        {
            return false;
        }
    }

    return true;
}


// * * * * * * * constantProperties Member Functions * * * * * * * * * * * * //

inline const Foam::Field<Foam::vector>&
Foam::polyatomic::constantProperties::siteReferencePositions() const
{
    return siteReferencePositions_;
}


inline const Foam::List<Foam::scalar>&
Foam::polyatomic::constantProperties::siteMasses() const
{
    return siteMasses_;
}


inline const Foam::List<Foam::scalar>&
Foam::polyatomic::constantProperties::siteCharges() const
{
    return siteCharges_;
}


inline const Foam::List<Foam::label>&
Foam::polyatomic::constantProperties::siteIds() const
{
    return siteIds_;
}


inline Foam::List<Foam::label>&
Foam::polyatomic::constantProperties::siteIds()
{
    return siteIds_;
}


inline const Foam::List<bool>&
Foam::polyatomic::constantProperties::pairPotentialSites() const
{
    return pairPotentialSites_;
}


inline bool Foam::polyatomic::constantProperties::pairPotentialSite
(
    label sId
) const
{
    label s = findIndex(siteIds_, sId);

    if (s == -1)
    {
        FatalErrorIn("polyatomicI.H") << nl
            << sId << " site not found."
            << nl << abort(FatalError);
    }

    return pairPotentialSites_[s];
}


inline const Foam::List<bool>&
Foam::polyatomic::constantProperties::electrostaticSites() const
{
    return electrostaticSites_;
}


inline bool Foam::polyatomic::constantProperties::electrostaticSite
(
    label sId
) const
{
    label s = findIndex(siteIds_, sId);

    if (s == -1)
    {
        FatalErrorIn
        (
            "polyatomic::constantProperties::electrostaticSite(label)"
        )   << sId << " site not found."
            << nl << abort(FatalError);
    }

    return electrostaticSites_[s];
}


inline const Foam::diagTensor&
Foam::polyatomic::constantProperties::momentOfInertia() const
{
    return momentOfInertia_;
}


inline bool Foam::polyatomic::constantProperties::linearMolecule() const
{
    return ((momentOfInertia_.xx() < 0) && (momentOfInertia_.yy() > 0));
}


inline bool Foam::polyatomic::constantProperties::pointMolecule() const
{
    return (momentOfInertia_.zz() < 0);
}


inline Foam::label Foam::polyatomic::constantProperties::degreesOfFreedom() const
{
    if (linearMolecule())
    {
        return 5;
    }
    else if (pointMolecule())
    {
        return 3;
    }
    else
    {
        return 6;
    }
}


inline Foam::scalar Foam::polyatomic::constantProperties::mass() const
{
    return mass_;
}


inline Foam::label Foam::polyatomic::constantProperties::nSites() const
{
    return siteIds_.size();

}


// * * * * * * * * * * * * polyatomic Member Functions * * * * * * * * * * * * //

inline const Foam::tensor& Foam::polyatomic::Q() const
{
    return Q_;
}


inline Foam::tensor& Foam::polyatomic::Q()
{
    return Q_;
}


inline const Foam::vector& Foam::polyatomic::v() const
{
    return v_;
}


inline Foam::vector& Foam::polyatomic::v()
{
    return v_;
}


inline const Foam::vector& Foam::polyatomic::a() const
{
    return a_;
}


inline Foam::vector& Foam::polyatomic::a()
{
    return a_;
}


inline const Foam::vector& Foam::polyatomic::pi() const
{
    return pi_;
}


inline Foam::vector& Foam::polyatomic::pi()
{
    return pi_;
}


inline const Foam::vector& Foam::polyatomic::tau() const
{
    return tau_;
}


inline Foam::vector& Foam::polyatomic::tau()
{
    return tau_;
}


inline const Foam::List<Foam::vector>& Foam::polyatomic::siteForces() const
{
    return siteForces_;
}


inline Foam::List<Foam::vector>& Foam::polyatomic::siteForces()
{
    return siteForces_;
}


inline const Foam::List<Foam::vector>& Foam::polyatomic::sitePositions() const
{
    return sitePositions_;
}


inline Foam::List<Foam::vector>& Foam::polyatomic::sitePositions()
{
    return sitePositions_;
}


inline const Foam::vector& Foam::polyatomic::specialPosition() const
{
    return specialPosition_;
}


inline Foam::vector& Foam::polyatomic::specialPosition()
{
    return specialPosition_;
}


inline Foam::scalar Foam::polyatomic::potentialEnergy() const
{
    return potentialEnergy_;
}


inline Foam::scalar& Foam::polyatomic::potentialEnergy()
{
    return potentialEnergy_;
}


inline const Foam::tensor& Foam::polyatomic::rf() const
{
    return rf_;
}


inline Foam::tensor& Foam::polyatomic::rf()
{
    return rf_;
}


inline Foam::label Foam::polyatomic::special() const
{
    return special_;
}


inline bool Foam::polyatomic::tethered() const
{
    return special_ == SPECIAL_TETHERED;
}


inline Foam::label Foam::polyatomic::id() const
{
    return id_;
}


// ************************************************************************* //

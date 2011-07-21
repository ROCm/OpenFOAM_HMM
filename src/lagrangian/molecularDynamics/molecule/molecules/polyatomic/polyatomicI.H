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
    sites_(),
    electrostaticSites_(),
    pairPotSites_(),
    momentOfInertia_(diagTensor(0, 0, 0)),
    mass_(0)
{}


inline Foam::polyatomic::constantProperties::constantProperties
(
    const dictionary& dict,
    const labelList& siteIds
)
:
    sites_(),
    electrostaticSites_(),
    pairPotSites_(),
    momentOfInertia_(),
    mass_(0)
{
    List<word> siteIdNames = dict.lookup("siteIds");

    List<vector> siteReferencePositions(dict.lookup("siteReferencePositions"));

    List<scalar> siteMasses(dict.lookup("siteMasses"));

    List<scalar> siteCharges(dict.lookup("siteCharges"));

    List<word> pairPotentialIds(dict.lookup("pairPotentialSiteIds"));

    sites_.setSize(siteReferencePositions.size());

    if
    (
        (siteIdNames.size() != sites_.size())
     || (siteCharges.size() != sites_.size())
    )
    {
            FatalErrorIn
            (
                "Foam::polyatomic::constantProperties::constantProperties"
                "("
                    "const dictionary& dict"
                ")"
            )
            << "Sizes of site id, charge and "
            << "referencePositions are not the same: " << nl
            << siteIdNames  << nl
            << siteReferencePositions << nl
            << siteCharges << nl
            << abort(FatalError);
    }

    electrostaticSites_.setSize(sites_.size(), -1);
    pairPotSites_.setSize(sites_.size(), -1);

    label pairI = 0;
    label electrostaticI = 0;

    forAll(sites_, sI)
    {
        sites_[sI] = constPropSite
        (
            siteReferencePositions[sI],
            siteMasses[sI],
            siteCharges[sI],
            siteIds[sI],
            siteIdNames[sI],
            (findIndex(pairPotentialIds, siteIdNames[sI]) != -1),   // pair
            (mag(siteCharges[sI]) > VSMALL)                         // charge
        );

        if (sites_[sI].pairPotentialSite())
        {
            pairPotSites_[pairI++] = sI;
        }

        if (sites_[sI].electrostaticSite())
        {
            electrostaticSites_[electrostaticI++] = sI;
        }

        if (!sites_[sI].pairPotentialSite() && !sites_[sI].electrostaticSite())
        {
            WarningIn
            (
                "Foam::polyatomic::constantProperties::constantProperties"
                "("
                    "const dictionary& dict"
                ")"
            )
                << siteIdNames[sI] << " is a non-interacting site." << endl;
        }
    }

    pairPotSites_.setSize(pairI);
    electrostaticSites_.setSize(electrostaticI);

    // Calculate the centre of mass of the body and subtract it from each
    // position

    vector centreOfMass(vector::zero);

    forAll(sites_, sI)
    {
        mass_ += sites_[sI].siteMass();

        centreOfMass +=
            sites_[sI].siteReferencePosition()*sites_[sI].siteMass();
    }

    centreOfMass /= mass_;

    forAll(sites_, sI)
    {
        sites_[sI].siteReferencePosition() -= centreOfMass;
    }

    if (sites_.size() == 1)
    {
        // Single site polyatomic - no rotational motion.

        sites_[0].siteReferencePosition() = vector::zero;

        momentOfInertia_ = diagTensor(-1, -1, -1);
    }
    else if (linearMoleculeTest())
    {
        // Linear polyatomic.

        Info<< nl << "Linear polyatomic." << endl;

        vector dir =
            sites_[1].siteReferencePosition()
          - sites_[0].siteReferencePosition();

        dir /= mag(dir);

        tensor Q = rotationTensor(dir, vector(1,0,0));

        forAll(sites_, sI)
        {
            sites_[sI].siteReferencePosition() =
                (Q & sites_[sI].siteReferencePosition());
        }

        // The rotation was around the centre of mass but remove any
        // components that have crept in due to floating point errors

        centreOfMass = vector::zero;

        forAll(sites_, sI)
        {
            centreOfMass +=
                sites_[sI].siteReferencePosition()*sites_[sI].siteMass();
        }

        centreOfMass /= mass_;

        forAll(sites_, sI)
        {
            sites_[sI].siteReferencePosition() -= centreOfMass;
        }

        diagTensor momOfInertia = diagTensor::zero;

        forAll(sites_, sI)
        {
            const vector& p(sites_[sI].siteReferencePosition());

            momOfInertia +=
                sites_[sI].siteMass()*diagTensor(0, p.x()*p.x(), p.x()*p.x());
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

        forAll(sites_, sI)
        {
            const vector& p(sites_[sI].siteReferencePosition());

            momOfInertia += sites_[sI].siteMass()*tensor
            (
                p.y()*p.y() + p.z()*p.z(), -p.x()*p.y(), -p.x()*p.z(),
                -p.y()*p.x(), p.x()*p.x() + p.z()*p.z(), -p.y()*p.z(),
                -p.z()*p.x(), -p.z()*p.y(), p.x()*p.x() + p.y()*p.y()
            );
        }

        if (eigenValues(momOfInertia).x() < VSMALL)
        {
            FatalErrorIn("polyatomic::constantProperties::constantProperties")
                << "An eigenvalue of the inertia tensor is zero, the "
                << "polyatomic " << dict.name() << " is not a valid 6DOF shape."
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

        forAll(sites_, sI)
        {
            sites_[sI].siteReferencePosition() =
                (Q & sites_[sI].siteReferencePosition());
        }

        // Recalculating the moment of inertia with the new site positions

        // The rotation was around the centre of mass but remove any
        // components that have crept in due to floating point errors

        centreOfMass = vector::zero;

        forAll(sites_, sI)
        {
            centreOfMass +=
                sites_[sI].siteReferencePosition()*sites_[sI].siteMass();
        }

        centreOfMass /= mass_;

        forAll(sites_, sI)
        {
            sites_[sI].siteReferencePosition() -= centreOfMass;
        }

        // Calculate the moment of inertia in the principle component
        // reference frame

        momOfInertia = tensor::zero;

        forAll(sites_, sI)
        {
            const vector& p(sites_[sI].siteReferencePosition());

            momOfInertia += sites_[sI].siteMass()*tensor
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
    const scalar temperature,
    const vector& bulkVelocity,
    const vector& specialPosition,
    const constantProperties& cP,
    polyatomic::trackingData& td,
    const label special,
    const label id
)
:
    particle(mesh, position, cellI, tetFaceI, tetPtI),
    Q_(I),
    v_(bulkVelocity),
    a_(vector::zero),
    pi_(vector::zero),
    tau_(vector::zero),
    specialPosition_(specialPosition),
    potentialEnergy_(0.0),
    rf_(tensor::zero),
    special_(special),
    id_(id),
    siteForces_(cP.nSites(), vector::zero),
    sitePositions_(cP.nSites())
{
    setSitePositions(cP);

    v_ += td.cloud().equipartitionLinearVelocity(temperature, cP.mass());

    if (!cP.pointMolecule())
    {
        pi_ = td.cloud().equipartitionAngularMomentum(temperature, cP);

        scalar phi(td.cloud().rndGen().scalar01()*twoPi);

        scalar theta(td.cloud().rndGen().scalar01()*twoPi);

        scalar psi(td.cloud().rndGen().scalar01()*twoPi);

        Q_ = tensor
        (
            cos(psi)*cos(phi) - cos(theta)*sin(phi)*sin(psi),
            cos(psi)*sin(phi) + cos(theta)*cos(phi)*sin(psi),
            sin(psi)*sin(theta),
            - sin(psi)*cos(phi) - cos(theta)*sin(phi)*cos(psi),
            - sin(psi)*sin(phi) + cos(theta)*cos(phi)*cos(psi),
            cos(psi)*sin(theta),
            sin(theta)*sin(phi),
            - sin(theta)*cos(phi),
            cos(theta)
        );
    }
}


Foam::polyatomic::polyatomic
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
    const constantProperties& cP,
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
    siteForces_(cP.nSites(), vector::zero),
    sitePositions_(cP.nSites())
{
    setSitePositions(cP);
}


// * * * * * * * constantProperties Private Member Functions * * * * * * * * //

inline bool Foam::polyatomic::constantProperties::linearMoleculeTest() const
{
    if (sites_.size() == 2)
    {
        return true;
    }

    vector refDir =
        sites_[1].siteReferencePosition()
      - sites_[0].siteReferencePosition();

    refDir /= mag(refDir);

    for
    (
        label i = 2;
        i < sites_.size();
        i++
    )
    {
        vector dir =
            sites_[i].siteReferencePosition()
          - sites_[i - 1].siteReferencePosition();

        dir /= mag(dir);

        if (mag(refDir & dir) < 1 - SMALL)
        {
            return false;
        }
    }

    return true;
}


// * * * * * * * constantProperties Member Functions * * * * * * * * * * * * //

inline const Foam::List<Foam::constPropSite>&
Foam::polyatomic::constantProperties::sites() const
{
    return sites_;
}


inline const Foam::List<Foam::label>&
Foam::polyatomic::constantProperties::pairPotSites() const
{
    return pairPotSites_;
}


inline const Foam::List<Foam::label>&
Foam::polyatomic::constantProperties::electrostaticSites() const
{
    return electrostaticSites_;
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


inline  bool Foam::polyatomic::constantProperties::pointMolecule() const
{
    return (momentOfInertia_.zz() < 0);
}


inline Foam::label
Foam::polyatomic::constantProperties::degreesOfFreedom() const
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
    return sites_.size();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

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

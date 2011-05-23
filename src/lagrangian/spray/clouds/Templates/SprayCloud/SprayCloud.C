/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2011 OpenCFD Ltd.
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

#include "SprayCloud.H"
#include "AtomizationModel.H"
#include "BreakupModel.H"
#include "CollisionModel.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class CloudType>
void Foam::SprayCloud<CloudType>::setModels()
{
    atomizationModel_.reset
    (
        AtomizationModel<SprayCloud<CloudType> >::New
        (
            this->subModelProperties(),
            *this
        ).ptr()
    );

    breakupModel_.reset
    (
        BreakupModel<SprayCloud<CloudType> >::New
        (
            this->subModelProperties(),
            *this
        ).ptr()
    );

    collisionModel_.reset
    (
        CollisionModel<SprayCloud<CloudType> >::New
        (
            this->subModelProperties(),
            *this
        ).ptr()
    );
}


template<class CloudType>
void Foam::SprayCloud<CloudType>::cloudReset
(
    SprayCloud<CloudType>& c
)
{
    CloudType::cloudReset(c);

    atomizationModel_.reset(c.atomizationModel_.ptr());
    breakupModel_.reset(c.breakupModel_.ptr());
    collisionModel_.reset(c.collisionModel_.ptr());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SprayCloud<CloudType>::SprayCloud
(
    const word& cloudName,
    const volScalarField& rho,
    const volVectorField& U,
    const dimensionedVector& g,
    const SLGThermo& thermo,
    bool readFields
)
:
    CloudType(cloudName, rho, U, g, thermo, false),
    sprayCloud(),
    cloudCopyPtr_(NULL),
    averageParcelMass_(this->injection().averageParcelMass()),
    atomizationModel_
    (
        AtomizationModel<SprayCloud<CloudType> >::New
        (
            this->particleProperties(),
            *this
        )
    ),
    breakupModel_
    (
        BreakupModel<SprayCloud<CloudType> >::New
        (
            this->particleProperties(),
            *this
        )
    ),
    collisionModel_
    (
        CollisionModel<SprayCloud<CloudType> >::New
        (
            this->particleProperties(),
            *this
        )
    )
{
    if (this->solution().active())
    {
        setModels();

        if (readFields)
        {
            parcelType::readFields(*this, this->composition());
        }
    }

    if (this->solution().resetSourcesOnStartup())
    {
        resetSourceTerms();
    }

    Info << "    Average parcel mass: " << averageParcelMass_ << endl;
}


template<class CloudType>
Foam::SprayCloud<CloudType>::SprayCloud
(
    SprayCloud<CloudType>& c,
    const word& name
)
:
    CloudType(c, name),
    sprayCloud(),
    cloudCopyPtr_(NULL),
    averageParcelMass_(c.averageParcelMass_),
    atomizationModel_(c.atomizationModel_->clone()),
    breakupModel_(c.breakupModel_->clone()),
    collisionModel_(c.collisionModel_->clone())
{}


template<class CloudType>
Foam::SprayCloud<CloudType>::SprayCloud
(
    const fvMesh& mesh,
    const word& name,
    const SprayCloud<CloudType>& c
)
:
    CloudType(mesh, name, c),
    sprayCloud(),
    cloudCopyPtr_(NULL),
    averageParcelMass_(0.0),
    atomizationModel_(NULL),
    breakupModel_(NULL),
    collisionModel_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SprayCloud<CloudType>::~SprayCloud()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::SprayCloud<CloudType>::checkParcelProperties
(
    parcelType& parcel,
    const scalar lagrangianDt,
    const bool fullyDescribed
)
{
    CloudType::checkParcelProperties
    (
        parcel,
        lagrangianDt,
        fullyDescribed
    );

    const scalarField& Y(parcel.Y());
    scalarField X(this->composition().liquids().X(Y));

    // override rho and cp from constantProperties
    parcel.Cp() = this->composition().liquids().Cp(parcel.pc(), parcel.T(), X);
    parcel.rho() = this->composition().liquids().rho(parcel.pc(), parcel.T(), X);

    // store the injection position and initial drop size
    parcel.position0() = parcel.position();
    parcel.d0() = parcel.d();

    parcel.y() = breakup().y0();
    parcel.yDot() = breakup().yDot0();
}


template<class CloudType>
void Foam::SprayCloud<CloudType>::storeState()
{
    cloudCopyPtr_.reset
    (
        static_cast<SprayCloud<CloudType>*>
        (
            clone(this->name() + "Copy").ptr()
        )
    );
}


template<class CloudType>
void Foam::SprayCloud<CloudType>::restoreState()
{
    cloudReset(cloudCopyPtr_());
    cloudCopyPtr_.clear();
}


template<class CloudType>
void Foam::SprayCloud<CloudType>::evolve()
{
    if (this->solution().canEvolve())
    {
        typename parcelType::template
            TrackingData<SprayCloud<CloudType> > td(*this);

        this->solve(td);
    }
}


template<class CloudType>
template<class TrackData>
void Foam::SprayCloud<CloudType>::motion(TrackData& td)
{
    const scalar dt = this->solution().trackTime();

    td.part() = TrackData::tpLinearTrack;
    CloudType::move(td, dt);

    this->updateCellOccupancy();


    if (collision().active())
    {
        label i = 0;
        forAllIter(typename SprayCloud<CloudType>, *this, iter)
        {
            label j = 0;
            forAllIter(typename SprayCloud<CloudType>, *this, jter)
            {
                if (j > i)
                {
                    parcelType& p = iter();
                    scalar Vi = this->mesh().V()[p.cell()];
                    scalarField X1(this->composition().liquids().X(p.Y()));
                    scalar sigma1 = this->composition().liquids().sigma(p.pc(), p.T(), X1);
                    scalar mp = p.mass()*p.nParticle();

                    parcelType& q = jter();
                    scalar Vj = this->mesh().V()[q.cell()];
                    scalarField X2(this->composition().liquids().X(q.Y()));
                    scalar sigma2 = this->composition().liquids().sigma(q.pc(), q.T(), X2);
                    scalar mq = q.mass()*q.nParticle();

                    bool updateProperties = collision().update
                    (
                        dt,
                        this->rndGen(),
                        p.position(),
                        mp,
                        p.d(),
                        p.nParticle(),
                        p.U(),
                        p.rho(),
                        p.T(),
                        p.Y(),
                        sigma1,
                        p.cell(),
                        Vi,
                        q.position(),
                        mq,
                        q.d(),
                        q.nParticle(),
                        q.U(),
                        q.rho(),
                        q.T(),
                        q.Y(),
                        sigma2,
                        q.cell(),
                        Vj
                    );

                    // for coalescence we need to update the density and
                    // the diameter cause of the temp/conc/mass-change
                    if (updateProperties)
                    {
                        if (mp > VSMALL)
                        {
                            scalarField Xp(this->composition().liquids().X(p.Y()));
                            p.rho() = this->composition().liquids().rho(p.pc(), p.T(), Xp);
                            p.Cp() = this->composition().liquids().Cp(p.pc(), p.T(), Xp);
                            scalar rhs = 6.0*mp/(p.nParticle()*p.rho()*constant::mathematical::pi);
                            p.d() = pow(rhs, 1.0/3.0);
                        }

                        if (mq > VSMALL)
                        {
                            scalarField Xq(this->composition().liquids().X(q.Y()));
                            q.rho() = this->composition().liquids().rho(q.pc(), q.T(), Xq);
                            q.Cp() = this->composition().liquids().Cp(q.pc(), q.T(), Xq);
                            scalar rhs = 6.0*mq/(q.nParticle()*q.rho()*constant::mathematical::pi);
                            q.d() = pow(rhs, 1.0/3.0);
                        }
                    }
                }
                j++;
            }
            i++;
        }

        // remove coalesced particles (diameter set to 0)
        forAllIter(typename SprayCloud<CloudType>, *this, iter)
        {
            parcelType& p = iter();
            if (p.mass() < VSMALL)
            {
                deleteParticle(p);
            }
        }
    }
}


template<class CloudType>
Foam::scalar Foam::SprayCloud<CloudType>::D
(
    const label i,
    const label j
) const
{
    scalar si = 0.0;
    scalar sj = 0.0;
    forAllConstIter(typename SprayCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();
        si += p.nParticle()*pow(p.d(), i);
        sj += p.nParticle()*pow(p.d(), j);
    }

    reduce(si, sumOp<scalar>());
    reduce(sj, sumOp<scalar>());
    sj = max(sj, VSMALL);

    return si/sj;
}


template<class CloudType>
Foam::scalar Foam::SprayCloud<CloudType>::liquidPenetration
(
    const scalar& prc
) const
{
    scalar distance = 0.0;
    scalar mTot = 0.0;

    label np = this->size();

    // arrays containing the parcels mass and
    // distance from injector in ascending order
    scalarField mass(np);
    scalarField dist(np);

    if (np > 0)
    {
        label n = 0;

        // first arrange the parcels in ascending order
        // the first parcel is closest to its injection position
        // and the last one is most far away.
        forAllConstIter(typename SprayCloud<CloudType>, *this, iter)
        {
            const parcelType& p = iter();
            scalar mi = p.nParticle()*p.mass();
            scalar di = mag(p.position() - p.position0());
            mTot += mi;

            // insert at the last place
            mass[n] = mi;
            dist[n] = di;

            label i = 0;
            bool found = false;

            // insert the parcel in the correct place
            // and move the others
            while (( i < n ) && (!found))
            {
                if (di < dist[i])
                {
                    found = true;
                    for (label j=n; j>i; j--)
                    {
                        mass[j] = mass[j-1];
                        dist[j] = dist[j-1];
                    }
                    mass[i] = mi;
                    dist[i] = di;
                }
                i++;
            }
            n++;
        }
    }

    reduce(mTot, sumOp<scalar>());

    if (np > 0)
    {
        scalar mLimit = prc*mTot;
        scalar mOff = (1.0 - prc)*mTot;

        if (np > 1)
        {
            // 'prc' is large enough that the parcel most far
            // away will be used, no need to loop...
            if (mLimit > mTot - mass[np-1])
            {
                distance = dist[np-1];
            }
            else
            {
                scalar mOffSum = 0.0;
                label i = np;

                while ((mOffSum < mOff) && (i>0))
                {
                    i--;
                    mOffSum += mass[i];
                }
                distance =
                    dist[i+1]
                  + (dist[i] - dist[i+1])*(mOffSum - mOff)
                   /mass[i+1] ;
            }
        }
        else
        {
            distance = dist[0];
        }
    }

    reduce(distance, maxOp<scalar>());

    return distance;
}


template<class CloudType>
void Foam::SprayCloud<CloudType>::info() const
{
    CloudType::info();
    scalar d32 = 1.0e+6*D(3, 2);
    scalar pen = liquidPenetration(0.95);

    Info << "    D32 (mu)                        = " << d32 << endl;
    Info << "    Liquid penetration 95% mass (m) = " << pen << endl;
}


// ************************************************************************* //

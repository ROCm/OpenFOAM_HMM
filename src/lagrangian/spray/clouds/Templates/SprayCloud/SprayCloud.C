/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2011 OpenCFD Ltd.
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
#include "PtrList.H"

template<class ParcelType>
void Foam::SprayCloud<ParcelType>::preEvolve()
{
    ReactingCloud<ParcelType>::preEvolve();
}


template<class ParcelType>
void Foam::SprayCloud<ParcelType>::evolveCloud()
{
    const volScalarField& T = this->carrierThermo().T();
    const volScalarField cp = this->carrierThermo().Cp();
    const volScalarField& p = this->carrierThermo().p();

    autoPtr<interpolation<scalar> > rhoInterp = interpolation<scalar>::New
    (
        this->interpolationSchemes(),
        this->rho()
    );

    autoPtr<interpolation<vector> > UInterp = interpolation<vector>::New
    (
        this->interpolationSchemes(),
        this->U()
    );

    autoPtr<interpolation<scalar> > muInterp = interpolation<scalar>::New
    (
        this->interpolationSchemes(),
        this->mu()
    );

    autoPtr<interpolation<scalar> > TInterp = interpolation<scalar>::New
    (
        this->interpolationSchemes(),
        T
    );

    autoPtr<interpolation<scalar> > cpInterp = interpolation<scalar>::New
    (
        this->interpolationSchemes(),
        cp
    );

    autoPtr<interpolation<scalar> > pInterp = interpolation<scalar>::New
    (
        this->interpolationSchemes(),
        p
    );

    typename ParcelType::trackData td
    (
        *this,
        constProps_,
        rhoInterp(),
        UInterp(),
        muInterp(),
        TInterp(),
        cpInterp(),
        pInterp(),
        this->g().value()
    );

    if (this->coupled())
    {
        resetSourceTerms();
    }

    if (collision().active())
    {

        label i = 0;
        scalar dt = this->db().time().deltaTValue();
        forAllIter(typename Cloud<ParcelType>, *this, iter)
        {

            label j = 0;
            forAllIter(typename Cloud<ParcelType>, *this, jter)
            {
                if (j > i)
                {
                    ParcelType& p = iter();
                    scalar Vi = this->mesh().V()[p.cell()];
                    scalarField X1(this->composition().liquids().X(p.Y()));
                    scalar sigma1 = this->composition().liquids().sigma(p.pc(), p.T(), X1);
                    scalar mp = p.mass()*p.nParticle();

                    ParcelType& q = jter();
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
                            p.cp() = this->composition().liquids().cp(p.pc(), p.T(), Xp);
                            scalar rhs = 6.0*mp/(p.nParticle()*p.rho()*mathematicalConstant::pi);
                            p.d() = pow(rhs, 1.0/3.0);
                        }

                        if (mq > VSMALL)
                        {
                            scalarField Xq(this->composition().liquids().X(q.Y()));
                            q.rho() = this->composition().liquids().rho(q.pc(), q.T(), Xq);
                            q.cp() = this->composition().liquids().cp(q.pc(), q.T(), Xq);
                            scalar rhs = 6.0*mq/(q.nParticle()*q.rho()*mathematicalConstant::pi);
                            q.d() = pow(rhs, 1.0/3.0);
                        }

                    }
                }
                j++;
            }

            i++;
        }

        // remove coalesced particles (diameter set to 0)
        forAllIter(typename Cloud<ParcelType>, *this, iter)
        {
            ParcelType& p = iter();
            if (p.mass() < VSMALL)
            {
                deleteParticle(p);
            }
        }
    }

    Cloud<ParcelType>::move(td);
    this->injection().inject(td);
}


template<class ParcelType>
void Foam::SprayCloud<ParcelType>::postEvolve()
{
    ReactingCloud<ParcelType>::postEvolve();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::SprayCloud<ParcelType>::SprayCloud
(
    const word& cloudName,
    const volScalarField& rho,
    const volVectorField& U,
    const dimensionedVector& g,
    basicThermo& thermo,
    bool readFields
)
:
    ReactingCloud<ParcelType>(cloudName, rho, U, g, thermo, false),
    sprayCloud(),
    averageParcelMass_(this->injection().averageParcelMass()),
    constProps_(this->particleProperties()),
    atomizationModel_
    (
        AtomizationModel<SprayCloud<ParcelType> >::New
        (
            this->particleProperties(),
            *this
        )
    ),
    breakupModel_
    (
        BreakupModel<SprayCloud<ParcelType> >::New
        (
            this->particleProperties(),
            *this
        )
    ),
    collisionModel_
    (
        CollisionModel<SprayCloud<ParcelType> >::New
        (
            this->particleProperties(),
            *this
        )
    )
{
    if (readFields)
    {
        ParcelType::readFields(*this);
    }

    Info << "    Average parcel mass: " << averageParcelMass_ << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::SprayCloud<ParcelType>::~SprayCloud()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParcelType>
void Foam::SprayCloud<ParcelType>::checkParcelProperties
(
    ParcelType& parcel,
    const scalar lagrangianDt,
    const bool fullyDescribed
)
{
    ReactingCloud<ParcelType>::checkParcelProperties
    (
        parcel,
        lagrangianDt,
        fullyDescribed
    );

    const scalarField& Y(parcel.Y());
    scalarField X(this->composition().liquids().X(Y));

    // override rho and cp from constantProperties
    parcel.cp() = this->composition().liquids().cp(parcel.pc(), parcel.T(), X);
    parcel.rho() = this->composition().liquids().rho(parcel.pc(), parcel.T(), X);

    // store the injection position and initial drop size
    parcel.position0() = parcel.position();
    parcel.d0() = parcel.d();

    parcel.y() = breakup().y0();
    parcel.yDot() = breakup().yDot0();
}

template<class ParcelType>
Foam::scalar Foam::SprayCloud<ParcelType>::D(const label i, const label j) const
{
    scalar si = 0.0;
    scalar sj = 0.0;
    forAllConstIter(typename Cloud<ParcelType>, *this, iter)
    {
        const ParcelType& p = iter();
        si += p.nParticle()*pow(p.d(), i);
        sj += p.nParticle()*pow(p.d(), j);
    }

    reduce(si, sumOp<scalar>());
    reduce(sj, sumOp<scalar>());
    sj = max(sj, VSMALL);

    return si/sj;
}


template<class ParcelType>
Foam::scalar Foam::SprayCloud<ParcelType>::liquidPenetration(const scalar& prc) const
{

    scalar distance = 0.0;
    scalar mTot = 0.0;

    label Np = this->size();

    // arrays containing the parcels mass and
    // distance from injector in ascending order
    scalarField mass(Np);
    scalarField dist(Np);

    if (Np > 0)
    {
        label n = 0;

        // first arrange the parcels in ascending order
        // the first parcel is closest to its injection position
        // and the last one is most far away.
        forAllConstIter(typename Cloud<ParcelType>, *this, iter)
        {
            const ParcelType& p = iter();
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
            while ( ( i < n ) && ( !found ) )
            {
                if (di < dist[i])
                {
                    found = true;
                    for(label j=n; j>i; j--)
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

    if (Np > 0)
    {

        scalar mLimit = prc*mTot;
        scalar mOff = (1.0 - prc)*mTot;

        if (Np > 1)
        {

            // 'prc' is large enough that the parcel most far
            // away will be used, no need to loop...
            if (mLimit > mTot - mass[Np-1])
            {
                distance = dist[Np-1];
            }
            else
            {
                scalar mOffSum = 0.0;
                label i = Np;

                while ((mOffSum < mOff) && (i>0))
                {
                    i--;
                    mOffSum += mass[i];
                }
                distance = dist[i+1] + (dist[i]-dist[i+1])*(mOffSum - mOff)/mass[i+1] ;
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

template<class ParcelType>
void Foam::SprayCloud<ParcelType>::resetSourceTerms()
{
    ReactingCloud<ParcelType>::resetSourceTerms();
}


template<class ParcelType>
void Foam::SprayCloud<ParcelType>::evolve()
{
    if (this->active())
    {
        preEvolve();

        evolveCloud();

        postEvolve();

        info();
        Info<< endl;
    }
}


template<class ParcelType>
void Foam::SprayCloud<ParcelType>::info() const
{
    ReactingCloud<ParcelType>::info();
    scalar d32 = 1.0e+6*D(3,2);
    scalar pen = liquidPenetration(0.95);

    Info << "    D32 (mu)                        = " << d32 << endl;
    Info << "    Liquid penetration 95% mass (m) = " << pen << endl;
}


// ************************************************************************* //

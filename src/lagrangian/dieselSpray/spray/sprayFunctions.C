/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "spray.H"
#include "mathematicalConstants.H"


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::spray::injectedMass(const scalar t) const
{
    scalar sum = 0.0;

    forAll(injectors_, i)
    {
        sum += injectors_[i].properties()->injectedMass(t);
    }

    return sum;
}


Foam::scalar Foam::spray::totalMassToInject() const
{
    scalar sum = 0.0;

    forAll(injectors_, i)
    {
        sum += injectors_[i].properties()->mass();
    }

    return sum;
}


Foam::scalar Foam::spray::injectedEnthalpy
(
    const scalar time
) const
{
    scalar sum = 0.0;
    label Nf = fuels_->components().size();

    forAll(injectors_, i)
    {
        scalar T = injectors_[i].properties()->T(time);
        scalarField X(injectors_[i].properties()->X());
        scalar pi = 1.0e+5;
        scalar hl = fuels_->hl(pi, T, X);
        scalar Wl = fuels_->W(X);
        scalar hg = 0.0;

        for (label j=0; j<Nf; j++)
        {
            label k = liquidToGasIndex_[j];
            hg += gasProperties()[k].H(T)*gasProperties()[k].W()*X[j]/Wl;
        }

        sum += injectors_[i].properties()->injectedMass(time)*(hg-hl);
    }

    return sum;
}


Foam::scalar Foam::spray::liquidMass() const
{
    scalar sum = 0.0;

    forAllConstIter(spray, *this, iter)
    {
        sum += iter().m();
    }

    if (twoD())
    {
        sum *= constant::mathematical::twoPi/angleOfWedge();
    }

    reduce(sum, sumOp<scalar>());

    return sum;
}


Foam::scalar Foam::spray::liquidEnthalpy() const
{
    scalar sum = 0.0;
    label Nf = fuels().components().size();

    forAllConstIter(spray, *this, iter)
    {
        scalar T = iter().T();
        scalar pc = p()[iter().cell()];
        scalar hlat = fuels().hl(pc, T, iter().X());
        scalar hg = 0.0;
        scalar Wl = fuels().W(iter().X());

        for (label j=0; j<Nf; j++)
        {
            label k = liquidToGasIndex_[j];

            hg +=
                gasProperties()[k].H(T)*gasProperties()[k].W()*iter().X()[j]
               /Wl;
        }

        scalar h = hg - hlat;
        sum += iter().m()*h;
    }

    if (twoD())
    {
        sum *= constant::mathematical::twoPi/angleOfWedge();
    }

    reduce(sum, sumOp<scalar>());

    return sum;
}


Foam::scalar Foam::spray::liquidTotalEnthalpy() const
{
    scalar sum = 0.0;
    label Nf = fuels().components().size();

    forAllConstIter(spray, *this, iter)
    {
        label cellI = iter().cell();
        scalar T = iter().T();
        scalar pc = p()[cellI];
        scalar rho = fuels().rho(pc, T, iter().X());
        scalar hlat = fuels().hl(pc, T, iter().X());
        scalar hg = 0.0;
        scalar Wl = fuels().W(iter().X());

        for (label j=0; j<Nf; j++)
        {
            label k = liquidToGasIndex_[j];
            hg +=
                gasProperties()[k].H(T)*gasProperties()[k].W()*iter().X()[j]
               /Wl;
        }

        scalar psat = fuels().pv(pc, T, iter().X());

        scalar h = hg - hlat + (pc - psat)/rho;
        sum += iter().m()*h;
    }

    if (twoD())
    {
        sum *= constant::mathematical::twoPi/angleOfWedge();
    }

    reduce(sum, sumOp<scalar>());

    return sum;
}


Foam::scalar Foam::spray::liquidKineticEnergy() const
{
    scalar sum = 0.0;

    forAllConstIter(spray, *this, iter)
    {
        const scalar ke = pow(mag(iter().U()), 2.0);
        sum += iter().m()*ke;
    }

    if (twoD())
    {
        sum *= constant::mathematical::twoPi/angleOfWedge();
    }

    reduce(sum, sumOp<scalar>());

    return 0.5*sum;

}


Foam::scalar Foam::spray::injectedLiquidKineticEnergy() const
{
    return injectedLiquidKE_;
}


Foam::scalar Foam::spray::liquidPenetration(const scalar prc) const
{
    return liquidPenetration(0, prc);
}


Foam::scalar Foam::spray::liquidPenetration
(
    const label nozzlei,
    const scalar prc
) const
{

    label nHoles = injectors_[nozzlei].properties()->nHoles();
    vector ip(vector::zero);
    if (nHoles > 1)
    {
        for (label i=0;i<nHoles;i++)
        {
            ip += injectors_[nozzlei].properties()->position(i);
        }
        ip /= nHoles;
    }
    else
    {
        ip = injectors_[nozzlei].properties()->position(0);
    }

//    vector ip = injectors_[nozzlei].properties()->position();
    scalar d = 0.0;
    scalar mTot = 0.0;

    label Np = size();

    // arrays containing the parcels mass and
    // distance from injector in ascending order
    scalarField m(Np);
    scalarField dist(Np);
    label n = 0;

    if (Np > 1)
    {
        // first arrange the parcels in ascending order
        // the first parcel is closest to injector
        // and the last one is most far away.
        spray::const_iterator first = begin();
        m[n] = first().m();
        dist[n] = mag(first().position() - ip);

        mTot += m[n];

        for
        (
            spray::const_iterator iter = ++first;
            iter != end();
            ++iter
        )
        {
            scalar de = mag(iter().position() - ip);
            scalar me = iter().m();
            mTot += me;

            n++;

            label i = 0;
            bool found = false;

            // insert the parcel in the correct place
            // and move the others
            while ( ( i < n-1 ) && ( !found ) )
            {
                if (de < dist[i])
                {
                    found = true;
                    for (label j=n; j>i; j--)
                    {
                        m[j]    = m[j-1];
                        dist[j] = dist[j-1];
                    }
                    m[i]    = me;
                    dist[i] = de;
                }
                i++;
            }

            if (!found)
            {
                m[n]    = me;
                dist[n] = de;
            }
        }
    }

    reduce(mTot, sumOp<scalar>());

    if (Np > 1)
    {
        scalar mLimit = prc*mTot;
        scalar mOff = (1.0 - prc)*mTot;

        // 'prc' is large enough that the parcel most far
        // away will be used, no need to loop...
        if (mLimit > mTot - m[Np-1])
        {
            d = dist[Np-1];
        }
        else
        {
            scalar mOffSum = 0.0;
            label i = Np;

            while ((mOffSum < mOff) && (i>0))
            {
                i--;
                mOffSum += m[i];
            }
            d = dist[i];
        }

    }
    else
    {
        if (Np > 0)
        {
            spray::const_iterator iter = begin();
            d = mag(iter().position() - ip);
        }
    }

    reduce(d, maxOp<scalar>());

    return d;
}


Foam::scalar Foam::spray::smd() const
{
    scalar numerator = 0.0, denominator = VSMALL;

    forAllConstIter(spray, *this, iter)
    {
        label cellI = iter().cell();
        scalar Pc = p()[cellI];
        scalar T = iter().T();
        scalar rho = fuels_->rho(Pc, T, iter().X());

        scalar tmp = iter().N(rho)*pow(iter().d(), 2.0);
        numerator += tmp*iter().d();
        denominator += tmp;
    }

    reduce(numerator, sumOp<scalar>());
    reduce(denominator, sumOp<scalar>());

    return numerator/denominator;
}


Foam::scalar Foam::spray::maxD() const
{
    scalar maxD = 0.0;

    forAllConstIter(spray, *this, iter)
    {
        maxD = max(maxD, iter().d());
    }

    reduce(maxD, maxOp<scalar>());

    return maxD;
}


void Foam::spray::calculateAmbientPressure()
{
    ambientPressure_ = p_.average().value();
}


void Foam::spray::calculateAmbientTemperature()
{
    ambientTemperature_ = T_.average().value();
}


// ************************************************************************* //

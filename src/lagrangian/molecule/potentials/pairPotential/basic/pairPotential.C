/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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

\*----------------------------------------------------------------------------*/

#include "pairPotential.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::pairPotential::setLookupTables()
{
    label N = label((rCut_ - rMin_)/dr_) + 1;

    forceLookup_.setSize(N);

    energyLookup_.setSize(N);

    forAll(forceLookup_, k)
    {
        forceLookup_[k] = force(k*dr_ + rMin_);

        energyLookup_[k] = energy(k*dr_ + rMin_);
    }

    forceLookup_[N-1] = 0.0;
    energyLookup_[N-1] = 0.0;
}


void Foam::pairPotential::setConstants()
{
    scalar nr = n(rCut_);

    u_at_rCut_ =
        epsilon_
       *(
            (6.0/(nr - 6.0))*Foam::pow( (rCut_/rm_), -nr)
          - (nr/(nr - 6.0))*Foam::pow( (rCut_/rm_), -6)
        );

    du_by_dr_at_rCut_ =
        -6.0 * epsilon_ * gamma_
       *(
            Foam::pow( (rCut_/rm_),-nr)
           *(
                (rCut_/rm_)
               *(1.0 / (nr - 6.0) + log(rCut_ / rm_) + 1.0)
              + (m_/gamma_)
              - 1.0
            )
          - Foam::pow((rCut_/rm_),-6.0)
            *(rCut_ / ( (nr - 6.0) * rm_) + (nr/gamma_))
        )
       /(rCut_*(nr - 6.0));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pairPotential::pairPotential()
{}


Foam::pairPotential::pairPotential(const dictionary& pPDict)
:
    forceLookup_(),
    energyLookup_(),
    m_(),
    gamma_(),
    rm_(),
    epsilon_(),
    rCut_(),
    rCutSqr_(),
    u_at_rCut_(),
    du_by_dr_at_rCut_(),
    rMin_(),
    dr_()
{

    m_ = readScalar(pPDict.lookup("m"));

    gamma_ = readScalar(pPDict.lookup("gamma"));

    rm_ = readScalar(pPDict.lookup("rm"));

    epsilon_ = readScalar(pPDict.lookup("epsilon"));

    rCut_ = readScalar(pPDict.lookup("rCut"));

    rCutSqr_ = rCut_ * rCut_;

    rMin_ = readScalar(pPDict.lookup("rMin"));

    dr_ = readScalar(pPDict.lookup("dr"));

    setConstants();

    setLookupTables();
}


Foam::pairPotential::pairPotential
(
    const scalar m,
    const scalar gamma,
    const scalar rm,
    const scalar epsilon,
    const scalar rCut,
    const scalar rMin,
    const scalar dr
)
:
    forceLookup_(),
    energyLookup_(),
    m_(m),
    gamma_(gamma),
    rm_(rm),
    epsilon_(epsilon),
    rCut_(rCut),
    rCutSqr_(rCut*rCut),
    u_at_rCut_(),
    du_by_dr_at_rCut_(),
    rMin_(rMin),
    dr_(dr)
{
    setConstants();

    setLookupTables();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pairPotential::~pairPotential()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::pairPotential::write(Ostream& os) const
{
    os<< "Shifted force Maitland Smith pair potential."
        << nl << tab << "m = " << m_
        << nl << tab << "gamma = " << gamma_
        << nl << tab << "rm = " << rm_
        << nl << tab << "epsilon = " << epsilon_
        << nl << tab << "rCut = " << rCut_
        << nl << tab << "rMin = " << rMin_
        << nl << tab << "dr = " << dr_
        << endl;
}


Foam::scalar Foam::pairPotential::n(const scalar rIJMag) const
{
    return (m_ + gamma_*(rIJMag/rm_ - 1.0));
}


Foam::scalar Foam::pairPotential::force(const scalar rIJMag) const
{
    scalar nr(n(rIJMag));

    scalar r = rIJMag/rm_;

    scalar rN = Foam::pow(r,-nr);

    scalar f =
      - epsilon_
       *(
           - 6.0/(nr - 6.0)/(nr - 6.0)*(gamma_/rm_)*rN
           + 6.0/(nr - 6.0)*rN
             *((gamma_ - m_)/rIJMag - (gamma_/rm_)*(log(r) + 1.0))
           + 6.0/(nr - 6.0)
             *(1.0/(nr - 6.0)*(gamma_/rm_) + nr/rIJMag)*Foam::pow(r, -6)
        )
      + du_by_dr_at_rCut_;

    return f;
}


Foam::scalar Foam::pairPotential::forceLookup(const scalar rIJMag) const
{
    scalar k_rIJ = (rIJMag - rMin_)/dr_;

    label k(k_rIJ);

    if (k < 0)
    {
        FatalErrorIn("pairPotential.C") << nl
        << "rIJMag less than rMin" << nl
        << abort(FatalError);
    }

    scalar f =
        (k_rIJ - k)*forceLookup_[k+1]
      + (k + 1 - k_rIJ)*forceLookup_[k];

    return f;
}


Foam::List< Foam::Pair< Foam::scalar > >
Foam::pairPotential::forceTable() const
{
    List< Pair<scalar> > forceTab(forceLookup_.size());

    forAll(forceLookup_,k)
    {
        forceTab[k].first() = rMin_ + k*dr_;

        forceTab[k].second() = forceLookup_[k];
    }

    return forceTab;
}


Foam::scalar Foam::pairPotential::energy(const scalar rIJMag) const
{
    scalar nr = n(rIJMag);

    scalar e =
        epsilon_
       *(
            (6.0 / (nr - 6.0))*Foam::pow( rIJMag/rm_, -nr)
          - (nr / (nr - 6.0))*Foam::pow( rIJMag/rm_, -6)
        )
      - u_at_rCut_
      - (rIJMag - rCut_)
       *du_by_dr_at_rCut_;

    return e;
}


Foam::scalar Foam::pairPotential::energyLookup(const scalar rIJMag) const
{
    scalar k_rIJ = (rIJMag - rMin_)/dr_;

    label k(k_rIJ);

    if (k < 0)
    {
        FatalErrorIn("pairPotential.C") << nl
        << "rIJMag less than rMin" << nl
        << abort(FatalError);
    }

    scalar e =
        (k_rIJ - k)*energyLookup_[k+1]
      + (k + 1 - k_rIJ)*energyLookup_[k];

    return e;
}


Foam::List< Foam::Pair< Foam::scalar > >
    Foam::pairPotential::energyTable() const
{
    List< Pair<scalar> > energyTab(energyLookup_.size());

    forAll(energyLookup_,k)
    {
        energyTab[k].first() = rMin_ + k*dr_;

        energyTab[k].second() = energyLookup_[k];
    }

    return energyTab;
}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const pairPotential& pot)
{
    pot.write(os);
    os.check
    (
        "Foam::Ostream& Foam::pperator<<(Ostream& f, const pairPotential& pot"
    );
    return os;
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2014-2016 OpenFOAM Foundation
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

#include "atmBoundaryLayer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

atmBoundaryLayer::atmBoundaryLayer(const Time& time, const polyPatch& pp)
:
    time_(time),
    patch_(pp),
    flowDir_(time, "flowDir"),
    zDir_(time, "zDir"),
    kappa_(0.41),
    Cmu_(0.09),
    Uref_(time, "Uref"),
    Zref_(time, "Zref"),
    z0_(),
    zGround_()
{}


atmBoundaryLayer::atmBoundaryLayer
(
    const Time& time,
    const polyPatch& pp,
    const dictionary& dict
)
:
    time_(time),
    patch_(pp),
    flowDir_(TimeFunction1<vector>(time, "flowDir", dict)),
    zDir_(TimeFunction1<vector>(time, "zDir", dict)),
    kappa_(dict.lookupOrDefault<scalar>("kappa", 0.41)),
    Cmu_(dict.lookupOrDefault<scalar>("Cmu", 0.09)),
    Uref_(TimeFunction1<scalar>(time, "Uref", dict)),
    Zref_(TimeFunction1<scalar>(time, "Zref", dict)),
    z0_(PatchFunction1<scalar>::New(pp, "z0", dict)),
    zGround_(PatchFunction1<scalar>::New(pp, "zGround", dict))
{}


atmBoundaryLayer::atmBoundaryLayer
(
    const atmBoundaryLayer& abl,
    const fvPatch& patch,
    const fvPatchFieldMapper& mapper
)
:
    time_(abl.time_),
    patch_(patch.patch()),
    flowDir_(abl.flowDir_),
    zDir_(abl.zDir_),
    kappa_(abl.kappa_),
    Cmu_(abl.Cmu_),
    Uref_(abl.Uref_),
    Zref_(abl.Zref_),
    z0_(abl.z0_.clone(patch_)),
    zGround_(abl.zGround_.clone(patch_))
{}


atmBoundaryLayer::atmBoundaryLayer(const atmBoundaryLayer& abl)
:
    time_(abl.time_),
    patch_(abl.patch_),
    flowDir_(abl.flowDir_),
    zDir_(abl.zDir_),
    kappa_(abl.kappa_),
    Cmu_(abl.Cmu_),
    Uref_(abl.Uref_),
    Zref_(abl.Zref_),
    z0_(abl.z0_.clone(patch_)),
    zGround_(abl.zGround_.clone(patch_))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

vector atmBoundaryLayer::flowDir() const
{
    const scalar t = time_.timeOutputValue();
    vector dir(flowDir_.value(t));
    const scalar magDir = mag(dir);

    if (magDir < SMALL)
    {
        FatalErrorInFunction
            << "magnitude of " << flowDir_.name()
            << " vector must be greater than zero"
            << abort(FatalError);
    }

    return dir/magDir;
}


vector atmBoundaryLayer::zDir() const
{
    const scalar t = time_.timeOutputValue();
    vector dir(zDir_.value(t));
    const scalar magDir = mag(dir);

    if (magDir < SMALL)
    {
        FatalErrorInFunction
            << "magnitude of " << zDir_.name()
            << " vector must be greater than zero"
            << abort(FatalError);
    }

    return dir/magDir;
}


tmp<scalarField> atmBoundaryLayer::Ustar(const scalarField& z0) const
{
    const scalar t = time_.timeOutputValue();
    const scalar Uref = Uref_.value(t);
    const scalar Zref = Zref_.value(t);

    return kappa_*Uref/(log((Zref + z0)/z0));
}


void atmBoundaryLayer::autoMap(const fvPatchFieldMapper& mapper)
{
    z0_->autoMap(mapper);
    zGround_->autoMap(mapper);
}


void atmBoundaryLayer::rmap
(
    const atmBoundaryLayer& abl,
    const labelList& addr
)
{
    z0_->rmap(abl.z0_(), addr);
    zGround_->rmap(abl.zGround_(), addr);
}


tmp<vectorField> atmBoundaryLayer::U(const vectorField& p) const
{
    const scalar t = time_.timeOutputValue();
    const scalarField zGround(zGround_->value(t));
    const scalarField z0(max(z0_->value(t), ROOTVSMALL));

    scalarField Un((Ustar(z0)/kappa_)*log(((zDir() & p) - zGround + z0)/z0));

    return flowDir()*Un;
}


tmp<scalarField> atmBoundaryLayer::k(const vectorField& p) const
{
    const scalar t = time_.timeOutputValue();
    const scalarField z0(max(z0_->value(t), ROOTVSMALL));

    return sqr(Ustar(z0))/sqrt(Cmu_);
}


tmp<scalarField> atmBoundaryLayer::epsilon(const vectorField& p) const
{
    const scalar t = time_.timeOutputValue();
    const scalarField zGround(zGround_->value(t));
    const scalarField z0(max(z0_->value(t), ROOTVSMALL));

    return pow3(Ustar(z0))/(kappa_*((zDir() & p) - zGround + z0));
}


void atmBoundaryLayer::write(Ostream& os) const
{
    z0_->writeData(os) ;
    flowDir_.writeData(os);
    zDir_.writeData(os);
    os.writeEntry("kappa", kappa_);
    os.writeEntry("Cmu", Cmu_);
    Uref_.writeData(os);
    Zref_.writeData(os);
    zGround_->writeData(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

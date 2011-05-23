/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2011 OpenCFD Ltd.
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

#include "BreakupModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::BreakupModel<CloudType>::BreakupModel
(
    CloudType& owner
)
:
    dict_(dictionary::null),
    owner_(owner),
    coeffDict_(dictionary::null),
    solveOscillationEq_(false),
    y0_(0.0),
    yDot0_(0.0),
    TABComega_(0.0),
    TABCmu_(0.0),
    TABWeCrit_(0.0)
{}


template<class CloudType>
Foam::BreakupModel<CloudType>::BreakupModel
(
    const dictionary& dict,
    CloudType& owner,
    const word& type
)
:
    dict_(dict),
    owner_(owner),
    coeffDict_(dict.subDict(type + "Coeffs")),
    solveOscillationEq_(dict_.lookup("solveOscillationEq")),
    y0_(0.0),
    yDot0_(0.0),
    TABComega_(0.0),
    TABCmu_(0.0),
    TABWeCrit_(0.0)
{
    if (solveOscillationEq_)
    {
        dictionary TABcoeffsDict(dict.subDict("TABCoeffs"));
        y0_ = readScalar(TABcoeffsDict.lookup("y0"));
        yDot0_ = readScalar(TABcoeffsDict.lookup("yDot0"));
        TABComega_ = readScalar(TABcoeffsDict.lookup("Comega"));
        TABCmu_ = readScalar(TABcoeffsDict.lookup("Cmu"));
        TABWeCrit_ = readScalar(TABcoeffsDict.lookup("WeCrit"));
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::BreakupModel<CloudType>::~BreakupModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
template<class CloudType>
const CloudType& Foam::BreakupModel<CloudType>::owner() const
{
    return owner_;
}


template<class CloudType>
const Foam::dictionary& Foam::BreakupModel<CloudType>::dict() const
{
    return dict_;
}


template<class CloudType>
const Foam::dictionary& Foam::BreakupModel<CloudType>::coeffDict() const
{
    return coeffDict_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "NewBreakupModel.C"

// ************************************************************************* //


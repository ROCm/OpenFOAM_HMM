/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2021 OpenCFD Ltd.
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

inline const Foam::primitivePatch& Foam::AMIInterpolation::srcPatch0() const
{
    if (!tsrcPatch0_.valid())
    {
        FatalErrorInFunction
            << "tsrcPatch0Ptr_ not set"
            << abort(FatalError);
    }

    return tsrcPatch0_();
}


inline const Foam::primitivePatch& Foam::AMIInterpolation::tgtPatch0() const
{

    if (!ttgtPatch0_.valid())
    {
        FatalErrorInFunction
            << "ttgtPatch0Ptr_ not set"
            << abort(FatalError);
    }

    return ttgtPatch0_();
}


inline bool Foam::AMIInterpolation::upToDate() const
{
    return upToDate_;
}


inline bool& Foam::AMIInterpolation::upToDate()
{
    return upToDate_;
}


inline bool Foam::AMIInterpolation::distributed() const
{
    return singlePatchProc_ == -1;
}


inline bool Foam::AMIInterpolation::requireMatch() const
{
    return requireMatch_;
}


inline bool Foam::AMIInterpolation::setRequireMatch(const bool flag)
{
    requireMatch_ = flag;
    return requireMatch_;
}


inline bool Foam::AMIInterpolation::mustMatchFaces() const
{
    return requireMatch_ && !applyLowWeightCorrection();
}


inline bool Foam::AMIInterpolation::reverseTarget() const
{
    return reverseTarget_;
}


inline Foam::scalar Foam::AMIInterpolation::lowWeightCorrection() const
{
    return lowWeightCorrection_;
}


inline bool Foam::AMIInterpolation::applyLowWeightCorrection() const
{
    return lowWeightCorrection_ > 0;
}


inline Foam::label Foam::AMIInterpolation::singlePatchProc() const
{
    return singlePatchProc_;
}


inline const Foam::List<Foam::scalar>& Foam::AMIInterpolation::srcMagSf() const
{
    return srcMagSf_;
}


inline Foam::List<Foam::scalar>& Foam::AMIInterpolation::srcMagSf()
{
    return srcMagSf_;
}


inline const Foam::labelListList& Foam::AMIInterpolation::srcAddress() const
{
    return srcAddress_;
}


inline Foam::labelListList& Foam::AMIInterpolation::srcAddress()
{
    return srcAddress_;
}


inline const Foam::scalarListList& Foam::AMIInterpolation::srcWeights() const
{
    return srcWeights_;
}


inline Foam::scalarListList& Foam::AMIInterpolation::srcWeights()
{
    return srcWeights_;
}


inline const Foam::scalarField& Foam::AMIInterpolation::srcWeightsSum() const
{
    return srcWeightsSum_;
}


inline Foam::scalarField& Foam::AMIInterpolation::srcWeightsSum()
{
    return srcWeightsSum_;
}


inline const Foam::pointListList& Foam::AMIInterpolation::srcCentroids() const
{
    return srcCentroids_;
}


inline Foam::pointListList& Foam::AMIInterpolation::srcCentroids()
{
    return srcCentroids_;
}


inline const Foam::mapDistribute& Foam::AMIInterpolation::srcMap() const
{
    return *srcMapPtr_;
}


inline const Foam::List<Foam::scalar>& Foam::AMIInterpolation::tgtMagSf() const
{
    return tgtMagSf_;
}


inline Foam::List<Foam::scalar>& Foam::AMIInterpolation::tgtMagSf()
{
    return tgtMagSf_;
}


inline const Foam::labelListList& Foam::AMIInterpolation::tgtAddress() const
{
    return tgtAddress_;
}


inline Foam::labelListList& Foam::AMIInterpolation::tgtAddress()
{
    return tgtAddress_;
}


inline const Foam::scalarListList& Foam::AMIInterpolation::tgtWeights() const
{
    return tgtWeights_;
}


inline Foam::scalarListList& Foam::AMIInterpolation::tgtWeights()
{
    return tgtWeights_;
}


inline const Foam::scalarField& Foam::AMIInterpolation::tgtWeightsSum() const
{
    return tgtWeightsSum_;
}


inline Foam::scalarField& Foam::AMIInterpolation::tgtWeightsSum()
{
    return tgtWeightsSum_;
}


inline const Foam::mapDistribute& Foam::AMIInterpolation::tgtMap() const
{
    return *tgtMapPtr_;
}

// ************************************************************************* //

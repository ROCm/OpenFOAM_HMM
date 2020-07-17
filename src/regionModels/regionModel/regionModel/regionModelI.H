/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "regionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline const Foam::fvMesh&
Foam::regionModels::regionModel::primaryMesh() const
{
    return primaryMesh_;
}


inline const Foam::Time& Foam::regionModels::regionModel::time() const
{
    return time_;
}


inline Foam::Switch Foam::regionModels::regionModel::active() const
{
    return active_;
}


inline Foam::Switch Foam::regionModels::regionModel::infoOutput() const
{
    return infoOutput_;
}


inline const Foam::word& Foam::regionModels::regionModel::modelName() const
{
    return modelName_;
}


inline const Foam::fvMesh& Foam::regionModels::regionModel::regionMesh() const
{
    const fvMesh* regionPtr = time_.getObjectPtr<fvMesh>(regionName_);

    if (!regionPtr)
    {
        FatalErrorInFunction
            << "Region mesh not available" << abort(FatalError);
    }

    return *regionPtr;
}


inline Foam::fvMesh& Foam::regionModels::regionModel::regionMesh()
{
    fvMesh* regionPtr = time_.getObjectPtr<fvMesh>(regionName_);

    if (!regionPtr)
    {
        FatalErrorInFunction
            << "Region mesh not available" << abort(FatalError);
    }

    return *regionPtr;
}


inline const Foam::dictionary& Foam::regionModels::regionModel::coeffs() const
{
    return coeffs_;
}


inline const Foam::dictionary&
Foam::regionModels::regionModel::solution() const
{
    return regionMesh().solutionDict();
}


inline const Foam::IOdictionary&
Foam::regionModels::regionModel::outputProperties() const
{
    if (!outputPropertiesPtr_)
    {
        FatalErrorInFunction
            << "outputProperties dictionary not available"
            << abort(FatalError);
    }

    return *outputPropertiesPtr_;
}


inline Foam::IOdictionary&
Foam::regionModels::regionModel::outputProperties()
{
    if (!outputPropertiesPtr_)
    {
        FatalErrorInFunction
            << "outputProperties dictionary not available"
            << abort(FatalError);
    }

    return *outputPropertiesPtr_;
}


inline bool Foam::regionModels::regionModel::isCoupledPatch
(
    const label regionPatchi
) const
{
    forAll(intCoupledPatchIDs_, i)
    {
        if (intCoupledPatchIDs_[i] == regionPatchi)
        {
            return true;
        }
    }

    return false;
}


inline bool Foam::regionModels::regionModel::isRegionPatch
(
    const label primaryPatchi
) const
{
    forAll(primaryPatchIDs_, i)
    {
        if (primaryPatchIDs_[i] == primaryPatchi)
        {
            return true;
        }
    }

    return false;
}


inline const Foam::labelList&
Foam::regionModels::regionModel::primaryPatchIDs() const
{
    return primaryPatchIDs_;
}


inline const Foam::labelList&
Foam::regionModels::regionModel::intCoupledPatchIDs() const
{
    return intCoupledPatchIDs_;
}


inline Foam::label Foam::regionModels::regionModel::regionPatchID
(
    const label primaryPatchID
) const
{
    forAll(primaryPatchIDs_, i)
    {
        if (primaryPatchIDs_[i] == primaryPatchID)
        {
            return intCoupledPatchIDs_[i];
        }
    }

    return -1;
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

#include "regionFaModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline const Foam::fvMesh&
Foam::regionModels::regionFaModel::primaryMesh() const
{
    return primaryMesh_;
}


inline const Foam::Time& Foam::regionModels::regionFaModel::time() const
{
    return time_;
}


inline const Foam::Switch& Foam::regionModels::regionFaModel::active() const
{
    return active_;
}


inline const Foam::Switch& Foam::regionModels::regionFaModel::infoOutput() const
{
    return infoOutput_;
}


inline const Foam::word& Foam::regionModels::regionFaModel::modelName() const
{
    return modelName_;
}


inline const Foam::faMesh& Foam::regionModels::regionFaModel::regionMesh() const
{
    const auto* regionPtr = time_.findObject<faMesh>(regionName_);

    if (regionPtr)
    {
        return *regionPtr;
    }
    else if (!regionMeshPtr_.valid())
    {
        FatalErrorInFunction
            << "Region mesh not available" << abort(FatalError);
    }

    return *regionMeshPtr_;
}


inline Foam::faMesh& Foam::regionModels::regionFaModel::regionMesh()
{
    auto* regionPtr = time_.getObjectPtr<faMesh>(regionName_);

    if (regionPtr)
    {
        return *regionPtr;
    }
    else if (!regionMeshPtr_.valid())
    {
        FatalErrorInFunction
            << "Region mesh not available" << abort(FatalError);
    }

    return *regionMeshPtr_;
}


inline const Foam::dictionary& Foam::regionModels::regionFaModel::coeffs() const
{
    return coeffs_;
}

inline const Foam::IOdictionary&
Foam::regionModels::regionFaModel::outputProperties() const
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
Foam::regionModels::regionFaModel::outputProperties()
{
    if (!outputPropertiesPtr_)
    {
        FatalErrorInFunction
            << "outputProperties dictionary not available"
            << abort(FatalError);
    }

    return *outputPropertiesPtr_;
}


inline const Foam::dictionary&
Foam::regionModels::regionFaModel::solution() const
{
    return regionMesh().solutionDict();
}


inline Foam::label Foam::regionModels::regionFaModel::patchID() const
{
    return patchID_;
}


// ************************************************************************* //

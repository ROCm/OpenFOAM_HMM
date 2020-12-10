/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2020 OpenCFD Ltd.
------------------------------------------------------------------------------
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::word& Foam::fa::option::name() const
{
    return name_;
}


inline const Foam::fvMesh& Foam::fa::option::mesh() const
{
    return mesh_;
}


inline const Foam::fvPatch& Foam::fa::option::patch() const
{
    return patch_;
}


inline const Foam::dictionary& Foam::fa::option::coeffs() const
{
    return coeffs_;
}


inline bool Foam::fa::option::active() const
{
    return active_;
}


inline void Foam::fa::option::setApplied(const label fieldi)
{
    applied_[fieldi] = true;
}


inline Foam::Switch& Foam::fa::option::active()
{
    return active_;
}


inline const Foam::word& Foam::fa::option::regionName() const
{
    return regionName_;
}


inline const Foam::faMesh& Foam::fa::option::regionMesh() const
{
    if (regionMeshPtr_.valid())
    {
        return regionMeshPtr_();
    }
    else
    {
        FatalErrorInFunction
            << "Region mesh not available" << abort(FatalError);
    }
    return *(new faMesh(mesh_));
}


inline const Foam::volSurfaceMapping& Foam::fa::option::vsm() const
{
    if (vsmPtr_.valid())
    {
        return vsmPtr_();
    }
    else
    {
        FatalErrorInFunction
            << "vsmPtr not available" << abort(FatalError);
    }
    return *(new volSurfaceMapping(regionMeshPtr_()));
}


// ************************************************************************* //

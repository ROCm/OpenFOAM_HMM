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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline const Foam::word& Foam::cyclicAMIPolyPatch::nbrPatchName() const
{
    return nbrPatchName_;
}


inline const Foam::cyclicAMIPolyPatch&
Foam::cyclicAMIPolyPatch::nbrPatch() const
{
    const polyPatch& pp = this->boundaryMesh()[nbrPatchID()];
    return refCast<const cyclicAMIPolyPatch>(pp);
}


inline Foam::coupledPolyPatch::transformType
Foam::cyclicAMIPolyPatch::transform() const
{
    return transform_;
}


inline const Foam::vector& Foam::cyclicAMIPolyPatch::rotationAxis() const
{
    return rotationAxis_;
}


inline const Foam::point& Foam::cyclicAMIPolyPatch::rotationCentre() const
{
    return rotationCentre_;
}


inline const Foam::vector& Foam::cyclicAMIPolyPatch::separationVector() const
{
    return separationVector_;
}


// ************************************************************************* //

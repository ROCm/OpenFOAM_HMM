/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013 OpenFOAM Foundation
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

template<class SourcePatch, class TargetPatch>
inline const SourcePatch&
Foam::AMIMethod<SourcePatch, TargetPatch>::srcPatch() const
{
    return srcPatch0_;
}


template<class SourcePatch, class TargetPatch>
inline const TargetPatch&
Foam::AMIMethod<SourcePatch, TargetPatch>::tgtPatch() const
{
    if (extendedTgtPatchPtr_)
    {
        return extendedTgtPatchPtr_();
    }

    return tgtPatch0_;
}


template<class SourcePatch, class TargetPatch>
bool Foam::AMIMethod<SourcePatch, TargetPatch>::distributed() const
{
    return singleMeshProc_ == -1;
}


template<class SourcePatch, class TargetPatch>
inline const Foam::labelList&
Foam::AMIMethod<SourcePatch, TargetPatch>::srcNonOverlap() const
{
    return srcNonOverlap_;
}


// ************************************************************************* //

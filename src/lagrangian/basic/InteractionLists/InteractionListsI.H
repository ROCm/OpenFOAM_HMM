/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2010 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParticleType>
const Foam::polyMesh& Foam::InteractionLists<ParticleType>::mesh() const
{
    return mesh_;
}


template<class ParticleType>
const Foam::mapDistribute&
Foam::InteractionLists<ParticleType>::map() const
{
    return mapPtr_();
}


template<class ParticleType>
const Foam::globalIndexAndTransform&
Foam::InteractionLists<ParticleType>::globalTransforms() const
{
    return globalTransforms_;
}


template<class ParticleType>
const Foam::labelListList& Foam::InteractionLists<ParticleType>::dil() const
{
    return dil_;
}


template<class ParticleType>
const Foam::labelListList&
Foam::InteractionLists<ParticleType>::directWallFaces() const
{
    return directWallFaces_;
}


template<class ParticleType>
const Foam::labelListList& Foam::InteractionLists<ParticleType>::ril() const
{
    return ril_;
}


template<class ParticleType>
const Foam::List<Foam::labelPair>&
Foam::InteractionLists<ParticleType>::cellIndexAndTransformToDistribute() const
{
    return cellIndexAndTransformToDistribute_;
}


template<class ParticleType>
const Foam::labelListList&
Foam::InteractionLists<ParticleType>::rilInverse() const
{
    return rilInverse_;
}


template<class ParticleType>
const Foam::List<Foam::IDLList<ParticleType> >&
Foam::InteractionLists<ParticleType>::referredParticles() const
{
    return referredParticles_;
}


template<class ParticleType>
Foam::List<Foam::IDLList<ParticleType> >&
Foam::InteractionLists<ParticleType>::referredParticles()
{
    return referredParticles_;
}




// ************************************************************************* //

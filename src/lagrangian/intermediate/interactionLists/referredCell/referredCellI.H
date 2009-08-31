/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2009 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParticleType>
inline Foam::label Foam::referredCell<ParticleType>::sourceProc() const
{
    return sourceProc_;
}


template<class ParticleType>
inline Foam::label Foam::referredCell<ParticleType>::sourceCell() const
{
    return sourceCell_;
}


template<class ParticleType>
inline const Foam::vector& Foam::referredCell<ParticleType>::offset() const
{
    return offset_;
}


template<class ParticleType>
inline const Foam::tensor& Foam::referredCell<ParticleType>::rotation() const
{
    return rotation_;
}


template<class ParticleType>
inline const Foam::vectorList&
Foam::referredCell<ParticleType>::vertexPositions() const
{
    return vertexPositions_;
}


template<class ParticleType>
inline const Foam::edgeList& Foam::referredCell<ParticleType>::edges() const
{
    return edges_;
}


template<class ParticleType>
inline const Foam::labelListList&
Foam::referredCell<ParticleType>::faces() const
{
    return faces_;
}


template<class ParticleType>
inline const Foam::vectorList&
Foam::referredCell<ParticleType>::faceCentres() const
{
    return faceCentres_;
}


template<class ParticleType>
inline const Foam::vectorList&
Foam::referredCell<ParticleType>::faceAreas() const
{
    return faceAreas_;
}


template<class ParticleType>
inline Foam::labelList& Foam::referredCell<ParticleType>::realCells()
{
    return realCellsForInteraction_;
}


template<class ParticleType>
inline const Foam::labelList&
Foam::referredCell<ParticleType>::realCellsForInteraction() const
{
    return realCellsForInteraction_;
}


// ************************************************************************* //

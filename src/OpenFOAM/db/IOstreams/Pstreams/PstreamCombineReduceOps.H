/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
    Copyright (C) 2022 OpenCFD Ltd.
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

InClass
    Foam

Description
    Combination-Reduction operation for a parallel run.  The
    information from all nodes is collected on the master node,
    combined using the given combination function and the result is
    broadcast to all nodes

\*---------------------------------------------------------------------------*/

#ifndef Foam_PstreamCombineReduceOps_H
#define Foam_PstreamCombineReduceOps_H

#include "UPstream.H"
#include "Pstream.H"
#include "ops.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class T, class CombineOp>
void combineReduce
(
    const List<UPstream::commsStruct>& comms,
    T& value,
    const CombineOp& cop,
    const int tag,
    const label comm
)
{
    Pstream::combineGather(comms, value, cop, tag, comm);
    Pstream::broadcast(value, comm);
}


template<class T, class CombineOp>
void combineReduce
(
    T& value,
    const CombineOp& cop,
    const int tag = UPstream::msgType(),
    const label comm = UPstream::worldComm
)
{
    if (UPstream::parRun() && UPstream::nProcs(comm) > 1)
    {
        const auto& comms = UPstream::whichCommunication(comm);

        Pstream::combineGather(comms, value, cop, tag, comm);
        Pstream::broadcast(value, comm);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

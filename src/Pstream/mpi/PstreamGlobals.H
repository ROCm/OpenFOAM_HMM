/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2015 OpenFOAM Foundation
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

Namespace
    Foam::PstreamGlobals

Description
    Global functions and variables for working with parallel streams,
    but principally for mpi

SourceFiles
    PstreamGlobals.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_PstreamGlobals_H
#define Foam_PstreamGlobals_H

#include "DynamicList.H"
#include <mpi.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace PstreamGlobals
{

//- Outstanding non-blocking operations.
extern DynamicList<MPI_Request> outstandingRequests_;
extern DynamicList<label> freedRequests_;

//- Max outstanding message tag operations.
extern int nTags_;

//- Free'd message tags
extern DynamicList<int> freedTags_;

// Current communicators, which may be allocated or predefined
// (eg, MPI_COMM_SELF, MPI_COMM_WORLD)
extern DynamicList<MPI_Comm> MPICommunicators_;

// Groups associated with the currrent communicators.
extern DynamicList<MPI_Group> MPIGroups_;


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

//- Fatal if comm is outside the allocated range
void checkCommunicator(const label comm, const label toProcNo);

//- Reuse previously freed request locations or push request onto list
//- of outstanding requests.
//
//  \return index of request within outstandingRequests_
inline label push_request(MPI_Request request)
{
    label index;

    if (freedRequests_.size())
    {
        index = freedRequests_.back();
        freedRequests_.pop_back();
        outstandingRequests_[index] = request;
    }
    else
    {
        index = outstandingRequests_.size();
        outstandingRequests_.push_back(request);
    }

    return index;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace PstreamGlobals
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

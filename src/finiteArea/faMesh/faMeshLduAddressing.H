/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2016-2017 Wikki Ltd
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

Class
    Foam::faMeshLduAddressing

Description
    lduAddressing wrapper for faMesh

SourceFiles
    faMeshLduAddressing.C

Author
    Zeljko Tukovic, FMENA
    Hrvoje Jasak, Wikki Ltd.

\*---------------------------------------------------------------------------*/

#ifndef faMeshLduAddressing_H
#define faMeshLduAddressing_H

#include "lduAddressing.H"
#include "faMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class faMeshLduAddressing Declaration
\*---------------------------------------------------------------------------*/

class faMeshLduAddressing
:
    public lduAddressing
{
    // Private data

        //- Lower as a subList of allOwner
        labelList::subList lowerAddr_;

        //- Upper as a reference to neighbour
        const labelList& upperAddr_;

        //- Patch addressing as a list of sublists
        List<const labelUList*> patchAddr_;

        //- Patch field evaluation schedule
        const lduSchedule& patchSchedule_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        faMeshLduAddressing(const faMeshLduAddressing&);

        //- Disallow default bitwise assignment
        void operator=(const faMeshLduAddressing&);


public:

    // Constructors

        //- Construct from components
        faMeshLduAddressing(const faMesh& mesh)
        :
            lduAddressing(mesh.nFaces()),
            lowerAddr_
            (
                labelList::subList
                (
                    mesh.edgeOwner(),
                    mesh.nInternalEdges()
                )
            ),
            upperAddr_(mesh.edgeNeighbour()),
            patchAddr_(mesh.boundary().size()),
            patchSchedule_(mesh.globalData().patchSchedule())
        {
            forAll(mesh.boundary(), patchI)
            {
                patchAddr_[patchI] = &mesh.boundary()[patchI].edgeFaces();
            }
        }

    //-Destructor
    virtual ~faMeshLduAddressing()
    {}


    // Member Functions

        //- Return number of interfaces
        virtual label nPatches() const
        {
            return patchAddr_.size();
        }

        //- Return lower addressing (i.e. lower label = upper triangle)
        virtual const labelUList& lowerAddr() const
        {
            return lowerAddr_;
        }

        //- Return upper addressing (i.e. upper label)
        virtual const labelUList& upperAddr() const
        {
            return upperAddr_;
        }

        //- Return patch addressing
        virtual const labelUList& patchAddr(const label i) const
        {
            return *patchAddr_[i];
        }

        // Return patch field evaluation schedule
        virtual const lduSchedule& patchSchedule() const
        {
            return patchSchedule_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

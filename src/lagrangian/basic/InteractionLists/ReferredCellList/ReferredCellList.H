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

Class
    Foam::ReferredCellList

Description

SourceFiles
    ReferredCellListI.H
    ReferredCellList.C

\*---------------------------------------------------------------------------*/

#ifndef ReferredCellList_H
#define ReferredCellList_H

#include "ReferredCell.H"
#include "sendingReferralList.H"
#include "receivingReferralList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
template<class ParticleType>
class InteractionLists;

/*---------------------------------------------------------------------------*\
                      Class ReferredCellList Declaration
\*---------------------------------------------------------------------------*/

template<class ParticleType>
class ReferredCellList
:
    public List<ReferredCell<ParticleType> >
{
    // Private data

        const InteractionLists<ParticleType>& il_;

        //- Dummy cloud to give to particles
        Cloud<ParticleType> cloud_;

        List<sendingReferralList> cellSendingReferralLists_;

        List<receivingReferralList> cellReceivingReferralLists_;



    // Private Member Functions

        //- Refer and rerefer cells to construct referred interactions
        void buildReferredCellList
        (
            bool pointPointListBuild
        );

        //- Build referralLists which define how to send information
        //  to ReferredCells to source cells
        void buildCellReferralLists();

        //- Store the referred particles in the referred cells
        void storeParticles
        (
            const receivingReferralList& rRL,
            const labelList& destinationReferredCell,
            IDLList<ParticleType>& particlesToReferIn
        );


public:

    // Constructors

        //- Construct lists by searching the mesh
        ReferredCellList
        (
            InteractionLists<ParticleType>& il,
            bool pointPointListBuild
        );

        //- Construct from file
        ReferredCellList (InteractionLists<ParticleType>& il);


    // Destructor

        ~ReferredCellList();


    // Member Functions

        void referParticles
        (
            const List<DynamicList<ParticleType*> >& cellOccupancy
        );

        inline const InteractionLists<ParticleType>& il() const;

        inline label nInteractingProcs() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ReferredCellListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "ReferredCellList.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

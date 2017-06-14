/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
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

Class
    Foam::oversetGAMGInterface

Description
    GAMG agglomerated overset interface.

SourceFiles
    oversetGAMGInterface.C

\*---------------------------------------------------------------------------*/

#ifndef oversetGAMGInterface_H
#define oversetGAMGInterface_H

#include "GAMGInterface.H"
#include "oversetLduInterface.H"
#include "mapDistribute.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class oversetGAMGInterface Declaration
\*---------------------------------------------------------------------------*/

class oversetGAMGInterface
:
    public GAMGInterface,
    virtual public oversetLduInterface
{
    // Private data

        //- Reference for the oversetLduInterface from which this is
        //  agglomerated
        const oversetLduInterface& fineOversetInterface_;

        //- Temporary copy of the fine level restrict map. Cleared upon
        //  calculating stencils below
        mutable labelList restrictMap_;

        labelListList stencil_;

        autoPtr<mapDistribute> mapPtr_;

        List<scalarList> cellInterpolationWeights_;

        scalarField normalisation_;

        labelList interpolationCells_;

        scalarList cellInterpolationWeight_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        oversetGAMGInterface(const oversetGAMGInterface&);

        //- Disallow default bitwise assignment
        void operator=(const oversetGAMGInterface&);


public:

    //- Runtime type information
    TypeName("overset");


    // Constructors

        //- Construct from fine level interface,
        //  local and neighbour restrict addressing
        oversetGAMGInterface
        (
            const label index,
            const lduInterfacePtrsList& coarseInterfaces,
            const lduInterface& fineInterface,
            const labelField& localRestrictAddressing,
            const labelField& neighbourRestrictAddressing,
            const label fineLevelIndex,
            const label coarseComm
        );


    //- Destructor
    virtual ~oversetGAMGInterface();


    // Member Functions

        // Interface transfer functions

            //- Transfer and return internal field adjacent to the interface
            virtual tmp<labelField> internalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const labelUList& iF
            ) const;


        //- Overset interface functions

            //- Name of lduInterface
            virtual const word& name() const
            {
                return fineOversetInterface_.name();
            }

            //- Am I the master interface
            virtual bool master() const
            {
                return fineOversetInterface_.master();
            }

            //- GAMG restriction (fine-to-coarse)
            virtual const labelList& restrictMap() const
            {
                return restrictMap_;
            }

            //- Donor stencil
            virtual const labelListList& stencil() const
            {
                return stencil_;
            }

            //- Map for obtaining data in stencil order
            virtual const mapDistribute& cellInterpolationMap() const
            {
                return mapPtr_();
            }

            //- Weights in stencil order
            virtual const List<scalarList>& cellInterpolationWeights() const
            {
                return cellInterpolationWeights_;
            }

            //- Normalisation ot matrix; for explicit contributions
            virtual const scalarField& normalisation() const
            {
                return normalisation_;
            }

            //- Acceptor cells
            virtual const labelList& interpolationCells() const
            {
                return interpolationCells_;
            }

            //- Underrelaxation for acceptor cells
            virtual const scalarList& cellInterpolationWeight() const
            {
                return cellInterpolationWeight_;
            }


        // I/O

            //- Write to stream
            virtual void write(Ostream&) const
            {
                //TBD.
                NotImplemented;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

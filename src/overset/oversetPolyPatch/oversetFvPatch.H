/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2017 OpenCFD Ltd.
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
    Foam::oversetFvPatch

Description
    Patch for indicating interpolated boundaries (in overset meshes).

SourceFiles
    oversetFvPatch.C

\*---------------------------------------------------------------------------*/

#ifndef oversetFvPatch_H
#define oversetFvPatch_H

#include "fvPatch.H"
#include "oversetPolyPatch.H"
#include "fvBoundaryMesh.H"
#include "lduInterface.H"
#include "oversetLduInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class oversetFvPatch Declaration
\*---------------------------------------------------------------------------*/

class oversetFvPatch
:
    public lduInterface,
    public oversetLduInterface,
    public fvPatch
{
    // Private data

        const oversetPolyPatch& oversetPolyPatch_;

        //- Temporary copy of the fine level restrict map. Cleared upon
        //  calculating stencils below
        mutable labelList restrictMap_;


public:

    //- Runtime type information
    TypeName(oversetPolyPatch::typeName_());


    // Constructors

        //- Construct from polyPatch
        oversetFvPatch(const polyPatch& patch, const fvBoundaryMesh& bm)
        :
            fvPatch(patch, bm),
            oversetPolyPatch_(refCast<const oversetPolyPatch>(patch))
        {}


    // Member functions

        // Access

            //- Return local reference cast into the overset patch
            const oversetPolyPatch& oversetPatch() const
            {
                return oversetPolyPatch_;
            }


        // Access

            //- Return faceCell addressing
            virtual const labelUList& faceCells() const
            {
                return oversetPolyPatch_.faceCells();
            }


        // Interface transfer functions

            //- Return the values of the given internal data adjacent to
            //  the interface as a field
            virtual tmp<labelField> interfaceInternalField
            (
                const labelUList& internalData
            ) const;

            //- Return neighbour field
            virtual tmp<labelField> internalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const labelUList& internalData
            ) const;


        //- Overset interface functions

            //- Name of interface (for debugging)
            virtual const word& name() const
            {
                return oversetPolyPatch_.name();
            }

            //- Am I the master interface
            virtual bool master() const
            {
                return oversetPolyPatch_.master();
            }

            //- GAMG restriction (fine-to-coarse)
            virtual const labelList& restrictMap() const
            {
                return restrictMap_;
            }

            //- Donor stencil
            virtual const labelListList& stencil() const;

            //- Map for obtaining data in stencil order
            virtual const mapDistribute& cellInterpolationMap() const;

            //- Weights in stencil order
            virtual const List<scalarList>& cellInterpolationWeights() const;

            //- Normalisation of matrix; for explicit contributions
            virtual const scalarField& normalisation() const;

            //- Acceptor cells
            virtual const labelList& interpolationCells() const;

            //- Underrelaxation for acceptor cells
            virtual const scalarList& cellInterpolationWeight() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    Foam::cyclicFaPatch

Description
    Cyclic-plane patch.

Author
    Zeljko Tukovic, FMENA
    Hrvoje Jasak, Wikki Ltd.

SourceFiles
    cyclicFaPatch.C

\*---------------------------------------------------------------------------*/

#ifndef cyclicFaPatch_H
#define cyclicFaPatch_H

#include "coupledFaPatch.H"
#include "cyclicLduInterface.H"
#include "cyclicPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class cyclicFaPatch Declaration
\*---------------------------------------------------------------------------*/

class cyclicFaPatch
:
    public coupledFaPatch,
    public cyclicLduInterface
{
    // Private member functions

        void calcTransforms();

protected:

    // Protected static data

        //- Relative tolerance (for geometric matching). Is factor of
        //  maximum edge length per face.
        static const scalar matchTol_;


    // Protected Member functions

        //- Make patch weighting factors
        void makeWeights(scalarField&) const;

        //- Make patch face - neighbour cell distances
        void makeDeltaCoeffs(scalarField&) const;


public:

    //- Runtime type information
    TypeName("cyclic");


    // Constructors

        //- Construct from dictionary
        cyclicFaPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const faBoundaryMesh& bm
        )
        :
            coupledFaPatch(name, dict, index, bm)
        {}


    //- Destructor
    virtual ~cyclicFaPatch()
    {}


    // Member functions

        // Access

            //- Is this the master side?  Yes: it contains both sets of faces
            virtual bool master() const
            {
                return true;
            }

            //- Return neighbour
            virtual label neighbPatchID() const
            {
                NotImplemented;
                return index();
            }

            virtual bool owner() const
            {
                return master();
            }

            //- Return processor number
            virtual const cyclicLduInterface& neighbPatch() const
            {
                NotImplemented;
                return *this;
            }

            //- Return face transformation tensor
            virtual const tensorField& forwardT() const
            {
                return coupledFaPatch::forwardT();
            }

            //- Return neighbour-cell transformation tensor
            virtual const tensorField& reverseT() const
            {
                return coupledFaPatch::reverseT();
            }

            //- Initialise the calculation of the patch geometry
            virtual void initGeometry();

            //- Calculate the patch geometry
            virtual void calcGeometry();

            //- Initialise the patches for moving points
            virtual void initMovePoints(const pointField&);

            //- Correct patches after moving points
            virtual void movePoints(const pointField&);

            //- Return delta (P to N) vectors across coupled patch
            virtual tmp<vectorField> delta() const;


        // Interface transfer functions

            //- Return the values of the given internal data adjacent to
            //  the interface as a field
            virtual tmp<labelField> interfaceInternalField
            (
                const labelUList& internalData
            ) const;

            //- Transfer and return neighbour field
            virtual tmp<labelField> transfer
            (
                const Pstream::commsTypes commsType,
                const labelUList& interfaceData
            ) const;

            //- Return neighbour field
            virtual tmp<labelField> internalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const labelUList& internalData
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

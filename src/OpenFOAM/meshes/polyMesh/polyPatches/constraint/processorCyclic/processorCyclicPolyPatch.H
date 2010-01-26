/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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
    Foam::processorCyclicPolyPatch

Description
    Neighbour processor patch.

    Note: morph patch face ordering is geometric.

SourceFiles
    processorCyclicPolyPatch.C

\*---------------------------------------------------------------------------*/

#ifndef processorCyclicPolyPatch_H
#define processorCyclicPolyPatch_H

#include "processorPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class processorCyclicPolyPatch Declaration
\*---------------------------------------------------------------------------*/

class processorCyclicPolyPatch
:
    public processorPolyPatch
{
    // Private data

        //- Name of originating patch
        const word referPatchName_;

        //- Index of originating patch
        mutable label referPatchID_;


    // Private member functions


protected:

    // Protected Member functions

        //- Initialise the calculation of the patch geometry
        void initGeometry(PstreamBuffers&);

//        //- Initialise the calculation of the patch geometry with externally
//        //  provided geometry
//        virtual void initGeometry
//        (
//            const primitivePatch& referPatch,
//            UList<point>&,
//            UList<point>&,
//            UList<point>&
//        )
//        {
//            notImplemented("processorCyclicPolyPatch::initGeometry(..)");
//        }

        //- Calculate the patch geometry
        void calcGeometry(PstreamBuffers&);

        //- Calculate the patch geometry with externally
        //  provided geometry
        virtual void calcGeometry
        (
            const primitivePatch& referPatch,
            const UList<point>& thisCtrs,
            const UList<point>& thisAreas,
            const UList<point>& thisCc,
            const UList<point>& nbrCtrs,
            const UList<point>& nbrAreas,
            const UList<point>& nbrCc
        )
        {
            notImplemented("processorCyclicPolyPatch::calcGeometry(..)");
        }

        //- Initialise the patches for moving points
        void initMovePoints(PstreamBuffers&, const pointField&);

        //- Correct patches after moving points
        void movePoints(PstreamBuffers&, const pointField&);

        //- Initialise the update of the patch topology
        virtual void initUpdateMesh(PstreamBuffers&);

        //- Update of the patch topology
        virtual void updateMesh(PstreamBuffers&);


public:

    //- Runtime type information
    TypeName("processorCyclic");


    // Constructors

        //- Construct from components
        processorCyclicPolyPatch
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const polyBoundaryMesh& bm,
            const int myProcNo,
            const int neighbProcNo,
            const word& referPatchName
        );

        //- Construct from dictionary
        processorCyclicPolyPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const polyBoundaryMesh&
        );

        //- Construct as copy, resetting the boundary mesh
        processorCyclicPolyPatch
        (
            const processorCyclicPolyPatch&,
            const polyBoundaryMesh&
        );

        //- Construct as given the original patch and resetting the
        //  face list and boundary mesh information
        processorCyclicPolyPatch
        (
            const processorCyclicPolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart,
            const word& referPatchName
        );

        //- Construct given the original patch and a map
        processorCyclicPolyPatch
        (
            const processorCyclicPolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const unallocLabelList& mapAddressing,
            const label newStart
        );


        //- Construct and return a clone, resetting the boundary mesh
        virtual autoPtr<polyPatch> clone(const polyBoundaryMesh& bm) const
        {
            return autoPtr<polyPatch>(new processorCyclicPolyPatch(*this, bm));
        }

        //- Construct and return a clone, resetting the face list
        //  and boundary mesh
        virtual autoPtr<polyPatch> clone
        (
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart,
            const word& referPatchName
        ) const
        {
            return autoPtr<polyPatch>
            (
                new processorCyclicPolyPatch
                (
                    *this,
                    bm,
                    index,
                    newSize,
                    newStart,
                    referPatchName
                )
            );
        }

        //- Construct and return a clone, resetting the face list
        //  and boundary mesh
        virtual autoPtr<polyPatch> clone
        (
            const polyBoundaryMesh& bm,
            const label index,
            const unallocLabelList& mapAddressing,
            const label newStart
        ) const
        {
            return autoPtr<polyPatch>
            (
                new processorCyclicPolyPatch
                (
                    *this,
                    bm,
                    index,
                    mapAddressing,
                    newStart
                )
            );
        }


    // Destructor

        virtual ~processorCyclicPolyPatch();


    // Member functions

        const word& referPatchName() const
        {
            return referPatchName_;
        }

        //- Referring patchID.
        label referPatchID() const
        {
            if (referPatchID_ == -1)
            {
                referPatchID_ = this->boundaryMesh().findPatchID
                (
                    referPatchName_
                );
                if (referPatchID_ == -1)
                {
                    FatalErrorIn
                    (
                        "processorCyclicPolyPatch::referPatchID() const"
                    )   << "Illegal referPatch name " << referPatchName_
                        << endl << "Valid patch names are "
                        << this->boundaryMesh().names()
                        << exit(FatalError);
                }
            }
            return referPatchID_;
        }

        const coupledPolyPatch& referPatch() const
        {
            const polyPatch& pp = this->boundaryMesh()[referPatchID()];
            return refCast<const processorCyclicPolyPatch>(pp);
        }

        //- Transform a patch-based field from other side to this side.
        virtual bool doTransform() const
        {
            return referPatch().doTransform();
        }
        virtual void transform(scalarField& l) const
        {
            referPatch().transform(l);
        }
        virtual void transform(vectorField& l) const
        {
            referPatch().transform(l);
        }
        virtual void transform(sphericalTensorField& l) const
        {
            referPatch().transform(l);
        }
        virtual void transform(diagTensorField& l) const
        {
            referPatch().transform(l);
        }
        virtual void transform(symmTensorField& l) const
        {
            referPatch().transform(l);
        }
        virtual void transform(tensorField& l) const
        {
            referPatch().transform(l);
        }

        //- Transform a patch-based position from other side to this side
        virtual void transformPosition(pointField& l) const
        {
            referPatch().transform(l);
        }

        //- Are the planes separated.
        virtual bool separated() const
        {
            return referPatch().separated();
        }

        //- If the planes are separated the separation vector.
        virtual const vector& separation() const
        {
            return referPatch().separation();
        }

        //- Are the cyclic planes parallel.
        virtual bool parallel() const
        {
            return referPatch().parallel();
        }

        //- Return face transformation tensor.
        virtual const tensor& forwardT() const
        {
            return referPatch().forwardT();
        }

        //- Return neighbour-cell transformation tensor.
        virtual const tensor& reverseT() const
        {
            return referPatch().reverseT();
        }

        //- Initialize ordering for primitivePatch. Does not
        //  refer to *this (except for name() and type() etc.)
        virtual void initOrder(PstreamBuffers&, const primitivePatch&) const;

        //- Return new ordering for primitivePatch.
        //  Ordering is -faceMap: for every face
        //  index of the new face -rotation:for every new face the clockwise
        //  shift of the original face. Return false if nothing changes
        //  (faceMap is identity, rotation is 0), true otherwise.
        virtual bool order
        (
            PstreamBuffers&,
            const primitivePatch&,
            labelList& faceMap,
            labelList& rotation
        ) const;


        //- Write the polyPatch data as a dictionary
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

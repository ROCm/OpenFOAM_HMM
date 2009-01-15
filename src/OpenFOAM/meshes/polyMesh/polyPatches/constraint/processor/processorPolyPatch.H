/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
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
    Foam::processorPolyPatch

Description
    Neighbour processor patch.

    Note: morph patch face ordering comes geometric or topological.
    Geometric: no cyclics allowed (assumes faces coincident)
    Topological: needs unmodified faces on both sides to correspond. Also
    needs at least one per connected patch area (so all patch faces can be
    visited from an unmodified face)

SourceFiles
    processorPolyPatch.C
    processorPolyPatchMorph.C

\*---------------------------------------------------------------------------*/

#ifndef processorPolyPatch_H
#define processorPolyPatch_H

#include "coupledPolyPatch.H"
#include "polyBoundaryMesh.H"
#include "faceListFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class processorPolyPatch Declaration
\*---------------------------------------------------------------------------*/

class processorPolyPatch
:
    public coupledPolyPatch
{
    // Private data

        int myProcNo_;
        int neighbProcNo_;

        // Per patch information

            //- Originating patches
            labelList patchIDs_;

            //- Slice starts
            labelList starts_;

        //- Processor-neighbbour patch face centres
        vectorField neighbFaceCentres_;

        //- Processor-neighbbour patch face areas
        vectorField neighbFaceAreas_;

        //- Processor-neighbbour patch neighbour cell centres
        vectorField neighbFaceCellCentres_;

        //- Corresponding neighbouring local point label for every local point
        //  (so localPoints()[i] == neighb.localPoints()[neighbPoints_[i]])
        mutable labelList* neighbPointsPtr_;

        //- Corresponding neighbouring local edge label for every local edge
        //  (so edges()[i] == neighb.edges()[neighbEdges_[i]])
        mutable labelList* neighbEdgesPtr_;


    // Private member functions

        void checkSubPatches() const;

protected:

    // Protected Member functions

        //- Initialise the calculation of the patch geometry
        void initGeometry();

        //- Initialise the calculation of the patch geometry with externally
        //  provided geometry
        virtual void initGeometry
        (
            const primitivePatch& referPatch,
            UList<point>&,
            UList<point>&,
            UList<point>&
        )
        {
            notImplemented("processorPolyPatch::initGeometry(..)");
        }

        //- Calculate the patch geometry
        void calcGeometry();

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
            notImplemented("processorPolyPatch::calcGeometry(..)");
        }

        //- Initialise the patches for moving points
        void initMovePoints(const pointField&);

        //- Correct patches after moving points
        void movePoints(const pointField&);

        //- Initialise the update of the patch topology
        virtual void initUpdateMesh();

        //- Update of the patch topology
        virtual void updateMesh();


public:

    //- Runtime type information
    TypeName("processor");


    // Constructors

        //- Construct from components
        processorPolyPatch
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const polyBoundaryMesh& bm,
            const int myProcNo,
            const int neighbProcNo,
            const labelList& patchIDs,
            const labelList& starts
        );

        //- Construct from dictionary
        processorPolyPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const polyBoundaryMesh&
        );

        //- Construct as copy, resetting the boundary mesh
        processorPolyPatch
        (
            const processorPolyPatch&,
            const polyBoundaryMesh&
        );

        //- Construct as given the original patch and resetting the
        //  face list and boundary mesh information
        processorPolyPatch
        (
            const processorPolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart,
            const labelList& patchIDs,
            const labelList& starts
        );

        //- Construct and return a clone, resetting the boundary mesh
        virtual autoPtr<polyPatch> clone(const polyBoundaryMesh& bm) const
        {
            return autoPtr<polyPatch>(new processorPolyPatch(*this, bm));
        }

        //- Construct and return a clone, resetting the face list
        //  and boundary mesh
        virtual autoPtr<polyPatch> clone
        (
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart,
            const labelList& patchIDs,
            const labelList& starts
        ) const
        {
            return autoPtr<polyPatch>
            (
                new processorPolyPatch
                (
                    refCast<const processorPolyPatch>(*this),
                    bm,
                    index,
                    newSize,
                    newStart,
                    patchIDs,
                    starts
                )
            );
        }


    // Destructor

        virtual ~processorPolyPatch();


    // Member functions

        //- Return processor number
        int myProcNo() const
        {
            return myProcNo_;
        }

        //- Return neigbour processor number
        int neighbProcNo() const
        {
            return neighbProcNo_;
        }

        //- Does the processor own the patch ?
        bool owner() const
        {
            return (myProcNo_ < neighbProcNo_);
        }

        //- Is the processor the patch neighbour ?
        bool neighbour() const
        {
            return !owner();
        }

        const labelList& patchIDs() const
        {
            return patchIDs_;
        }

        const labelList& starts() const
        {
            return starts_;
        }

        //- Return processor-neighbbour patch face centres
        const vectorField& neighbFaceCentres() const
        {
            return neighbFaceCentres_;
        }

        //- Return processor-neighbbour patch face areas
        const vectorField& neighbFaceAreas() const
        {
            return neighbFaceAreas_;
        }

        //- Return processor-neighbbour patch neighbour cell centres
        const vectorField& neighbFaceCellCentres() const
        {
            return neighbFaceCellCentres_;
        }

        //- Are the planes separated.
        virtual bool separated() const
        {
            notImplemented("processorPolyPatch::separated(..)");
            return false;
        }

        //- If the planes are separated the separation vector.
        virtual const vector& separation() const
        {
            notImplemented("processorPolyPatch::separation(..)");
            return vector::zero;
        }

        //- Are the cyclic planes parallel.
        virtual bool parallel() const
        {
            forAll(patchIDs_, i)
            {
                label patchI = patchIDs_[i];

                if (patchI != -1)
                {
                    if
                    (
                       !refCast<const coupledPolyPatch>
                        (
                            boundaryMesh()[patchI]
                        ).parallel()
                    )
                    {
                        return false;
                    }
                }
            }
            return true;
        }

        //- Return face transformation tensor.
        virtual const tensor& forwardT() const
        {
            notImplemented("processorPolyPatch::forwardT(..)");
            return tensor::zero;
        }

        //- Return neighbour-cell transformation tensor.
        virtual const tensor& reverseT() const
        {
            notImplemented("processorPolyPatch::reverseT(..)");
            return tensor::zero;
        }

        //- Find sub patch for local face. -1 for internal faces or label of
        //  cyclic patch.
        label whichSubPatch(const label facei) const;

        //- Slice patch
        autoPtr<primitivePatch> subPatch(const label subI) const
        {
            const faceList& thisFaces = static_cast<const faceList&>(*this);

            return autoPtr<primitivePatch>
            (
                new primitivePatch
                (
                    faceSubList
                    (
                        thisFaces,
                        this->starts_[subI+1]-this->starts_[subI],
                        this->starts_[subI]
                    ),
                    this->points()
                )
            );
        }

        //- Slice patchlist to subpatch
        template<class T>
        const typename List<T>::subList subSlice
        (
            const UList<T>& l,
            const label subI
        ) const
        {
            return typename List<T>::subList
            (
                l,
                this->starts_[subI+1]-this->starts_[subI],
                this->starts_[subI]
            );
        }

        //- Slice patchlist to subpatch
        template<class T>
        typename List<T>::subList subSlice
        (
            UList<T>& l,
            const label subI
        )
        {
            return typename List<T>::subList
            (
                l,
                this->starts_[subI+1]-this->starts_[subI],
                this->starts_[subI]
            );
        }

        //- Slice patchField to subpatch
        template<class T>
        const typename Field<T>::subField subSlice
        (
            const Field<T>& l,
            const label subI
        ) const
        {
            return typename Field<T>::subList
            (
                l,
                this->starts_[subI+1]-this->starts_[subI],
                this->starts_[subI]
            );
        }

        //- Slice patchField to subpatch
        template<class T>
        typename Field<T>::subField subSlice
        (
            Field<T>& l,
            const label subI
        )
        {
            return typename Field<T>::subList
            (
                l,
                this->starts_[subI+1]-this->starts_[subI],
                this->starts_[subI]
            );
        }

        //- Return neighbour point labels. This is for my local point (-1 or)
        //  the corresponding local point on the other side. It is -1 if
        //  there are multiple corresponding points on this or the other side
        //  (can happen for cyclics being converted into proc patches)
        const labelList& neighbPoints() const;

        //- Return neighbour edge labels. This is for my local edge (-1 or) the
        //  corresponding local edge on the other side. See above for -1 cause.
        const labelList& neighbEdges() const;

        ////- Transform a position
        //virtual void transformPosition(point&, const label facei) const;

        //- Initialize ordering for primitivePatch. Does not
        //  refer to *this (except for name() and type() etc.)
        virtual void initOrder(const primitivePatch&) const;

        //- Return new ordering for primitivePatch.
        //  Ordering is -faceMap: for every face
        //  index of the new face -rotation:for every new face the clockwise
        //  shift of the original face. Return false if nothing changes
        //  (faceMap is identity, rotation is 0), true otherwise.
        virtual bool order
        (
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

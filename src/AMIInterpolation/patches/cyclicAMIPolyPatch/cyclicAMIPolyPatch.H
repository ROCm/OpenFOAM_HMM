/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2011 OpenCFD Ltd.
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
    Foam::cyclicAMIPolyPatch

Description
    Cyclic patch for Arbitrary Mesh Interface (AMI)

SourceFiles
    cyclicAMIPolyPatch.C

\*---------------------------------------------------------------------------*/

#ifndef cyclicAMIPolyPatch_H
#define cyclicAMIPolyPatch_H

#include "coupledPolyPatch.H"
#include "AMIPatchToPatchInterpolation.H"
#include "polyBoundaryMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class cyclicAMIPolyPatch Declaration
\*---------------------------------------------------------------------------*/

class cyclicAMIPolyPatch
:
    public coupledPolyPatch
{

private:

    // Private data

        //- Name of other half
        const word nbrPatchName_;

        //- Index of other half
        mutable label nbrPatchID_;

        //- Neighbour delta field to *this patch
        vectorField nbrPatchDelta_;


        // Transformations

            //- Type of transformation - rotational or translational
            transformType transform_;


            // For rotation

                //- Axis of rotation for rotational cyclics
                vector rotationAxis_;

                //- point on axis of rotation for rotational cyclics
                point rotationCentre_;


            // For translation

                //- Translation vector
                vector separationVector_;


        //- AMI interpolation class
        AMIPatchToPatchInterpolation* AMIPtr_;

        //- Projection surface
        autoPtr<searchableSurface> surfPtr_;


    // Private Member Functions

        //- Return face ID of face at max distance from rotation axis
        label findFaceMaxRadius(const pointField& faceCentres) const;

        void calcTransforms
        (
            const primitivePatch& half0,
            const pointField& half0Ctrs,
            const vectorField& half0Areas,
            const pointField& half1Ctrs,
            const vectorField& half1Areas
        );


protected:

    // Protected Member Functions

        //- Recalculate the transformation tensors
        virtual void calcTransforms();

        //- Initialise the calculation of the patch geometry
        virtual void initGeometry(PstreamBuffers&);

        //- Calculate the patch geometry
        virtual void calcGeometry(PstreamBuffers&);

        //- Initialise the patches for moving points
        virtual void initMovePoints(PstreamBuffers& pBufs, const pointField&);

        //- Correct patches after moving points
        virtual void movePoints(PstreamBuffers& pBufs, const pointField&);

        //- Initialise the update of the patch topology
        virtual void initUpdateMesh(PstreamBuffers&);

        //- Update of the patch topology
        virtual void updateMesh(PstreamBuffers&);


public:

//    //- Declare friendship with processorCyclicPolyPatch
//    friend class processorCyclicPolyPatch;

    //- Runtime type information
    TypeName("cyclicAMI");


    // Constructors

        //- Construct from (base couped patch) components
        cyclicAMIPolyPatch
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const polyBoundaryMesh& bm
        );

        //- Construct from dictionary
        cyclicAMIPolyPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const polyBoundaryMesh& bm
        );

        //- Construct as copy, resetting the boundary mesh
        cyclicAMIPolyPatch(const cyclicAMIPolyPatch&, const polyBoundaryMesh&);

        //- Construct given the original patch and resetting the
        //  face list and boundary mesh information
        cyclicAMIPolyPatch
        (
            const cyclicAMIPolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart,
            const word& nbrPatchName
        );

        //- Construct given the original patch and a map
        cyclicAMIPolyPatch
        (
            const cyclicAMIPolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const labelUList& mapAddressing,
            const label newStart
        );


        //- Construct and return a clone, resetting the boundary mesh
        virtual autoPtr<polyPatch> clone(const polyBoundaryMesh& bm) const
        {
            return autoPtr<polyPatch>(new cyclicAMIPolyPatch(*this, bm));
        }

        //- Construct and return a clone, resetting the face list
        //  and boundary mesh
        virtual autoPtr<polyPatch> clone
        (
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart
        ) const
        {
            return autoPtr<polyPatch>
            (
                new cyclicAMIPolyPatch
                (
                    *this,
                    bm,
                    index,
                    newSize,
                    newStart,
                    nbrPatchName_
                )
            );
        }

        //- Construct and return a clone, resetting the face list
        //  and boundary mesh
        virtual autoPtr<polyPatch> clone
        (
            const polyBoundaryMesh& bm,
            const label index,
            const labelUList& mapAddressing,
            const label newStart
        ) const
        {
            return autoPtr<polyPatch>
            (
                new cyclicAMIPolyPatch
                (
                    *this,
                    bm,
                    index,
                    mapAddressing,
                    newStart
                )
            );
        }


    //- Destructor
    virtual ~cyclicAMIPolyPatch();


    // Member Functions

        // Access

            //- Neighbour patch name
            const word& nbrPatchName() const
            {
                return nbrPatchName_;
            }

            //- Neighbour patch ID
            virtual label nbrPatchID() const;

            //- Does this side own the patch ?
            virtual bool owner() const
            {
                return index() < nbrPatchID();
            }

            const cyclicAMIPolyPatch& nbrPatch() const
            {
                const polyPatch& pp = this->boundaryMesh()[nbrPatchID()];
                return refCast<const cyclicAMIPolyPatch>(pp);
            }


            const vectorField& nbrPatchDelta() const
            {
                return nbrPatchDelta_;
            }

            const AMIPatchToPatchInterpolation& AMI();


            // Transformations

                //- Type of transform
                transformType transform() const
                {
                    return transform_;
                }

                //- Axis of rotation for rotational cyclic AMI
                const vector& rotationAxis() const
                {
                    return rotationAxis_;
                }

                //- point on axis of rotation for rotational cyclic AMI
                const point& rotationCentre() const
                {
                    return rotationCentre_;
                }

                //- Translation vector for translational cyclic AMI
                const vector& separationVector() const
                {
                    return separationVector_;
                }

                //- Transform patch-based positions from nbr side to this side
                virtual void transformPosition(pointField&) const;

                //- Transform a patch-based position from nbr side to this side
                virtual void transformPosition
                (
                    point& l,
                    const label faceI
                ) const;

                template<class Type>
                tmp<Field<Type> > interpolateToSource(const Field<Type>& fld) const
                {
                    return AMIPtr_->interpolateToSource(fld);
                }

                template<class Type>
                tmp<Field<Type> > interpolateToSource(const tmp<Field<Type> >& tFld) const
                {
                    return AMIPtr_->interpolateToSource(tFld);
                }


        //- Calculate the patch geometry
        virtual void calcGeometry
        (
            const primitivePatch& referPatch,
            const pointField& thisCtrs,
            const vectorField& thisAreas,
            const pointField& thisCc,
            const pointField& nbrCtrs,
            const vectorField& nbrAreas,
            const pointField& nbrCc
        );

        //- Initialize ordering for primitivePatch. Does not
        //  refer to *this (except for name() and type() etc.)
        virtual void initOrder
        (
            PstreamBuffers&,
            const primitivePatch&
        ) const;

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

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
    Foam::cyclicPolyPatch

Description
    Cyclic plane patch.

    Note: morph patch face ordering uses geometric matching so with the
    following restrictions:
        -coupled patches should be flat planes.
        -no rotation in patch plane

    Uses coupledPolyPatch::calcFaceTol to calculate
    tolerance per face which might need tweaking.

    Switch on 'cyclicPolyPatch' debug flag to write .obj files to show
    the matching.

SourceFiles
    cyclicPolyPatch.C

\*---------------------------------------------------------------------------*/

#ifndef cyclicPolyPatch_H
#define cyclicPolyPatch_H

#include "coupledPolyPatch.H"
//#include "SubField.H"
//#include "FixedList.H"
#include "edgeList.H"
//#include "transform.H"
#include "polyBoundaryMesh.H"
#include "diagTensorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class cyclicPolyPatch Declaration
\*---------------------------------------------------------------------------*/

class cyclicPolyPatch
:
    public coupledPolyPatch
{
public:

    enum transformType
    {
        UNKNOWN,
        ROTATIONAL,
        TRANSLATIONAL
    };
    static const NamedEnum<transformType, 3> transformTypeNames;


private:

    // Private data

        //- Name of other half
        const word neighbPatchName_;

        //- Index of other half
        mutable label neighbPatchID_;

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


        //- List of edges formed from connected points. e[0] is the point on
        //  the first half of the patch, e[1] the corresponding point on the
        //  second half.
        mutable edgeList* coupledPointsPtr_;

        //- List of connected edges. e[0] is the edge on the first half of the
        //  patch, e[1] the corresponding edge on the second half.
        mutable edgeList* coupledEdgesPtr_;

        //- Temporary storage of owner side patch during ordering.
        mutable autoPtr<primitivePatch> ownerPatchPtr_;

        // Transformation

            //- Has non-zero separation
            mutable bool separated_;

            //- Offset (distance) vector from one side of the couple to the other
            mutable vector separation_;

            //- No rotation
            mutable bool parallel_;

            //- Face transformation tensor
            mutable tensor forwardT_;

            //- Neighbour-cell transformation tensor
            mutable tensor reverseT_;


    // Private member functions

        //- Find amongst selected faces the one with the largest area
        static label findMaxArea(const pointField&, const faceList&);

        void calcTransforms();

        void calcTransforms
        (
            const primitivePatch& half0,
            const UList<point>& half0Ctrs,
            const UList<point>& half0Areas,
            const UList<point>& half1Ctrs,
            const UList<point>& half1Areas
        );


        // Face ordering

            //- Calculate geometric factors of the two halves.
            void getCentresAndAnchors
            (
                const primitivePatch& pp0,
                const primitivePatch& pp1,

                pointField& half0Ctrs,
                pointField& half1Ctrs,
                pointField& anchors0,
                scalarField& tols
            ) const;

            //- For rotational cases, try to find a unique face on each side
            //  of the cyclic.
            label getConsistentRotationFace(const pointField&) const;


protected:

    // Protected Member functions

        //- Initialise the calculation of the patch geometry
        virtual void initGeometry();

        //- Initialise the calculation of the patch geometry
        virtual void initGeometry
        (
            const primitivePatch& referPatch,
            UList<point>& nbrCtrs,
            UList<point>& nbrAreas,
            UList<point>& nbrCc
        );

        //- Calculate the patch geometry
        virtual void calcGeometry();

        //- Calculate the patch geometry
        virtual void calcGeometry
        (
            const primitivePatch& referPatch,
            const UList<point>& thisCtrs,
            const UList<point>& thisAreas,
            const UList<point>& thisCc,
            const UList<point>& nbrCtrs,
            const UList<point>& nbrAreas,
            const UList<point>& nbrCc
        );

        //- Initialise the patches for moving points
        virtual void initMovePoints(const pointField&);

        //- Correct patches after moving points
        virtual void movePoints(const pointField&);

        //- Initialise the update of the patch topology
        virtual void initUpdateMesh();

        //- Update of the patch topology
        virtual void updateMesh();

public:

    //- Runtime type information
    TypeName("cyclic");


    // Constructors

        //- Construct from components
        cyclicPolyPatch
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const polyBoundaryMesh& bm
        );

        //- Construct from dictionary
        cyclicPolyPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const polyBoundaryMesh& bm
        );

        //- Construct as copy, resetting the boundary mesh
        cyclicPolyPatch(const cyclicPolyPatch&, const polyBoundaryMesh&);

        //- Construct given the original patch and resetting the
        //  face list and boundary mesh information
        cyclicPolyPatch
        (
            const cyclicPolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart,
            const word& neighbPatchName
        );

        //- Construct and return a clone, resetting the boundary mesh
        virtual autoPtr<polyPatch> clone(const polyBoundaryMesh& bm) const
        {
            return autoPtr<polyPatch>(new cyclicPolyPatch(*this, bm));
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
                new cyclicPolyPatch
                (
                    *this,
                    bm,
                    index,
                    newSize,
                    newStart,
                    neighbPatchName_
                )
            );
        }


    // Destructor

        virtual ~cyclicPolyPatch();


    // Member Functions

        const word& neighbPatchName() const
        {
            return neighbPatchName_;
        }

        //- Neighbour patchID.
        label neighbPatchID() const
        {
            if (neighbPatchID_ == -1)
            {
                neighbPatchID_ = this->boundaryMesh().findPatchID
                (
                    neighbPatchName_
                );
                if (neighbPatchID_ == -1)
                {
                    FatalErrorIn("cyclicPolyPatch::neighbPatchID() const")
                        << "Illegal neighbourPatch name " << neighbPatchName_
                        << endl << "Valid patch names are "
                        << this->boundaryMesh().names()
                        << exit(FatalError);
                }
            }
            return neighbPatchID_;
        }

        bool owner() const
        {
            return index() < neighbPatchID();
        }

        bool neighbour() const
        {
            return !owner();
        }

        const cyclicPolyPatch& neighbPatch() const
        {
            const polyPatch& pp = this->boundaryMesh()[neighbPatchID()];
            return refCast<const cyclicPolyPatch>(pp);
        }

        //- Return connected points (from patch local to neighbour patch local)
        //  Demand driven calculation. Does primitivePatch::clearOut after
        //  calculation!
        const edgeList& coupledPoints() const;

        //- Return connected edges (from patch local to neighbour patch local).
        //  Demand driven calculation. Does primitivePatch::clearOut after
        //  calculation!
        const edgeList& coupledEdges() const;

        //- Transform a patch-based field on this side to a field on the
        //  other side.
        template<class T>
        void transform(Field<T>& l) const
        {
            if (!parallel())
            {
                transform(reverseT_, l);
            }
        }
        virtual void transform(scalarField& l) const
        {
            transform(l);
        }
        virtual void transform(vectorField& l) const
        {
            transform(l);
        }
        virtual void transform(sphericalTensorField& l) const
        {
            transform(l);
        }
        virtual void transform(diagTensorField& l) const
        {
            transform(l);
        }
        virtual void transform(symmTensorField& l) const
        {
            transform(l);
        }
        virtual void transform(tensorField& l) const
        {
            transform(l);
        }

        //- Transform a patch-based position on this side to a position
        //  the other side.
        virtual void transformPosition(pointField& l) const;


        // Low level geometric information

            //- Are the coupled planes separated
            virtual bool separated() const
            {
                return separated_;
            }

            virtual const vector& separation() const
            {
                return separation_;
            }

            //- Are the cyclic planes parallel
            virtual bool parallel() const
            {
                return parallel_;
            }

            //- Return face transformation tensor
            virtual const tensor& forwardT() const
            {
                return forwardT_;
            }

            //- Return neighbour-cell transformation tensor
            virtual const tensor& reverseT() const
            {
                return reverseT_;
            }

        label transformGlobalFace(const label facei) const
        {
            label offset = facei-start();
            label neighbStart = neighbPatch().start();

            if (offset >= 0 && offset < size())
            {
                return neighbStart+offset;
            }
            else
            {
                FatalErrorIn
                (
                    "cyclicPolyPatch::transformGlobalFace(const label) const"
                )   << "Face " << facei << " not in patch " << name()
                    << exit(FatalError);
                return -1;
            }
        }

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

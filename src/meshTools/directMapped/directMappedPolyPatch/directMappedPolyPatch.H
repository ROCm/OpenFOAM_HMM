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
    Foam::directMappedPolyPatch

Description
    Determines a mapping between patch face centres and mesh cell or face
    centres and processors they're on.

Note
    Storage is not optimal. It stores all face centres and cells on all
    processors to keep the addressing calculation simple.

SourceFiles
    directMappedPolyPatch.C

\*---------------------------------------------------------------------------*/

#ifndef directMappedPolyPatch_H
#define directMappedPolyPatch_H

#include "polyPatch.H"
#include "labelPair.H"
#include "Tuple2.H"
#include "pointIndexHit.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class directMappedPolyPatch Declaration
\*---------------------------------------------------------------------------*/

class directMappedPolyPatch
:
    public polyPatch
{

public:

        //- Mesh items to sample
        enum sampleMode
        {
            NEARESTCELL,
            NEARESTPATCHFACE,
            NEARESTFACE
        };

private:

    // Private data

        static const NamedEnum<sampleMode, 3> sampleModeNames_;

        //- What to sample
        sampleMode mode_;

        //- Patch (only if NEARESTBOUNDARY)
        const word samplePatch_;

        //- Offset vector
        const vector offset_;



        // Derived information

            //- Communication schedule.
            mutable autoPtr<List<labelPair> > schedulePtr_;

            //- Cells/faces to sample per processor
            mutable autoPtr<labelListList> sendLabelsPtr_;

            //- Patch faces to receive per processor
            mutable autoPtr<labelListList> receiveFaceLabelsPtr_;


    // Private Member Functions

        //- Collect single list of samples and originating processor+face.
        void collectSamples
        (
            pointField&,
            labelList& patchFaceProcs,
            labelList& patchFaces,
            pointField& patchFc
        ) const;

        //- Find cells/faces containing samples
        void findSamples
        (
            const pointField&,
            labelList& sampleProcs,     // processor containing sample
            labelList& sampleIndices,   // local index of cell/face
            pointField& sampleLocations // actual representative location
        ) const;

        //- Calculate matching
        void calcMapping() const;


    // Private class

        //- Private class for finding nearest
        //  - point+local index
        //  - sqr(distance)
        //  - processor
        typedef Tuple2<pointIndexHit, Tuple2<scalar, label> > nearInfo;

        class nearestEqOp
        {

        public:

            void operator()(nearInfo& x, const nearInfo& y) const
            {
                if (y.first().hit())
                {
                    if (!x.first().hit())
                    {
                        x = y;
                    }
                    else if (y.second().first() < x.second().first())
                    {
                        x = y;
                    }
                }
            }
        };



protected:

        void clearOut();

        //- Initialise the calculation of the patch geometry
        virtual void initGeometry()
        {
            clearOut();
        }

        //- Calculate the patch geometry
        virtual void calcGeometry()
        {
            clearOut();
        }

        //- Initialise the patches for moving points
        virtual void initMovePoints(const pointField&)
        {
            clearOut();
        }

        //- Correct patches after moving points
        virtual void movePoints(const pointField& p)
        {
            clearOut();
        }

        //- Initialise the update of the patch topology
        virtual void initUpdateMesh()
        {
            clearOut();
        }

        //- Update of the patch topology
        virtual void updateMesh()
        {
            clearOut();
        }


public:

    //- Runtime type information
    TypeName("directMappedPatch");


    // Constructors

        //- Construct from components
        directMappedPolyPatch
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const polyBoundaryMesh& bm
        );

        //- Construct from dictionary
        directMappedPolyPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const polyBoundaryMesh& bm
        );

        //- Construct as copy, resetting the boundary mesh
        directMappedPolyPatch
        (
            const directMappedPolyPatch&,
            const polyBoundaryMesh&
        );

        //- Construct given the original patch and resetting the
        //  face list and boundary mesh information
        directMappedPolyPatch
        (
            const directMappedPolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart
        );

        //- Construct and return a clone, resetting the boundary mesh
        virtual autoPtr<polyPatch> clone(const polyBoundaryMesh& bm) const
        {
            return autoPtr<polyPatch>(new directMappedPolyPatch(*this, bm));
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
                new directMappedPolyPatch(*this, bm, index, newSize, newStart)
            );
        }


    // Destructor

        ~directMappedPolyPatch();


    // Member functions

        //- What to sample
        const sampleMode& mode() const
        {
            return mode_;
        }

        //- Patch (only if NEARESTBOUNDARY)
        const word& samplePatch() const
        {
            return samplePatch_;
        }

        //- Offset vector (from patch faces to destination mesh objects)
        const vector& offset() const
        {
            return offset_;
        }

        //- Access to communication.
        const List<labelPair>& schedule() const
        {
            if (!schedulePtr_.valid())
            {
                calcMapping();
            }
            return schedulePtr_();
        }

        //- Cells/faces to sample per processor
        const labelListList& sendLabels() const
        {
            if (debug)
            {
                Pout<< "Asking for sendLabels." << endl;
            }

            if (!sendLabelsPtr_.valid())
            {
                if (debug)
                {
                    Pout<< "Calculating mapping." << endl;
                    calcMapping();
                }
            }
            return sendLabelsPtr_();
        }

        //- Patch faces to receive per processor
        const labelListList& receiveFaceLabels() const
        {
            if (debug)
            {
                Pout<< "Asking for receiveFaceLabels." << endl;
            }

            if (!receiveFaceLabelsPtr_.valid())
            {
                if (debug)
                {
                    Pout<< "Calculating mapping." << endl;
                    calcMapping();
                }
            }
            return receiveFaceLabelsPtr_();
        }

        //- Write the polyPatch data as a dictionary
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

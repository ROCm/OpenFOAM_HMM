/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::directMappedPatchBase

Description
    Determines a mapping between patch face centres and mesh cell or face
    centres and processors they're on.

    If constructed from dictionary:
        // Region to sample (default is region0)
        sampleRegion region0;

        // What to sample:
        // - nearestCell : sample nearest cell
        // - nearestPatchFace : nearest face on selected patch
        // - nearestFace : nearest boundary face on any patch
        sampleMode nearestCell;

        // If sampleMod is nearestPatchFace : patch to find faces of
        samplePatch movingWall;

        // How to supply offset (w.r.t. my patch face centres):
        // - uniform : single offset vector
        // - nonuniform : per-face offset vector
        // - normal : using supplied distance and face normal
        offsetMode uniform;

        // According to offsetMode (see above) supply one of
        // offset, offsets or distance
        offset  (1 0 0);

    Note: if offsetMode is 'normal' it uses outwards pointing normals. So
    supply a negative distance if sampling inside the domain.


Note
    Storage is not optimal. It temporary collects all (patch)face centres
    on all processors to keep the addressing calculation simple.

SourceFiles
    directMappedPatchBase.C

\*---------------------------------------------------------------------------*/

#ifndef directMappedPatchBase_H
#define directMappedPatchBase_H

#include "pointField.H"
#include "Tuple2.H"
#include "pointIndexHit.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyPatch;
class polyMesh;
class mapDistribute;

/*---------------------------------------------------------------------------*\
                      Class directMappedPatchBase Declaration
\*---------------------------------------------------------------------------*/

class directMappedPatchBase
{

public:

        //- Mesh items to sample
        enum sampleMode
        {
            NEARESTCELL,        // nearest cell
            NEARESTPATCHFACE,   // faces on selected patch
            NEARESTFACE         // nearest face
        };

        //- How to project face centres
        enum offsetMode
        {
            UNIFORM,            // single offset vector
            NONUNIFORM,         // per-face offset vector
            NORMAL              // use face normal + distance
        };


    //- Helper class for finding nearest
    // Nearest:
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

private:

    // Private data

        static const NamedEnum<sampleMode, 3> sampleModeNames_;

        static const NamedEnum<offsetMode, 3> offsetModeNames_;

        //- Patch to sample
        const polyPatch& patch_;

        //- Region to sample
        const word sampleRegion_;

        //- What to sample
        const sampleMode mode_;

        //- Patch (only if NEARESTPATCHFACE)
        const word samplePatch_;

        //- How to obtain samples
        offsetMode offsetMode_;

        //- Offset vector (uniform)
        vector offset_;

        //- Offset vector (nonuniform)
        vectorField offsets_;

        //- Offset distance (normal)
        scalar distance_;

        //- Same region
        const bool sameRegion_;


        // Derived information

            //- Communication schedule:
            //  - Cells/faces to sample per processor
            //  - Patch faces to receive per processor
            //  - schedule
            mutable autoPtr<mapDistribute> mapPtr_;


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


public:

    //- Runtime type information
    TypeName("directMappedPatchBase");


    // Constructors

        //- Construct from patch
        directMappedPatchBase(const polyPatch&);

        //- Construct with offsetMode=non-uniform
        directMappedPatchBase
        (
            const polyPatch& pp,
            const word& sampleRegion,
            const sampleMode sampleMode,
            const word& samplePatch,
            const vectorField& offsets
        );

        //- Construct from offsetMode=uniform
        directMappedPatchBase
        (
            const polyPatch& pp,
            const word& sampleRegion,
            const sampleMode sampleMode,
            const word& samplePatch,
            const vector& offset
        );

        //- Construct from offsetMode=normal and distance
        directMappedPatchBase
        (
            const polyPatch& pp,
            const word& sampleRegion,
            const sampleMode sampleMode,
            const word& samplePatch,
            const scalar distance
        );

        //- Construct from dictionary
        directMappedPatchBase(const polyPatch&, const dictionary&);

        //- Construct as copy, resetting patch
        directMappedPatchBase(const polyPatch&, const directMappedPatchBase&);

        //- Construct as copy, resetting patch, map original data
        directMappedPatchBase
        (
            const polyPatch&,
            const directMappedPatchBase&,
            const labelUList& mapAddressing
        );


    //- Destructor
    virtual ~directMappedPatchBase();


    // Member functions

        void clearOut();

        //- What to sample
        const sampleMode& mode() const
        {
            return mode_;
        }

        //- Region to sample
        const word& sampleRegion() const
        {
            return sampleRegion_;
        }

        //- Patch (only if NEARESTPATCHFACE)
        const word& samplePatch() const
        {
            return samplePatch_;
        }

        //- Offset vector (from patch faces to destination mesh objects)
        const vector& offset() const
        {
            return offset_;
        }

        //- Offset vector (from patch faces to destination mesh objects)
        const vectorField& offsets() const
        {
            return offsets_;
        }

        //- Return reference to the parallel distribution map
        const mapDistribute& map() const
        {
            if (mapPtr_.empty())
            {
                calcMapping();
            }
            return mapPtr_();
        }

        //- Cached sampleRegion != mesh.name()
        bool sameRegion() const
        {
            return sameRegion_;
        }

        //- Get the region mesh
        const polyMesh& sampleMesh() const;

        //- Get the patch on the region
        const polyPatch& samplePolyPatch() const;

        //- Get the sample points
        tmp<pointField> samplePoints() const;

        //- Write as a dictionary
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

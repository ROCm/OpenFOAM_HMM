/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::wallBoundedStreamLineParticle

Description
    Particle class that samples fields as it passes through. Used in streamline
    calculation.

SourceFiles
    wallBoundedStreamLineParticle.C

\*---------------------------------------------------------------------------*/

#ifndef wallBoundedStreamLineParticle_H
#define wallBoundedStreamLineParticle_H

#include "particle.H"
#include "autoPtr.H"
#include "interpolation.H"
#include "vectorList.H"
#include "InfoProxy.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class wallBoundedStreamLineParticleCloud;

/*---------------------------------------------------------------------------*\
                     Class wallBoundedStreamLineParticle Declaration
\*---------------------------------------------------------------------------*/

class wallBoundedStreamLineParticle
:
    public particle
{

public:

    //- Class used to pass tracking data to the trackToFace function
    class trackingData
    :
        public particle::TrackingData<Cloud<wallBoundedStreamLineParticle> >
    {

    public:


        const PtrList<interpolation<scalar> >& vsInterp_;
        const PtrList<interpolation<vector> >& vvInterp_;
        const label UIndex_;
        const bool trackForward_;
        const scalar trackLength_;

        const PackedBoolList& isWallPatch_;

        DynamicList<vectorList>& allPositions_;
        List<DynamicList<scalarList> >& allScalars_;
        List<DynamicList<vectorList> >& allVectors_;


        // Constructors

            trackingData
            (
                Cloud<wallBoundedStreamLineParticle>& cloud,
                const PtrList<interpolation<scalar> >& vsInterp,
                const PtrList<interpolation<vector> >& vvInterp,
                const label UIndex,
                const bool trackForward,
                const scalar trackLength,
                const PackedBoolList& isWallPatch,

                DynamicList<List<point> >& allPositions,
                List<DynamicList<scalarList> >& allScalars,
                List<DynamicList<vectorList> >& allVectors
            )
            :
                particle::TrackingData<Cloud<wallBoundedStreamLineParticle> >
                (
                    cloud
                ),
                vsInterp_(vsInterp),
                vvInterp_(vvInterp),
                UIndex_(UIndex),
                trackForward_(trackForward),
                trackLength_(trackLength),
                isWallPatch_(isWallPatch),

                allPositions_(allPositions),
                allScalars_(allScalars),
                allVectors_(allVectors)
            {}
    };


private:

    // Private data

        //- Particle is on mesh edge:
        //      const face& f = mesh.faces()[tetFace()]
        //      const edge e(f[meshEdgeStart_], f.nextLabel(meshEdgeStart_));
        //  Note that this real edge
        //  is also one of the edges of the face-triangle (from
        //  tetFace()+tetPt()).
        label meshEdgeStart_;

        //- Particle is on diagonal edge:
        //      const face& f = mesh.faces()[tetFace()]
        //      label faceBasePtI = mesh.tetBasePtIs()[faceI];
        //      label diagPtI = (faceBasePtI+diagEdge_)%f.size();
        //      const edge e(f[faceBasePtI], f[diagPtI]);
        label diagEdge_;

        //- Lifetime of particle. Particle dies when reaches 0.
        label lifeTime_;

        //- sampled positions
        DynamicList<point> sampledPositions_;

        //- sampled scalars
        List<DynamicList<scalar> > sampledScalars_;

        //- sampled vectors
        List<DynamicList<vector> > sampledVectors_;


    // Private Member Functions

        //- Construct current edge
        edge currentEdge() const;

        //- Check if inside current tet
        //void checkInside() const;

        //- Check if on current edge
        //void checkOnEdge() const;

        //- Check if point on triangle
        //void checkOnTriangle(const point&) const;

        //- Cross mesh edge into different face on same cell
        void crossEdgeConnectedFace(const edge& meshEdge);

        //- Cross diagonal edge into different triangle on same face,cell
        void crossDiagonalEdge();

        //- Track through single triangle
        scalar trackFaceTri(const vector& endPosition, label& minEdgeI);

        //- Do all patch interaction
        void patchInteraction(trackingData& td, const scalar trackFraction);

        //- Is current triangle in the track direction
        bool isTriAlongTrack(const point& endPosition) const;

        //- Equivalent of trackToFace
        scalar trackToEdge
        (
            trackingData& td,
            const vector& endPosition
        );

        //- Interpolate all quantities; return interpolated velocity.
        vector interpolateFields
        (
            const trackingData&,
            const point&,
            const label cellI,
            const label faceI
        );

        //- Interpolate all quantities, return updated velocity.
        vector sample(trackingData& td);

public:

    // Constructors

        //- Construct from components
        wallBoundedStreamLineParticle
        (
            const polyMesh& c,
            const vector& position,
            const label cellI,
            const label tetFaceI,
            const label tetPtI,
            const label meshEdgeStart,
            const label diagEdge,
            const label lifeTime
        );

        //- Construct from Istream
        wallBoundedStreamLineParticle
        (
            const polyMesh& c,
            Istream& is,
            bool readFields = true
        );

        //- Construct copy
        wallBoundedStreamLineParticle(const wallBoundedStreamLineParticle& p);

        //- Construct and return a clone
        autoPtr<particle> clone() const
        {
            return autoPtr<particle>(new wallBoundedStreamLineParticle(*this));
        }

        //- Factory class to read-construct particles used for
        //  parallel transfer
        class iNew
        {
            const polyMesh& mesh_;

        public:

            iNew(const polyMesh& mesh)
            :
                mesh_(mesh)
            {}

            autoPtr<wallBoundedStreamLineParticle> operator()
            (
                Istream& is
            ) const
            {
                return autoPtr<wallBoundedStreamLineParticle>
                (
                    new wallBoundedStreamLineParticle(mesh_, is, true)
                );
            }
        };


    // Member Functions

        // Tracking

            //- Track all particles to their end point
            bool move(trackingData&, const scalar trackTime);


            //- Overridable function to handle the particle hitting a patch
            //  Executed before other patch-hitting functions
            bool hitPatch
            (
                const polyPatch&,
                trackingData& td,
                const label patchI,
                const scalar trackFraction,
                const tetIndices& tetIs
            );

            //- Overridable function to handle the particle hitting a wedge
            void hitWedgePatch
            (
                const wedgePolyPatch&,
                trackingData& td
            );

            //- Overridable function to handle the particle hitting a
            //  symmetryPlane
            void hitSymmetryPatch
            (
                const symmetryPolyPatch&,
                trackingData& td
            );

            //- Overridable function to handle the particle hitting a cyclic
            void hitCyclicPatch
            (
                const cyclicPolyPatch&,
                trackingData& td
            );

            //- Overridable function to handle the particle hitting a
            //- processorPatch
            void hitProcessorPatch
            (
                const processorPolyPatch&,
                trackingData& td
            );

            //- Overridable function to handle the particle hitting a wallPatch
            void hitWallPatch
            (
                const wallPolyPatch&,
                trackingData& td,
                const tetIndices&
            );

            //- Overridable function to handle the particle hitting a polyPatch
            void hitPatch
            (
                const polyPatch&,
                trackingData& td
            );


        // Info

            //- Return info proxy.
            //  Used to print particle information to a stream
            InfoProxy<wallBoundedStreamLineParticle> info() const
            {
                return *this;
            }


        // I-O

            //- Read
            static void readFields(Cloud<wallBoundedStreamLineParticle>&);

            //- Write
            static void writeFields
            (
                const Cloud<wallBoundedStreamLineParticle>&
            );


    // Ostream Operator

        friend Ostream& operator<<
        (
            Ostream&,
            const wallBoundedStreamLineParticle&
        );

        friend Ostream& operator<<
        (
            Ostream&,
            const InfoProxy<wallBoundedStreamLineParticle>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    Foam::streamLineParticle

Description
    Particle class that samples fields as it passes through. Used in streamline
    calculation.

SourceFiles
    streamLineParticle.C

\*---------------------------------------------------------------------------*/

#ifndef streamLineParticle_H
#define streamLineParticle_H

#include "Particle.H"
#include "autoPtr.H"
#include "interpolationCellPoint.H"
#include "vectorList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class streamLineParticleCloud;

/*---------------------------------------------------------------------------*\
                           Class streamLineParticle Declaration
\*---------------------------------------------------------------------------*/

class streamLineParticle
:
    public Particle<streamLineParticle>
{

public:

    //- Class used to pass tracking data to the trackToFace function
    class trackData
    :
        public Particle<streamLineParticle>::trackData
    {

    public:


        const PtrList<interpolationCellPoint<scalar> >& vsInterp_;
        const PtrList<interpolationCellPoint<vector> >& vvInterp_;
        const label UIndex_;
        const bool trackForward_;

        DynamicList<vectorList>& allPositions_;
        List<DynamicList<scalarList> >& allScalars_;
        List<DynamicList<vectorList> >& allVectors_;


        // Constructors

            trackData
            (
                Cloud<streamLineParticle>& cloud,
                const PtrList<interpolationCellPoint<scalar> >& vsInterp,
                const PtrList<interpolationCellPoint<vector> >& vvInterp,
                const label UIndex,
                const bool trackForward,
                DynamicList<List<point> >& allPositions,
                List<DynamicList<scalarList> >& allScalars,
                List<DynamicList<vectorList> >& allVectors
            )
            :
                Particle<streamLineParticle>::trackData(cloud),
                vsInterp_(vsInterp),
                vvInterp_(vvInterp),
                UIndex_(UIndex),
                trackForward_(trackForward),
                allPositions_(allPositions),
                allScalars_(allScalars),
                allVectors_(allVectors)
            {}
    };


private:

    // Private data

        //- Lifetime of particle. Particle dies when reaches 0.
        label lifeTime_;

        //- sampled positions
        DynamicList<point> sampledPositions_;

        //- sampled scalars
        List<DynamicList<scalar> > sampledScalars_;

        //- sampled vectors
        List<DynamicList<vector> > sampledVectors_;


    // Private Member Functions

        //- Interpolate all quantities; return interpolated velocity.
        vector interpolateFields
        (
            const streamLineParticle::trackData&,
            const point&,
            const label cellI
        );


public:

    //- Run-time type information
    TypeName("streamLineParticle");



    // Constructors

        //- Construct from components
        streamLineParticle
        (
            const Cloud<streamLineParticle>& c,
            const vector& position,
            const label celli,
            const label lifeTime
        );

        //- Construct from Istream
        streamLineParticle
        (
            const Cloud<streamLineParticle>& c,
            Istream& is,
            bool readFields = true
        );

        //- Construct copy
        streamLineParticle(const streamLineParticle& c);

        //- Construct and return a clone
        autoPtr<streamLineParticle> clone() const
        {
            return autoPtr<streamLineParticle>
            (
                new streamLineParticle(*this)
            );
        }


    // Member Functions



        // Tracking

            //- Track all particles to their end point
            bool move(trackData&);


            //- Overridable function to handle the particle hitting a patch
            //  Executed before other patch-hitting functions
            bool hitPatch
            (
                const polyPatch&,
                streamLineParticle::trackData& td,
                const label patchI
            );
            bool hitPatch
            (
                const polyPatch&,
                int&,
                const label patchI
           );

            //- Overridable function to handle the particle hitting a wedge
            void hitWedgePatch
            (
                const wedgePolyPatch&,
                streamLineParticle::trackData& td
            );
            void hitWedgePatch
            (
                const wedgePolyPatch&,
                int&
            );

            //- Overridable function to handle the particle hitting a
            //  symmetryPlane
            void hitSymmetryPatch
            (
                const symmetryPolyPatch&,
                streamLineParticle::trackData& td
            );
            void hitSymmetryPatch
            (
                const symmetryPolyPatch&,
                int&
            );

            //- Overridable function to handle the particle hitting a cyclic
            void hitCyclicPatch
            (
                const cyclicPolyPatch&,
                streamLineParticle::trackData& td
            );
            void hitCyclicPatch
            (
                const cyclicPolyPatch&,
                int&
            );

            //- Overridable function to handle the particle hitting a
            //- processorPatch
            void hitProcessorPatch
            (
                const processorPolyPatch&,
                streamLineParticle::trackData& td
            );
            void hitProcessorPatch
            (
                const processorPolyPatch&,
                int&
            );

            //- Overridable function to handle the particle hitting a wallPatch
            void hitWallPatch
            (
                const wallPolyPatch&,
                streamLineParticle::trackData& td
            );
            void hitWallPatch
            (
                const wallPolyPatch&,
                int&
            );

            //- Overridable function to handle the particle hitting a polyPatch
            void hitPatch
            (
                const polyPatch&,
                streamLineParticle::trackData& td
            );
            void hitPatch
            (
                const polyPatch&,
                int&
            );


        // I-O

            //- Read
            static void readFields(Cloud<streamLineParticle>&);

            //- Write
            static void writeFields(const Cloud<streamLineParticle>&);


    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const streamLineParticle&);


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

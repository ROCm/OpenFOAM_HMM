/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2011 OpenCFD Ltd.
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
    Foam::monoatomic

Description
    Foam::monoatomic

SourceFiles
    monoatomicI.H
    monoatomic.C
    monoatomicIO.C

\*---------------------------------------------------------------------------*/

#ifndef monoatomic_H
#define monoatomic_H

#include "particle.H"
#include "IOstream.H"
#include "autoPtr.H"
#include "diagTensor.H"
#include "constPropSite.H"
#include "MoleculeCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class monoatomic Declaration
\*---------------------------------------------------------------------------*/

class monoatomic
:
    public particle
{

public:

    // Values of special that are less than zero are for built-in functionality.
    // Values greater than zero are user specifiable/expandable (i.e. test
    // special_ >= SPECIAL_USER)

    enum specialTypes
    {
        SPECIAL_TETHERED = -1,
        SPECIAL_FROZEN   = -2,
        NOT_SPECIAL      = 0,
        SPECIAL_USER     = 1
    };

    //- Class to hold monoatomic constant properties
    class constantProperties
    {
        // Private data

            //- Sites of mass, charge or interaction
            FixedList<constPropSite, 1> sites_;

            //- Which sites require electrostatic interactions
            FixedList<label, 1> electrostaticSites_;

            //- Which sites require pair interactions
            FixedList<label, 1> pairPotSites_;

            //-
            scalar mass_;


    public:

        //-
        inline constantProperties();

        //- Construct from dictionary
        inline constantProperties
        (
            const dictionary& dict,
            const List<label>& siteIds
        );

        // Member functions

        //-
        inline const FixedList<constPropSite, 1>& sites() const;

        //-
        inline const FixedList<label, 1>& pairPotSites() const;

        //-
        inline const FixedList<label, 1>& electrostaticSites() const;

        //-
        inline label degreesOfFreedom() const;

        //-
        inline scalar mass() const;

        //-
        inline label nSites() const;
    };


    //- Class used to pass tracking data to the trackToFace function
    class trackingData
    :
        public particle::TrackingData<MoleculeCloud<monoatomic> >
    {
    public:

        enum trackPart
        {
            tpFirstVelocityHalfStep,
            tpLinearTrack,
            tpRotationalTrack,
            tpSecondVelocityHalfStep,
            tpAccess
        };

    private:

        // label specifying which part of the integration algorithm is taking
        label part_;


    public:

        // Constructors

            trackingData(MoleculeCloud<monoatomic>& cloud, trackPart part)
            :
                particle::TrackingData<MoleculeCloud<monoatomic> >(cloud),
                part_(part)
            {}

        // Member functions

            inline label part() const
            {
                return part_;
            }
    };


private:

    // Private data

          //- Linear velocity of monoatomic
        vector v_;

        //- Total linear acceleration of monoatomic
        vector a_;

            //-
        vector specialPosition_;

        //-
        scalar potentialEnergy_;

        // - r_ij f_ij, stress dyad
        tensor rf_;

//         // - r_ij outer product f_ij: virial contribution
//         tensor rDotf_;

        //-
        label special_;

        //-
        label id_;

        //-
        FixedList<vector, 1> siteForces_;

        //-
        FixedList<vector, 1> sitePositions_;


public:

    //- Runtime type information
    TypeName("monoatomic");

    friend class Cloud<monoatomic>;

    // Constructors

        //- Construct with macroscopic description
        inline monoatomic
        (
            const polyMesh& mesh,
            const vector& position,
            const label cellI,
            const label tetFaceI,
            const label tetPtI,
            const scalar temperature,
            const vector& bulkVelocity,
            const vector& specialPosition,
            const constantProperties& constProps,
            trackingData& td,
            const label special,
            const label id
        );

        //- Construct from all components
        inline monoatomic
        (
            const polyMesh& mesh,
            const vector& position,
            const label cellI,
            const label tetFaceI,
            const label tetPtI,
            const vector& v,
            const vector& a,
            const vector& specialPosition,
            const constantProperties& constProps,
            const label special,
            const label id
        );

        //- Construct from Istream
        monoatomic
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields = true
        );

        //- Construct and return a clone
        autoPtr<particle> clone() const
        {
            return autoPtr<particle>(new monoatomic(*this));
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

            autoPtr<monoatomic> operator()(Istream& is) const
            {
                return autoPtr<monoatomic>(new monoatomic(mesh_, is, true));
            }
        };


    // Member Functions

        // Tracking

            //-
            bool move(trackingData&, const scalar trackTime);

            //-
            virtual void transformProperties(const tensor& T);

            //-
            virtual void transformProperties(const vector& separation);

            //-
            void setSitePositions(const constantProperties& constProps);

            //-
            void setSiteSizes(label size);


        // Access

            //-
            inline const vector& v() const;

            //-
            inline vector& v();

            //-
            inline const vector& a() const;

            //-
            inline vector& a();

            //-
            inline const FixedList<vector, 1>& siteForces() const;

            //-
            inline FixedList<vector, 1>& siteForces();

            //-
            inline const FixedList<vector, 1>& sitePositions() const;

            //-
            inline FixedList<vector, 1>& sitePositions();

            //-
            inline const vector& specialPosition() const;

            //-
            inline vector& specialPosition();

            //-
            inline scalar potentialEnergy() const;

            //-
            inline scalar& potentialEnergy();

            //-
            inline const tensor& rf() const;

            //-
            inline tensor& rf();

            //-
            inline label special() const;

            //-
            inline bool tethered() const;

            //-
            inline label id() const;


    // Member Operators

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

        //- Overridable function to handle the particle hitting a processorPatch
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


    // I-O

        //- Read
        static void readFields(Cloud<monoatomic>& mC);

        //- Write
        static void writeFields(const Cloud<monoatomic>& mC);

        //- Show info
        static void info(trackingData& td);


    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const monoatomic&);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "monoatomicI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    Foam::polyatomic

Description
    Foam::polyatomic

SourceFiles
    polyatomicI.H
    polyatomic.C
    polyatomicIO.C

\*---------------------------------------------------------------------------*/

#ifndef polyatomic_H
#define polyatomic_H

#include "particle.H"
#include "IOstream.H"
#include "autoPtr.H"
#include "diagTensor.H"
#include "constPropSite.H"
#include "MoleculeCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class ParcelType>
class DsmcParcel;

/*---------------------------------------------------------------------------*\
                           Class polyatomic Declaration
\*---------------------------------------------------------------------------*/

class polyatomic
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

    //- Class to hold polyatomic constant properties
    class constantProperties
    {
        // Private data

            //- Sites of mass, charge or interaction
            List<constPropSite> sites_;

            //- Which sites require electrostatic interactions
            List<label> electrostaticSites_;

            //- Which sites require pair interactions
            List<label> pairPotSites_;

            //- Moment of intertia (in principal axis configiration)
            diagTensor momentOfInertia_;

            //-
            scalar mass_;


        // Private Member Functions

            //-
            bool linearMoleculeTest() const;


    public:

        //-
        inline constantProperties();

        //- Construct from dictionary
        inline constantProperties
        (
            const dictionary& dict,
            const labelList& siteIds
        );

        // Member functions

        //-
        inline const List<constPropSite>& sites() const;

        //-
        inline const List<label>& pairPotSites() const;

        //-
        inline const List<label>& electrostaticSites() const;

        //-
        inline const diagTensor& momentOfInertia() const;

        //-
        inline bool linearMolecule() const;

        //-
        inline bool pointMolecule() const;

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
        public particle::TrackingData<MoleculeCloud<polyatomic> >
    {
        // label specifying which part of the integration algorithm is taking
        label part_;


    public:

        // Constructors

            trackingData(MoleculeCloud<polyatomic>& cloud, label part)
            :
                particle::TrackingData<MoleculeCloud<polyatomic> >(cloud),
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

        //- Orientation, stored as the rotation tensor to transform
        //  from the polyatomic to the global reference frame, i.e.:
        //  globalVector = Q_ & bodyLocalVector
        //  bodyLocalVector = Q_.T() & globalVector
        tensor Q_;

        //- Linear velocity of polyatomic
        vector v_;

        //- Total linear acceleration of polyatomic
        vector a_;

        //- Angular momentum of polyatomic, in body local reference frame
        vector pi_;

        //- Total torque on polyatomic, in body local reference frame
        vector tau_;

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
        List<vector> siteForces_;

        //-
        List<vector> sitePositions_;


    // Private Member Functions


        //-
        tensor rotationTensorX(scalar deltaT) const;

        //-
        tensor rotationTensorY(scalar deltaT) const;

        //-
        tensor rotationTensorZ(scalar deltaT) const;


public:

    //- Runtime type information
    TypeName("polyatomic");

    friend class Cloud<polyatomic>;

    // Constructors

        //- Construct from all components
        inline polyatomic
        (
            const polyMesh& mesh,
            const vector& position,
            const label cellI,
            const label tetFaceI,
            const label tetPtI,
            const tensor& Q,
            const vector& v,
            const vector& a,
            const vector& pi,
            const vector& tau,
            const vector& specialPosition,
            const constantProperties& constProps,
            const label special,
            const label id
        );

        //- Construct from Istream
        polyatomic
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields = true
        );

        //- Construct and return a clone
        autoPtr<particle> clone() const
        {
            return autoPtr<particle>(new polyatomic(*this));
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

            autoPtr<polyatomic> operator()(Istream& is) const
            {
                return autoPtr<polyatomic>(new polyatomic(mesh_, is, true));
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
            inline const tensor& Q() const;

            //-
            inline tensor& Q();

            //-
            inline const vector& v() const;

            //-
            inline vector& v();

            //-
            inline const vector& a() const;

            //-
            inline vector& a();

            //-
            inline const vector& pi() const;

            //-
            inline vector& pi();

            //-
            inline const vector& tau() const;

            //-
            inline vector& tau();

            //-
            inline const List<vector>& siteForces() const;

            //-
            inline List<vector>& siteForces();

            //-
            inline const List<vector>& sitePositions() const;

            //-
            inline List<vector>& sitePositions();

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
        static void readFields(Cloud<polyatomic>& mC);

        //- Write
        static void writeFields(const Cloud<polyatomic>& mC);


    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const polyatomic&);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "polyatomicI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    Foam::molecule

Description
    Foam::molecule

SourceFiles
    moleculeI.H
    molecule.C
    moleculeIO.C

\*---------------------------------------------------------------------------*/

#ifndef molecule_H
#define molecule_H

#include "Particle.H"
#include "IOstream.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class moleculeCloud;

/*---------------------------------------------------------------------------*\
                           Class molecule Declaration
\*---------------------------------------------------------------------------*/

class molecule
:
    public Particle<molecule>
{

public:

    //- Class to hold molecule constant properties
    class constantProperties
    {

        // Private data

        List<vector> siteReferencePositions_;

        List<scalar> siteCharges_;

        List<label> siteIds_;

        diagTensor momentOfInertia_;

        scalar mass_;


    public:

        //- Constructor
        inline constantProperties();

        // Member functions

        inline const List<vector>& siteReferencePositions() const;

        inline const List<scalar>& siteCharges() const;

        inline const List<label>& siteIds() const;

        inline const diagTensor& momentOfInertia() const;

        inline scalar mass() const;

    };


    //- Class used to pass tracking data to the trackToFace function
    class trackData
    :
        public Particle<molecule>::trackData
    {
        moleculeCloud& molCloud_;

        // label specifying which part of the integration algorithm is taking
        // place (i.e. leapfrog 1 or leapfrog 2. Predictor or Corrector)
        label part_;


    public:

        // Constructors

            inline trackData
            (
                moleculeCloud& molCloud,
                label part
            );

        // Member functions

            inline moleculeCloud& molCloud();

            inline label part() const;
    };


private:

    // Private data

        //- Be careful with the ordering of data.
        //  It has an impact on binary transfer:
        //    -# Put the largest data members 1st
        //    -# Pair up labels,
        //    -# Don't go scalar-label, scalar-label, because in 64bit mode,
        //       the labels will be padded by 4bytes.

        // Centre of mass velocities and accelerations

        vector u_;

        vector a_;

        vector omega_;

        vector alpha_;

        List<vector> siteForces_;

        List<vector> sitePostions_;

        vector specialPosition_;

        scalar potentialEnergy_;

        // - r_ij f_ij, stress dyad
        tensor rf_;

        const constantProperties& constProps_;

        label special_;

        label id_;

public:

    friend class Cloud<molecule>;

    // Constructors

        //- Construct from components
        inline molecule
        (
            const Cloud<molecule>& c,
            const vector& position,
            const label celli,
            const scalar mass,
            const vector& U,
            const vector& A,
            const vector& specialPosition,
            const label special,
            const label id
        );

        //- Construct from Istream
        molecule
        (
            const Cloud<molecule>& c,
            Istream& is,
            bool readFields = true
        );

        //- Construct and return a clone
        autoPtr<molecule> clone() const
        {
            return autoPtr<molecule>(new molecule(*this));
        }


    // Member Functions

        void transformProperties(const tensor& T);

        void transformProperties(const vector& separation);

        // Access

            inline label id() const;

            inline scalar mass() const;

            inline const vector& U() const;
            inline vector& U();

            inline const vector& A() const;
            inline vector& A();

            inline scalar potentialEnergy() const;
            inline scalar& potentialEnergy();

            inline const tensor& rf() const;
            inline tensor& rf();

            inline label special() const;

            inline const vector& specialPosition() const;
            inline vector& specialPosition();


        //- Tracking
        bool move(trackData&);


    // Member Operators

        //- Overridable function to handle the particle hitting a processorPatch
        void hitProcessorPatch
        (
            const processorPolyPatch&,
            molecule::trackData& td
        );

        //- Overridable function to handle the particle hitting a processorPatch
        //  without trackData
        void hitProcessorPatch
        (
            const processorPolyPatch&,
            int&
        );

        //- Overridable function to handle the particle hitting a wallPatch
        void hitWallPatch
        (
            const wallPolyPatch&,
            molecule::trackData& td
        );

        //- Overridable function to handle the particle hitting a wallPatch
        //  without trackData
        void hitWallPatch
        (
            const wallPolyPatch&,
            int&
        );

        //- Overridable function to handle the particle hitting a polyPatch
        void hitPatch
        (
            const polyPatch&,
            molecule::trackData& td
        );

        //- Overridable function to handle the particle hitting a polyPatch
        //  without trackData
        void hitPatch
        (
            const polyPatch&,
            int&
        );

        static void readFields(moleculeCloud& mC);

        static void writeFields(const moleculeCloud& mC);


    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const molecule&);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "moleculeI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

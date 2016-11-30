/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenCFD Ltd.
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
    Foam::injectedParticle

Description

SourceFiles
    injectedParticle.C
    injectedParticleIO.C

\*---------------------------------------------------------------------------*/

#ifndef injectedParticle_H
#define injectedParticle_H

#include "particle.H"
#include "IOstream.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class injectedParticle;

// Forward declaration of friend functions

Ostream& operator<<
(
    Ostream&,
    const injectedParticle&
);

/*---------------------------------------------------------------------------*\
                         Class injectedParticle Declaration
\*---------------------------------------------------------------------------*/

class injectedParticle
:
    public particle
{
public:

    //- Size in bytes of the fields
    static const std::size_t sizeofFields;


protected:

    // Protected data

        // Particle properties

            //- Start of injection [s]
            scalar soi_;

            //- Diameter [m]
            scalar d_;

            //- Velocity [m/s]
            vector U_;


public:

    // Static data members

        //- Runtime type information
        TypeName("injectedParticle");

        //- String representation of properties
        AddToPropertyList
        (
            particle,
            " soi"
          + " d"
          + " (Ux Uy Uz)";
        );

        //- String representation of property types
        AddToPropertyTypes
        (
            particle,
            "{scalar"
          + " scalar"
          + " vector}"
        );


    // Constructors

        //- Construct from components
        inline injectedParticle
        (
            const polyMesh& mesh,
            const vector& position,
            const scalar soi,
            const scalar d,
            const vector& U
        );

        //- Construct from Istream
        injectedParticle
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields = true
        );

        //- Construct as a copy
        injectedParticle(const injectedParticle& p);

        //- Construct as a copy
        injectedParticle(const injectedParticle& p, const polyMesh& mesh);

        //- Construct and return a (basic particle) clone
        virtual autoPtr<particle> clone() const
        {
            return autoPtr<particle>(new injectedParticle(*this));
        }

        //- Construct and return a (basic particle) clone
        virtual autoPtr<particle> clone(const polyMesh& mesh) const
        {
            return autoPtr<particle>(new injectedParticle(*this, mesh));
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

            autoPtr<injectedParticle> operator()(Istream& is) const
            {
                return autoPtr<injectedParticle>
                (
                    new injectedParticle(mesh_, is, true)
                );
            }
        };


    // Member Functions

        // Access

            //- Return const access to the start of injection
            inline scalar soi() const;

            //- Return const access to diameter
            inline scalar d() const;

            //- Return const access to velocity
            inline const vector& U() const;


        // Edit

            //- Return start of injection
            inline scalar& soi();

            //- Return access to diameter
            inline scalar& d();

            //- Return access to velocity
            inline vector& U();


        // I-O

            //- Read
            static void readFields(Cloud<injectedParticle>& c);

            //- Write
            static void writeFields(const Cloud<injectedParticle>& c);


    // Ostream Operator

        friend Ostream& operator<<
        (
            Ostream&,
            const injectedParticle&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "injectedParticleI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

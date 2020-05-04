/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2019 OpenCFD Ltd.
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
    Primarily stores particle properties so that it can be injected at a later
    time.  Note that this stores its own local position as opposed to the
    base particle class barycentric coordinates since the particle is not
    (usually) attached to a mesh, and instead used for post-processing.

SourceFiles
    injectedParticle.C
    injectedParticleIO.C

SeeAlso
    Foam::functionObjects::extractEulerianParticles

\*---------------------------------------------------------------------------*/

#ifndef injectedParticle_H
#define injectedParticle_H

#include "particle.H"
#include "IOstream.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class injectedParticle;
Ostream& operator<<(Ostream&, const injectedParticle&);

/*---------------------------------------------------------------------------*\
                         Class injectedParticle Declaration
\*---------------------------------------------------------------------------*/

class injectedParticle
:
    public particle
{
protected:

    // Protected data

        // Particle properties

            //- Position
            point position_;

            //- Tag
            label tag_;

            //- Start of injection [s]
            scalar soi_;

            //- Diameter [m]
            scalar d_;

            //- Velocity [m/s]
            vector U_;


public:

    // Static Data Members

        //- Size in bytes of the fields
        static const std::size_t sizeofFields;

        //- Runtime type information
        TypeName("injectedParticle");

        //- String representation of properties
        AddToPropertyList
        (
            particle,
            " tag"
          + " soi"
          + " d"
          + " (Ux Uy Uz)";
        );


    // Constructors

        //- Construct from a position and a cell.
        //  Searches for the rest of the required topology.
        //  Other properties are zero initialised.
        inline injectedParticle
        (
            const polyMesh& mesh,
            const vector& position,
            const label celli = -1
        );

        //- Construct from components
        inline injectedParticle
        (
            const polyMesh& mesh,
            const vector& position,
            const label tag,
            const scalar soi,
            const scalar d,
            const vector& U,
            const bool doLocate = true
        );

        //- Construct from Istream
        injectedParticle
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields = true,
            bool newFormat = true
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

            //- Return const access to the tag
            inline label tag() const;

            //- Return const access to the start of injection
            inline scalar soi() const;

            //- Return const access to diameter
            inline scalar d() const;

            //- Return const access to velocity
            inline const vector& U() const;


        // Edit

            //- Return the tag
            inline label& tag();

            //- Return the start of injection
            inline scalar& soi();

            //- Return access to diameter
            inline scalar& d();

            //- Return access to velocity
            inline vector& U();


        // I-O

            //- Read fields
            static void readFields(Cloud<injectedParticle>& c);

            //- Write individual parcel properties to stream
            void writeProperties
            (
                Ostream& os,
                const wordRes& filters,
                const word& delim,
                const bool namesOnly
            ) const;

            //- Write fields
            static void writeFields(const Cloud<injectedParticle>& c);

            //- Read particle fields as objects from the obr registry
            static void readObjects
            (
                Cloud<injectedParticle>& c,
                const objectRegistry& obr
            );

            //- Write particle fields as objects into the obr registry
            static void writeObjects
            (
                const Cloud<injectedParticle>& c,
                objectRegistry& obr
            );

            //- Write the particle position and cell
            //  Note: This uses the local particle position, and bypasses the
            //  barycentric description
            virtual void writePosition(Ostream&) const;


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

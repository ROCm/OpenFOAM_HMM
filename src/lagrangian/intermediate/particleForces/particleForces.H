/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2010 OpenCFD Ltd.
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
    Foam::particleForces

Description
    Provides a mechanism to calculate particle forces
    Note: forces are force per unit mass (accelerations)

SourceFiles
    particleForces.C

\*---------------------------------------------------------------------------*/

#ifndef particleForces_H
#define particleForces_H

#include "dictionary.H"
#include "Switch.H"
#include "vector.H"
#include "volFieldsFwd.H"
#include "interpolation.H"
#include "tetIndices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward class declarations
class fvMesh;

/*---------------------------------------------------------------------------*\
                      Class particleForces Declaration
\*---------------------------------------------------------------------------*/

class particleForces
{
    // Private data

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- The particleForces dictionary
        const dictionary& dict_;

        //- Gravity
        const vector g_;

        //- Velocity gradient field
        const volTensorField* gradUPtr_;

        //- HdotGradH interpolator
        const interpolation<vector>* HdotGradHInterPtr_;


        // Forces to include in particle motion evaluation

            //- Gravity
            Switch gravity_;

            //- Virtual mass
            Switch virtualMass_;

            //- Virtual mass force coefficient
            scalar Cvm_;

            //- Pressure gradient
            Switch pressureGradient_;

            //- Paramagnetic force
            Switch paramagnetic_;

            //- Magnetic susceptibility of particle
            scalar magneticSusceptibility_;


        // Additional info

            //- Name of velocity field - default = "U"
            const word UName_;

            //- Name of paramagnetic field strength field - default =
            //  "HdotGradH"
            const word HdotGradHName_;


    // Private Member Functions

        //- Delete cached carrier fields
        void deleteFields();


public:

    // Constructors

        //- Construct null from mesh reference
        particleForces(const fvMesh& mesh);

        //- Construct from mesh, dictionary and gravity
        particleForces
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const vector& g
        );

        //- Construct copy
        particleForces(const particleForces& f);


    //- Destructor
    ~particleForces();


    // Member Functions

        // Access

            //- Return the particleForces dictionary
            const dictionary& dict() const;

            //- Return the gravity vector
            const vector& g() const;

            //- Return gravity force activate switch
            Switch gravity() const;

            //- Return virtual mass force activate switch
            Switch virtualMass() const;

            //- Return virtual mass force coefficient
            scalar Cvm() const;

            //- Return pressure gradient force activate switch
            Switch pressureGradient() const;

            //- Return name of velocity field
            const word& UName() const;

            //- Return paramagnetic force activate switch
            Switch paramagnetic() const;

            //- Return magnetic susceptibility
            scalar magneticSusceptibility() const;

            //- Return name of paramagnetic field strength field
            const word& HdotGradHName() const;


       // Evaluation

            //- Cache carrier fields
            void cacheFields
            (
                const bool store,
                const dictionary& interpolationSchemes
            );

            //- Calculate action/reaction forces between carrier and particles
            vector calcCoupled
            (
                const vector& position,
                const tetIndices& tetIs,
                const scalar dt,
                const scalar rhoc,
                const scalar rho,
                const vector& Uc,
                const vector& U,
                const scalar d
            ) const;

            //- Calculate external forces applied to the particles
            vector calcNonCoupled
            (
                const vector& position,
                const tetIndices& tetIs,
                const scalar dt,
                const scalar rhoc,
                const scalar rho,
                const vector& Uc,
                const vector& U,
                const scalar d
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

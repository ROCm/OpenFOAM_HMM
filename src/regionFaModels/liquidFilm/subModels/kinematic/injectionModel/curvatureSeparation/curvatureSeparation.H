/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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
    Foam::regionModels::areaSurfaceFilmModels::curvatureSeparation

Description
    Curvature film separation model

    Assesses film curvature via the mesh geometry and calculates a force
    balance of the form:

        F_sum = F_inertial + F_body + F_surface

    If F_sum < 0, the film separates. Similarly, if F_sum > 0 the film will
    remain attached.

    Based on description given by
        Owen and D. J. Ryley. The flow of thin liquid films around corners.
        International Journal of Multiphase Flow, 11(1):51-62, 1985.

SourceFiles
    curvatureSeparation.C

\*---------------------------------------------------------------------------*/

#ifndef curvatureSeparation_H
#define curvatureSeparation_H

#include "injectionModel.H"
#include "faCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                   Class curvatureSeparation Declaration
\*---------------------------------------------------------------------------*/

class curvatureSeparation
:
    public injectionModel
{
    // Private Member Functions

        //- No copy construct
        curvatureSeparation(const curvatureSeparation&) = delete;

        //- No copy assignment
        void operator=(const curvatureSeparation&) = delete;


protected:

    // Protected Data

        //- Gradient of surface normals
        areaTensorField gradNHat_;

        //- Minimum gravity driven film thickness (non-dimensionalised delta/R1)
        scalar deltaByR1Min_;

        //- List of radii for patches - if patch not defined, radius
        //- calculated based on mesh geometry
        scalar definedPatchRadii_;

        //- Magnitude of gravity vector
        scalar magG_;

        //- Direction of gravity vector
        vector gHat_;

        //- Threshold force for separation
        scalar fThreshold_;

        //- Minimum inv R1 for separation
        scalar minInvR1_;


    // Protected Member Functions

        //- Calculate local (inverse) radius of curvature
        tmp<areaScalarField> calcInvR1
        (
            const areaVectorField& U,
            const scalarField& calcCosAngle
        ) const;

        //- Calculate the cosine of the angle between gravity vector and
        //- cell out flow direction
        tmp<scalarField> calcCosAngle(const edgeScalarField& phi) const;


public:

    //- Runtime type information
    TypeName("curvatureSeparation");


    // Constructors

        //- Construct from surface film model
        curvatureSeparation
        (
            liquidFilmBase& film,
            const dictionary& dict
        );


    //- Destructor
    virtual ~curvatureSeparation() = default;


    // Member Functions

        // Evolution

            //- Correct
            virtual void correct
            (
                scalarField& availableMass,
                scalarField& massToInject,
                scalarField& diameterToInject
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

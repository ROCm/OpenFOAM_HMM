/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2020 PCOpt/NTUA
    Copyright (C) 2013-2020 FOSS GP
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::incompressible::FIBase

Description
    Base class for Field Integral-based sensitivity derivatives

SourceFiles
    FIBase.C

\*---------------------------------------------------------------------------*/

#ifndef FIBaseIncompressible_H
#define FIBaseIncompressible_H

#include "adjointSensitivityIncompressible.H"
#include "shapeSensitivitiesIncompressible.H"
#include "adjointEikonalSolverIncompressible.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace incompressible
{

/*---------------------------------------------------------------------------*\
                            Class FIBase Declaration
\*---------------------------------------------------------------------------*/

class FIBase
:
    public shapeSensitivities
{
protected:

    // Protected data

        //- grad(dx/db) multiplier
        volTensorField gradDxDbMult_;

        //- div(dx/db) multiplier
        scalarField divDxDbMult_;

        //- dx/db multiplier coming from fvOptions
        vectorField optionsDxDbMult_;

        //- Include distance variation in sens computation
        bool includeDistance_;

        //- Adjoint eikonal equation solver
        autoPtr<adjointEikonalSolver> eikonalSolver_;


    // Protected Member Functions

        //- Read options and update solver pointers if necessary
        void read();


private:

    // Private Member Functions

        //- No copy construct
        FIBase(const FIBase&) = delete;

        //- No copy assignment
        void operator=(const FIBase&) = delete;


public:

    //- Runtime type information
    TypeName("volumetricBSplinesFI");


    // Constructors

        //- Construct from components
        FIBase
        (
            const fvMesh& mesh,
            const dictionary& dict,
            incompressibleVars& primalVars,
            incompressibleAdjointVars& adjointVars,
            objectiveManager& objectiveManager
        );


    //- Destructor
    virtual ~FIBase() = default;


    // Member Functions

        //- Read dict if changed
        virtual bool readDict(const dictionary& dict);

        //- Accumulate sensitivity integrands
        virtual void accumulateIntegrand(const scalar dt);

        //- Assemble sensitivities
        virtual void assembleSensitivities() = 0;

        //- Zero sensitivity fields and their constituents
        virtual void clearSensitivities();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

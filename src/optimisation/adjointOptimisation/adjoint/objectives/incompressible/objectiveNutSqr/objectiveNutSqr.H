/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2020 PCOpt/NTUA
    Copyright (C) 2013-2020 FOSS GP
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
    Foam::objectives::objectiveNutSqr

Description
    Objective qualitatively quantifying noise through the integral of the
    squared turbulent viscosity over specified cellZones. Requires the adjoint
    to the turbulence model to be incorporated into the optimisation loop.

    Reference:
    \verbatim
        Objective function initially presented in

            Papoutsis-Kiachagias, E. M., Magoulas, N., Mueller, J., Othmer, C.,
            & Giannakoglou, K. C. (2015).
            Noise reduction in car aerodynamics using a surrogate objective
            function and the continuous adjoint method with wall functions.
            Computers & Fluids, 122(20), 223-232.
            https://doi.org/10.1016/j.compfluid.2015.09.002
    \endverbatim

SourceFiles
    objectiveNutSqrNoise.C

\*---------------------------------------------------------------------------*/

#ifndef objectiveNutSqr_H
#define objectiveNutSqr_H

#include "objectiveIncompressible.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace objectives
{

/*---------------------------------------------------------------------------*\
                       Class objectiveNutSqr Declaration
\*---------------------------------------------------------------------------*/

class objectiveNutSqr
:
    public objectiveIncompressible
{
    // Private data

        //- Where to define the objective
        labelList zones_;


public:

    //- Runtime type information
    TypeName("nutSqr");


    // Constructors

        //- From components
        objectiveNutSqr
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const word& adjointSolverName,
            const word& primalSolverName
        );


    //- Destructor
    virtual ~objectiveNutSqr() = default;


    // Member Functions

        //- Return the objective function value
        scalar J();

        //- Update field to be added to the adjoint turbulence model PDE
        //  The equivalent for the second turbulence model variable should
        //  be implemented for turbulence models other than SA
        void update_dJdTMvar1();

        //- Update field to be added to be added to volume-based
        //- sensitivity derivatives, emerging from delta ( dV ) / delta b
        void update_divDxDbMultiplier();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace objectiveNutSqrNoise
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

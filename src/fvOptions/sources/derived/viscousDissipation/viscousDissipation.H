/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd
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
    Foam::fv::viscousDissipation

Group
    grpFvOptionsSources

Description
    Calculates and applies the viscous dissipation energy source to the energy
    equation.

Usage
    Example usage:
    \verbatim
    fields          (h);                    // Name of energy field
    \endverbatim

SourceFiles
    viscousDissipation.C

\*---------------------------------------------------------------------------*/

#ifndef viscousDissipation_H
#define viscousDissipation_H

#include "fvOption.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                    Class viscousDissipation Declaration
\*---------------------------------------------------------------------------*/

class viscousDissipation
:
    public option
{
    // Private data

        //- Name of velocity field; default = U
        word UName_;

        //- Name of the rho field for incompressible solvers
        word rhoName_;

        //- Density for single phase flows
        dimensionedScalar rho_;


     // Private Member Functions

        //- Return the viscosity field
        tmp<volSymmTensorField> devRhoReff() const;

        //- Disallow default bitwise copy construct
        viscousDissipation(const viscousDissipation&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const viscousDissipation&) = delete;

private:

    // Private member functions


        //- Return rho field
        tmp<volScalarField> rho() const;


public:

    //- Runtime type information
    TypeName("viscousDissipation");


    // Constructors

        //- Construct from explicit source name and mesh
        viscousDissipation
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    // Member Functions

        // Evaluate

            //- Add explicit contribution to compressible energy equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<scalar>& eqn,
                const label fieldi
            );


        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict)
            {
                return true;
            }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2011 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License

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
    Foam::combustionModels::infinitelyFastChemistry

Description
    Simple infinitely fast chemistry combustion model based on the principle
    mixed is burnt. Additional parameter C is used to distribute the heat
    release rate.in time

SourceFiles
    infinitelyFastChemistry.C

\*---------------------------------------------------------------------------*/

#ifndef infinitelyFastChemistry_H
#define infinitelyFastChemistry_H

#include "combustionModel.H"
#include "singleStepReactingMixture.H"
#include "thermoPhysicsTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{

/*---------------------------------------------------------------------------*\
                  Class infinitelyFastChemistry Declaration
\*---------------------------------------------------------------------------*/

class infinitelyFastChemistry
:
    public combustionModel
{
    // Private data

        //- Model constant
        scalar C_;

        //- Reference to singleStepReactingMixture mixture
        singleStepReactingMixture<gasThermoPhysics>& singleMixture_;

        //- Normalised consumption rate of (fu - fres)
        volScalarField wFuelNorm_;


    // Private Member Functions

        //- Disallow copy construct
        infinitelyFastChemistry(const infinitelyFastChemistry&);

        //- Disallow default bitwise assignment
        void operator=(const infinitelyFastChemistry&);


public:

    //- Runtime type information
    TypeName("infinitelyFastChemistry");


    // Constructors

        //- Construct from components
        infinitelyFastChemistry
        (
            const dictionary& combustionProps,
            hsCombustionThermo& thermo,
            const compressible::turbulenceModel& turbulence,
            const surfaceScalarField& phi,
            const volScalarField& rho
        );


    //- Destructor
    virtual ~infinitelyFastChemistry();


    // Member Functions

        // Evolution

            //- Correct combustion rate
            virtual void correct();

            //- Fuel consumption rate matrix, i.e. source term for fuel equation
            virtual tmp<fvScalarMatrix> R(volScalarField& Y) const;

            //- Heat release rate calculated from fuel consumption rate matrix
            virtual tmp<volScalarField> dQ() const;

            //- Return normalised consumption rate of (fu - fres)
            virtual tmp<volScalarField> wFuelNorm() const;


    // I-O

            //- Update properties from given dictionary
            virtual bool read(const dictionary& combustionProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

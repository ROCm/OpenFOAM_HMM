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
    Foam::combustionModel

Description
    Base class for all non-premixed combustion models based on single step
    chemistry

SourceFiles
    combustionModel.C

\*---------------------------------------------------------------------------*/

#ifndef combustionModel_H
#define combustionModel_H

#include "IOdictionary.H"
#include "hsCombustionThermo.H"
#include "turbulenceModel.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class combustionModel Declaration
\*---------------------------------------------------------------------------*/

class combustionModel
{

protected:

    // Protected data

        //- Dictionary of coefficients for the particular model
        dictionary coeffs_;

        //- Reference to the thermodynamics
        hsCombustionThermo& thermo_;

        //- Reference to the turbulence model
        const compressible::turbulenceModel& turbulence_;

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Reference to mass-flux field
        const surfaceScalarField& phi_;

        //- Reference to the density field
        const volScalarField& rho_;


private:

    // Private Member Functions

        //- Disallow copy construct
        combustionModel(const combustionModel&);

        //- Disallow default bitwise assignment
        void operator=(const combustionModel&);


public:

    //- Runtime type information
    TypeName("combustionModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            combustionModel,
            dictionary,
            (
                const dictionary& combustionProperties,
                hsCombustionThermo& thermo,
                const compressible::turbulenceModel& turbulence,
                const surfaceScalarField& phi,
                const volScalarField& rho
            ),
            (
                combustionProperties,
                thermo,
                turbulence,
                phi,
                rho
            )
        );


    // Selectors

        //- Return a reference to the selected combustion model
        static autoPtr<combustionModel> New
        (
            const dictionary& combustionProperties,
            hsCombustionThermo& thermo,
            const compressible::turbulenceModel& turbulence,
            const surfaceScalarField& phi,
            const volScalarField& rho
        );


    // Constructors

        //- Construct null from components
        combustionModel
        (
            const dictionary& combustionProps,
            hsCombustionThermo& thermo,
            const compressible::turbulenceModel& turbulence,
            const surfaceScalarField& phi,
            const volScalarField& rho
        );

        //- Construct from components
        combustionModel
        (
            const word& modelType,
            const dictionary& combustionProperties,
            hsCombustionThermo& thermo,
            const compressible::turbulenceModel& turbulence,
            const surfaceScalarField& phi,
            const volScalarField& rho
        );


    //- Destructor
    virtual ~combustionModel();


    // Member Functions

        // Access functions

            //- Access combustion dictionary
            const dictionary coeffs() const
            {
                return coeffs_;
            }


    // Evolution

        //- Correct combustion rate
        virtual void correct();

        //- Fuel consumption rate matrix, i.e. source term for fuel equation
        virtual tmp<fvScalarMatrix> R(volScalarField& Y) const;

        //- Heat release rate calculated from fuel consumption rate matrix
        virtual tmp<volScalarField> dQ() const;

        //- Return normalised consumption rate of (fu - fres)
        virtual tmp<Foam::volScalarField> wFuelNorm() const;


    // I-O

        //- Update properties from given dictionary
        virtual bool read(const dictionary& combustionProps);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

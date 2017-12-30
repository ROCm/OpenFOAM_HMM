/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
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
    Foam::functionObjects::energyTransport

Group
    grpSolversFunctionObjects

Description
    Evolves a simplified energy transport equation for incompressible flows.
    It takes into account the inertia, conduction and convection terms plus
    a source.

    - The field name must be temperature and its BC's specified in the time
      directory.
    - The turbulence model should be incompressible
    - In order to use in a incompressible multi phase a list of thermal
      properties are needed. See bellow


Usage
    Example of function object specification to solve a energy transport
    equation for a single phase flow plus a source term
    \verbatim
    functions
    {
        energy
        {
            type            energyTransport;
            libs             ("libenergyTransportFunctionObjects.so");

            enabled         true;
            writeControl    outputTime;
            writeInterval   1;

            field           T;

            // volumetric Flux
            phi             phi;

            // Thermal properties
            Cp              Cp    [J/kg/K]  1e3;
            kappa           kappa [W/m/K]   0.0257;
            rhoInf          rho   [kg/m^3]  1.2;

            write           true;

            fvOptions
            {
                viscousDissipation
                {
                    type            viscousDissipation;
                    enabled         true;

                    viscousDissipationCoeffs
                    {
                        fields          (T);
                        rhoInf          $....rhoInf;
                    }
                }
            }
        }
    }
    \endverbatim

    Example of function object specification to solve a energy transport
    equation for a multiphase phase flow plus a source term

    equation:
    \verbatim
    functions
    {
        energy
        {
            type            energyTransport;
            libs            ("libenergyTransportFunctionObjects.so");

            enabled         true;
            writeControl    outputTime;
            writeInterval   1;

            field           T;

            // rho field name
            rho             rho;
            // mass flux for multiphase
            phi             rhoPhi;

            write           true;

            // Thermal properties of the phases
            phaseThermos
            {
                alpha.air
                {
                    Cp          1e3;
                    kappa       0.0243;
                }
                alpha.mercury
                {
                    Cp          140;
                    kappa       8.2;
                }
                alpha.oil
                {
                    Cp          2e3;
                    kappa       0.2;
                }
                alpha.water
                {
                    Cp          4e3;
                    kappa       0.6;
                }
            }


            fvOptions
            {
                viscousDissipation
                {
                    type            viscousDissipation;
                    enabled         true;

                    viscousDissipationCoeffs
                    {
                        fields          (T);
                        rho             rho; //rho Field
                    }
                }
            }
        }
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property     | Description             | Required    | Default value
        type         | Type name: energyTransport | yes      |
        field        | Name of the scalar field | no         | T
        phi          | Name of flux field      | no          | phi
        rho          | Name of density field   | no          | rho
        nCorr        | Number of correctors    | no          | 0
        schemesField | Name of field to specify schemes | no | field name
        fvOptions    | List of scalar sources  | no          |
        Cp           | Heat capacity for single phase | no   | 0
        rhoInf       | Density for single phase | no         | 0
        kappa        | Thermal conductivity for single phase | no   | 0
        Prt          | Turbulent Prandt number | no          | 1.0
        phaseThermos | Dictionary for multi-phase thermo |no | null
        fvOptions    | Opotional extra sources | no          | null
    \endtable

See also
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    energyTransport.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_energyTransport_H
#define functionObjects_energyTransport_H

#include "fvMeshFunctionObject.H"
#include "volFields.H"
#include "fvOptionList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class energyTransport Declaration
\*---------------------------------------------------------------------------*/

class energyTransport
:
    public fvMeshFunctionObject
{
    // Private data

        //- Name of the transport field.
        word fieldName_;

        //- Name of flux field
        word phiName_;

        //- Name of density field
        word rhoName_;

        //- Number of corrector iterations (optional)
        label nCorr_;

        //- Name of field whose schemes are used (optional)
        word schemesField_;

        //- Run-time selectable finite volume options, e.g. sources, constraints
        fv::optionList fvOptions_;

        //- Dictionary for multiphase thermos
        dictionary multiphaseThermo_;

        //- List of phase names
        wordList phaseNames_;

        //- List of phase heat capacities
        PtrList<dimensionedScalar> Cps_;

        //- List of phase thermal diffusivity for temperature [J/m/s/K]
        PtrList<dimensionedScalar> kappas_;

        //- Unallocated phase list
        UPtrList<volScalarField> phases_;

        //- Heat capacity for single phase flows
        dimensionedScalar Cp_;

        //- Thermal diffusivity for temperature for single phase flows
        dimensionedScalar kappa_;

        //- Density for single phase flows
        dimensionedScalar rho_;

        //- Turbulent Prandt number
        dimensionedScalar Prt_;

        //- rhoCp
        volScalarField rhoCp_;


    // Private Member Functions

        //- Return reference to registered transported field
        volScalarField& transportedField();

        //- Return the diffusivity field
        tmp<volScalarField> kappaEff() const;

        //- Return rho field
        tmp<volScalarField> rho() const;

        //- Return Cp
        tmp<volScalarField> Cp() const;

        //- Return kappa
        tmp<volScalarField> kappa() const;

        //- Disallow default bitwise copy construct
        energyTransport(const energyTransport&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const energyTransport&) = delete;


public:

    //- Runtime type information
    TypeName("energyTransport");


    // Constructors

        //- Construct from Time and dictionary
        energyTransport
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~energyTransport();


    // Member Functions

        //- Read the energyTransport data
        virtual bool read(const dictionary&);

        //- Calculate the energyTransport
        virtual bool execute();

        //- Do nothing.
        //  The volScalarField is registered and written automatically
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
